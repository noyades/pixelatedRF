# -*- coding: utf-8 -*-
import os
import time, sys
import random
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.interpolate import interp1d
from vt_rrfc import *

# Load constants and design choices. Assumes 2 layer PCB with 1 oz copper and 19 mil total thickness
fc = 5.0e9 # Operating center frequency for electrical length calculations 
           # (make this smaller than or equal to the desired operating frequency
z0 = 50 # Desired characteristic impedance of launches
EL = 90 # Desired unit of electrical length in degrees
t = 2.8 # Thickness of conductor in mils
cond = 5.88e7 # Conductivity of the conductors
h = 30 # height of conductor above substrate
t_air = 2*(t+h) # thickness of the air above the conductor layer
er = 4.5 # relative permittivity of the substrate material

sub1 = microstrip_sub(t, cond, h, er, fc)

simulator = 'ADS' # This controls the simulation to be used. Right now there are two valid values 'ADS' or 'EMX'
libName = 'dbsEmSim' # This is the name of the ADS workspace that was created to run sims in
sim = True # This controls whether a simulation is run or not.
view = False # This controls if the GDS is viewed after each creation 
write = True # Control whether output files are written or not
ports = 2 # For now, code makes either 2, 3 or 4 ports
sides = 2 # For now, code can put ports on 2, 3 or 4 sides, with constraints that are spelled out in rrfc
corner = 'normal' # Options are 'overlap' or 'noverlap'=non-overlap. Use overlap to predict what will happen
                  # in corners due to manufacturing, or noverlap to guarantee non-overlap at corners
pixelSize = 18 # the size of the randomized pixel in mils. Typically contrained by a PCB manufacturer.
minPixel = 6 # the size of the minimum possible pixel in mils. Typically contrained by a PCB manufacturer.
layoutRes = 1 # If wanted a sub-pixel grid, a factor greater than 1 can be set here.
layoutUnit = 25.4e-6/layoutRes # Set the layout unit to mils
pathName = '/home/jswalling/pythonWork/rrfc/dbs/' # Base path for file creation

w_l, l_l = microstrip_calc.synthMicrostrip(sub1, z0, 30);

designName = "twoPort_random_pixelSize=" + str(pixelSize)
outFile = pathName + 'data/' + designName + "_init"

os.chdir(pathName)
connectMap = [0, 0, 0, 0, 0, 0]
sym = 'xy-axis'
x = datetime.now()
rrfc1 = rrfc(unit=layoutUnit,ports=ports,sides=sides,corner=corner,connect=connectMap,\
          minPix=minPixel,pixelSize=pixelSize,layoutRes=layoutRes,scale='',launchLen=30,seed=x,\
          sim=simulator,view=view,write=write,outF=outFile,sym=sym,portPosition='')

## Reference Design 
## Uncomment the following line if you have a pixelMap you would like to start with.
refPixelSize = 18
csv_file = '/home/jswalling/pythonWork/rrfc/dbs/data/twoPort_random_pixelSize=' + str(refPixelSize) + '_start.csv'
gds_file = '/home/jswalling/pythonWork/rrfc/dbs/data/twoPort_random_pixelSize=' + str(pixelSize) + '_init.gds'

launch_pixels = round(l_l/refPixelSize)

if csv_file != 0:
  refPixMap = np.loadtxt(csv_file, delimiter=',')
  if refPixelSize >= pixelSize:
    reSizePixMap = np.repeat(refPixMap,int(refPixelSize/pixelSize),axis=0)
    reSizePixMap = np.repeat(reSizePixMap,int(refPixelSize/pixelSize),axis=1)
    yProto = np.size(reSizePixMap,0)
    xProto = np.size(reSizePixMap,1)
    print(xProto,yProto)
    launch_l_pixels = launch_pixels*int(refPixelSize/pixelSize)
    xProto = xProto - 2*launch_l_pixels #temporarily adjust for launch will add patch 
    portPosition, _, _, _, _, _, _ = randomGDS_dim(sub1, rrfc1, xProto*pixelSize, yProto*pixelSize, z0)
  else:
    a = np.zeros((int(np.size(refPixMap,0)/int(pixelSize/refPixelSize)),int(np.size(refPixMap,1)/int(pixelSize/refPixelSize)),int(pixelSize/refPixelSize)),dtype=int)
    for p in range(1,int(pixelSize/refPixelSize)):
      a[:,:,p] = refPixMap[p::int(pixelSize/refPixelSize),p::int(pixelSize/refPixelSize)]
    reSizePixMap = np.mean(a, axis=2)
    reSizePixMap = np.where(reSizePixMap >= 0.5, 1, 0)
    yProto = np.size(reSizePixMap,0)
    xProto = np.size(reSizePixMap,1)
    launch_l_pixels = int(launch_pixels*refPixelSize/pixelSize)
    xProto = xProto - 2*launch_l_pixels #temporarily adjust for launch will add patch 
    portPosition, _, _, _, _, _, _ = randomGDS_dim(sub1, rrfc1, xProto*pixelSize, yProto*pixelSize, z0)
else:
  reSizePixMap = False
  launch_l_pixels = launch_pixels
  xProto = 37*int(refPixelSize/pixelSize) - 2*launch_l_pixels #temporarily adjust for launch will add patch 
  yProto = 19*int(refPixelSize/pixelSize)
  portPosition, xBoard, yBoard, csv_file, gds_file, cell, _ = randomGDS_dim(sub1, rrfc1, xProto*pixelSize, yProto*pixelSize, z0)

## Define optimization parameters
rows = int(yProto)
cols = int(xProto)+2*launch_l_pixels
pixels = rows * (cols-2*launch_l_pixels)
list_counter = 0

# Initialize pixel map and reshape to flat array
if np.any(reSizePixMap) != 0: 
  flatInitPixMap = reSizePixMap[:,launch_l_pixels:cols-launch_l_pixels]
  flatInitPixMap = flatInitPixMap.reshape(1,rows*(cols-2*launch_l_pixels))
else:
  initPixMap = np.loadtxt(csv_file, delimiter=',')
  flatInitPixMap = initPixMap[:,launch_l_pixels:cols-launch_l_pixels]
  flatInitPixMap = flatInitPixMap.reshape(1,rows*(cols-2*launch_l_pixels))

# Define optimization parameters
max_iterations = 3
simulPositions = 1	 # This sets how many pixels will be flipped on each iteration in the DBS algorithm. Must bean integer >= 1
pixels = int(xProto * yProto)
rw1_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw2_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw3_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw4_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw5_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw6_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw7_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw8_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw9_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))

# Passband #1
fc1 = 2.44e9
fbw1 = 120e6
fp1s = fc1 - fbw1/2
fp1e = fc1 + fbw1/2
spIntFreq = np.arange(0,10.005e9,5e6) # frequency vector to interpolate data to
# ideally will ensure that fs and fe land on grid
spIntFreq = spIntFreq.tolist()

fp1s_index = spIntFreq.index(fp1s)
fp1e_index = spIntFreq.index(fp1e)

# Stopband #1
fc2 = 3.675e9
fbw2 = 1850e6    
fp2s = fc2 - fbw2/2
fp2e = fc2 + fbw2/2

fp2s_index = spIntFreq.index(fp2s)
fp2e_index = spIntFreq.index(fp2e)

# Passband #2
fc3 = 6.15e9
fbw3 = 2100e6
fp3s = fc3 - fbw3/2
fp3e = fc3 + fbw3/2

fp3s_index = spIntFreq.index(fp3s)
fp3e_index = spIntFreq.index(fp3e)

# Stopband #2
fc4 = 8.850e9
fbw4 = 2300e6
fp4s = fc4 - fbw4/2
fp4e = fc4 + fbw4/2

fp4s_index = spIntFreq.index(fp4s)
fp4e_index = spIntFreq.index(fp4e)

# Stopband #3
fc5 = 0.1e9
fbw5 = 50e6
fp5s = fc5 - fbw5/2
fp5e = fc5 + fbw5/2

fp5s_index = spIntFreq.index(fp5s)
fp5e_index = spIntFreq.index(fp5e)

def cost_function(pixMap):
  pOps = pixOps(pixelSize = pixelSize)
  newPixMap = pixOps.updatePixels(pOps,pixMap)
  dbsPixMap = newPixMap.reshape(rows,cols-2*launch_l_pixels)
  if np.any(reSizePixMap) != 0:
    dbsPixMap = np.concatenate((reSizePixMap[:,0:launch_l_pixels],dbsPixMap,reSizePixMap[:,cols-launch_l_pixels:cols]),1)
  else:
    dbsPixMap = np.concatenate((initPixMap[:,0:launch_l_pixels],dbsPixMap,initPixMap[:,cols-launch_l_pixels:cols]),1)

  global list_counter

  outFile = pathName + 'data/' + designName + '_' + str(list_counter)
  csv_file2 = outFile + ".csv"
  # Export Pixel Map file
  np.savetxt(csv_file2, dbsPixMap, fmt = '%d', delimiter = ",")  

  #rfc2 = rfc(unit=layoutUnit,pixelSize=pixelSize,sim=simulator,\
  #      view=False,write=True,outF=csv_file2)
  rfc2 = rrfc(unit=layoutUnit,ports=ports,sides=sides,corner=corner,connect=connectMap,\
          minPix=minPixel,pixelSize=pixelSize,layoutRes=layoutRes,scale='',launchLen=30,seed=x,\
          sim=simulator,view=False,write=True,outF=csv_file2,sym=sym,portPosition='')
  cell = recreateGDS_file(rfc2,_,_)
  
  gds_file2 = gds_file.replace('_init.gds', '_' + str(list_counter) + '.gds')
  print(gds_file2)
  em1 = emSim(workingPath = pathName, adsLibName = libName, gdsFile = gds_file2,\
              csvFile = csv_file2, numPorts = ports, portPositions = portPosition,\
              gdsCellName = cell, dataFile = outFile)
  emSim.momRun(em1)

  # Read the init s-parameter data from the file and pixelMap
  dataPath = pathName + 'data/spfiles/afs/' + designName + '_' + str(list_counter) + '.afs'
  initFreq, initS = readCiti(dataPath)
  reFreq = np.arange(0,10.005e9,5e6)
  s11interp = interp1d(initFreq,initS[:,0])
  s12interp = interp1d(initFreq,initS[:,1])
  s21interp = interp1d(initFreq,initS[:,2])
  s22interp = interp1d(initFreq,initS[:,3])
  reS11 = s11interp(spIntFreq)
  reS12 = s12interp(spIntFreq)
  reS21 = s21interp(spIntFreq)
  reS22 = s22interp(spIntFreq)

  aS11 = abs(reS11) 
  aS21 = abs(reS21)
  aS22 = abs(reS22)
  
  ilMax = 10**(-2/20)
  rmse = np.zeros((fp1e_index-fp1s_index),dtype=float)
  rmse[aS21[fp1s_index:fp1e_index] < ilMax] = 1
  rmse = np.sqrt(np.sum( ((1/(fp1e_index-fp1s_index))*rmse*(ilMax - aS21[fp1s_index:fp1e_index]))**2,0 ))
  rw1 = rmse

  rmse = np.zeros((fp3e_index-fp3s_index),dtype=float)
  rmse[aS21[fp3s_index:fp3e_index] < ilMax] = 1
  rmse = np.sqrt(np.sum( ((1/(fp3e_index-fp3s_index))*rmse*(ilMax - aS21[fp3s_index:fp3e_index]))**2,0 ))
  rw3 = rmse

  rlMin = 10**(-12/20)
  rmse = np.zeros((fp1e_index-fp1s_index),dtype=float)
  rmse[aS11[fp1s_index:fp1e_index] > rlMin] = 1
  rmse = np.sqrt(np.sum( ((1/(fp1e_index-fp1s_index))*rmse*(aS11[fp1s_index:fp1e_index] - rlMin))**2,0 ))
  rw4 = rmse

  rmse = np.zeros((fp1e_index-fp1s_index),dtype=float)
  rmse[aS22[fp1s_index:fp1e_index] > rlMin] = 1
  rmse = np.sqrt(np.sum( ((1/(fp1e_index-fp1s_index))*rmse*(aS22[fp1s_index:fp1e_index] - rlMin))**2,0 ))
  rw5 = rmse

  rmse = np.zeros((fp3e_index-fp3s_index),dtype=float)
  rmse[aS11[fp3s_index:fp3e_index] > rlMin] = 1
  rmse = np.sqrt(np.sum( ((1/(fp3e_index-fp3s_index))*rmse*(aS11[fp3s_index:fp3e_index] - rlMin))**2,0 ))
  rw6 = rmse

  rmse = np.zeros((fp3e_index-fp3s_index),dtype=float)
  rmse[aS22[fp3s_index:fp3e_index] > rlMin] = 1
  rmse = np.sqrt(np.sum( ((1/(fp3e_index-fp3s_index))*rmse*(aS22[fp3s_index:fp3e_index] - rlMin))**2,0 ))
  rw7 = rmse

  ilMin = 10**(-40/20)
  rmse = np.zeros((fp2e_index-fp2s_index),dtype=float)
  rmse[aS21[fp2s_index:fp2e_index] > ilMin] = 1
  rmse = np.sqrt(np.sum( ((1/(fp2e_index-fp2s_index))*rmse*(aS21[fp2s_index:fp2e_index] - ilMin))**2,0 ))
  rw2 = rmse

  rmse = np.zeros((fp4e_index-fp4s_index),dtype=float)
  rmse[aS21[fp4s_index:fp4e_index] > ilMin] = 1
  rmse = np.sqrt(np.sum( ((1/(fp4e_index-fp4s_index))*rmse*(aS21[fp4s_index:fp4e_index] - ilMin))**2,0 ))
  rw8 = rmse

  rmse = np.zeros((fp5e_index-fp5s_index),dtype=float)
  rmse[aS21[fp5s_index:fp5e_index] > ilMin] = 1
  rmse = np.sqrt(np.sum( ((1/(fp5e_index-fp5s_index))*rmse*(aS21[fp5s_index:fp5e_index] - ilMin))**2,0 ))
  rw9 = rmse

  rw1_list[list_counter] = rw1
  rw2_list[list_counter] = rw2
  rw3_list[list_counter] = rw3
  rw4_list[list_counter] = rw4
  rw5_list[list_counter] = rw5
  rw6_list[list_counter] = rw6
  rw7_list[list_counter] = rw7
  rw8_list[list_counter] = rw8
  rw9_list[list_counter] = rw9
  list_counter += 1
  fomUpdate = "Current RMSE Errors:" + "\nIL (ISM):" + str(rw1) + "\nIL (UNII):" + str(rw3) + \
                                       "\nATT (LOW):" + str(rw9) + "\nATT (MID):" + str(rw2) + \
                                       "\nATT (HIGH):" + str(rw8) + "\nRL (ISM):" + str(rw4+rw6) + \
                                       "\nRL (UNII):" + str(rw5+rw7) + "\nTotal Error:" + \
                                       str(rw1 + rw2 + rw3 + rw4 + rw5 + rw6 + rw7 + rw8 + rw9) + "\n"
  #fomUpdate = "Current RMSE Errors:" + "\nIL (ISM):" + str(rw1) + "\nIL (UNII):" + str(rw3) + \
  #                                     "\nATT (LOW):" + str(rw9) + "\nATT (MID):" + str(rw2) + \
  #                                     "\nATT (HIGH):" + str(rw8) + "\nTotal Error:" + \
  #                                     str(rw1 + rw2 + rw3 + rw8 + rw9) + "\n"
  
  print(fomUpdate)
  return rw1 + rw2 + rw3 + rw4 + rw5 + rw6 + rw7 + rw8 + rw9
  #return rw1 + rw2 + rw3 + rw8 + rw9
  #return max(rw1, rw2, rw3, rw4, rw5, rw6, rw7, rw8, rw9) #minimax


def call_back():
  print("Size of remained:", DBS.get_remained_size())
  print("Number of iteration:", DBS.get_iteration_number())
  print("The minimum cost:", DBS.get_cost())
  print("Best Solution:\n", DBS.get_best_solution())
  ## save results into files
  np.savetxt(pathName + "data/cg_curve", DBS.cg_curve, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
  csvFinal = pathName + 'data/' + designName + '_finalDesign.csv'
  if np.any(reSizePixMap) != 0:
    np.savetxt(csvFinal, np.concatenate((reSizePixMap[:,0:launch_l_pixels], DBS.best_solution.reshape(rows,cols-2*launch_l_pixels),\
       reSizePixMap[:,cols-launch_l_pixels:cols]),1), fmt = '%d', delimiter = ",")
  else:
    np.savetxt(csvFinal, np.concatenate((initPixMap[:,0:launch_l_pixels], DBS.best_solution.reshape(rows,cols-2*launch_l_pixels),\
       initPixMap[:,cols-launch_l_pixels:cols]),1), fmt = '%d', delimiter = ",") 

DBS = dbsAlgo(pixels, cost_function, rows, 'none', simulPositions, max_iterations, call_back, initial_solution=flatInitPixMap)

## save initial solution into the file
np.save("init_solution", DBS.best_solution)

## save initial solution into the file
csvInit = pathName + 'data/' + designName + '_initialDesign.csv'
if np.any(reSizePixMap) != 0:
  np.savetxt(csvInit, np.concatenate((reSizePixMap[:,0:launch_l_pixels], DBS.best_solution.reshape(rows,cols-2*launch_l_pixels),\
     reSizePixMap[:,cols-launch_l_pixels:cols]),1), fmt = '%d', delimiter = ",")
else:
  np.savetxt(csvInit, np.concatenate((initPixMap[:,0:launch_l_pixels], DBS.best_solution.reshape(rows,cols-2*launch_l_pixels),\
     initPixMap[:,cols-launch_l_pixels:cols]),1), fmt = '%d', delimiter = ",") 

## start DBS optimization
DBS.run()
 
## save results into files
np.save("solution", DBS.best_solution.reshape(rows,cols-2*launch_l_pixels))
np.savetxt(pathName + "data/cg_curve", DBS.cg_curve, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt(pathName + "data/RW1", np.array(rw1_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt(pathName + "data/RW2", np.array(rw2_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt(pathName + "data/RW3", np.array(rw3_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
#np.savetxt(pathName + "data/RW4", np.array(rw4_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
#np.savetxt(pathName + "data/RW5", np.array(rw4_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
#np.savetxt(pathName + "data/RW6", np.array(rw4_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
#np.savetxt(pathName + "data/RW7", np.array(rw4_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt(pathName + "data/RW8", np.array(rw4_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt(pathName + "data/RW9", np.array(rw4_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)

## plot for the cost
x = range(0, DBS.cg_curve.shape[0])
y = DBS.cg_curve
xlabel = "Operation Times"
ylabel = "Cost (Lower Better)"
plt.figure()
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.plot(x, y)
plt.savefig("PBS", dpi=600)
plt.show()
plt.close()
plt.clf()
