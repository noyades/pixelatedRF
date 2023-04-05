# -*- coding: utf-8 -*-
import os
import time, sys
import random
import numpy as np
import cmath as cm
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
libName = 'dbsEmSimDiplexer' # This is the name of the ADS workspace that was created to run sims in
sim = True # This controls whether a simulation is run or not.
view = False # This controls if the GDS is viewed after each creation 
write = True # Control whether output files are written or not
ports = 3 # For now, code makes either 2, 3 or 4 ports
sides = 2 # For now, code can put ports on 2, 3 or 4 sides, with constraints that are spelled out in rrfc
corner = 'normal' # Options are 'overlap' or 'noverlap'=non-overlap. Use overlap to predict what will happen
                  # in corners due to manufacturing, or noverlap to guarantee non-overlap at corners
pixelSize = 9 # the size of the randomized pixel in mils. Typically contrained by a PCB manufacturer.
minPixel = 6 # the size of the minimum possible pixel in mils. Typically contrained by a PCB manufacturer.
layoutRes = 1 # If wanted a sub-pixel grid, a factor greater than 1 can be set here.
layoutUnit = 25.4e-6/layoutRes # Set the layout unit to mils
pathName = '/home/jswalling/pythonWork/rrfc/dbsDiplexer/' # Base path for file creation

w_l, l_l = microstrip_calc.synthMicrostrip(sub1, z0, 30);

designName = "threePort_random_pixelSize=" + str(pixelSize)
outFile = pathName + 'data/' + designName + "_init"

os.chdir(pathName)

connectMap = 0#[0, 0, 0, 0, 0, 0]
sym = 'xy-axis'
x = datetime.now()
rrfc1 = rrfc(unit=layoutUnit,ports=ports,sides=sides,corner=corner,connect=connectMap,\
             minPix=minPixel,pixelSize=pixelSize,layoutRes=layoutRes,scale='',launchLen=30,seed=x,sim=simulator,view=view,\
             write=write,outF=outFile,sym=sym,portPosition='')
## Reference Design 
## Uncomment the following line if you have a pixelMap you would like to start with.

refPixelSize = 9 #Starting Pixelsize (If your optimization started with larger pixels and is being split to smaller pixels
csv_file = '/home/jswalling/pythonWork/rrfc/dbsDiplexer/data/threePort_random_pixelSize=' + str(refPixelSize) + '_start.csv'
gds_file = '/home/jswalling/pythonWork/rrfc/dbsDiplexer/data/threePort_random_pixelSize=' + str(pixelSize) + '_init.gds'

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
    for p in range(0,int(pixelSize/refPixelSize)):
      a[:,:,p] = refPixMap[p::int(pixelSize/refPixelSize),p::int(pixelSize/refPixelSize)]
    reSizePixMap = np.mean(a, axis=2)
    reSizePixMap = np.where(reSizePixMap >= 0.5, 1, 0)
    print(reSizePixMap)
    yProto = np.size(reSizePixMap,0)
    xProto = np.size(reSizePixMap,1)
    print(xProto,yProto)
    launch_l_pixels = int(launch_pixels*refPixelSize/pixelSize)
    print('ttt' + str(launch_l_pixels))
    xProto = xProto - 2*launch_l_pixels #temporarily adjust for launch will add patch 
    portPosition, _, _, _, _, _, _ = randomGDS_dim(sub1, rrfc1, xProto*pixelSize, yProto*pixelSize, z0)
else:
  reSizePixMap = False
  launch_l_pixels = launch_pixels
  xProto = 40 - 2*launch_l_pixels #temporarily adjust for launch will add patch 
  yProto = 27
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
simulPositions = 1 # This sets how many pixels will be flipped on each iteration in the DBS algorithm. Must bean integer >= 1
pixels = int(xProto * yProto)
rw1_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw2_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw3_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw4_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw5_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw6_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw7_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))

# Setup frequency vector for interpolation
spIntFreq = np.arange(0,10.005e9,5e6) # frequency vector to interpolate data to
# ideally will ensure that fs and fe land on grid
spIntFreq = spIntFreq.tolist()

# Starting and ending index of the frequency vector
startIndex = spIntFreq.index(0)
endIndex = spIntFreq.index(10e9)

# Define UNII (1-4) Passband
fc1 = 5.535e9
fbw1 = 780e6
fp1s = fc1 - fbw1/2
fp1e = fc1 + fbw1/2
fp1s_index = spIntFreq.index(fp1s) 
fp1e_index = spIntFreq.index(fp1e)

# Define UNII (5-8) Passband
fc2 = 6.525e9
fbw2 = 1200e6
fp2s = fc2 - fbw2/2
fp2e = fc2 + fbw2/2
fp2s_index = spIntFreq.index(fp2s)
fp2e_index = spIntFreq.index(fp2e)

# Define extra stop frequency
fc3 = 4.2e9
fbw3 = 400e6
fp3s = fc3 - fbw3/2
fp3e = fc3 + fbw3/2
fp3s_index = spIntFreq.index(fp3s)
fp3e_index = spIntFreq.index(fp3e)

# Define low frequency stop-band for UNII (1-4)
deltaBand = 250e6
fp1ls_index = startIndex
fp1le = fp1s - deltaBand
fp1le_index = spIntFreq.index(fp1le) 

# Define high frequency stop-band for UNII (1-4)
fp1he_index = endIndex
fp1hs = fp1e + deltaBand
fp1hs_index = spIntFreq.index(fp1hs)

# Define low frequency stop-band for UNII (5-8)
fp2ls_index = startIndex
fp2le = fp2s - deltaBand
fp2le_index = spIntFreq.index(fp2le) 

# Define high frequency stop-band for UNII (5-8)
fp2he_index = endIndex
fp2hs = fp2e + deltaBand
fp2hs_index = spIntFreq.index(fp2hs)

def cost_function(pixMap):
  pOps = pixOps(pixelSize = pixelSize)
  newPixMap = pixOps.updatePixels(pOps,pixMap)
  dbsPixMap = newPixMap.reshape(rows,cols-2*launch_l_pixels)
  if np.any(reSizePixMap) != 0:
    dbsPixMap = np.concatenate((reSizePixMap[:,0:launch_l_pixels],dbsPixMap,reSizePixMap[:,cols-launch_l_pixels:cols]),1)
  else:
    dbsPixMap = np.concatenate((initPixMap[:,0:launch_l_pixels],dbsPixMap,initPixMap[:,cols-launch_l_pixels:cols]),1)

  #print(dbsPixMap)
  global list_counter

  outFile = pathName + 'data/' + designName + '_' + str(list_counter)
  csv_file2 = outFile + ".csv"
  # Export Pixel Map file
  np.savetxt(csv_file2, dbsPixMap, fmt = '%d', delimiter = ",")  

  #rfc2 = rfc(unit=layoutUnit,pixelSize=pixelSize,sim=simulator,\
  #      view=False,write=True,outF=csv_file2)
  rfc2 = rrfc(unit=layoutUnit,ports=ports,sides=sides,corner=corner,connect=connectMap,\
              minPix=minPixel,pixelSize=pixelSize,layoutRes=layoutRes,scale='',launchLen=30,seed=x,sim=simulator,\
              view=False,write=True,outF=csv_file2,sym=sym,portPosition='')
  cell = recreateGDS_file(rfc2,_,_)
  
  gds_file2 = gds_file.replace('_init.gds', '_' + str(list_counter) + '.gds')
  
  em1 = emSim(workingPath = pathName, adsLibName = libName, gdsFile = gds_file2,\
              csvFile = csv_file2, numPorts = ports, portPositions = portPosition,\
              gdsCellName = cell, dataFile = outFile)
  emSim.momRun(em1)

  # Read the init s-parameter data from the file and pixelMap
  dataPath = pathName + 'data/spfiles/afs/' + designName + '_' + str(list_counter) + '.afs'
  initFreq, initS = readCiti(dataPath)
  s11interp = interp1d(initFreq,initS[:,0])
  s12interp = interp1d(initFreq,initS[:,1])
  s13interp = interp1d(initFreq,initS[:,2])
  s21interp = interp1d(initFreq,initS[:,3])
  s22interp = interp1d(initFreq,initS[:,4])
  s23interp = interp1d(initFreq,initS[:,5])
  s31interp = interp1d(initFreq,initS[:,6])
  s32interp = interp1d(initFreq,initS[:,7])
  s33interp = interp1d(initFreq,initS[:,8])
  reS11 = s11interp(spIntFreq)
  reS12 = s12interp(spIntFreq)
  reS13 = s13interp(spIntFreq)
  reS21 = s21interp(spIntFreq)
  reS22 = s22interp(spIntFreq)
  reS23 = s23interp(spIntFreq)
  reS31 = s31interp(spIntFreq)
  reS32 = s32interp(spIntFreq)
  reS33 = s33interp(spIntFreq)

  aS11 = abs(reS11)
  pS11 = np.angle(reS11, deg=True)
  aS12 = abs(reS12)
  pS12 = np.angle(reS12, deg=True)
  aS13 = abs(reS13)
  pS13 = np.angle(reS13, deg=True)
  aS21 = abs(reS21)
  pS21 = np.angle(reS21, deg=True)
  aS22 = abs(reS22)
  pS22 = np.angle(reS22, deg=True)
  aS23 = abs(reS23)
  pS23 = np.angle(reS23, deg=True)
  aS31 = abs(reS31)
  pS31 = np.angle(reS31, deg=True)
  aS32 = abs(reS32)
  pS32 = np.angle(reS32, deg=True)
  aS33 = abs(reS33)
  pS33 = np.angle(reS33, deg=True)
  
  ilMax = 10**(-3/20)
  rmse = np.zeros((fp1e_index-fp1s_index),dtype=float)
  rmse[aS32[fp1s_index:fp1e_index] < ilMax] = 1
  rmse = np.sqrt(np.sum( ((1/(fp1e_index-fp1s_index))*rmse*(ilMax - aS32[fp1s_index:fp1e_index]))**2,0 ))
  rw1 = rmse

  rmse = np.zeros((fp2e_index-fp2s_index),dtype=float)
  rmse[aS31[fp2s_index:fp2e_index] < ilMax] = 1
  rmse = np.sqrt(np.sum( ((1/(fp2e_index-fp2s_index))*rmse*(ilMax - aS31[fp2s_index:fp2e_index]))**2,0 ))
  rw2 = rmse

  #Error for Low Frequency Stop-band for UNII (1-4) Filter
  ilMin = 10**(-40/20)
  rmse = np.zeros((fp1le_index-fp1ls_index),dtype=float)
  rmse[aS32[fp1ls_index:fp1le_index] > ilMin] = 1
  rmse = np.sqrt(np.sum( ((1/(fp1le_index-fp1ls_index))*rmse*(aS32[fp1ls_index:fp1le_index]-ilMin))**2,0 ))
  rw3 = rmse

  #Error for High Frequency Stop-band for UNII (1-4) Filter
  ilMin = 10**(-40/20)
  rmse = np.zeros((fp1he_index-fp1hs_index),dtype=float)
  rmse[aS32[fp1hs_index:fp1he_index] > ilMin] = 1
  rmse = np.sqrt(np.sum( ((1/(fp1he_index-fp1hs_index))*rmse*(aS32[fp1hs_index:fp1he_index]-ilMin))**2,0 ))
  rw4 = rmse

  #Error for Low Frequency Stop-band for UNII (5-8) Filter
  ilMin = 10**(-40/20)
  rmse = np.zeros((fp2le_index-fp2ls_index),dtype=float)
  rmse[aS31[fp2ls_index:fp2le_index] > ilMin] = 1
  rmse = np.sqrt(np.sum( ((1/(fp2le_index-fp2ls_index))*rmse*(aS31[fp2ls_index:fp2le_index]-ilMin))**2,0 ))
  rw5 = rmse

  #Error for High Frequency Stop-band for UNII (5-8) Filter
  ilMin = 10**(-40/20)
  rmse = np.zeros((fp2he_index-fp2hs_index),dtype=float)
  rmse[aS31[fp2hs_index:fp1he_index] > ilMin] = 1
  rmse = np.sqrt(np.sum( ((1/(fp1he_index-fp2hs_index))*rmse*(aS31[fp2hs_index:fp2he_index]-ilMin))**2,0 ))
  rw6 = rmse

  #Error for Problem band for UNII (5-8) Filter
  ilMin = 10**(-40/20)
  rmse = np.zeros((fp3e_index-fp3s_index),dtype=float)
  rmse[aS31[fp3s_index:fp3e_index] > ilMin] = 1
  rmse = np.sqrt(np.sum( ((1/(fp3e_index-fp3s_index))*rmse*(aS31[fp3s_index:fp3e_index]-ilMin))**2,0 ))
  rw7 = rmse
  
  """
  ilMin = 10**(-10/20)
  rmse = np.zeros((endIndex-startIndex),dtype=float)
  rmse[aS21[startIndex:endIndex] > ilMin] = 1
  rmse = np.sqrt(np.sum( ((1/(endIndex-startIndex))*rmse*(aS21[startIndex:endIndex]-ilMin))**2,0 ))
  rw7 = rmse
  """

  rw1_list[list_counter] = rw1
  rw2_list[list_counter] = rw2
  rw3_list[list_counter] = rw3
  rw4_list[list_counter] = rw4
  rw5_list[list_counter] = rw5
  rw6_list[list_counter] = rw6
  rw7_list[list_counter] = rw7
  list_counter += 1
  #fomUpdate = "Current RMSE Errors:" + "\nIL UNII(1-4):" + str(rw1) + "\nIL UNII(5-8):" + str(rw2) + \
  #                                       "\nATT UNII(1-4):" + str(rw3+rw4) + "\nATT UNII(5-8):" + str(rw5+rw6) + \
  #                                       "\nInput Port Isolation (S12):" + str(rw7) + "\nTotal Error:" + \
  #                                       str(rw1 + rw2 + rw3 + rw4 + rw5 + rw6 + rw7) + "\n"
  fomUpdate = "Current RMSE Errors:" + "\nIL UNII(1-4):" + str(rw1) + "\nIL UNII(5-8):" + str(rw2) + \
                                         "\nATT UNII(1-4):" + str(rw3+rw4) + "\nATT UNII(5-8):" + str(rw5+rw6) + \
                                         "\nTotal Error:" + \
                                         str(rw1 + rw2 + rw3 + rw4 + rw5 + rw6 + rw7) + "\n"
  print(fomUpdate)
  return rw1 + rw2 + rw3 + rw4 + rw5 + rw6 + rw7 #sum of errors
  #return max(rw1, rw2, rw3, rw4, rw5, rw6, rw7) #minimax

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
csvInit = pathName + 'data/' + designName + '_initialDesign.csv'
if np.any(reSizePixMap) != 0:
  np.savetxt(csvInit, np.concatenate((reSizePixMap[:,0:launch_l_pixels], DBS.best_solution.reshape(rows,cols-2*launch_l_pixels),\
     reSizePixMap[:,cols-launch_l_pixels:cols]),1), fmt = '%d', delimiter = ",")
else:
  np.savetxt(csvInit, np.concatenate((initPixMap[:,0:launch_l_pixels], DBS.best_solution.reshape(rows,cols-2*launch_l_pixels),\
     initPixMap[:,cols-launch_l_pixels:cols]),1), fmt = '%d', delimiter = ",") 

## start BPS optimization
DBS.run()
 
## save results into files
np.savetxt(pathName + "data/cg_curve", DBS.cg_curve, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt(pathName + "data/RW1", np.array(rw1_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt(pathName + "data/RW2", np.array(rw2_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt(pathName + "data/RW3", np.array(rw3_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt(pathName + "data/RW4", np.array(rw4_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt(pathName + "data/RW5", np.array(rw5_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt(pathName + "data/RW6", np.array(rw6_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt(pathName + "data/RW7", np.array(rw7_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)

## plot for the cost
x = range(0, DBS.cg_curve.shape[0])
y = DBS.cg_curve
xlabel = "Operation Times"
ylabel = "Cost (Lower Better)"
plt.figure()
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.plot(x, y)
plt.savefig("DBS", dpi=600)
plt.show()
plt.close()
plt.clf()
