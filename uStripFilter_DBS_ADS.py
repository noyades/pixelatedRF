# -*- coding: utf-8 -*-
import os
import time, sys
import random
import numpy as np
import matplotlib as plt
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
corner = 'overlap'
pixelSize = 9 # the size of the randomized pixel in mils. Typically contrained by a PCB manufacturer.
layoutUnit = 25.4e-6 # Set the layout unit to mils
pathName = '/home/jswalling/pythonWork/rrfc/dbs/' # Base path for file creation

# This code will generate a prototype filter using the classical design method for the stepped impedance filter. 
# The prototype has a fixed order and filter type that can be defined, and ultimately this filter will be used 
# to define the boundaries of the random structures that will be designed and simulated to train an Neural Network.
filtOrder = 3
filtType = 'butter'
max_w = 20
min_w = 1

Zo_h, ereff_h = microstrip_calc.calcMicrostrip(sub1, min_w*pixelSize, 1000);
Zo_l, ereff_l = microstrip_calc.calcMicrostrip(sub1, max_w*pixelSize, 1000);
w_h = min_w*pixelSize
w_l = max_w*pixelSize
print(Zo_h,Zo_l)

designName = "twoPort_random_pixelSize=" + str(pixelSize)
outFile = pathName + 'data/' + designName + "_init"

os.chdir(pathName)
connectMap = [0, 0, 0, 0, 0, 0]
sym = 'xy-axis'
x = 1
rrfc1 = rrfc(unit=layoutUnit,ports=ports,sides=sides,corner=corner,connect=connectMap,\
          pixelSize=pixelSize,seed=x,sim=simulator,view=view,write=write,\
          outF=outFile,sym=sym)

## Reference Design 
## Uncomment the following line if you have a pixelMap you would like to start with.
refPixelSize = 9
csv_file = '/home/jswalling/pythonWork/rrfc/dbs/data/twoPort_random_pixelSize=' + str(refPixelSize) + '_start.csv'
gds_file = '/home/jswalling/pythonWork/rrfc/dbs/data/twoPort_random_pixelSize=' + str(pixelSize) + '_init.gds'

if csv_file != 0:
  refPixMap = np.loadtxt(csv_file, delimiter=',')
  reSizePixMap = np.repeat(refPixMap,int(refPixelSize/pixelSize),axis=0)
  reSizePixMap = np.repeat(reSizePixMap,int(refPixelSize/pixelSize),axis=1)
  yProto = np.size(reSizePixMap,0)
  xProto = np.size(reSizePixMap,1)
  print(xProto,yProto)
  launch_l_pixels = 12*int(refPixelSize/pixelSize)
  xProto = xProto - 2*launch_l_pixels #temporarily adjust for launch will add patch 
  portPosition, _, _, _, _, _, _ = randomGDS_dim(sub1, rrfc1, xProto*pixelSize, yProto*pixelSize, z0)
else:
  reSizePixMap = False
  launch_l_pixels = 6*int(refPixelSize/pixelSize)
  xProto = 37*int(refPixelSize/pixelSize) - 2*launch_l_pixels #tempora#temporarily adjust for launch will add patch 
  yProto = 19*int(refPixelSize/pixelSize)
  portPosition, xBoard, yBoard, csv_file, gds_file, cell, _ = randomGDS_dim(sub1, rrfc1, xProto*pixelSize, yProto*pixelSize, z0)
"""
#designName = "twoPort_order=" + str(filtOrder) + "_type=" +\
#              filtType + "_pixelSize=" + str(pixelSize)
designName = "twoPort_random_pixelSize=" + str(pixelSize)
outFile = pathName + 'data/' + designName + "_init"
launch_l_pixels = 6*2
rfc1 = rfc(unit=layoutUnit,pixelSize=pixelSize,sim=simulator,\
        view=True,write=True,outF=outFile)
portPosition, xProto, yProto, csv_file, gds_file, cell, launch_l_pixels = uStripSteppedImpFilterGDS(sub1, rfc1, filtType, filtOrder, \
                                                     w_h, w_l, Zo_h, Zo_l, z0)

xProto = xProto - 2*launch_l_pixels*pixelSize #temporarily adjust for launch will add patch 
print(xProto, yProto)
"""

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
pixels = int(xProto * yProto)
rw1_list = np.zeros((max_iterations * pixels + 1))
rw2_list = np.zeros((max_iterations * pixels + 1))
rw3_list = np.zeros((max_iterations * pixels + 1))
rw4_list = np.zeros((max_iterations * pixels + 1))
rw5_list = np.zeros((max_iterations * pixels + 1))
rw6_list = np.zeros((max_iterations * pixels + 1))
rw7_list = np.zeros((max_iterations * pixels + 1))
rw8_list = np.zeros((max_iterations * pixels + 1))

fc1 = 2.44e9
fbw1 = 120e6
fp1s = fc1 - fbw1/2
fp1e = fc1 + fbw1/2
spIntFreq = np.arange(0,10.005e9,5e6) # frequency vector to interpolate data to
# ideally will ensure that fs and fe land on grid
spIntFreq = spIntFreq.tolist()

fp1s_index = spIntFreq.index(fp1s)
fp1e_index = spIntFreq.index(fp1e)

fc2 = 3.8e9
fbw2 = 2000e6    
fp2s = fc2 - fbw2/2
fp2e = fc2 + fbw2/2

fp2s_index = spIntFreq.index(fp2s)
fp2e_index = spIntFreq.index(fp2e)

fc3 = 6.125e9
fbw3 = 2000e6
fp3s = fc3 - fbw3/2
fp3e = fc3 + fbw3/2

fp3s_index = spIntFreq.index(fp3s)
fp3e_index = spIntFreq.index(fp3e)

fc4 = 8.625e9
fbw4 = 2750e6
fp4s = fc4 - fbw4/2
fp4e = fc4 + fbw4/2

fp4s_index = spIntFreq.index(fp4s)
fp4e_index = spIntFreq.index(fp4e)

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
          pixelSize=pixelSize,seed=x,sim=simulator,view=False,write=True,\
          outF=csv_file2,sym=sym)
  cell = recreateGDS_file(rfc2)
  
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
  
  ilMax = 10**(-1/20)
  rmse = np.zeros((fp1e_index-fp1s_index),dtype=float)
  rmse[aS21[fp1s_index:fp1e_index] < ilMax] = 1
  rmse = np.sqrt(np.sum( ((1/(fp1e_index-fp1s_index))*rmse*(ilMax - aS21[fp1s_index:fp1e_index]))**2,0 ))
  rw1 = rmse

  rmse = np.zeros((fp3e_index-fp3s_index),dtype=float)
  rmse[aS21[fp3s_index:fp3e_index] < ilMax] = 1
  rmse = np.sqrt(np.sum( ((1/(fp3e_index-fp3s_index))*rmse*(ilMax - aS21[fp3s_index:fp3e_index]))**2,0 ))
  rw3 = rmse

  ilMin = 10**(-40/20)
  rmse = np.zeros((fp2e_index-fp2s_index),dtype=float)
  rmse[aS21[fp2s_index:fp2e_index] > ilMin] = 1
  rmse = np.sqrt(np.sum( ((1/(fp2e_index-fp2s_index))*rmse*(aS21[fp2s_index:fp2e_index] - ilMin))**2,0 ))
  rw2 = rmse

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
  rmse = np.zeros((fp4e_index-fp4s_index),dtype=float)
  rmse[aS21[fp4s_index:fp4e_index] > ilMin] = 1
  rmse = np.sqrt(np.sum( ((1/(fp4e_index-fp4s_index))*rmse*(aS21[fp4s_index:fp4e_index] - ilMin))**2,0 ))
  rw8 = rmse

  rw1_list[list_counter] = rw1
  rw2_list[list_counter] = rw2
  rw3_list[list_counter] = rw3
  rw4_list[list_counter] = rw4
  rw5_list[list_counter] = rw5
  rw6_list[list_counter] = rw6
  rw7_list[list_counter] = rw7
  rw8_list[list_counter] = rw8
  list_counter += 1
  print(rw1, rw2, rw3, rw4, rw5, rw6, rw7, rw8)
  return rw1 + rw2 + rw3 + rw4 + rw5 + rw6 + rw7 + rw8


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

simulPositions = 5
DBS = dbsAlgo(pixels, cost_function, simulPositions, max_iterations, call_back, initial_solution=flatInitPixMap)

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
np.savetxt(pathName + "data/RW4", np.array(rw4_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt(pathName + "data/RW5", np.array(rw4_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt(pathName + "data/RW6", np.array(rw4_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt(pathName + "data/RW7", np.array(rw4_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
np.savetxt(pathName + "data/RW8", np.array(rw4_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)

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
