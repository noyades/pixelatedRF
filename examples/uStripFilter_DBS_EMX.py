# -*- coding: utf-8 -*-
import skrf as rf
import os
import time, sys
import random
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.interpolate import interp1d
from vt_rrfc import *

# Load constants and design choices. Assumes 2 layer PCB with 1 oz copper and 19 mil total thickness
fc = 25.0e9 # Operating center frequency for electrical length calculations 
           # (make this smaller than or equal to the desired operating frequency
z0 = 50 # Desired characteristic impedance of launches
EL = 90 # Desired unit of electrical length in degrees
t = 0.129921 # Thickness of conductor (OI, 22nm) in mils
cond = 5.71e7 # Conductivity of the conductors (S/m)
h = 0.143511 # height of conductor above substrate (assumes dielectric between M1 and OI in 22nm Stack)
t_air = 2*(t+h) # thickness of the air above the conductor layer
er = 3.93 # relative permittivity of the substrate material (assumes dielectric between M1 and OI in 22nm Stack)

sub1 = microstrip_sub(t, cond, h, er, fc)

simulator = 'EMX' # This controls the simulation to be used. Right now there are two valid values 'ADS' or 'EMX'
sim = True # This controls whether a simulation is run or not.
view = False # This controls if the GDS is viewed after each creation 
write = True # Control whether output files are written or not
ports = 2 # For now, code makes either 2, 3 or 4 ports
sides = 2 # For now, code can put ports on 2, 3 or 4 sides, with constraints that are spelled out in rrfc
corner = 'standard' # Options are 'overlap' or 'noverlap'=non-overlap. Use overlap to predict what will happen
                  # in corners due to manufacturing, or noverlap to guarantee non-overlap at corners
pixel = 6
scale = 10
pixelSize = scale*pixel # the size of the randomized pixel in mils. Typically contrained by a PCB manufacturer.
minPixel = scale*1.0 # the size of the minimum possible pixel in mils. Typically contrained by a PCB manufacturer.
layoutRes = scale*25.4 # If wanted a sub-pixel grid, a factor greater than 1 can be set here.
layoutUnit = 25.4e-6/layoutRes # Set the layout unit to mils
procFile = '/software/RFIC/PDK/globalFoundries/22FDX-EXT/release/Emagnetic/EMX/9M_2Mx_5Cx_1Jx_1Ox_LBthick/22fdsoi_9M_2Mx_5Cx_1Jx_1Ox_LBthick_nominal_detailed.encrypted.proc'
pathName = '/home/jswalling/pythonWork/rrfc/dbsEmx/' # Base path for file creation

w_l, l_l = microstrip_calc.synthMicrostrip(sub1, z0, 2)
print('quarter wave line=',25.4*w_l,25.4*l_l)
w_90, l_90 = microstrip_calc.synthMicrostrip(sub1, z0, 30)

designName = "twoPort_random_pixelSize=" + str(pixelSize)
outFile = pathName + 'data/' + designName + "_init"

os.chdir(pathName)
connectMap = [0, 0, 0, 0, 0, 0]
sym = 'asym'
x = datetime.now()
rrfc1 = rrfc(unit=layoutUnit,ports=ports,sides=sides,corner=corner,connect=connectMap,\
          minPix=minPixel,pixelSize=pixelSize,layoutRes=layoutRes,scale=scale,launchLen=2,seed=x,\
          sim=simulator,view=view,write=write,outF=outFile,sym=sym,portPosition='')

## Reference Design 
## Uncomment the following line if you have a pixelMap you would like to start with.
refPixelSize = scale*pixel
csv_file = '/home/jswalling/pythonWork/rrfc/dbsEmx/data/twoPort_random_pixelSize=' + str(refPixelSize) + '_start.csv'
gds_file = '/home/jswalling/pythonWork/rrfc/dbsEmx/data/twoPort_random_pixelSize=' + str(pixelSize) + '_init.gds'

launch_pixels = round(25.4*l_l/pixel)
if csv_file != 0:
  refPixMap = np.loadtxt(csv_file, delimiter=',')
  reSizePixMap = np.repeat(refPixMap,int(refPixelSize/pixelSize),axis=0)
  reSizePixMap = np.repeat(reSizePixMap,int(refPixelSize/pixelSize),axis=1)
  yProto = np.size(reSizePixMap,0)
  xProto = np.size(reSizePixMap,1)
  print(xProto,yProto)
  launch_l_pixels = launch_pixels*int(refPixelSize/pixelSize)
  xProto = xProto - 2*launch_l_pixels #temporarily adjust for launch will add patch 
  portPosition, _, _, _, _ = randomPixMap_dim(rrfc1, sub1, xProto*pixelSize, yProto*pixelSize, z0)
else:
  reSizePixMap = False
  launch_l_pixels = launch_pixels
  xProto = int(25.4*l_90/(pixel*1.5)) - 2*launch_l_pixels #temporarily adjust for launch will add patch 
  yProto = int(25.4*l_90/(pixel*3)) #Needs to be odd number...working on bug fix
  print(xProto,yProto)
  portPosition, xBoard, yBoard, csv_file, _ = randomPixMap_dim(rrfc1, sub1, xProto*pixelSize, yProto*pixelSize, z0)

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
simulPositions = 5 # This sets how many pixels will be flipped on each iteration in the DBS algorithm. Must bean integer >= 1
pixels = int(xProto * yProto)
rw1_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw2_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw3_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw8_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw9_list = np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))

# Passband #1
fc1 = 5*2.44e9
fbw1 = 5*120e6
fp1s = fc1 - fbw1/2
fp1e = fc1 + fbw1/2
spIntFreq = np.arange(0,100.005e9,5e6) # frequency vector to interpolate data to
# ideally will ensure that fs and fe land on grid
spIntFreq = spIntFreq.tolist()

fp1s_index = spIntFreq.index(fp1s)
fp1e_index = spIntFreq.index(fp1e)

# Stopband #1
fc2 = 5*3.8e9
fbw2 = 5*1600e6    
fp2s = fc2 - fbw2/2
fp2e = fc2 + fbw2/2

fp2s_index = spIntFreq.index(fp2s)
fp2e_index = spIntFreq.index(fp2e)

# Passband #2
fc3 = 5*6.125e9
fbw3 = 5*2000e6
fp3s = fc3 - fbw3/2
fp3e = fc3 + fbw3/2

fp3s_index = spIntFreq.index(fp3s)
fp3e_index = spIntFreq.index(fp3e)

# Stopband #2
fc4 = 5*8.625e9
fbw4 = 5*2750e6
fp4s = fc4 - fbw4/2
fp4e = fc4 + fbw4/2

fp4s_index = spIntFreq.index(fp4s)
fp4e_index = spIntFreq.index(fp4e)

# Stopband #3
fc5 = 5*0.1e9
fbw5 = 5*50e6
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
          minPix=minPixel,pixelSize=pixelSize,layoutRes=layoutRes,scale=scale,launchLen=30,seed=x,\
          sim=simulator,view=view,write=True,outF=csv_file2,sym=sym,portPosition=portPosition)
  cell = createGDSFromPixMap_file(rfc2,sub1,50,'none')  

  gds_file2 = gds_file.replace('_init.gds', '_' + str(list_counter) + '.gds')
  print(gds_file2)
  em1 = emSim(workingPath = pathName, adsLibName = '', gdsFile = gds_file2,\
              csvFile = csv_file2, numPorts = ports, portPositions = portPosition,\
              gdsCellName = cell, dataFile = outFile)
  emSim.emxRun(em1,procFile)

  # Read the init s-parameter data from the file and pixelMap
  dataPath = pathName + 'data/spfiles/' + designName + '_' + str(list_counter) + '.s2p'
  se_ntwk = rf.Network(dataPath)
  initFreq = se_ntwk.frequency._f

  s11interp = interp1d(initFreq,np.transpose(se_ntwk.s11._s))
  s12interp = interp1d(initFreq,np.transpose(se_ntwk.s12._s))
  s21interp = interp1d(initFreq,np.transpose(se_ntwk.s21._s))
  s22interp = interp1d(initFreq,np.transpose(se_ntwk.s22._s))
  reS11 = s11interp(spIntFreq)
  reS12 = s12interp(spIntFreq)
  reS21 = s21interp(spIntFreq)
  reS22 = s22interp(spIntFreq)

  aS11 = abs(reS11).flatten()
  pS11 = np.angle(reS11, deg=True).flatten()
  aS12 = abs(reS12).flatten()
  pS12 = np.angle(reS12, deg=True).flatten()
  aS21 = abs(reS21).flatten()
  pS21 = np.angle(reS21, deg=True).flatten()
  aS22 = abs(reS22).flatten()
  pS22 = np.angle(reS22, deg=True).flatten()

  ilMax = 10**(-2/20)
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
  rw8_list[list_counter] = rw8
  rw9_list[list_counter] = rw9
  list_counter += 1
  fomUpdate = "Current RMSE Errors:" + "\nIL (ISM):" + str(rw1) + "\nIL (UNII):" + str(rw3) + \
                                       "\nATT (LOW):" + str(rw9) + "\nATT (MID):" + str(rw2) + \
                                       "\nATT (HIGH):" + str(rw8) + "\nTotal Error:" + \
                                       str(rw1 + rw2 + rw3 + rw8 + rw9) + "\n"
  
  print(fomUpdate)
  #return rw1 + rw2 + rw3 + rw4 + rw5 + rw6 + rw7 + rw8 + rw9
  return rw1 + rw2 + rw3 + rw8 + rw9
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
