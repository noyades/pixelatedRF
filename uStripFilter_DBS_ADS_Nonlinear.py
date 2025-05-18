# -*- coding: utf-8 -*-
import os
import time, sys
import random
import numpy as np
import math
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.interpolate import interp1d
from vt_rrfc import *

# Load constants and design choices. Assumes 2 layer PCB with 1 oz copper and 19 mil total thickness
fc = 8.0e9 # Operating center frequency for electrical length calculations 
           # (make this smaller than or equal to the desired operating frequency
z0 = 50 # Desired characteristic impedance of launches
EL = 90 # Desired unit of electrical length in degrees
t = 1.4 # Thickness of conductor in mils
cond = 5.88e7 # Conductivity of the conductors
h = 30 # height of conductor above substrate
t_air = 2*(t+h) # thickness of the air above the conductor layer
er = 3.66 # relative permittivity of the substrate material

sub1 = MicrostripSub(t, cond, h, er, fc)

simulator = 'ADS' # This controls the simulation to be used. Right now there are two valid values 'ADS' or 'EMX'
libName = 'dbsEmSim' # This is the name of the ADS workspace that was created to run sims in
sim = True # This controls whether a simulation is run or not.
view = False # This controls if the GDS is viewed after each creation 
write = True # Control whether output files are written or not
ports = 2 # For now, code makes either 2, 3 or 4 ports
sides = 2 # For now, code can put ports on 2, 3 or 4 sides, with constraints that are spelled out in rrfc
corner = 'normal' # Options are 'overlap' or 'noverlap'=non-overlap. Use overlap to predict what will happen
                  # in corners due to manufacturing, or noverlap to guarantee non-overlap at corners
scale = 1
pixelSize = scale*8 # the size of the randomized pixel in mils. Typically contrained by a PCB manufacturer.
minPixel = scale*6 # the size of the minimum possible pixel in mils. Typically contrained by a PCB manufacturer.
layoutRes = scale*1 # If wanted a sub-pixel grid, a factor greater than 1 can be set here.
layoutUnit = 25.4e-6/layoutRes # Set the layout unit to mils
pathName = '/home/jswalling/pythonWork/rrfc/dbs/' # Base path for file creation

w_l, l_l = MicrostripCalc.synth_microstrip(sub1, z0, 30);

designName = "twoPort_30x30_pixelSize=" + str(pixelSize)
outFile = pathName + 'data/' + designName + "_init"

os.chdir(pathName)
connectMap = [0, 0, 0, 0, 0, 0]
sym = 'xy-axis'
x = datetime.now()
rrfc1 = RandomComponent(unit=layoutUnit,
  ports=ports,sides=sides,
  corner=corner,
  connect=connectMap,
  minPix=minPixel,
  pixelSize=pixelSize,
  layoutRes=layoutRes,
  scale=scale,
  launchLen=30,
  seed=x,
  sim=simulator,
  view=view,
  write=write,
  outF=outFile,
  sym=sym,
  portPosition='')

## Reference Design 
## Uncomment the following line if you have a pixelMap you would like to start with.
refPixelSize = 8
csv_file = 0#'/home/jswalling/pythonWork/rrfc/dbs/data/twoPort_28x40_pixelSize=' + str(refPixelSize) + '_start.csv'
gds_file = '/home/jswalling/pythonWork/rrfc/dbs/data/twoPort_28x40_pixelSize=' + str(pixelSize) + '_init.gds'

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
    portPosition, _, _, _, _, _, _ = rrfc1.random_gds_dim(sub1, xProto*pixelSize, yProto*pixelSize, z0)
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
    portPosition, _, _, _, _, _, _ = rrfc1.random_gds_dim(sub1, xProto*pixelSize, yProto*pixelSize, z0)
else:
  reSizePixMap = False
  launch_l_pixels = launch_pixels
  xProto = 48*int(refPixelSize/pixelSize) - 2*launch_l_pixels #temporarily adjust for launch will add patch 
  yProto = 29*int(refPixelSize/pixelSize)
  portPosition, xBoard, yBoard, csv_file, gds_file, cell, _ = rrfc1.random_gds_dim(sub1, xProto*pixelSize, yProto*pixelSize, z0)

## Define optimization parameters
rows = int(yProto)
cols = int(xProto)+2*launch_l_pixels
pixels = rows * (cols-2*launch_l_pixels)
list_counter = 0

print('Rows=',rows,' Cols=',cols)
# Initialize pixel map and reshape to flat array
if np.any(reSizePixMap) != 0: 
  flatInitPixMap = reSizePixMap[:,launch_l_pixels:cols-launch_l_pixels]
  flatInitPixMap = flatInitPixMap.reshape(1,rows*(cols-2*launch_l_pixels))
else:
  print(csv_file)
  initPixMap = np.loadtxt(csv_file, delimiter=',')
  print('Pixmap Shape=',initPixMap.shape)
  flatInitPixMap = initPixMap[:,launch_l_pixels:cols-launch_l_pixels]
  print('Pixmap ReShape=',flatInitPixMap.shape)
  flatInitPixMap = flatInitPixMap.reshape(1,rows*(cols-2*launch_l_pixels))

# Define optimization parameters
max_iterations = 3
simulPositions = 1	 # This sets how many pixels will be flipped on each iteration in the DBS algorithm. Must bean integer >= 1
pixels = int(xProto * yProto)
rw1_list = np.zeros(900)#np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw2_list = np.zeros(900)#np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))
rw3_list = np.zeros(900)#np.zeros(int(max_iterations * np.ceil(pixels/simulPositions) + 1))

# Define frequencies for optimization
spIntFreq = np.arange(0,20.002e9,2e6) # frequency vector to interpolate data to
# ideally will ensure that fs and fe land on grid
spIntFreq = spIntFreq.tolist()

fs1s = 3.562e9
fs1s_index = spIntFreq.index(fs1s)
fs1e = 5.624e9
fs1e_index = spIntFreq.index(fs1e)

fp1s = 7.124e9
fp1s_index = spIntFreq.index(fp1s)
fp1e = 9.124e9
fp1e_index = spIntFreq.index(fp1e)

fs2s = 10.624e9
fs2s_index = spIntFreq.index(fs2s)
fs2e = 20.00e9
fs2e_index = spIntFreq.index(fs2e)

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
  rfc2 = RandomUstrip(unit=layoutUnit,ports=ports,sides=sides,corner=corner,connect=connectMap,\
          minPix=minPixel,pixelSize=pixelSize,layoutRes=layoutRes,scale=scale,launchLen=30,seed=x,\
          sim=simulator,view=False,write=True,outF=csv_file2,sym=sym,portPosition='')
  cell = rfc2.create_gds_from_pixmap_file(sub1,z0,_)
  
  gds_file2 = gds_file.replace('_init.gds', '_' + str(list_counter) + '.gds')
  print(gds_file2)
  em1 = EmSim(workingPath = pathName, adsLibName = libName, gdsFile = gds_file2,\
              csvFile = csv_file2, numPorts = ports, portPositions = portPosition,\
              gdsCellName = cell, dataFile = outFile)
  em1.mom_run()

  # Read the init s-parameter data from the file and pixelMap
  dataPath = pathName + 'data/spfiles/afs/' + designName + '_' + str(list_counter) + '.afs'
  initFreq, initS = readCiti(dataPath)
  reFreq = np.arange(0,20.002e9,2e6)
  s11interp = interp1d(initFreq,initS[:,0])
  s12interp = interp1d(initFreq,initS[:,1])
  s21interp = interp1d(initFreq,initS[:,2])
  s22interp = interp1d(initFreq,initS[:,3])
  reS11 = s11interp(spIntFreq)
  reS12 = s12interp(spIntFreq)
  reS21 = s21interp(spIntFreq)
  reS22 = s22interp(spIntFreq)

  aS11 = 20*np.log10(abs(reS11)) 
  aS21 = 20*np.log10(abs(reS21))
  aS22 = 20*np.log10(abs(reS22))
  
  #AS11_pass = np.log10(EmSim.sigmoid(abs(aS11) - 10))
  #AS22_pass = np.log10(EmSim.sigmoid(abs(aS22) - 10))
  #AS21_pass = np.log10(EmSim.sigmoid(abs(aS21) - 0))
  #AS11_stop = np.log10(EmSim.sigmoid(abs(aS11) - 0))
  #AS22_stop = np.log10(EmSim.sigmoid(abs(aS22) - 0))
  #AS21_stop = np.log10(EmSim.sigmoid(abs(aS21) - 20))

  fs1_s11 = 20*np.log10(np.sqrt((1/(fs1e_index-fs1s_index))*np.sum((10**((aS11[fs1s_index:fs1e_index])/20))**2)))
  fs1_s21 = 20*np.log10(np.sqrt((1/(fs1e_index-fs1s_index))*np.sum((10**((aS21[fs1s_index:fs1e_index])/20))**2)))
  fs1_s22 = 20*np.log10(np.sqrt((1/(fs1e_index-fs1s_index))*np.sum((10**((aS22[fs1s_index:fs1e_index])/20))**2)))

  fp1_s11 = 20*np.log10(np.sqrt((1/(fp1e_index-fp1s_index))*np.sum((10**((aS11[fp1s_index:fp1e_index])/20))**2)))
  fp1_s21 = 20*np.log10(np.sqrt((1/(fp1e_index-fp1s_index))*np.sum((10**((aS21[fp1s_index:fp1e_index])/20))**2)))
  fp1_s22 = 20*np.log10(np.sqrt((1/(fp1e_index-fp1s_index))*np.sum((10**((aS22[fp1s_index:fp1e_index])/20))**2)))

  fs2_s11 = 20*np.log10(np.sqrt((1/(fs2e_index-fs2s_index))*np.sum((10**((aS11[fs2s_index:fs2e_index])/20))**2)))
  fs2_s21 = 20*np.log10(np.sqrt((1/(fs2e_index-fs2s_index))*np.sum((10**((aS21[fs2s_index:fs2e_index])/20))**2)))
  fs2_s22 = 20*np.log10(np.sqrt((1/(fs2e_index-fs2s_index))*np.sum((10**((aS22[fs2s_index:fs2e_index])/20))**2)))

  print('RMS S-params:\nS11(<7.124G)=' + str(fs1_s11) + \
                     '\nS21(<7.124G)=' + str(fs1_s21) + \
                     '\nS22(<7.124G)=' + str(fs1_s22) + \
                     '\nS11(7.124-9.124G)=' + str(fp1_s11) + \
                     '\nS21(7.124-9.124G)=' + str(fp1_s21) + \
                     '\nS22(7.124-9.124G)=' + str(fp1_s22) + \
                     '\nS11(>9.124G)=' + str(fs2_s11) + \
                     '\nS21(>9.124G)=' + str(fs2_s21) + \
                     '\nS22(>9.124G)=' + str(fs2_s22))

  # Passband rewards. Goal is to maximize S21 in the passband while minimizing S11 and S22. Added non-linear functions to 
  # reduce overfitting to a single parameter. The gain (S21) is weighted more heavily than the return loss (factor of 5 here)
  rw1 = abs(fp1_s11) / (1 + np.exp(-2 + abs(fp1_s21))) + abs(fp1_s22) / (1 + np.exp(-2 + abs(fp1_s21))) + 5*(10**((fp1_s21 + 2)/20) / (1 + np.exp(15-abs(fp1_s11)) + np.exp(15-abs(fp1_s22))))

  # Stopband rewards. Goal is to minimize S21 in the stopband while maximizing S11 and S22. Added non-linear functions to 
  # reduce overfitting to a single parameter. The gain (S21) is weighted more heavily than the return loss (factor of 5 here)
  rw2 = 10**((-35-(fs1_s21))/20) #/ (1 + np.exp(3-abs(aS11[fs1_index])) + np.exp(3-abs(aS22[fs1_index])))
  rw3 = 10**((-35-(fs2_s21))/20) #/ (1 + np.exp(3-abs(aS11[fs1_index])) + np.exp(3-abs(aS22[fs1_index])))


  rw1_list[list_counter] = rw1
  rw2_list[list_counter] = rw2
  rw3_list[list_counter] = rw3

  list_counter += 1
  #fomUpdate = "Current RMSE Errors:" + "\nIL (ISM):" + str(rw1) + "\nIL (UNII):" + str(rw3) + \
  #                                     "\nATT (LOW):" + str(rw9) + "\nATT (MID):" + str(rw2) + \
  #                                     "\nATT (HIGH):" + str(rw8) + "\nRL (ISM):" + str(rw4+rw6) + \
  #                                     "\nRL (UNII):" + str(rw5+rw7) + "\nTotal Error:" + \
  #                                     str(rw1 + rw2 + rw3 + rw4 + rw5 + rw6 + rw7 + rw8 + rw9) + "\n"
  fomUpdate = "Current RMSE Errors:" + "\nfp1:" + str(rw1) + "\nfp2:" + str(rw2) + "\nfs1:" + str(rw3) + "\nTotal Error:" + str(1 / (1/rw1 + 1/rw2 + 1/rw3)) + "\n"
  print(fomUpdate)
  return 1 / (1/rw1 + 1/rw2 + 1/rw3)
  #return rw1 + rw2 + rw3 + rw8 + rw9
  #return max(rw1, rw2, rw3, rw4, rw5, rw6, rw7, rw8, rw9) #minimax


def call_back():
  print("Size of remained:", DBS.get_remained_size())
  print("Number of iteration:", DBS.get_iteration_number())
  print("The minimum cost:", DBS.get_cost())
  print("Best Solution:\n", DBS.get_best_solution())
  ## save results into files
  np.savetxt(pathName + "data/cg_curve", DBS.cg_curve, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
  np.savetxt(pathName + "data/RW1", np.array(rw1_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
  np.savetxt(pathName + "data/RW2", np.array(rw2_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
  np.savetxt(pathName + "data/RW3", np.array(rw3_list), fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
  csvFinal = pathName + 'data/' + designName + '_finalDesign.csv'
  if np.any(reSizePixMap) != 0:
    np.savetxt(csvFinal, np.concatenate((reSizePixMap[:,0:launch_l_pixels], DBS.best_solution.reshape(rows,cols-2*launch_l_pixels),\
       reSizePixMap[:,cols-launch_l_pixels:cols]),1), fmt = '%d', delimiter = ",")
  else:
    np.savetxt(csvFinal, np.concatenate((initPixMap[:,0:launch_l_pixels], DBS.best_solution.reshape(rows,cols-2*launch_l_pixels),\
       initPixMap[:,cols-launch_l_pixels:cols]),1), fmt = '%d', delimiter = ",") 

DBS = dbsAlgo(pixels, cost_function, rows, 'none', 'max', simulPositions, max_iterations, call_back, initial_solution=flatInitPixMap)

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