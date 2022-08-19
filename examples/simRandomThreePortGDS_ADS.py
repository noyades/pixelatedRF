import os
import time, sys
import random
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
libName = 'threePortEmSim' # This is the name of the ADS workspace that was created to run sims in
sym = 'xy-axis' # Do you want the random pixels symmetric about 'x-axis' or 'y-axis'
sim = False # This controls whether a simulation is run or not.
view = False # This controls if the GDS is viewed after each creation 
write = True # Control whether output files are written or not
ports = 3 # For now, code makes either 2, 3 or 4 ports
sides = 2 # For now, code can put ports on 2, 3 or 4 sides, with constraints that are spelled out in rrfc
pixelSize = 18 # the size of the randomized pixel in mils. Typically contrained by a PCB manufacturer.
layoutUnit = 25.4e-6 # Set the layout unit to mils

w_l, l_l = microstrip_calc.synthMicrostrip(sub1, z0, 120);

print(w_l,l_l)

pathName = '/home/jswalling/pythonWork/rrfc/3port/' # Base path for file creation

xRect, yRect = 25*pixelSize, 15*pixelSize
"""
connectMap is a map for connections to be enforced: 1_2 1_3 1_4 2_3 2_4 3_4
if any position in array is a 1, files will only be printed if that 
connectivity is true. 1_2 means that port 1 and 2 are connected 1_3 = port 1 
and 3, etc. This essentially forces a DC connection to exist between the 
ports and more than one connection can be enforced at a time. connectMap = 
[1, 1, 0, 1, 0, 0] would enforce connections between ports 1 and 2, ports 1 
and 3 and ports 2 and 3 as an example
"""
connectMap = [1, 1, 0, 1, 0, 0]

y = 1294
for x in range(73580,100000): # Run 100 iterations of file generation and simulation.
  random.seed(x)
  symSelect = random.randint(0, 3)
  if symSelect == 0:
    sym = 'x-axis'
  elif symSelect == 1:
    sym = 'y-axis'
  elif symSelect == 2:
    sym = 'xy-axis'
  else:
    sym = 'asym'
  """
  conSelect = random.randint(0,7)
  if conSelect == 0:
    connectMap = [0, 0, 0, 0, 0, 0]
  elif conSelect == 1:
    connectMap = [0, 0, 1, 0, 0, 0]
  elif conSelect == 2:
    connectMap = [0, 1, 0, 0, 0, 0]
  elif conSelect == 3:
    connectMap = [0, 1, 1, 0, 0, 0]
  elif conSelect == 4:
    connectMap = [1, 0, 0, 0, 0, 0]
  elif conSelect == 5:
    connectMap = [1, 0, 1, 0, 0, 0]
  elif conSelect == 6:
    connectMap = [1, 1, 0, 0, 0, 0]
  else:
    connectMap = [1, 1, 1, 0, 0, 0]
  """
  data_file = pathName + 'data/' + "randomGDSThreePort_450x270" + "_pixelSize=" +\
              str(pixelSize) + "_sim=" + str(y)
  rrfc1 = rrfc(unit=layoutUnit,ports=ports,sides=sides,connect=connectMap,\
          pixelSize=pixelSize,seed=x,sim=simulator,view=view,write=write,\
          outF=data_file,sym=sym)
  portPosition, xBoard, yBoard, csv_file, gds_file, cell = randomGDS_dim(sub1, \
                                                    rrfc1, xRect, yRect, z0)

  # checking if files were written. When connectivity is enforced, files are only written for
  # the random seeds that meet the connection criterion. If no file is written, no simulation 
  # will be run. Also, a seedList is created that gives he list of seeds that generated DC 
  # connections for the forced connectivity ports.
  if csv_file == '':
    sim = False
  else:
    y += 1
    print('Seed=' + str(x) + ' File is not empty')
    sim = True
    fileList = pathName + 'seedList_threePort_450x270.txt'
    f = open(fileList, 'a+')
    f.write(str(x) + ' ' + sym + ' connect=0' + '\n')
    f.close

  if sim == True:
    # Import GDS into ADS environment and setup environment for simulation
    em1 = emSim(workingPath = pathName, adsLibName = libName, gdsFile = gds_file,\
              csvFile = csv_file, numPorts = ports, portPositions = portPosition,\
              gdsCellName = cell, dataFile = data_file)
    emSim.momRun(em1)

