import os
import time, sys, random
from vt_rrfc import MicrostripSub, MicrostripCalc, UStrip, RandomUstrip, RandomComponent, EmSim

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

sub1 = MicrostripSub(t, cond, h, er, fc)

simulator = 'ADS' # This controls the simulation to be used. Right now there are two valid values 'ADS' or 'EMX'
libName = 'MyFirstWorkspace' # This is the name of the ADS workspace that was created to run sims in
sym = 'y-axis' # Do you want the random pixels symmetric about 'x-axis' or 'y-axis'
sim = False # This controls whether a simulation is run or not.
view = True # This controls if the GDS is viewed after each creation 
write = True # Control whether output files are written or not
ports = 2 # For now, code makes either 2, 3 or 4 ports
sides = 2 # For now, code can put ports on 2, 3 or 4 sides, with constraints that are spelled out in rrfc
scale = 1
pixelSize = scale*18 # the size of the randomized pixel in mils. Typically contrained by a PCB manufacturer.
minPixel = scale*6
layoutRes = scale*1 # Transmission line calculators are in mils, so since the target here is um, 25.4 converts to um
layoutUnit = 25.4e-6/layoutRes # Set the layout unit to microns/layout resolution

# Define Zo,h and Zo,l for the filter sections: Assume that the minimum width is 1 line of pixels; this
# defines the width for Zo,h. Assume that the maximum width is 15 pixels; this defines the width for Zo,l.
filtType = 'butter'
filtOrder = 3
max_w = 9	
min_w = 1

Zo_h, ereff_h = MicrostripCalc.calc_microstrip(sub1, min_w*pixelSize, 1000);
Zo_l, ereff_l = MicrostripCalc.calc_microstrip(sub1, max_w*pixelSize, 1000);
w_h = min_w*pixelSize
w_l = max_w*pixelSize
print(Zo_h,Zo_l)

pathName = 'H:\\Shared drives\\RRFC\\SimFiles\\randomGDS\\' # Base path for file creation

# This code will generate a prototype filter using the classical design method for the stepped impedance filter. 
# The prototype has a fixed order and filter type that can be defined, and ultimately this filter will be used 
# to define the boundaries of the random structures that will be designed and simulated to train an Neural Network.
os.chdir(pathName)
outFile = pathName + 'data\\' + "steppedImpFilter_pixelSize=" + str(pixelSize) + "_order=" + str(filtOrder) + '_' + filtType
ustrip1 = UStrip(
    unit=layoutUnit,
    corner='standard',
    pixelSize=pixelSize,
    layoutRes=layoutRes,
    scale=scale,
    sim=simulator,
    view=False,
    write=False,
    outF=outFile)
_, xProto, yProto, _, _, _, launch_l_pixels = ustrip1.u_strip_stepped_imp_filter_gds(
    sub1, filtType, filtOrder, w_h, w_l, Zo_h, Zo_l, z0)
yProto = 19 * pixelSize
xProto = 37 * pixelSize
xProto = xProto - 2 * launch_l_pixels * pixelSize #temporarily adjust for launch will add patch 
print(xProto, yProto)
"""
connectMap is a map for connections to be enforced: 1_2 1_3 1_4 2_3 2_4 3_4
if any position in array is a 1, files will only be printed if that 
connectivity is true. 1_2 means that port 1 and 2 are connected 1_3 = port 1 
and 3, etc. This essentially forces a DC connection to exist between the 
ports and more than one connection can be enforced at a time. connectMap = 
[1, 1, 0, 1, 0, 0] would enforce connections between ports 1 and 2, ports 1 
and 3 and ports 2 and 3 as an example
"""
connectMap = [1, 0, 0, 0, 0, 0]

y = 92795
for x in range(230372,300000): # Run 100 iterations of file generation and simulation.
  random.seed(x)
  symSelect = random.randint(0, 3)
  conSelect = random.randint(0, 1)
  if symSelect == 0:
    sym = 'x-axis'
    if conSelect == 0:
      connectMap = [0, 0, 0, 0, 0, 0]
    else:
      connectMap = [1, 0, 0, 0, 0, 0]
  elif symSelect == 1:
    sym = 'y-axis'
    if conSelect == 0:
      connectMap = [0, 0, 0, 0, 0, 0]
    else:
      connectMap = [1, 0, 0, 0, 0, 0]
  elif symSelect == 2:
    sym = 'xy-axis'
    if conSelect == 0:
      connectMap = [0, 0, 0, 0, 0, 0]
    else:
      connectMap = [1, 0, 0, 0, 0, 0]
  else:
    sym = 'asym'
    if conSelect == 0:
      connectMap = [0, 0, 0, 0, 0, 0]
    else:
      connectMap = [1, 0, 0, 0, 0, 0]

  data_file = pathName + 'data\\' + "randomGDSSteppedImpFilter_Type=" + filtType + "_order=" \
              + str(filtOrder) + "_pixelSize=" + str(pixelSize) + "_sim=" + str(y)

  rComp1 = RandomComponent(
      unit=layoutUnit,
      ports=ports,
      sides=sides,
      corner='standard',
      connect=connectMap,
      minPix=minPixel,
      pixelSize=pixelSize,
      layoutRes=layoutRes,
      scale=1,
      launchLen=30,
      seed=x,
      sim=simulator,
      view=view,
      write=write,
      outF=data_file,
      sym=sym,
      portPosition='')
  portPosition, xBoard, yBoard, csv_file, gds_file, cell, _ = rComp1.random_gds_dim(sub1, xProto, yProto, z0)

  # checking if files were written. When connectivity is enforced, files are only written for
  # the random seeds that meet the connection criterion. If no file is written, no simulation 
  # will be run. Also, a seedList is created that gives he list of seeds that generated DC 
  # connections for the forced connectivity ports.
  if csv_file == '':
    sim = False
    print('Seed=' + str(x) + ' File is empty')
  else:
    y += 1
    print('Seed=' + str(x) + ' File is not empty')
    sim = True
    fileList = pathName + 'seedList.txt'
    f = open(fileList, 'a+')
    f.write(str(x) + ' ' + sym + ' connect=' + str(conSelect) + '\n')
    f.close

  if sim == True:
    # Import GDS into ADS environment and setup environment for simulation
    em1 = EmSim(workingPath = pathName, 
                adsLibName = libName, 
                gdsFile = gds_file,
                csvFile = csv_file, 
                numPorts = ports, 
                portPositions = portPosition,
                gdsCellName = cell, 
                dataFile = data_file)
    EmSim.mom_run(em1)

