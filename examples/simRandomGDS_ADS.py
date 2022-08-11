import os
import time, sys
from vt_rrfc import microstrip as ustrip
from vt_rrfc import ustripComponents as usc
from vt_rrfc import ustripRandomComponents as usrc
from vt_rrfc import adsAelGeneration as ael

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

sub1 = ustrip.microstrip_sub(t, cond, h, er, fc)

simulator = 'ADS' # This controls the simulation to be used. Right now there are two valid values 'ADS' or 'EMX'
sym = 'x-axis' # Do you want the random pixels symmetric about 'x-axis' or 'y-axis'
sim = False # This controls whether a simulation is run or not.
view = False # This controls if the GDS is viewed after each creation 
write = True # Control whether output files are written or not
ports = 2 # For now, code makes either 2, 3 or 4 ports
sides = 2 # For now, code can put ports on 2, 3 or 4 sides, with constraints that are spelled out in rrfc
pixelSize = 18 # the size of the randomized pixel in mils. Typically contrained by a PCB manufacturer.
layoutUnit = 25.4e-6 # Set the layout unit to mils

# Define Zo,h and Zo,l for the filter sections: Assume that the minimum width is 1 line of pixels; this
# defines the width for Zo,h. Assume that the maximum width is 15 pixels; this defines the width for Zo,l.
filtType = 'butter'
filtOrder = 3
max_w = 9	
min_w = 1

Zo_h, ereff_h = ustrip.microstrip_calc.calcMicrostrip(sub1, min_w*pixelSize, 1000);
Zo_l, ereff_l = ustrip.microstrip_calc.calcMicrostrip(sub1, max_w*pixelSize, 1000);
w_h = min_w*pixelSize
w_l = max_w*pixelSize
print(Zo_h,Zo_l)

pathName = '/home/jswalling/pythonWork/rrfc/' # Base path for file creation

# This code will generate a prototype filter using the classical design method for the stepped impedance filter. 
# The prototype has a fixed order and filter type that can be defined, and ultimately this filter will be used 
# to define the boundaries of the random structures that will be designed and simulated to train an Neural Network.
os.chdir(pathName)
outFile = pathName + 'data/' + "steppedImpFilter_pixelSize=" + str(pixelSize) + "_order=" + str(filtOrder) + '_' + filtType
_, xBoard, yBoard, _, _, _ = usc.uStripSteppedImpFilterGDS(sub1, filtType, filtOrder, \
                                                     w_h, w_l, Zo_h, Zo_l, pixelSize, z0, \
                                                     simulator, view, False, outFile)
print(xBoard, yBoard)
"""
connectMap is a map for connections to be enforced: 1_2 1_3 1_4 2_3 2_4 3_4
if any position in array is a 1, files will only be printed if that 
connectivity is true. 1_2 means that port 1 and 2 are connected 1_3 = port 1 
and 3, etc. This essentially forces a DC connection to exist between the 
ports and more than one connection can be enforced at a time. connectMap = 
[1, 1, 0, 1, 0, 0] would enforce connections between ports 1 and 2, ports 1 
and 3 and ports 2 and 3 as an example
"""
connectMap = [0, 0, 0, 0, 0, 0]

y = 2019
for x in range(0,100): # Run 100 iterations of file generation and simulation.
  data_file = pathName + 'data/' + "randomGDSSteppedImpFilter_Type=" + filtType + "_order=" \
              + str(filtOrder) + "_pixelSize=" + str(pixelSize) + "_sim=" + str(y)
  rrfc1 = usrc.rrfc(unit=layoutUnit,ports=ports,sides=sides,connect=connectMap,\
          pixelSize=pixelSize,seed=x,sim=simulator,view=view,write=write,\
          outF=data_file,sym=sym)
  portPosition, xBoard, yBoard, csv_file, gds_file, cell = usrc.randomGDS_dim(sub1, \
                                                      rrfc1, xBoard, yBoard, z0)

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
    fileList = pathName + 'seedList.txt'
    f = open(fileList, 'a+')
    f.write(str(x)+' '+ sym+ '\n')
    f.close

  if sim == True:
    # Import GDS into ADS environment and setup environment for simulation
    libName = 'MyFirstWorkspace'	
    aelName = 'autoloadEMSim.dem'
    ael.createOpenAel(pathName, libName, gds_file, ports, portPosition, aelName, cell)

    # An example setup-tools is included in the repository. This is an example of the tool setup
    # script for CAD tools used by the RFIC group in MICS at Virginia Tech. All that is needed 
    # for this script is the ADS setup (HPEESOF) which will put ads and adsMomWrapper on the
    # path and setup the licensing. Future variants of the script will add hooks for EMX simulation    
    command = 'source /software/RFIC/cadtools/cadence/setups/setup-tools; ads -m ' \
              + pathName + libName + '_wrk/' + aelName + ' &'
    os.system(command)

    time.sleep(30)
    print('We are still working')

    dataSet = "randomGDSSteppedImpFilter_Type=" + filtType + "_order=" + str(filtOrder) + "_pixelSize=" + \
              str(pixelSize) + "_sim=" + str(y-1)
    # Run Momentum Simulation
    os.chdir(pathName + libName + '_wrk/simulation/' + libName + '_lib/' + cell + '/layout/emSetup_MoM/')
    commands = ['source /software/RFIC/cadtools/cadence/setups/setup-tools',
                'echo $HPEESOF_DIR',
                'which adsMomWrapper',
	        'adsMomWrapper --dsName=' + dataSet + ' --dsDir=' + pathName + 'data/ -O -3D proj proj']
    command = 'source /software/RFIC/cadtools/cadence/setups/setup-tools; \
              adsMomWrapper --dsName=' + dataSet + ' --dsDir=' + pathName + 'data/ -O -3D proj proj'
    os.system(command)
    
    # Clean up after Momentum Simulation to prepare for next simulation
    aelCloseName = 'autoCloseEMSim.dem'
    ael.createCloseAel(pathName,libName,aelCloseName, cell)

    os.chdir(pathName)
    commands = ['mv ' + csv_file + ' ' + pathName + '/data/pixelMaps/.',
	        'mv ' + gds_file + ' ' + pathName + '/data/gds/.',
	        'source /software/RFIC/cadtools/cadence/setups/setup-tools',
	        'echo $HPEESOF_DIR',
	        'ads -m ' + pathName + libName + '_wrk/' + aelCloseName + ' &',
	        'rm -rf ' + pathName + libName + '_wrk/simulation/*']
    command = 'mv ' + csv_file + ' ' + pathName + '/data/pixelMaps/.; mv ' + gds_file + ' ' + pathName + '/data/gds/.'
    os.system(command)
    command = 'source /software/RFIC/cadtools/cadence/setups/setup-tools; ads -m ' + \
              pathName + libName + '_wrk/' + aelCloseName + ' &'
    os.system(command)
    command = 'rm -rf ' + pathName + libName + '_wrk/simulation/*'
    os.system(command)

