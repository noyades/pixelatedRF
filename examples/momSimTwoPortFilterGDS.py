# Python imports
import os
import time, sys, os
from vt_rrfc import microstrip as ustrip
from vt_rrfc import ustripComponents as usc
from vt_rrfc import ustripRandomComponents as usrc
from vt_rrfc import adsAelGeneration as ael

# Load constants and design choices. Assumes 2 layer PCB with 2 oz copper and 30 mil total thickness
fc = 5.0e9 # Operating center frequency for electrical length calculations (make this smaller than or equal to the desired operating frequency
z0 = 50 # Desired characteristic impedance of launches
EL = 90 # Desired unit of electrical length in degrees
t = 2.8 # Thickness of conductor in mils
cond = 5.88e7 # Conductivity of the conductors
h = 30 # height of conductor above substrate
t_air = 2*(t+h) # thickness of the air above the conductor layer
er = 4.5 # relative permittivity of the substrate material

sub1 = ustrip.microstrip_sub(t_metal=t, metalCond=cond, t_sub=h, relPerm=er, freq=fc)

simulator = 'ADS' # This controls the simulation to be used. Right now there are two valid values 'ADS' or 'EMX'
sim = True # This controls whether a simulation is run or not.
view = False # This controls if the GDS is viewed after each creation 
write = True # Control whether output files are written or not
ports = 2 # For now, code makes either 2, 3 or 4 ports
sides = 2 # For now, code can put ports on 2, 3 or 4 sides, with constraints that are spelled out in rrfc
pixelSize = 18 # the size of the randomized pixel in mils. Typically contrained by a PCB manufacturer.

# Define Zo,h and Zo,l for the filter sections: Assume that the minimum width is 1 line of pixels; this
# defines the width for Zo,h. Assume that the maximum width is 15 pixels; this defines the width for Zo,l.
filtType = 'cheby3'
filterOrder = 3
max_w = 9
min_w = 1

Zo_h, ereff_h = ustrip.microstrip_calc.calcMicrostrip(sub1, min_w*pixelSize, 1000);
Zo_l, ereff_l = ustrip.microstrip_calc.calcMicrostrip(sub1, max_w*pixelSize, 1000);
w_h = min_w*pixelSize
w_l = max_w*pixelSize
print(Zo_h,Zo_l)

pathName = '/home/jswalling/pythonWork/rrfc/' # Base path for file creation
os.chdir(pathName)
#gds.uStripSteppedImpFilterGDS(filtType, filterOrder, w_h, w_l, Zo_h, Zo_l, pixelSize, fc, z0, t, h, er, simulator, view)

for x in range(3,4): # Run 100 iterations of file generation and simulation.
  data_file = pathName + 'data/' + "steppedImpFilter_pixelSize=" + str(pixelSize) + "_order=" + str(x) \
              + '_' + filtType

  portPosition, xBoard, yBoard, csv_file, gds_file, cell = usc.uStripSteppedImpFilterGDS(sub1, filtType, x, \
                                                     w_h, w_l, Zo_h, Zo_l, pixelSize, z0, simulator, view, \
                                                     write, data_file)
  
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

    dataSet = "steppedImpFilter_pixelSize=" + str(pixelSize) + "_order=" + str(x) \
              + '_' + filtType
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
