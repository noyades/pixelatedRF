import os
import time, sys
from vt_rrfc import *

# Load constants and design choices. Assumes 2 layer PCB with 1 oz copper and 19 mil total thickness
fc = 15.0e9 # Operating center frequency for electrical length calculations 
           # (make this smaller than or equal to the desired operating frequency
z0 = 50 # Desired characteristic impedance of launches
EL = 90 # Desired unit of electrical length in degrees
t = 0.1114 # Thickness of conductor in mils
cond = 3.21e7 # Conductivity of the conductors (S/m)
h = 27.823 # height of conductor above substrate
t_air = 2*(t+h) # thickness of the air above the conductor layer
er = 11.83 # relative permittivity of the substrate material

sub1 = microstrip_sub(t, cond, h, er, fc)

simulator = 'EMX' # This controls the simulation to be used. Right now there are two valid values 'ADS' or 'EMX'
sym = 'x-axis' # Do you want the random pixels symmetric about 'x-axis' or 'y-axis'
sim = False # This controls whether a simulation is run or not.
view = True # This controls if the GDS is viewed after each creation 
write = True # Control whether output files are written or not
ports = 2 # For now, code makes either 2, 3 or 4 ports
sides = 2 # For now, code can put ports on 2, 3 or 4 sides, with constraints that are spelled out in rrfc
pixelSize = 6 # the size of the randomized pixel in microns. Typically contrained by a PCB manufacturer.
layoutUnit = 1e-6 # Set the layout unit to microns

# Define Zo,h and Zo,l for the filter sections: Assume that the minimum width is 1 line of pixels; this
# defines the width for Zo,h. Assume that the maximum width is 15 pixels; this defines the width for Zo,l.
filtType = 'butter'
filtOrder = 3
max_w = 4
min_w = 1

Zo_h, ereff_h = ustrip.microstrip_calc.calcMicrostrip(sub1, min_w*pixelSize, 1000);
Zo_l, ereff_l = ustrip.microstrip_calc.calcMicrostrip(sub1, max_w*pixelSize, 1000);
w_h = min_w*pixelSize
w_l = max_w*pixelSize
print(Zo_h,Zo_l)
"""
procFile = '/software/RFIC/PDK/globalFoundries/22FDX-EXT/release/Emagnetic/EMX/10M_2Mx_4Cx_2Bx_2Jx_LBthick/22fdsoi_10M_2Mx_4Cx_2Bx_2Jx_LBthick_nominal_detailed.encrypted.proc'
pathName = '/home/jswalling/pythonWork/rrfc/' # Base path for file creation

# This code will generate a prototype filter using the classical design method for the stepped impedance filter. 
# The prototype has a fixed order and filter type that can be defined, and ultimately this filter will be used 
# to define the boundaries of the random structures that will be designed and simulated to train an Neural Network.
os.chdir(pathName)
outFile = pathName + 'data/' + "steppedImpFilter_pixelSize=" + str(pixelSize) + "_order=" + str(filtOrder) + '_' + filtType + '_EMX_22nm'
rrfc1 = usc.rrfc(unit=layoutUnit,pixelSize=pixelSize,sim=simulator,\
        view=view,write=write,outF=outFile)
_, xBoard, yBoard, _, _, _ = usc.uStripSteppedImpFilterGDS(sub1, rrfc1, \
                                          filtType, filtOrder, w_h, w_l, Zo_h, Zo_l, z0)
print(xBoard, yBoard)
"""
"""
connectMap is a map for connections to be enforced: 1_2 1_3 1_4 2_3 2_4 3_4
if any position in array is a 1, files will only be printed if that 
connectivity is true. 1_2 means that port 1 and 2 are connected 1_3 = port 1 
and 3, etc. This essentially forces a DC connection to exist between the 
ports and more than one connection can be enforced at a time. connectMap = 
[1, 1, 0, 1, 0, 0] would enforce connections between ports 1 and 2, ports 1 
and 3 and ports 2 and 3 as an example
"""
"""
connectMap = [1, 0, 0, 0, 0, 0]

y = 0
for x in range(0,10000): # Run 100 iterations of file generation and simulation.
  data_file = pathName + 'data/' + "randomGDSSteppedImpFilter_Type=" + filtType + "_order=" \
              + str(filtOrder) + "_pixelSize=" + str(pixelSize) + "_sim=" + str(y) + '_EMX_22nm'
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
    fileList = pathName + 'seedList_EMX.txt'
    f = open(fileList, 'a+')
    f.write(str(x)+' '+ sym+ '\n')
    f.close

  if sim == True:
    # Call EMX Simulation
    
    # An example setup-tools is included in the repository. This is an example of the tool setup
    # script for CAD tools used by the RFIC group in MICS at Virginia Tech. All that is needed 
    # for this script is the EMX setup which will put emx on the path. You also need a
    # pointer to the proc file you are using.     
    command = 'source /software/RFIC/cadtools/cadence/setups/setup-tools; emx --edge-width=2 ' \
              + gds_file + ' RANDOM ' + procFile + ' --verbose=3 --touchstone -s ' \
              + data_file + '.s2p --sweep 0 30e9 --key=EMXkey --max-memory=60' \
              + ' --include-command-line --include-port-order'
    os.system(command)

    # Clean up after sim
    os.chdir(pathName)
    command = 'mv ' + csv_file + ' ' + pathName + '/data/pixelMaps/.; mv ' + gds_file + ' ' + \
              pathName + '/data/gds/.; mv ' + data_file + '.s2p ' + pathName + '/data/spfiles/.'
    os.system(command)
"""
