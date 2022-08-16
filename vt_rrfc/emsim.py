import os
import time, sys
from vt_rrfc import *

class emSim:
  def __init__(self, 
               workingPath: str,
               adsLibName: str,
               gdsFile: str,
               csvFile: str,
               numPorts: int,
               portPositions,
               gdsCellName: str,
               dataFile: str):
    self.pathName = workingPath
    self.libName = adsLibName
    self.gds_file = gdsFile
    self.csv_file = csvFile
    self.ports = numPorts
    self.portPosition = portPositions
    self.cell = gdsCellName
    self.dataF = dataFile
  """
  This code assumes that in your path, you will have an ADS workspace created amd that at 
  the same level of hierarchy you will have a data folder to collect simulation results 
  and design results in. The structure should be:
    ./YOUR_ADS_Workspace_wrk
    ./data
    ./data/gds
    ./data/pixelMaps
    ./data/spfiles
    ./data/spfiles/afs
  This is necessary because of how the simulation environment is created in ADS, and is not
  necessary for EMX simulations.

  WorkingPath --> base path where files will be placed
  adsLibName --> name of the ADS library that the simulation will be loaded. This library 
                 must exist before simulations can be called and run. Working to see if it
                 is possible to create a library using AEL, but no luck so far
  gdsFile --> This is the artwork file that will be loaded into ADS to be simulated
  csvFile --> This is a binary pixel map that shows the pixels that are pop'd (1) and 
              unpop'd (0)
  numPorts --> The number of ports that ADS will add to the drawing. 
  portPositions --> Vector that has the location of up to four ports...will work to expand 
                    possible port positions
  cell --> The cell name of the layout in the GDS file. For now, there needs to be a cell 
           name in the library with an emSetup located under the cell for the simulation 
           to be launched. Working to see if this can be automated, but no luck so far.
  dataFile --> This is the name of the data file that will be created at the end of the
               simulation. 
  """
  def momRun(self):
    # Import GDS into ADS environment and setup environment for simulation
    aelName = 'autoloadEMSim.dem'
    createOpenAel(self.pathName, self.libName, self.gds_file, self.ports, \
                      self.portPosition, aelName, self.cell)

    # This is must setup ADS tools to be on the path. It points to a setup script for the ECE
    # department at Virginia Tech. The script for CAD tools is used by the RFIC group in MICS 
    # at Virginia Tech. From a basic perspective, all that is needed for this script is to
    # setup $HPEESOF according to Keysights instructions and to add ADS binaries (e.g., ads, 
    # adsMomWrapper, etc.) onto the path. It must also provide a path to the license server, 
    # if the tools and license are setup by default in your system, this command can be
    # bypassed. Once the tools are on the path, ads is launched with
    # the AEL script that is created above the ADS setup (HPEESOF) which will put ads and 
    # adsMomWrapper on the path and setup the licensing. Future variants of the script will 
    # add hooks for EMX simulation    
    command = 'source /software/RFIC/cadtools/cadence/setups/setup-tools; ads -m ' \
              + self.pathName + self.libName + '_wrk/' + aelName + ' &'
    os.system(command)

    # Import of the gds and automation of port addition takes a few seconds. I add a pause 
    # to ensure that the import is fully done before starting the sim (can probably be 
    # shortened. Will later look to just add feature to only proceed when the execution 
    # above is complete
    time.sleep(30)
    print('We are still working')

    dataSet = self.dataF.replace(self.pathName + 'data/','') 
    print(dataSet)
    # Run Momentum Simulation
    os.chdir(self.pathName + self.libName + '_wrk/simulation/' + self.libName + '_lib/' + \
             self.cell + '/layout/emSetup_MoM/')
    command = 'source /software/RFIC/cadtools/cadence/setups/setup-tools; \
              adsMomWrapper --dsName=' + dataSet + ' --dsDir=' + self.pathName + 'data/ -O -3D proj proj'
              #adsMomWrapper -O -3D proj proj'
    os.system(command)
    
    # Create the ADS Dataset
    os.chdir(self.pathName + self.libName + '_wrk/simulation/' + self.libName + '_lib/' + \
             self.cell + '/layout/emSetup_MoM/')
    command = 'source /software/RFIC/cadtools/cadence/setups/setup-tools; \
              adsMomWrapper -CD proj proj'
    os.system(command)
    
    # Clean up after Momentum Simulation to prepare for next simulation
    aelCloseName = 'autoCloseEMSim.dem'
    createCloseAel(self.pathName,self.libName,aelCloseName, self.cell)
  
    os.chdir(self.pathName)
    # Eveything in its right place :)
    command = 'mv ' + self.libName + '_wrk/simulation/' + \
              self.libName + '_lib/' + self.cell + '/layout/emSetup_MoM/proj.afs ' + \
              'data/spfiles/afs/' + dataSet + \
              '.afs; mv ' + self.csv_file + ' ' + self.pathName + \
              '/data/pixelMaps/.; mv ' + self.gds_file + ' ' + self.pathName + \
              '/data/gds/.'
    os.system(command)
    # This will delete the layout view, leaving the emSetup so that the environment
    # is preparted for the next simulation
    command = 'source /software/RFIC/cadtools/cadence/setups/setup-tools; ads -m ' + \
              self.pathName + self.libName + '_wrk/' + aelCloseName + ' &'
    os.system(command)
    # This removes the simulation path so the environment can create a fresh one for 
    # the next simulation. All data should have been moved in steps above.
    command = 'rm -rf ' + self.pathName + self.libName + '_wrk/simulation/*'
    os.system(command)

    time.sleep(10)
    print('Cleaning up!')
  
  def emxRun(self, procFile):
    # Call EMX Simulation
    
    # An example setup-tools is included in the repository. This is an example of the tool setup
    # script for CAD tools used by the RFIC group in MICS at Virginia Tech. All that is needed 
    # for this script is the EMX setup which will put emx on the path. You also need a
    # pointer to the proc file you are using.     
    command = 'source /software/RFIC/cadtools/cadence/setups/setup-tools; emx --edge-width=2 ' \
              + self.gds_file + ' RANDOM ' + procFile + ' --verbose=3 --touchstone -s ' \
              + self.dataF + '.s2p --sweep 0 30e9 --key=EMXkey --max-memory=60' \
              + ' --include-command-line --include-port-order'
    os.system(command)

    # Clean up after sim
    os.chdir(self.pathName)
    command = 'mv ' + self.csv_file + ' ' + self.pathName + '/data/pixelMaps/.; mv ' + \
              self.gds_file + ' ' + self.pathName + '/data/gds/.; mv ' + self.data_file + \
              '.s2p ' + self.pathName + '/data/spfiles/.'
    os.system(command)
