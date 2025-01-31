import math
import meep as mp
#import meep.adjoint as mpa
#import autograd.numpy as npa
#import meep_materials
#from meep.materials import Cu

#from autograd import tensor_jacobian_product, grad
import gdspy
import os
import time, sys
from vt_rrfc import *
import matplotlib.pyplot as plt
import subprocess 
import psutil

class EmSim:
  def __init__(self, 
               workingPath: str,
               adsLibName: str,
               gdsFile: str,
               csvFile: str,
               numPorts: int,
               portPositions: list,
               gdsCellName: str,
               dataFile: str,
               Substrate: object = None):
    self.pathName = workingPath
    self.libName = adsLibName
    self.gds_file = gdsFile
    self.csv_file = csvFile
    self.ports = numPorts
    self.portPosition = portPositions
    self.cell = gdsCellName
    self.dataF = dataFile
    self.sub = Substrate
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
  def mom_run(self):
    """
    Import GDS into ADS environment, set up the environment for simulation, and run the Momentum simulation.
    """
    aelName = 'autoloadEMSim.dem'
    AelTools.create_open_ael(self.pathName, self.libName, self.gds_file, self.ports, self.portPosition, aelName, self.cell)

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
    ##command = "source /software/RFIC/cadtools/cadence/setups/setup-tools; ads -m " \
    ##          + self.pathName + self.libName + "_wrk/" + aelName + " &"
    command = f"source /software/RFIC/cadtools/cadence/setups/setup-tools; ads -m {self.pathName}{self.libName}_wrk/{aelName} &"
    #os.system(command)
    ##adsPrep = subprocess.run(command, shell=True)
    subprocess.run(command, shell=True)
    time.sleep(30)

    # Import of the gds and automation of port addition takes a few seconds. I add a pause 
    # to ensure that the import is fully done before starting the sim (can probably be 
    # shortened. Will later look to just add feature to only proceed when the execution 
    # above is complete
    #time.sleep(90)
    print('We are still working')

    dataSet = self.dataF.replace(self.pathName + 'data/','') 
    print(dataSet)
    # Run Momentum Simulation
    os.chdir(self.pathName + self.libName + '_wrk/simulation/' + self.libName + '_lib/' + \
             self.cell + '/layout/emSetup_MoM/')
    command = 'source /software/RFIC/cadtools/cadence/setups/setup-tools; \
              adsMomWrapper --dsName=' + dataSet + ' --dsDir=' + self.pathName + 'data/ -O -3D proj proj'
              #adsMomWrapper -O -3D proj proj'
    #os.system(command)
    subprocess.run(command, shell=True , check=True)
    
    # Create the ADS Dataset
    os.chdir(self.pathName + self.libName + '_wrk/simulation/' + self.libName + '_lib/' + \
             self.cell + '/layout/emSetup_MoM/')
    command = 'source /software/RFIC/cadtools/cadence/setups/setup-tools; \
              adsMomWrapper -CD proj proj'
    #os.system(command)
    subprocess.run(command, shell=True, check=True)
    
    # Clean up after Momentum Simulation to prepare for next simulation
    aelCloseName = 'autoCloseEMSim.dem'
    AelTools.create_close_ael(self.pathName,self.libName,aelCloseName, self.cell)
  
    os.chdir(self.pathName)
    # Eveything in its right place :)
    command = 'mv ' + self.libName + '_wrk/simulation/' + \
              self.libName + '_lib/' + self.cell + '/layout/emSetup_MoM/proj.afs ' + \
              'data/spfiles/afs/' + dataSet + \
              '.afs; mv ' + self.csv_file + ' ' + self.pathName + \
              '/data/pixelMaps/.; mv ' + self.gds_file + ' ' + self.pathName + \
              '/data/gds/.'
    #os.system(command)
    subprocess.run(command, shell=True, check=True)
    # This will delete the layout view, leaving the emSetup so that the environment
    # is preparted for the next simulation
    command = 'source /software/RFIC/cadtools/cadence/setups/setup-tools; ads -m ' + \
              self.pathName + self.libName + '_wrk/' + aelCloseName + ' &'
    #os.system(command)
    subprocess.run(command, shell=True, check=True)
    # This removes the simulation path so the environment can create a fresh one for 
    # the next simulation. All data should have been moved in steps above.
    command = 'rm -rf ' + self.pathName + self.libName + '_wrk/simulation/*'
    #os.system(command)
    subprocess.run(command, shell=True, check=True)

    time.sleep(10)
    print('Cleaning up!')

  def sigmoid(x):
    return 1 / (1 + npa.exp(-x))

  def emx_run(self, procFile):
    # Call EMX Simulation

    if self.ports == 1:
      emsPorts = '-p P000=p1 -p P001=p2 -i P000 '
      sports = '.s1p'
    elif self.ports == 2:
      emsPorts = '-p P000=p1 -p P001=p2 -p P002=p3 -p P003=p4 -i P000 -i P001 '
      sports = '.s2p'
    elif self.ports == 3:
      emsPorts = '-p P000=p1 -p P001=p2 -p P002=p3 -p P003=p4 -p P004=p5 ' + \
                 '-p P005=p6 -i P000 -i P001 -i P002 '
      sports = '.s3p'
    elif self.ports == 4:
      emsPorts = '-p P000=p1 -p P001=p2 -p P002=p3 -p P003=p4 -p P004=p5 ' + \
                 '-p P005=p6 -p P006=p7 -p P007=p8 -i P000 -i P001 -i P002 ' + \
                 '-i P003 '
      sports = '.s4p'
    else:
      print('Only 1, 2, 3, or 4 ports are currently supported')

    # An example setup-tools is included in the repository. This is an example of the tool setup
    # script for CAD tools used by the RFIC group in MICS at Virginia Tech. All that is needed 
    # for this script is the EMX setup which will put emx on the path. You also need a
    # pointer to the proc file you are using.     
    #+ self.gds_file + ' RANDOM ' + procFile + ' -e 0.2 -t 0.2 -v 0.2 --3d=* ' \
    command = 'source /software/RFIC/cadtools/cadence/setups/setup-tools; emx ' \
              + self.gds_file + ' ' + self.cell + ' ' + procFile + ' -e 1 -t 1 -v 0.5 --3d=* ' \
              + emsPorts + '--sweep 0 1e+11 --sweep-stepsize 1e+08 --verbose=3 --print-command-line -l ' \
              + '2 --dump-connectivity --quasistatic --dump-connectivity ' \
              + '--parallel=0 --simultaneous-frequencies=0 --recommended-memory ' \
              + '--key=EMXkey --format=touchstone -s ' + self.dataF + sports

    os.system(command)

    # Clean up after sim
    # '.s2p ' + self.pathName + '/data/spfiles/.'
    os.chdir(self.pathName)
    command = 'mv ' + self.csv_file + ' ' + self.pathName + '/data/pixelMaps/.; mv ' + \
              self.gds_file + ' ' + self.pathName + '/data/gds/.; mv ' + self.dataF + \
              sports + ' ' + self.pathName + '/data/spfiles/.'
    os.system(command)

  def meep_run(self):
    
    # Define specific boundary conditions for meep
    res = 50 # number of pixels per mil
    three_d = False # Do a full 3D calculation or no
    gds = gdspy.GdsLibrary(infile=self.gds_file) # load the GDS file
    pml_size = 1.0
    CELL_LAYER = 0
    PORT1_LAYER = 1
    PORT2_LAYER = 2
    PORT3_LAYER = 3
    PORT4_LAYER = 4
    SOURCE_LAYER = 5
    METAL_LAYER = 11
    
    # Define materials and frequencies
    fr4 = mp.Medium(epsilon=4.5)
    freq = 5e9
    dpml = 0
    cell_thickness = dpml + self.sub.t_metal + self.sub.t_sub + self.sub.t_metal + 2*(self.sub.t_metal+self.sub.t_sub) + dpml
    cell_zmin = 0
    cell_zmax = cell_zmin + cell_thickness
    cu_zmax = 0.5*self.sub.t_metal
    cu_zmin = -0.5*self.sub.t_metal

    # Read the cell size and volumes for the sources and monitors from the GDS file
    rrc = mp.get_GDSII_prisms(Cu, self.gds_file, METAL_LAYER, cu_zmin, cu_zmax)
    gnd = mp.get_GDSII_prisms(Cu, self.gds_file, CELL_LAYER, cu_zmin-self.sub.t_sub, cu_zmax-self.sub.t_sub-self.sub.t_metal)
    cell = mp.GDSII_vol(self.gds_file, CELL_LAYER, 0, 0)
    src_vol = mp.GDSII_vol(self.gds_file, SOURCE_LAYER, cu_zmax-self.sub.t_sub-self.sub.t_metal, cu_zmin)
    p1 = mp.GDSII_vol(self.gds_file, PORT1_LAYER, zmin=cu_zmin, zmax=cu_zmax)
    p2 = mp.GDSII_vol(self.gds_file, PORT2_LAYER, zmin=cu_zmin, zmax=cu_zmax)

    for np in range(len(rrc)):
      rrc[np].center -= gnd.center
      for nv in range(len(rrc[np].vertices)):
        rrc[np].vertices[nv] -= gnd.center
    src_vol.center -= gnd.center
    geometry = rrc + gnd
    sources = [mp.Source(mp.GaussianSource(freq, fwidth=0.5*freq),
                         component=mp.Ey,
                         center = src_vol.center,
                         size = src_vol.size)]
    sim = mp.Simulation(cell_size=gnd.size,
                        geometry=geometry,
                        sources=sources,
                        resolution=res)
    freqs = np.linspace(0,10e9,101)
    s_params = sim.get_S_parameters(freqs, 1)
    plt.plot(freqs/1e9, np.abs(s_params[0, 0, :]), label='S11')
    print("hello world")
    #p3 = mp.GDSII_vol(gds_file, PORT3_LAYER, zmin=cu_zmin, zmax=cu_zmax)
    #p4 = mp.GDSII_vol(gds_file, PORT4_LAYER, zmin=cu_zmin, zmax=cu_zmax)