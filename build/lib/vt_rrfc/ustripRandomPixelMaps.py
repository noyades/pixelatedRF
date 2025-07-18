# -*- coding: utf-8 -*-
import csv
import math
import random
import numpy as np
import gdspy
import pya
#from svglib.svglib import svg2rlg
#from reportlab.graphics import renderPM
#from vt_rrfc import * 
from vt_rrfc.randomDesigns import *
from vt_rrfc.microstrip import MicrostripCalc
from vt_rrfc.connectivity import Connectivity

#from skimage import measure

class RandomUstrip:
  def __init__(self, 
               unit,
               ports: int, 
               sides: int,
               corner: str, 
               connect,
               minPix: int,
               pixelSize: int,
               layoutRes: int,
               scale: int,
               launchLen,
               seed: int, 
               sim: str, 
               view: bool, 
               write: bool, 
               outF: str,
               sym: str,
               portPosition):
    """ define the microstrip substrate
    Args:
    unit = grid unit of the layout. Set to 1e-6 for um or 25.4e-6 for mil for example
    ports = number of ports
    sides = the number of sides ports are placed on
    corner = how the corner connection is realize for simulation 
    sym = whether to force symmetry about 'x-axis' or 'y-axis'. Leave empty for no symmetry 
    minPix = minimum pixel size in the process
    pixelSize = size of the randomized pixel
    layoutRes = integer dividend on pixel size to set the layout resolution 
    launchLen = electrical length of the launches at the design frequency in degrees
    seed = random seed number
    sim = boolean flag on what simulator to assume (this is mainly for port placement and setup)
    view = boolean flag on whether to view gds (do not use in sim mode)
    write = boolean flag on whether to write files
    outF = Output file string
    Port definitions must take the following forms, presently no others are acceptable and using any wrong combination will cause
    program to exit with an error message:
    ports = 2, sides = 2: This will create a rectangle with one port (1) on the west side and one port (2) on the west side of 
                          the structure
    ports = 3, sides = 2: This will create a rectangle with one port (1) on the southwest side, one port (2) on the northwest 
                          side and one port (3) on the east side
                          of the structure
    ports = 3, sides = 3: This will create a rectangle with one port (1) on the west side, one port (2) on the east side and one
                          port (3) on the north side of the structure
    ports = 4, sides = 2: This will create a rectangle with one port (1) on the southwest side, one port (2) on the northwest 
                          side, one port (3) on the southeast side, and one port (4) on the northeast side of the structure
    ports = 4, sides = 4: This will create a rectangle with one port (1) on the west side, one port (2) on the north side,
                          one port (3) on the east side, and one port (4) on the south side of the structure
    
    connect = array that dictates the connections between different ports and whether to enforce strict DC connections (e.g.,
              continuity between ports
    connectMap is a map for connections to be enforced: 1_2 1_3 1_4 2_3 2_4 3_4
    if any position in array is a 1, files will only be printed if that 
    connectivity is true. 1_2 means that port 1 and 2 are connected 1_3 = port 1 
    and 3, etc. This essentially forces a DC connection to exist between the 
    ports and more than one connection can be enforced at a time. connectMap = 
    [1, 1, 0, 1, 0, 0] would enforce connections between ports 1 and 2, ports 1 
    and 3 and ports 2 and 3 as an example
    portPosition: For now is a 8 element array with the xy coords for the port positions [x1,y1,x2,y2,x3,y3,x4,y4]; Can be generated locally
    """    
    self.unit = unit
    self.ports = ports
    self.corner = corner 
    self.sides = sides
    self.sym = sym
    self.minPix = minPix
    self.pixelSize = pixelSize
    self.layoutRes = layoutRes
    self.scale = scale
    self.launchLen = launchLen
    self.seed = seed
    self.write = write
    self.sim = sim
    self.view = view
    self.outF = outF
    self.connect = connect
    self.portPos = portPosition
    self.mil2um = self.layoutRes / self.scale
    if self.scale == 0 or self.layoutRes == 0 or self.pixelSize == 0:
      raise ValueError("Scale, layout resolution and pixel size must be non-zero values.")

  def _calc_launch_pixels(self, sub, imp: float):
    """
    Calculate the launch pixels for the microstrip.

    Args:
    sub: Microstrip substrate initialization.
    imp (float): Characteristic impedance of the launch paths.

    Returns:
    tuple: A tuple containing the number of launch length pixels, 
           number of launch width pixels, width of the launch, 
           and length of the launch.
    """  
    
    width_launch, length_launch = MicrostripCalc.synth_microstrip(sub, imp, self.launchLen)
    width_launch = width_launch * self.scale * self.mil2um # Convert from mils to layout resolution. layoutRes defaults to 1 for mils
    length_launch = length_launch * self.scale * self.mil2um # Convert from mils to layout resolution. layoutRes defaults to 1 for mils
    if width_launch < self.pixelSize: # Make sure that launch is at least 1 pixel wide
      width_launch = self.pixelSize
    if length_launch < self.pixelSize: # Make sure that launch is at least 1 pixel long
      length_launch = self.pixelSize
    launch_x_pixels = round(length_launch / self.pixelSize) # length of the line to connect to structure in number of pixels
    launch_y_pixels = round(width_launch / self.pixelSize) # width of the line to connect to structure in number of pixels
    length_launch = launch_x_pixels * self.pixelSize
    width_launch = launch_y_pixels * self.pixelSize
    return launch_x_pixels, launch_y_pixels, width_launch, length_launch

  def random_pixmap_dim(self, sub, x_dim: int, y_dim: int, imp: float):
    """
    This script will create a rectangular pixel grid that is x_dim wide by y_dim tall 
    according to the following rules:
    sub = initialization of the microstrip substrate that includes important material properties
    rrfc = initialization of important design parameters (e.g., port numbers, sides, seed number, etc.
    y_dim, x_dim = y and x dimensions of the design. These should be integer (ideally odd integer) multiples of the pixelSize
    imp = characteristic impedance of the launch paths from the ports

    The output of the script is a pixel map, which is a CSV file with ones and zeros in the positions of the pixels in the grid, 
    and a GDS file that can be imported into an EM simulator. It will also output an array with port positions that can be used 
    for placing ports in an ADS simulation. 
    """

    # Pixel Size
    launch_x_pixels, launch_y_pixels, width_launch, length_launch = self._calc_launch_pixels(sub, imp)

    # Set horizontal and vertical pixel limits
    y_pixels = math.ceil(y_dim / self.pixelSize) 
    x_pixels = math.ceil(x_dim / self.pixelSize)
    if y_pixels % 2 == 0 and launch_x_pixels % 2 !=0:
      y_pixels += 1
    #if x_pixels % 2 == 0:
    #  x_pixels += 1

    rd1 = RandomDesign(self.unit, 
                       x_pixels, 
                       y_pixels, 
                       self.ports,
                       self.sides, 
                       launch_x_pixels, 
                       launch_y_pixels, 
                       self.pixelSize, 
                       self.scale, 
                       self.sim,
                       self.sym,
                       self.seed)
    x_total, y_total, portPos = RandomDesign.gen_port_pos(rd1)
    pixelMap = np.zeros((int(self.scale * x_total / self.pixelSize), int(self.scale * y_total/self.pixelSize)),dtype=int)
    pixelMap = RandomDesign.add_launch(rd1, x_total, y_total, pixelMap)
    pixelMap = RandomDesign.random_pix_rect(rd1, x_dim, y_dim, pixelMap)

    pixMap = (np.transpose(pixelMap))

    rows = np.size(pixMap, 0)
    cols = np.size(pixMap, 1)

    portPosScale = [x / self.scale for x in portPos]
    c12, c13, c14, c23, c24, c34 = Connectivity.find_connectivity(pixelMap, self.pixelSize,\
                                   self.ports, self.sides, portPosScale, connectivity=1)
    connections = [c12, c13, c14, c23, c24, c34]

    if (self.connect == 0 or self.connect == connections) and self.write == True:
      csvFile = self.outF + ".csv"
      # Export Pixel Map file
      np.savetxt(csvFile, np.flipud(np.transpose(pixelMap)), fmt = '%d', delimiter = ",")
    else:
        csvFile = ''

    return portPos, x_total, y_total, csvFile, launch_x_pixels

  def create_gds_from_pixmap_file(self,sub,imp,gnd):

    pixMap = np.flipud(np.loadtxt(self.outF, delimiter=','))
    rows = np.size(pixMap, 0)
    cols = np.size(pixMap, 1)

    if self.sim == 'EMX':

      if self.portPos == '':
        # In EMX sims, ports are defined on metal layers. They are added here, so port position needs to be computed
        launch_x_pixels, launch_y_pixels, width_launch, length_launch = self._calc_launch_pixels(sub, imp)

        rd1 = RandomDesign(
          self.unit, 
          cols-2*launch_x_pixels, 
          rows, 
          self.ports,
          self.sides, 
          launch_x_pixels, 
          launch_y_pixels, 
          self.pixelSize,
          self.sim,
          self.sym,
          self.seed)
        
        x_total, y_total, portPos = rd1.gen_port_pos(rd1)
      else:

        portPos = self.portPos
        x_total = rows*self.pixelSize/self.scale
        y_total = cols*self.pixelSize/self.scale
        
      ly = pya.Layout()

      # Set the database unit to 1 um. Generally for EMX it is easier to work in um.
      ly.dbu = self.unit*1e6

      # Create Cell obj
      cellName = 'INVDESIGN'
      UNIT = ly.create_cell(cellName)

      # Create layer #'s
      l_top = ly.layer(59, 0) # layer for signal metal, corresponds to OI in FDX process
      l_bottom = ly.layer(15, 0) # Ground Layer corresponds to M1 in FDX process

      if gnd == 'mesh':
        # Draw outline
        # Make a mesh shield
        # Assume ground hash is 0.5um strip followed by 1um space, place it on the 
        # bottom metal layer
        for x in range(0, int(self.scale*x_total/self.pixelSize)):
          rect = UNIT.shapes(l_bottom).insert( pya.Box(0, 0, self.scale*y_total, self.scale).moved(0,(x+0.5)*self.pixelSize-0.5*self.scale))
        for x in range(0, int(self.scale*y_total/self.pixelSize)):
          rect = UNIT.shapes(l_bottom).insert( pya.Box(0, 0, self.scale, self.scale*x_total).moved((x+0.5)*self.pixelSize-0.5*self.scale,0))
      else: #flood
        for x in range(0, int(self.scale*x_total/self.pixelSize)):
          rect = UNIT.shapes(l_bottom).insert( pya.Box(0, 0, self.scale*y_total, self.pixelSize/3).moved(0,(x+0.5)*self.pixelSize-0.5*self.pixelSize/3))
        for x in range(0, int(self.scale*y_total/self.pixelSize)-1):
          rect = UNIT.shapes(l_bottom).insert( pya.Box(0, 0, self.pixelSize, self.scale*x_total).moved((x+0.5)*self.pixelSize-0.0*self.scale,0))
    
      if self.corner == 'overlap':
        for x in range(0, cols):
          for y in range(0, rows):
            if pixMap[y,x]  == 1:
              rect = pya.Box(0, 0, self.pixelSize, self.pixelSize).moved(x*self.pixelSize,\
                     y*self.pixelSize)
              UNIT.shapes(l_top).insert(rect)
              if gnd == 'cell':
                UNIT.shapes(l_bottom).insert(rect)
   
        # This defines a diamond shaped polygon that will be added to the corners of adjoining polygons to
        # guarantee overlap in the structure. This simulates under-etch of the metal in the corners so that
        # connectivity in the corners and its impact can be estimated.
        points = [pya.DPoint(-self.pixelSize/(10*self.scale), 0),
                  pya.DPoint(0, self.pixelSize/(10*self.scale)),
                  pya.DPoint(self.pixelSize/(10*self.scale), 0),
                  pya.DPoint(0,-self.pixelSize/(10*self.scale))]
        for x in range(0,cols-1):
          for y in range(0,rows-1):
            if pixMap[y,x] == 1 and pixMap[y+1,x+1] == 1:
              diam = pya.DPolygon(points).moved((x+1)*self.pixelSize/self.scale,(y+1)*self.pixelSize/self.scale)
              UNIT.shapes(l_top).insert(diam)
        
        for x in range(1,cols):
          for y in range(0,rows-1):
            if pixMap[y,x] == 1 and pixMap[y+1,x-1] == 1:
              diam = pya.DPolygon(points).moved((x)*self.pixelSize/self.scale,(y+1)*self.pixelSize/self.scale)
              UNIT.shapes(l_top).insert(diam)

      # When gds is fabricated it naturally creates overlap in the corners because of underetch of the metal. This routine 
      # creates non-overlapping polygons with enough space for manufacture. Space should be 6mil for PCB and typically 2um
      # in thick metal for ICs (A 0 in the pixel map corresponds to no metal and a 1 corresponds to metal
      elif self.corner == 'noverlap': 
        caseMap = pixMap
        for x in range(0,cols):
          for y in range(0,rows):
            if y > 0 and x > 0 and y < rows-1 and x < cols-1: #Do inner portions of map first, edges will be done last
              # First, locate all 1 pixels that are surrounded on all edges by 0 pixels
              # D0D
              # 010
              # D0D
              if pixMap[y,x] == 1 and pixMap[y,x-1] == 0 and pixMap[y,x+1] == 0 and pixMap[y-1,x] == 0 and pixMap[y+1,x] == 0:
                # First, find all pixels that have 0's adjacent and 1's on all diagonals
                if pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] != 0:
                  caseMap[y,x] = 2
                # First, find all pixels that have 0's adjacent and 1's on 3/4 diagonals
                elif pixMap[y+1,x+1] == 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] != 0:
                  caseMap[y,x] = 3
                elif pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] == 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] != 0:
                  caseMap[y,x] = 4
                elif pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] == 0:
                  caseMap[y,x] = 5
                elif pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] == 0 and pixMap[y-1,x-1] != 0:
                  caseMap[y,x] = 6
                # First, find all pixels that have 0's adjacent and 1's on 2/4 diagonals
                elif pixMap[y+1,x+1] == 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] == 0:
                  caseMap[y,x] = 7
                elif pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] == 0 and pixMap[y-1,x+1] == 0 and pixMap[y-1,x-1] != 0:
                  caseMap[y,x] = 8
              # First, locate all 1 pixels that are surrounded on right edge by 0 pixels
              # DDD
              # D10
              # DDD
              if pixMap[y,x] == 1 and pixMap[y,x+1] == 0:
                if pixMap[y+1,x+1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y+1,x] == 0 and pixMap[y-1,x] == 0:
                  caseMap[y,x] = 15
              # First, locate all 1 pixels that are surrounded on left edge by 0 pixels
              # DDD
              # 01D
              # DDD
              if pixMap[y,x] == 1 and pixMap[y,x-1] == 0:
                if pixMap[y+1,x-1] != 0 and pixMap[y-1,x-1] != 0 and pixMap[y+1,x] == 0 and pixMap[y-1,x] == 0:
                  caseMap[y,x] = 16
              # First, locate all 1 pixels that are surrounded on top edge by 0 pixels
              # D0D
              # D1D
              # DDD
              if pixMap[y,x] == 1 and pixMap[y+1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                # First, find all pixels that have 0's adjacent and 1's on top 2 diagonals
                if pixMap[y+1,x+1] != 0 and pixMap[y+1, x-1] != 0 and pixMap[y,x+1] == 0 and pixMap[y,x-1] == 0:
                  caseMap[y,x] = 9
                # First, find all pixels that have 0 right adjacent and 1's on top right diagonal
                elif pixMap[y+1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 10
                # First, find all pixels that have 0 left adjacent and 1's on top left diagonal
                elif pixMap[y+1,x-1] != 0 and pixMap[y,x-1] == 0: 
                  caseMap[y,x] = 11
              # First, locate all 1 pixels that are surrounded on bottom edge by 0 pixels
              # DDD
              # D1D
              # D0D
              if pixMap[y,x] == 1 and pixMap[y-1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                # First, find all pixels that have 0's adjacent and 1's on bottom 2 diagonals
                if pixMap[y-1,x+1] != 0 and pixMap[y-1, x-1] != 0 and pixMap[y,x+1] == 0 and pixMap[y,x-1] == 0:
                  caseMap[y,x] = 12
                # First, find all pixels that have 0 right adjacent and 1's on bottom right diagonal
                elif pixMap[y-1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 13
                # First, find all pixels that have 0 left adjacent and 1's on bottom left diagonal
                elif pixMap[y-1,x-1] != 0 and pixMap[y,x-1] == 0: 
                  caseMap[y,x] = 14
            if y == 0 and x > 0 and x < cols-1: #Do bottom row
              # First, locate all 1 pixels that are surrounded on top edge by 0 pixels
              # D0D
              # D1D
              # DDD
              if pixMap[y,x] == 1 and pixMap[y+1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                # First, find all pixels that have 0's adjacent and 1's on top 2 diagonals
                if pixMap[y+1,x+1] != 0 and pixMap[y+1, x-1] != 0 and pixMap[y,x+1] == 0 and pixMap[y,x-1] == 0:
                  caseMap[y,x] = 9
                # First, find all pixels that have 0 right adjacent and 1's on top right diagonal
                elif pixMap[y+1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 10
                # First, find all pixels that have 0 left adjacent and 1's on top left diagonal
                elif pixMap[y+1,x-1] != 0 and pixMap[y,x-1] == 0: 
                  caseMap[y,x] = 11
            if y == rows-1 and x > 0 and x < cols-1: #Do top row
              if pixMap[y,x] == 1 and pixMap[y-1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                # First, find all pixels that have 0's adjacent and 1's on bottom 2 diagonals
                if pixMap[y-1,x+1] != 0 and pixMap[y-1, x-1] != 0 and pixMap[y,x+1] == 0 and pixMap[y,x-1] == 0:
                  caseMap[y,x] = 12
                # First, find all pixels that have 0 right adjacent and 1's on bottom right diagonal
                elif pixMap[y-1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 13
                # First, find all pixels that have 0 left adjacent and 1's on bottom left diagonal
                elif pixMap[y-1,x-1] != 0 and pixMap[y,x-1] == 0: 
                  caseMap[y,x] = 14
            if y > 0 and x == 0 and y < rows-1: #Do first column
              if pixMap[y,x] == 1 and pixMap[y+1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                # First, find all pixels that have 0 right adjacent and 1's on top right diagonal
                if pixMap[y+1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 10
                # First, find all pixels that have 0 right adjacent and 1's on bottom right diagonal
            if pixMap[y,x] == 1 and pixMap[y-1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
              if pixMap[y-1,x+1] != 0 and pixMap[y,x+1] == 0: 
                caseMap[y,x] = 13
            if pixMap[y,x] == 1 and pixMap[y,x+1] == 0:
              if pixMap[y+1,x+1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y+1,x] == 0 and pixMap[y-1,x] == 0:
                caseMap[y,x] = 15

        for x in range(0,cols):
          for y in range(0,rows):
            if caseMap[y,x] == 1:
              rect = pya.Box(0, 0, self.pixelSize, self.pixelSize).moved(x*self.pixelSize,\
                     y*self.pixelSize)
              UNIT.shapes(l_top).insert(rect)
              if gnd == 'cell':
                UNIT.shapes(l_bottom).insert(rect)
            if caseMap[y,x] == 2:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
	                (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
	                ((self.minPix/2)*np.sqrt(2),self.pixelSize),
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = pya.DPolygon(points).moved(x*self.pixelSize,y*self.pixelSize)
              UNIT.shapes(l_top).insert(poly)
              if gnd == 'cell':
                UNIT.shapes(l_bottom).insert(rect)
            if caseMap[y,x] == 3:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2),self.pixelSize),
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = pya.DPolygon(points).moved(x*self.pixelSize,y*self.pixelSize)
              UNIT.shapes(l_top).insert(poly)
              if gnd == 'cell':
                UNIT.shapes(l_bottom).insert(rect)
            if caseMap[y,x] == 4:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),self.pixelSize),
                        (0, self.pixelSize), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = pya.DPolygon(points).moved(x*self.pixelSize,y*self.pixelSize)
              UNIT.shapes(l_top).insert(poly)
              if gnd == 'cell':
                UNIT.shapes(l_bottom).insert(rect)
            if caseMap[y,x] == 5:
              points = [(0,0),                                            #  _ 7-sided polygon 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0),                  # / \
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)),                   #|   |
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)),    #|___/
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),self.pixelSize),
                        ((self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2))]
              poly = pya.DPolygon(points).moved(x*self.pixelSize,y*self.pixelSize)
              UNIT.shapes(l_top).insert(poly)
              if gnd == 'cell':
                UNIT.shapes(l_bottom).insert(rect)
            if caseMap[y,x] == 6:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize,0), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2),self.pixelSize),
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = pya.DPolygon(points).moved(x*self.pixelSize,y*self.pixelSize)
              UNIT.shapes(l_top).insert(poly)
              if gnd == 'cell':
                UNIT.shapes(l_bottom).insert(rect)
            if caseMap[y,x] == 7:
              points = [(0,0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2), self.pixelSize),
                        (0,self.pixelSize-(self.minPix/2)*np.sqrt(2))]
              poly = pya.DPolygon(points).moved(x*self.pixelSize,y*self.pixelSize)
              UNIT.shapes(l_top).insert(poly)
              if gnd == 'cell':
                UNIT.shapes(l_bottom).insert(rect)
          if caseMap[y,x] == 8:
            points = [((self.minPix/2)*np.sqrt(2),0), 
                      (self.pixelSize,0), 
                      (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                      (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                      (0, self.pixelSize),
                      (0,(self.minPix/2)*np.sqrt(2))]
            poly = pya.DPolygon(points).moved(x*self.pixelSize,y*self.pixelSize)
            UNIT.shapes(l_top).insert(poly)
            if gnd == 'cell':
              UNIT.shapes(l_bottom).insert(rect)
          if caseMap[y,x] == 9:
            points = [(0,0), 
                      (self.pixelSize,0), 
                      (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                      (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                      ((self.minPix/2)*np.sqrt(2), self.pixelSize), 
                      (0, self.pixelSize-(self.minPix/2)*np.sqrt(2))]
            poly = pya.DPolygon(points).moved(x*self.pixelSize,y*self.pixelSize)
            UNIT.shapes(l_top).insert(poly)
            if gnd == 'cell':
              UNIT.shapes(l_bottom).insert(rect)
          if caseMap[y,x] == 10:
            points = [(0,0), 
              (self.pixelSize,0), 
              (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
              (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
              (0, self.pixelSize)]
            poly = pya.DPolygon(points).moved(x*self.pixelSize,y*self.pixelSize)
            UNIT.shapes(l_top).insert(poly)
            if gnd == 'cell':
              UNIT.shapes(l_bottom).insert(rect)
          if caseMap[y,x] == 11:
            points = [(0,0), 
              (self.pixelSize,0), 
              (self.pixelSize, self.pixelSize), 
              ((self.minPix/2)*np.sqrt(2), self.pixelSize), 
              (0, self.pixelSize-(self.minPix/2)*np.sqrt(2))]
            poly = pya.DPolygon(points).moved(x*self.pixelSize,y*self.pixelSize)
            UNIT.shapes(l_top).insert(poly)
            if gnd == 'cell':
              UNIT.shapes(l_bottom).insert(rect)
          if caseMap[y,x] == 12:
            points = [((self.minPix/2)*np.sqrt(2),0), 
              (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
              (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
              (self.pixelSize, self.pixelSize), 
              (0, self.pixelSize), 
              (0, (self.minPix/2)*np.sqrt(2))]
            poly = pya.DPolygon(points).moved(x*self.pixelSize,y*self.pixelSize)
            UNIT.shapes(l_top).insert(poly)
            if gnd == 'cell':
              UNIT.shapes(l_bottom).insert(rect)
          if caseMap[y,x] == 13:
            points = [(0,0), 
              (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
              (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
              (self.pixelSize, self.pixelSize), 
              (0, self.pixelSize)]
            poly = pya.DPolygon(points).moved(x*self.pixelSize,y*self.pixelSize)
            UNIT.shapes(l_top).insert(poly)
            if gnd == 'cell':
              UNIT.shapes(l_bottom).insert(rect)
          if caseMap[y,x] == 14:
            points = [((self.minPix/2)*np.sqrt(2),0), 
              (self.pixelSize,0), 
              (self.pixelSize, self.pixelSize), 
              (0, self.pixelSize), 
              (0, (self.minPix/2)*np.sqrt(2))]
            poly = pya.DPolygon(points).moved(x*self.pixelSize,y*self.pixelSize)
            UNIT.shapes(l_top).insert(poly)
            if gnd == 'cell':
              UNIT.shapes(l_bottom).insert(rect)
          if caseMap[y,x] == 15:
            points = [(0,0), 
              (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
              (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
              (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
              (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
              (0, self.pixelSize)]
            poly = pya.DPolygon(points).moved(x*self.pixelSize,y*self.pixelSize)
            UNIT.shapes(l_top).insert(poly)
            if gnd == 'cell':
              UNIT.shapes(l_bottom).insert(rect)
          if caseMap[y,x] == 16:
            points = [((self.minPix/2)*np.sqrt(2),0), 
              (self.pixelSize,0), 
              (self.pixelSize, self.pixelSize), 
              ((self.minPix/2)*np.sqrt(2), self.pixelSize), 
              (0, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
              (0, (self.minPix/2)*np.sqrt(2))]
            poly = pya.DPolygon(points).moved(x*self.pixelSize,y*self.pixelSize)
            UNIT.shapes(l_top).insert(poly)
            if gnd == 'cell':
              UNIT.shapes(l_bottom).insert(rect)

      # This draws a gds of a pixelmap with no overlap or rounding in the corners. 
      else:
        for x in range(0, cols):
          for y in range(0, rows):
            if pixMap[y,x]  == 1:
              rect = pya.Box(0, 0, self.pixelSize, self.pixelSize).moved(x*self.pixelSize,\
                     y*self.pixelSize)
              UNIT.shapes(l_top).insert(rect)
              if gnd == 'cell':
                UNIT.shapes(l_bottom).insert(rect)

      if self.ports == 1:
        port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', portPos[0]/self.scale, portPos[1]/self.scale) )
        port_2 = UNIT.shapes(l_bottom).insert( pya.Text('p2', portPos[0]/self.scale, portPos[1]/self.scale) ) # Ground Port
      elif self.ports == 2:
        port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', portPos[0]/self.scale, portPos[1]/self.scale) )
        port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', portPos[2]/self.scale, portPos[3]/self.scale) )
        port_3 = UNIT.shapes(l_bottom).insert( pya.Text('p3', portPos[0]/self.scale, portPos[1]/self.scale) ) # Ground Port
        port_4 = UNIT.shapes(l_bottom).insert( pya.Text('p4', portPos[2]/self.scale, portPos[3]/self.scale) ) # Ground Port
      ### 7/10/2023 Still need to fix positioning issue for cases with 2, 3, 4 sides
      elif self.ports == 3:
        port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', portPos[0]/self.scale, portPos[1]/self.scale) )
        port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', portPos[2]/self.scale, portPos[3]/self.scale) )
        port_3 = UNIT.shapes(l_top).insert( pya.Text('p3', portPos[4]/self.scale, portPos[5]/self.scale) )
        port_4 = UNIT.shapes(l_bottom).insert( pya.Text('p4', portPos[0]/self.scale, portPos[1]/self.scale) ) # Ground Port
        port_5 = UNIT.shapes(l_bottom).insert( pya.Text('p5', portPos[2]/self.scale, portPos[3]/self.scale) ) # Ground Port
        port_6 = UNIT.shapes(l_bottom).insert( pya.Text('p6', portPos[4]/self.scale, portPos[5]/self.scale) ) # Ground Port
      elif self.ports == 4:
        port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', portPos[0]/self.scale, portPos[1]/self.scale) )
        port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', portPos[2]/self.scale, portPos[3]/self.scale) )
        port_3 = UNIT.shapes(l_top).insert( pya.Text('p3', portPos[4]/self.scale, portPos[5]/self.scale) )
        port_4 = UNIT.shapes(l_top).insert( pya.Text('p4', portPos[6]/self.scale, portPos[7]/self.scale) )
        port_5 = UNIT.shapes(l_bottom).insert( pya.Text('p5', portPos[0]/self.scale, portPos[1]/self.scale) ) # Ground Port
        port_6 = UNIT.shapes(l_bottom).insert( pya.Text('p6', portPos[2]/self.scale, portPos[3]/self.scale) ) # Ground Port
        port_7 = UNIT.shapes(l_bottom).insert( pya.Text('p7', portPos[4]/self.scale, portPos[5]/self.scale) ) # Ground Port
        port_8 = UNIT.shapes(l_bottom).insert( pya.Text('p8', portPos[6]/self.scale, portPos[7]/self.scale) ) # Ground Port

      gdsFile = self.outF.replace('csv', 'gds')
      # Export GDS

      ly.write(gdsFile)
      # View GDS      
      
    else:
      lib = gdspy.GdsLibrary()

      # Set the database unit to 1 mil
      lib.unit = self.unit

      # Create Cell obj
      cellName = 'INVDESIGN'
      gdspy.current_library = gdspy.GdsLibrary() # This line of code has to be here to reset the GDS library on every loop
      UNIT = lib.new_cell(cellName)

      l_bottom = {"layer": 10, "datatype": 0}
      l_top = {"layer": 11, "datatype": 0}
      l_sources = {"layer": 5, "datatype": 0}

      outline = gdspy.Rectangle((0, 0), (cols*self.pixelSize, rows*self.pixelSize), **l_bottom)
      UNIT.add(outline) 

      # When gds is fabricated it naturally creates overlap in the corners because of underetch of the metal. This routine 
      # creates non-overlapping polygons with enough space for manufacture. Space should be 6mil for PCB and typically 2um
      # in thick metal for ICs (A 0 in the pixel map corresponds to no metal and a 1 corresponds to metal
      if self.corner == 'noverlap': 
          caseMap = pixMap
          for x in range(0,cols):
            for y in range(0,rows):
              if y > 0 and x > 0 and y < rows-1 and x < cols-1: #Do inner portions of map first, edges will be done last
                # First, locate all 1 pixels that are surrounded on all edges by 0 pixels
                # D0D
                # 010
                # D0D
                if pixMap[y,x] == 1 and pixMap[y,x-1] == 0 and pixMap[y,x+1] == 0 and pixMap[y-1,x] == 0 and pixMap[y+1,x] == 0:
                  # First, find all pixels that have 0's adjacent and 1's on all diagonals
                  if pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] != 0:
                    caseMap[y,x] = 2
                  # First, find all pixels that have 0's adjacent and 1's on 3/4 diagonals
                  elif pixMap[y+1,x+1] == 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] != 0:
                    caseMap[y,x] = 3
                  elif pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] == 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] != 0:
                    caseMap[y,x] = 4
                  elif pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] == 0:
                    caseMap[y,x] = 5
                  elif pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] == 0 and pixMap[y-1,x-1] != 0:
                    caseMap[y,x] = 6
                  # First, find all pixels that have 0's adjacent and 1's on 2/4 diagonals
                  elif pixMap[y+1,x+1] == 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] == 0:
                    caseMap[y,x] = 7
                  elif pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] == 0 and pixMap[y-1,x+1] == 0 and pixMap[y-1,x-1] != 0:
                    caseMap[y,x] = 8
                # First, locate all 1 pixels that are surrounded on right edge by 0 pixels
                # DDD
                # D10
                # DDD
                if pixMap[y,x] == 1 and pixMap[y,x+1] == 0:
                  if pixMap[y+1,x+1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y+1,x] == 0 and pixMap[y-1,x] == 0:
                    caseMap[y,x] = 15
                # First, locate all 1 pixels that are surrounded on left edge by 0 pixels
                # DDD
                # 01D
                # DDD
                if pixMap[y,x] == 1 and pixMap[y,x-1] == 0:
                  if pixMap[y+1,x-1] != 0 and pixMap[y-1,x-1] != 0 and pixMap[y+1,x] == 0 and pixMap[y-1,x] == 0:
                    caseMap[y,x] = 16
                # First, locate all 1 pixels that are surrounded on top edge by 0 pixels
                # D0D
                # D1D
                # DDD
                if pixMap[y,x] == 1 and pixMap[y+1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                  # First, find all pixels that have 0's adjacent and 1's on top 2 diagonals
                  if pixMap[y+1,x+1] != 0 and pixMap[y+1, x-1] != 0 and pixMap[y,x+1] == 0 and pixMap[y,x-1] == 0:
                    caseMap[y,x] = 9
                  # First, find all pixels that have 0 right adjacent and 1's on top right diagonal
                  elif pixMap[y+1,x+1] != 0 and pixMap[y,x+1] == 0: 
                    caseMap[y,x] = 10
                  # First, find all pixels that have 0 left adjacent and 1's on top left diagonal
                  elif pixMap[y+1,x-1] != 0 and pixMap[y,x-1] == 0: 
                    caseMap[y,x] = 11
                # First, locate all 1 pixels that are surrounded on bottom edge by 0 pixels
                # DDD
                # D1D
                # D0D
                if pixMap[y,x] == 1 and pixMap[y-1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                  # First, find all pixels that have 0's adjacent and 1's on bottom 2 diagonals
                  if pixMap[y-1,x+1] != 0 and pixMap[y-1, x-1] != 0 and pixMap[y,x+1] == 0 and pixMap[y,x-1] == 0:
                    caseMap[y,x] = 12
                  # First, find all pixels that have 0 right adjacent and 1's on bottom right diagonal
                  elif pixMap[y-1,x+1] != 0 and pixMap[y,x+1] == 0: 
                    caseMap[y,x] = 13
                  # First, find all pixels that have 0 left adjacent and 1's on bottom left diagonal
                  elif pixMap[y-1,x-1] != 0 and pixMap[y,x-1] == 0: 
                    caseMap[y,x] = 14
              if y == 0 and x > 0 and x < cols-1: #Do bottom row
                # First, locate all 1 pixels that are surrounded on top edge by 0 pixels
                # D0D
                # D1D
                # DDD
                if pixMap[y,x] == 1 and pixMap[y+1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                  # First, find all pixels that have 0's adjacent and 1's on top 2 diagonals
                  if pixMap[y+1,x+1] != 0 and pixMap[y+1, x-1] != 0 and pixMap[y,x+1] == 0 and pixMap[y,x-1] == 0:
                    caseMap[y,x] = 9
                  # First, find all pixels that have 0 right adjacent and 1's on top right diagonal
                  elif pixMap[y+1,x+1] != 0 and pixMap[y,x+1] == 0: 
                    caseMap[y,x] = 10
                  # First, find all pixels that have 0 left adjacent and 1's on top left diagonal
                  elif pixMap[y+1,x-1] != 0 and pixMap[y,x-1] == 0: 
                    caseMap[y,x] = 11
              if y == rows-1 and x > 0 and x < cols-1: #Do top row
                if pixMap[y,x] == 1 and pixMap[y-1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                  # First, find all pixels that have 0's adjacent and 1's on bottom 2 diagonals
                  if pixMap[y-1,x+1] != 0 and pixMap[y-1, x-1] != 0 and pixMap[y,x+1] == 0 and pixMap[y,x-1] == 0:
                    caseMap[y,x] = 12
                  # First, find all pixels that have 0 right adjacent and 1's on bottom right diagonal
                  elif pixMap[y-1,x+1] != 0 and pixMap[y,x+1] == 0: 
                    caseMap[y,x] = 13
                  # First, find all pixels that have 0 left adjacent and 1's on bottom left diagonal
                  elif pixMap[y-1,x-1] != 0 and pixMap[y,x-1] == 0: 
                    caseMap[y,x] = 14
              if y > 0 and x == 0 and y < rows-1: #Do first column
                if pixMap[y,x] == 1 and pixMap[y+1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                  # First, find all pixels that have 0 right adjacent and 1's on top right diagonal
                  if pixMap[y+1,x+1] != 0 and pixMap[y,x+1] == 0: 
                    caseMap[y,x] = 10
                  # First, find all pixels that have 0 right adjacent and 1's on bottom right diagonal
                if pixMap[y,x] == 1 and pixMap[y-1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                  if pixMap[y-1,x+1] != 0 and pixMap[y,x+1] == 0: 
                    caseMap[y,x] = 13
                if pixMap[y,x] == 1 and pixMap[y,x+1] == 0:
                  if pixMap[y+1,x+1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y+1,x] == 0 and pixMap[y-1,x] == 0:
                    caseMap[y,x] = 15
          for x in range(0,cols):
            for y in range(0,rows):
              if caseMap[y,x] == 1:
                rect = gdspy.Rectangle((0, 0), (self.pixelSize, self.pixelSize), **l_top).translate(x*self.pixelSize,\
                       y*self.pixelSize)
                UNIT.add(rect)
              if caseMap[y,x] == 2:
                points = [((self.minPix/2)*np.sqrt(2),0), 
                          (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                          (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                          (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                          (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                          ((self.minPix/2)*np.sqrt(2),self.pixelSize),
                          (0, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                          (0, (self.minPix/2)*np.sqrt(2))]
                poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
                UNIT.add(poly)
              if caseMap[y,x] == 3:
                points = [((self.minPix/2)*np.sqrt(2),0), 
                          (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                          (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                          (self.pixelSize, self.pixelSize), 
                          ((self.minPix/2)*np.sqrt(2),self.pixelSize),
                          (0, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                          (0, (self.minPix/2)*np.sqrt(2))]
                poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
                UNIT.add(poly)
              if caseMap[y,x] == 4:
                points = [((self.minPix/2)*np.sqrt(2),0), 
                          (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                          (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                          (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                          (self.pixelSize-(self.minPix/2)*np.sqrt(2),self.pixelSize),
                          (0, self.pixelSize), 
                          (0, (self.minPix/2)*np.sqrt(2))]
                poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
                UNIT.add(poly)
              if caseMap[y,x] == 5:
                points = [(0,0),                                            #  _ 7-sided polygon 
                          (self.pixelSize-(self.minPix/2)*np.sqrt(2),0),                  # / \
                          (self.pixelSize, (self.minPix/2)*np.sqrt(2)),                   #|   |
                          (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)),    #|___/
                          (self.pixelSize-(self.minPix/2)*np.sqrt(2),self.pixelSize),
                          ((self.minPix/2)*np.sqrt(2), self.pixelSize), 
                          (0, self.pixelSize-(self.minPix/2)*np.sqrt(2))]
                poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
                UNIT.add(poly)
              if caseMap[y,x] == 6:
                points = [((self.minPix/2)*np.sqrt(2),0), 
                          (self.pixelSize,0), 
                          (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                          (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                          ((self.minPix/2)*np.sqrt(2),self.pixelSize),
                          (0, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                          (0, (self.minPix/2)*np.sqrt(2))]
                poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
                UNIT.add(poly)
              if caseMap[y,x] == 7:
                points = [(0,0), 
                          (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                          (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                          (self.pixelSize, self.pixelSize), 
                          ((self.minPix/2)*np.sqrt(2), self.pixelSize),
                          (0,self.pixelSize-(self.minPix/2)*np.sqrt(2))]
                poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
                UNIT.add(poly)
              if caseMap[y,x] == 8:
                points = [((self.minPix/2)*np.sqrt(2),0), 
                          (self.pixelSize,0), 
                          (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                          (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                          (0, self.pixelSize),
                          (0,(self.minPix/2)*np.sqrt(2))]
                poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
                UNIT.add(poly)
              if caseMap[y,x] == 9:
                points = [(0,0), 
                          (self.pixelSize,0), 
                          (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                          (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                          ((self.minPix/2)*np.sqrt(2), self.pixelSize), 
                          (0, self.pixelSize-(self.minPix/2)*np.sqrt(2))]
                poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
                UNIT.add(poly)
              if caseMap[y,x] == 10:
                points = [(0,0), 
                          (self.pixelSize,0), 
                          (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                          (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                          (0, self.pixelSize)]
                poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
                UNIT.add(poly)
              if caseMap[y,x] == 11:
                points = [(0,0), 
                          (self.pixelSize,0), 
                          (self.pixelSize, self.pixelSize), 
                          ((self.minPix/2)*np.sqrt(2), self.pixelSize), 
                          (0, self.pixelSize-(self.minPix/2)*np.sqrt(2))]
                poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
                UNIT.add(poly)
              if caseMap[y,x] == 12:
                points = [((self.minPix/2)*np.sqrt(2),0), 
                          (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                          (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                          (self.pixelSize, self.pixelSize), 
                          (0, self.pixelSize), 
                          (0, (self.minPix/2)*np.sqrt(2))]
                poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
                UNIT.add(poly)
              if caseMap[y,x] == 13:
                points = [(0,0), 
                          (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                          (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                          (self.pixelSize, self.pixelSize), 
                          (0, self.pixelSize)]
                poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
                UNIT.add(poly)
              if caseMap[y,x] == 14:
                points = [((self.minPix/2)*np.sqrt(2),0), 
                          (self.pixelSize,0), 
                          (self.pixelSize, self.pixelSize), 
                          (0, self.pixelSize), 
                          (0, (self.minPix/2)*np.sqrt(2))]
                poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
                UNIT.add(poly)
              if caseMap[y,x] == 15:
                points = [(0,0), 
                          (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                          (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                          (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                          (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                          (0, self.pixelSize)]
                poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
                UNIT.add(poly)
              if caseMap[y,x] == 16:
                points = [((self.minPix/2)*np.sqrt(2),0), 
                          (self.pixelSize,0), 
                          (self.pixelSize, self.pixelSize), 
                          ((self.minPix/2)*np.sqrt(2), self.pixelSize), 
                          (0, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                          (0, (self.minPix/2)*np.sqrt(2))]
                poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
                UNIT.add(poly)
        
      elif self.corner == 'overlap':
        for x in range(0, cols):
          for y in range(0, rows):
            if pixMap[y,x]  == 1:
              rect = gdspy.Rectangle((0, 0), (self.pixelSize, self.pixelSize), **l_top).translate(x*self.pixelSize,\
                     y*self.pixelSize)
              UNIT.add(rect)

        for x in range(0,cols-1):
          for y in range(0,rows-1):
            if pixMap[y,x] == 1 and pixMap[y+1,x+1] == 1:
              diam = gdspy.Rectangle((-0.707*self.pixelSize/10, -0.707*self.pixelSize/10), \
                                     (0.707*self.pixelSize/10, 0.707*self.pixelSize/10), **l_top).rotate(np.pi/4).\
                                     translate((x+1)*self.pixelSize,(y+1)*self.pixelSize)
              UNIT.add(diam)
        
        for x in range(1,cols):
          for y in range(0,rows-1):
            if pixMap[y,x] == 1 and pixMap[y+1,x-1] == 1:
              diam = gdspy.Rectangle((-0.707*self.pixelSize/10, -0.707*self.pixelSize/10), \
                                     (0.707*self.pixelSize/10, 0.707*self.pixelSize/10), **l_top).rotate(np.pi/4).\
                                     translate((x)*self.pixelSize,(y+1)*self.pixelSize)
              UNIT.add(diam)
      else: # Mo Geometric Modifications
        for x in range(0, cols):
          for y in range(0, rows):
            if pixMap[y,x]  == 1:
              rect = gdspy.Rectangle((0, 0), (self.pixelSize, self.pixelSize), **l_top).translate(x*self.pixelSize,\
                     y*self.pixelSize)
              UNIT.add(rect)

      gdsFile = self.outF.replace('csv', 'gds') 
      ##print(gdsFile) #debug 
      #if self.view == True:
      #  gdspy.LayoutViewer(lib) 
      ## Export GDS
      lib.write_gds(gdsFile)
      
    if self.view == True:
      # Load a GDSII file into a new library
      gdsii = gdspy.GdsLibrary(infile=gdsFile)
      gdspy.LayoutViewer(gdsii)    
    return cellName

  def recreate_gds_array(self,pixMap):

    lib = gdspy.GdsLibrary()

    # Set the database unit to 1 mil
    lib.unit = self.unit

    # Create Cell obj
    cellName = 'INVDESIGN'
    gdspy.current_library = gdspy.GdsLibrary() # This line of code has to be here to reset the GDS library on every loop
    UNIT = lib.new_cell(cellName)

    l_bottom = {"layer": 10, "datatype": 0}
    l_top = {"layer": 11, "datatype": 0}
    l_sources = {"layer": 5, "datatype": 0}

    rows = np.size(pixMap, 0)
    cols = np.size(pixMap, 1)
    print(rows, cols)

    outline = gdspy.Rectangle((0, 0), (cols*self.pixelSize, rows*self.pixelSize), **l_bottom)
    UNIT.add(outline) 

    # When gds is fabricated it naturally creates overlap in the corners because of underetch of the metal. This routine 
    # creates non-overlapping polygons with enough space for manufacture. Space should be 6mil for PCB and typically 2um
    # in thick metal for ICs (A 0 in the pixel map corresponds to no metal and a 1 corresponds to metal
    if self.corner == 'noverlap': 
        caseMap = pixMap
        for x in range(0,cols):
          for y in range(0,rows):
            if y > 0 and x > 0 and y < rows-1 and x < cols-1: #Do inner portions of map first, edges will be done last
              # First, locate all 1 pixels that are surrounded on all edges by 0 pixels
              # D0D
              # 010
              # D0D
              if pixMap[y,x] == 1 and pixMap[y,x-1] == 0 and pixMap[y,x+1] == 0 and pixMap[y-1,x] == 0 and pixMap[y+1,x] == 0:
                # First, find all pixels that have 0's adjacent and 1's on all diagonals
                if pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] != 0:
                  caseMap[y,x] = 2
                # First, find all pixels that have 0's adjacent and 1's on 3/4 diagonals
                elif pixMap[y+1,x+1] == 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] != 0:
                  caseMap[y,x] = 3
                elif pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] == 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] != 0:
                  caseMap[y,x] = 4
                elif pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] == 0:
                  caseMap[y,x] = 5
                elif pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] == 0 and pixMap[y-1,x-1] != 0:
                  caseMap[y,x] = 6
                # First, find all pixels that have 0's adjacent and 1's on 2/4 diagonals
                elif pixMap[y+1,x+1] == 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] == 0:
                  caseMap[y,x] = 7
                elif pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] == 0 and pixMap[y-1,x+1] == 0 and pixMap[y-1,x-1] != 0:
                  caseMap[y,x] = 8
              # First, locate all 1 pixels that are surrounded on right edge by 0 pixels
              # DDD
              # D10
              # DDD
              if pixMap[y,x] == 1 and pixMap[y,x+1] == 0:
                if pixMap[y+1,x+1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y+1,x] == 0 and pixMap[y-1,x] == 0:
                  caseMap[y,x] = 15
              # First, locate all 1 pixels that are surrounded on left edge by 0 pixels
              # DDD
              # 01D
              # DDD
              if pixMap[y,x] == 1 and pixMap[y,x-1] == 0:
                if pixMap[y+1,x-1] != 0 and pixMap[y-1,x-1] != 0 and pixMap[y+1,x] == 0 and pixMap[y-1,x] == 0:
                  caseMap[y,x] = 16
              # First, locate all 1 pixels that are surrounded on top edge by 0 pixels
              # D0D
              # D1D
              # DDD
              if pixMap[y,x] == 1 and pixMap[y+1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                # First, find all pixels that have 0's adjacent and 1's on top 2 diagonals
                if pixMap[y+1,x+1] != 0 and pixMap[y+1, x-1] != 0 and pixMap[y,x+1] == 0 and pixMap[y,x-1] == 0:
                  caseMap[y,x] = 9
                # First, find all pixels that have 0 right adjacent and 1's on top right diagonal
                elif pixMap[y+1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 10
                # First, find all pixels that have 0 left adjacent and 1's on top left diagonal
                elif pixMap[y+1,x-1] != 0 and pixMap[y,x-1] == 0: 
                  caseMap[y,x] = 11
              # First, locate all 1 pixels that are surrounded on bottom edge by 0 pixels
              # DDD
              # D1D
              # D0D
              if pixMap[y,x] == 1 and pixMap[y-1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                # First, find all pixels that have 0's adjacent and 1's on bottom 2 diagonals
                if pixMap[y-1,x+1] != 0 and pixMap[y-1, x-1] != 0 and pixMap[y,x+1] == 0 and pixMap[y,x-1] == 0:
                  caseMap[y,x] = 12
                # First, find all pixels that have 0 right adjacent and 1's on bottom right diagonal
                elif pixMap[y-1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 13
                # First, find all pixels that have 0 left adjacent and 1's on bottom left diagonal
                elif pixMap[y-1,x-1] != 0 and pixMap[y,x-1] == 0: 
                  caseMap[y,x] = 14
            if y == 0 and x > 0 and x < cols-1: #Do bottom row
              # First, locate all 1 pixels that are surrounded on top edge by 0 pixels
              # D0D
              # D1D
              # DDD
              if pixMap[y,x] == 1 and pixMap[y+1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                # First, find all pixels that have 0's adjacent and 1's on top 2 diagonals
                if pixMap[y+1,x+1] != 0 and pixMap[y+1, x-1] != 0 and pixMap[y,x+1] == 0 and pixMap[y,x-1] == 0:
                  caseMap[y,x] = 9
                # First, find all pixels that have 0 right adjacent and 1's on top right diagonal
                elif pixMap[y+1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 10
                # First, find all pixels that have 0 left adjacent and 1's on top left diagonal
                elif pixMap[y+1,x-1] != 0 and pixMap[y,x-1] == 0: 
                  caseMap[y,x] = 11
            if y == rows-1 and x > 0 and x < cols-1: #Do top row
              if pixMap[y,x] == 1 and pixMap[y-1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                # First, find all pixels that have 0's adjacent and 1's on bottom 2 diagonals
                if pixMap[y-1,x+1] != 0 and pixMap[y-1, x-1] != 0 and pixMap[y,x+1] == 0 and pixMap[y,x-1] == 0:
                  caseMap[y,x] = 12
                # First, find all pixels that have 0 right adjacent and 1's on bottom right diagonal
                elif pixMap[y-1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 13
                # First, find all pixels that have 0 left adjacent and 1's on bottom left diagonal
                elif pixMap[y-1,x-1] != 0 and pixMap[y,x-1] == 0: 
                  caseMap[y,x] = 14
            if y > 0 and x == 0 and y < rows-1: #Do first column
              if pixMap[y,x] == 1 and pixMap[y+1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                # First, find all pixels that have 0 right adjacent and 1's on top right diagonal
                if pixMap[y+1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 10
                # First, find all pixels that have 0 right adjacent and 1's on bottom right diagonal
              if pixMap[y,x] == 1 and pixMap[y-1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                if pixMap[y-1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 13
              if pixMap[y,x] == 1 and pixMap[y,x+1] == 0:
                if pixMap[y+1,x+1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y+1,x] == 0 and pixMap[y-1,x] == 0:
                  caseMap[y,x] = 15
        for x in range(0,cols):
          for y in range(0,rows):
            if caseMap[y,x] == 1:
              rect = gdspy.Rectangle((0, 0), (self.pixelSize, self.pixelSize), **l_top).translate(x*self.pixelSize,\
                     y*self.pixelSize)
              UNIT.add(rect)
            if caseMap[y,x] == 2:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2),self.pixelSize),
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 3:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2),self.pixelSize),
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 4:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),self.pixelSize),
                        (0, self.pixelSize), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 5:
              points = [(0,0),                                            #  _ 7-sided polygon 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0),                  # / \
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)),                   #|   |
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)),    #|___/
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),self.pixelSize),
                        ((self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 6:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize,0), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2),self.pixelSize),
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 7:
              points = [(0,0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2), self.pixelSize),
                        (0,self.pixelSize-(self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 8:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize,0), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        (0, self.pixelSize),
                        (0,(self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 9:
              points = [(0,0), 
                        (self.pixelSize,0), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 10:
              points = [(0,0), 
                        (self.pixelSize,0), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        (0, self.pixelSize)]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 11:
              points = [(0,0), 
                        (self.pixelSize,0), 
                        (self.pixelSize, self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 12:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize), 
                        (0, self.pixelSize), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 13:
              points = [(0,0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize), 
                        (0, self.pixelSize)]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 14:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize,0), 
                        (self.pixelSize, self.pixelSize), 
                        (0, self.pixelSize), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 15:
              points = [(0,0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        (0, self.pixelSize)]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 16:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize,0), 
                        (self.pixelSize, self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)

    elif self.corner == 'overlap':
      for x in range(0, cols):
        for y in range(0, rows):
          if pixMap[y,x]  == 1:
            rect = gdspy.Rectangle((0, 0), (self.pixelSize, self.pixelSize), **l_top).translate(x*self.pixelSize,\
                   y*self.pixelSize)
            UNIT.add(rect)

      for x in range(0,cols-1):
        for y in range(0,rows-1):
          if pixMap[y,x] == 1 and pixMap[y+1,x+1] == 1:
            diam = gdspy.Rectangle((-0.707*self.pixelSize/10, -0.707*self.pixelSize/10), \
                                   (0.707*self.pixelSize/10, 0.707*self.pixelSize/10), **l_top).rotate(np.pi/4).\
                                   translate((x+1)*self.pixelSize,(y+1)*self.pixelSize)
            UNIT.add(diam)
        
      for x in range(1,cols):
        for y in range(0,rows-1):
          if pixMap[y,x] == 1 and pixMap[y+1,x-1] == 1:
            diam = gdspy.Rectangle((-0.707*self.pixelSize/10, -0.707*self.pixelSize/10), \
                                   (0.707*self.pixelSize/10, 0.707*self.pixelSize/10), **l_top).rotate(np.pi/4).\
                                   translate((x)*self.pixelSize,(y+1)*self.pixelSize)
            UNIT.add(diam)
    else: # Mo Geometric Modifications
      for x in range(0, cols):
        for y in range(0, rows):
          if pixMap[y,x]  == 1:
            rect = gdspy.Rectangle((0, 0), (self.pixelSize, self.pixelSize), **l_top).translate(x*self.pixelSize,\
                   y*self.pixelSize)
            UNIT.add(rect)

    gdsFile = self.outF.replace('csv', 'gds') 
    gdspy.LayoutViewer(lib) 
    # Export GDS
    lib.write_gds(gdsFile)

  def print_image(self):

    lib = gdspy.GdsLibrary()
    
    # Set the database unit to 1 mil
    lib.unit = self.unit

    # Create Cell obj
    cellName = 'INVDESIGN'
    gdspy.current_library = gdspy.GdsLibrary() # This line of code has to be here to reset the GDS library on every loop
    UNIT = lib.new_cell(cellName)

    l_bottom = {"layer": 10, "datatype": 0}
    l_top = {"layer": 11, "datatype": 0}
    l_sources = {"layer": 5, "datatype": 0}

    pixMap = np.flipud(np.loadtxt(self.outF+'.csv', delimiter=','))
    print(pixMap)
    rows = np.size(pixMap, 0)
    cols = np.size(pixMap, 1)
    print(rows, cols)

    outline = gdspy.Rectangle((0, 0), (cols*self.pixelSize, rows*self.pixelSize), **l_bottom)
    UNIT.add(outline) 

    # When gds is fabricated it naturally creates overlap in the corners because of underetch of the metal. This routine 
    # creates non-overlapping polygons with enough space for manufacture. Space should be 6mil for PCB and typically 2um
    # in thick metal for ICs (A 0 in the pixel map corresponds to no metal and a 1 corresponds to metal
    if self.corner == 'noverlap': 
        caseMap = pixMap
        for x in range(0,cols):
          for y in range(0,rows):
            if y > 0 and x > 0 and y < rows-1 and x < cols-1: #Do inner portions of map first, edges will be done last
              # First, locate all 1 pixels that are surrounded on all edges by 0 pixels
              # D0D
              # 010
              # D0D
              if pixMap[y,x] == 1 and pixMap[y,x-1] == 0 and pixMap[y,x+1] == 0 and pixMap[y-1,x] == 0 and pixMap[y+1,x] == 0:
                # First, find all pixels that have 0's adjacent and 1's on all diagonals
                if pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] != 0:
                  caseMap[y,x] = 2
                # First, find all pixels that have 0's adjacent and 1's on 3/4 diagonals
                elif pixMap[y+1,x+1] == 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] != 0:
                  caseMap[y,x] = 3
                elif pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] == 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] != 0:
                  caseMap[y,x] = 4
                elif pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] == 0:
                  caseMap[y,x] = 5
                elif pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] == 0 and pixMap[y-1,x-1] != 0:
                  caseMap[y,x] = 6
                # First, find all pixels that have 0's adjacent and 1's on 2/4 diagonals
                elif pixMap[y+1,x+1] == 0 and pixMap[y+1,x-1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y-1,x-1] == 0:
                  caseMap[y,x] = 7
                elif pixMap[y+1,x+1] != 0 and pixMap[y+1,x-1] == 0 and pixMap[y-1,x+1] == 0 and pixMap[y-1,x-1] != 0:
                  caseMap[y,x] = 8
              # First, locate all 1 pixels that are surrounded on right edge by 0 pixels
              # DDD
              # D10
              # DDD
              if pixMap[y,x] == 1 and pixMap[y,x+1] == 0:
                if pixMap[y+1,x+1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y+1,x] == 0 and pixMap[y-1,x] == 0:
                  caseMap[y,x] = 15
              # First, locate all 1 pixels that are surrounded on left edge by 0 pixels
              # DDD
              # 01D
              # DDD
              if pixMap[y,x] == 1 and pixMap[y,x-1] == 0:
                if pixMap[y+1,x-1] != 0 and pixMap[y-1,x-1] != 0 and pixMap[y+1,x] == 0 and pixMap[y-1,x] == 0:
                  caseMap[y,x] = 16
              # First, locate all 1 pixels that are surrounded on top edge by 0 pixels
              # D0D
              # D1D
              # DDD
              if pixMap[y,x] == 1 and pixMap[y+1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                # First, find all pixels that have 0's adjacent and 1's on top 2 diagonals
                if pixMap[y+1,x+1] != 0 and pixMap[y+1, x-1] != 0 and pixMap[y,x+1] == 0 and pixMap[y,x-1] == 0:
                  caseMap[y,x] = 9
                # First, find all pixels that have 0 right adjacent and 1's on top right diagonal
                elif pixMap[y+1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 10
                # First, find all pixels that have 0 left adjacent and 1's on top left diagonal
                elif pixMap[y+1,x-1] != 0 and pixMap[y,x-1] == 0: 
                  caseMap[y,x] = 11
              # First, locate all 1 pixels that are surrounded on bottom edge by 0 pixels
              # DDD
              # D1D
              # D0D
              if pixMap[y,x] == 1 and pixMap[y-1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                # First, find all pixels that have 0's adjacent and 1's on bottom 2 diagonals
                if pixMap[y-1,x+1] != 0 and pixMap[y-1, x-1] != 0 and pixMap[y,x+1] == 0 and pixMap[y,x-1] == 0:
                  caseMap[y,x] = 12
                # First, find all pixels that have 0 right adjacent and 1's on bottom right diagonal
                elif pixMap[y-1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 13
                # First, find all pixels that have 0 left adjacent and 1's on bottom left diagonal
                elif pixMap[y-1,x-1] != 0 and pixMap[y,x-1] == 0: 
                  caseMap[y,x] = 14
            if y == 0 and x > 0 and x < cols-1: #Do bottom row
              # First, locate all 1 pixels that are surrounded on top edge by 0 pixels
              # D0D
              # D1D
              # DDD
              if pixMap[y,x] == 1 and pixMap[y+1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                # First, find all pixels that have 0's adjacent and 1's on top 2 diagonals
                if pixMap[y+1,x+1] != 0 and pixMap[y+1, x-1] != 0 and pixMap[y,x+1] == 0 and pixMap[y,x-1] == 0:
                  caseMap[y,x] = 9
                # First, find all pixels that have 0 right adjacent and 1's on top right diagonal
                elif pixMap[y+1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 10
                # First, find all pixels that have 0 left adjacent and 1's on top left diagonal
                elif pixMap[y+1,x-1] != 0 and pixMap[y,x-1] == 0: 
                  caseMap[y,x] = 11
            if y == rows-1 and x > 0 and x < cols-1: #Do top row
              if pixMap[y,x] == 1 and pixMap[y-1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                # First, find all pixels that have 0's adjacent and 1's on bottom 2 diagonals
                if pixMap[y-1,x+1] != 0 and pixMap[y-1, x-1] != 0 and pixMap[y,x+1] == 0 and pixMap[y,x-1] == 0:
                  caseMap[y,x] = 12
                # First, find all pixels that have 0 right adjacent and 1's on bottom right diagonal
                elif pixMap[y-1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 13
                # First, find all pixels that have 0 left adjacent and 1's on bottom left diagonal
                elif pixMap[y-1,x-1] != 0 and pixMap[y,x-1] == 0: 
                  caseMap[y,x] = 14
            if y > 0 and x == 0 and y < rows-1: #Do first column
              if pixMap[y,x] == 1 and pixMap[y+1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                # First, find all pixels that have 0 right adjacent and 1's on top right diagonal
                if pixMap[y+1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 10
                # First, find all pixels that have 0 right adjacent and 1's on bottom right diagonal
              if pixMap[y,x] == 1 and pixMap[y-1,x] == 0:# and (pixMap[y,x+1] == 0 or pixMap[y,x-1] == 0):
                if pixMap[y-1,x+1] != 0 and pixMap[y,x+1] == 0: 
                  caseMap[y,x] = 13
              if pixMap[y,x] == 1 and pixMap[y,x+1] == 0:
                if pixMap[y+1,x+1] != 0 and pixMap[y-1,x+1] != 0 and pixMap[y+1,x] == 0 and pixMap[y-1,x] == 0:
                  caseMap[y,x] = 15
        for x in range(0,cols):
          for y in range(0,rows):
            if caseMap[y,x] == 1:
              rect = gdspy.Rectangle((0, 0), (self.pixelSize, self.pixelSize), **l_top).translate(x*self.pixelSize,\
                     y*self.pixelSize)
              UNIT.add(rect)
            if caseMap[y,x] == 2:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2),self.pixelSize),
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 3:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2),self.pixelSize),
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 4:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),self.pixelSize),
                        (0, self.pixelSize), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 5:
              points = [(0,0),                                            #  _ 7-sided polygon 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0),                  # / \
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)),                   #|   |
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)),    #|___/
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),self.pixelSize),
                        ((self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 6:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize,0), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2),self.pixelSize),
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 7:
              points = [(0,0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2), self.pixelSize),
                        (0,self.pixelSize-(self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 8:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize,0), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        (0, self.pixelSize),
                        (0,(self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 9:
              points = [(0,0), 
                        (self.pixelSize,0), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 10:
              points = [(0,0), 
                        (self.pixelSize,0), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        (0, self.pixelSize)]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 11:
              points = [(0,0), 
                        (self.pixelSize,0), 
                        (self.pixelSize, self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 12:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize), 
                        (0, self.pixelSize), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 13:
              points = [(0,0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize), 
                        (0, self.pixelSize)]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 14:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize,0), 
                        (self.pixelSize, self.pixelSize), 
                        (0, self.pixelSize), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 15:
              points = [(0,0), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize, (self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (self.pixelSize-(self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        (0, self.pixelSize)]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 16:
              points = [((self.minPix/2)*np.sqrt(2),0), 
                        (self.pixelSize,0), 
                        (self.pixelSize, self.pixelSize), 
                        ((self.minPix/2)*np.sqrt(2), self.pixelSize), 
                        (0, self.pixelSize-(self.minPix/2)*np.sqrt(2)), 
                        (0, (self.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*self.pixelSize,y*self.pixelSize)
              UNIT.add(poly)

    elif self.corner == 'overlap':
      for x in range(0, cols):
        for y in range(0, rows):
          if pixMap[y,x]  == 1:
            rect = gdspy.Rectangle((0, 0), (self.pixelSize, self.pixelSize), **l_top).translate(x*self.pixelSize,\
                   y*self.pixelSize)
            UNIT.add(rect)

      for x in range(0,cols-1):
        for y in range(0,rows-1):
          if pixMap[y,x] == 1 and pixMap[y+1,x+1] == 1:
            diam = gdspy.Rectangle((-0.707*self.pixelSize/10, -0.707*self.pixelSize/10), \
                                   (0.707*self.pixelSize/10, 0.707*self.pixelSize/10), **l_top).rotate(np.pi/4).\
                                   translate((x+1)*self.pixelSize,(y+1)*self.pixelSize)
            UNIT.add(diam)
        
      for x in range(1,cols):
        for y in range(0,rows-1):
          if pixMap[y,x] == 1 and pixMap[y+1,x-1] == 1:
            diam = gdspy.Rectangle((-0.707*self.pixelSize/10, -0.707*self.pixelSize/10), \
                                   (0.707*self.pixelSize/10, 0.707*self.pixelSize/10), **l_top).rotate(np.pi/4).\
                                   translate((x)*self.pixelSize,(y+1)*self.pixelSize)
            UNIT.add(diam)
    else: # Mo Geometric Modifications
      for x in range(0, cols):
        for y in range(0, rows):
          if pixMap[y,x]  == 1:
            rect = gdspy.Rectangle((0, 0), (self.pixelSize, self.pixelSize), **l_top).translate(x*self.pixelSize,\
                   y*self.pixelSize)
            UNIT.add(rect)

    svgFile = self.outF.replace('csv', 'svg')  
    UNIT.write_svg(svgFile)
    pngFile = self.outF.replace('csv', 'png')  

    # read svg -> write png
    renderPM.drawToFile(svg2rlg(svgFile), pngFile, fmt='PNG')