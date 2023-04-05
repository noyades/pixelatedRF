# -*- coding: utf-8 -*-
import csv
import math
import random
import numpy as np
import gdspy
import pya
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPM
from vt_rrfc import * 
from vt_rrfc.randomDesigns import *
#from skimage import measure

class rrfc:
  def __init__(rrfc, 
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
    """    
    rrfc.unit = unit
    rrfc.ports = ports
    rrfc.corner = corner 
    rrfc.sides = sides
    rrfc.sym = sym
    rrfc.minPix = minPix
    rrfc.pixelSize = pixelSize
    rrfc.layoutRes = layoutRes
    rrfc.scale = scale
    rrfc.launchLen = launchLen
    rrfc.seed = seed
    rrfc.write = write
    rrfc.sim = sim
    rrfc.view = view
    rrfc.outF = outF
    rrfc.connect = connect
    rrfc.portPos = portPosition

def randomGDS_dim(
              sub,
              rrfc,
              x_dim: int,
              y_dim: int,
              imp: float):
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
  width_launch, length_launch = microstrip_calc.synthMicrostrip(sub, imp, rrfc.launchLen)
  width_launch = rrfc.layoutRes*width_launch # Convert from mils to layout resolution. layoutRes defaults to 1 for mils
  length_launch = rrfc.layoutRes*length_launch # Convert from mils to layout resolution. layoutRes defaults to 1 for mils
  if width_launch < rrfc.pixelSize: # Make sure that launch is at least 1 pixel wide
    width_launch = rrfc.pixelSize
  if length_launch < rrfc.pixelSize: # Make sure that launch is at least 1 pixel long
    length_launch = rrfc.pixelSize
  launch_l_pixels = round(length_launch/rrfc.pixelSize) # length of the line to connect to structure in number of pixels
  launch_w_pixels = round(width_launch/rrfc.pixelSize) # width of the line to connect to structure in number of pixels
  length_launch = launch_l_pixels*rrfc.pixelSize
  width_launch = launch_w_pixels*rrfc.pixelSize

  # Set horizontal and vertical pixel limits
  y_pixels = math.ceil(y_dim/rrfc.pixelSize) 
  x_pixels = math.ceil(x_dim/rrfc.pixelSize)
  if y_pixels % 2 == 0 and launch_l_pixels % 2 !=0:
    y_pixels += 1
  #if x_pixels % 2 == 0:
  #  x_pixels += 1

  rd1 = randomDesign(rrfc.unit, x_pixels, y_pixels, rrfc.ports,\
                    rrfc.sides, launch_l_pixels, launch_w_pixels, rrfc.pixelSize,rrfc.sim,rrfc.sym,rrfc.seed)
  x_total, y_total, portPos = randomDesign.genPortPos(rd1)
  pixelMap = np.zeros((int(x_total/rrfc.pixelSize),int(y_total/rrfc.pixelSize)),dtype=int)

  if rrfc.sim == 'EMX':
    ly = pya.Layout()

    # Set the database unit to 1 um. Generally for EMX it is easier to work in um.
    ly.dbu = rrfc.unit*1e6

    # Create Cell obj
    cellName = 'RANDOM'
    UNIT = ly.create_cell(cellName)

    # Create layer #'s
    l_top = ly.layer(59, 0) # layer for signal metal, corresponds to OI in FDX process
    l_bottom = ly.layer(15, 0) # Ground Layer corresponds to M1 in FDX process

    # Draw outline
    # pixelMap = np.zeros((int(x_total/rrfc.pixelSize),int(y_total/rrfc.pixelSize)),dtype=int)
    # Make a mesh shield
    # Assume ground hash is 0.5um strip followed by 1um space
    for x in range(0, int(y_total/rrfc.pixelSize)):
      rect = UNIT.shapes(l_bottom).insert( pya.Box(0, 0, x_total, 10).moved(0,(x+0.5)*rrfc.pixelSize-5))
    for x in range(0, int(x_total/rrfc.pixelSize)):
      rect = UNIT.shapes(l_bottom).insert( pya.Box(0, 0, 10, y_total).moved((x+0.5)*rrfc.pixelSize-5,0))

    #outline = UNIT.shapes(l_bottom).insert( pya.Box(0, 0, x_total, y_total) ) 

    pixelMap = randomDesign.addLaunch(rd1,x_total,y_total,pixelMap)

    pixelMap = randomDesign.randomPixRect(rd1, x_dim, y_dim, pixelMap)

    pixMap = (np.transpose(pixelMap))
    rows = np.size(pixMap, 0)
    cols = np.size(pixMap, 1)
    #print(rows, cols)

    if rrfc.corner == 'overlap':
      for x in range(0, cols):
        for y in range(0, rows):
          if pixMap[y,x]  == 1:
            rect = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize,\
                   y*rrfc.pixelSize)
            UNIT.shapes(l_top).insert(rect)
      
      points = [pya.DPoint(-rrfc.pixelSize/100, 0),
                pya.DPoint(0, rrfc.pixelSize/100),
                pya.DPoint(rrfc.pixelSize/100, 0),
                pya.DPoint(0,-rrfc.pixelSize/100)]
      for x in range(0,cols-1):
        for y in range(0,rows-1):
          if pixMap[y,x] == 1 and pixMap[y+1,x+1] == 1:
            diam = pya.DPolygon(points).moved((x+1)*rrfc.pixelSize/10,(y+1)*rrfc.pixelSize/10)
            UNIT.shapes(l_top).insert(diam)
      
      for x in range(1,cols):
        for y in range(0,rows-1):
          if pixMap[y,x] == 1 and pixMap[y+1,x-1] == 1:
            diam = pya.DPolygon(points).moved((x)*rrfc.pixelSize/10,(y+1)*rrfc.pixelSize/10)
            UNIT.shapes(l_top).insert(diam)

    if rrfc.ports == 1:
      port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', portPos[0]+rrfc.pixelSize/10, portPos[1]) )
      port_2 = UNIT.shapes(l_bottom).insert( pya.Text('p2', portPos[0]+rrfc.pixelSize/10, portPos[1]) ) # Ground Port
    elif rrfc.ports == 2:
      port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', portPos[0]+rrfc.pixelSize/10, portPos[1]) )
      port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', portPos[2]-rrfc.pixelSize/10, portPos[3]) )
      port_3 = UNIT.shapes(l_bottom).insert( pya.Text('p3', portPos[0]+rrfc.pixelSize/10, portPos[1]) ) # Ground Port
      port_4 = UNIT.shapes(l_bottom).insert( pya.Text('p4', portPos[2]-rrfc.pixelSize/10, portPos[3]) ) # Ground Port
    elif rrfc.ports == 3:
      port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', portPos[0]+rrfc.pixelSize/10, portPos[1]) )
      port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', portPos[2]+rrfc.pixelSize/10, portPos[3]) )
      port_3 = UNIT.shapes(l_top).insert( pya.Text('p3', portPos[4]-rrfc.pixelSize/10, portPos[5]) )
      port_4 = UNIT.shapes(l_bottom).insert( pya.Text('p4', portPos[0]+rrfc.pixelSize/10, portPos[1]) ) # Ground Port
      port_5 = UNIT.shapes(l_bottom).insert( pya.Text('p5', portPos[2]+rrfc.pixelSize/10, portPos[3]) ) # Ground Port
      port_6 = UNIT.shapes(l_bottom).insert( pya.Text('p6', portPos[4]-rrfc.pixelSize/10, portPos[5]) ) # Ground Port
    elif rrfc.ports == 4:
      port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', portPos[0]+rrfc.pixelSize/10, portPos[1]) )
      port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', portPos[2]+rrfc.pixelSize/10, portPos[3]) )
      port_3 = UNIT.shapes(l_top).insert( pya.Text('p3', portPos[4]-rrfc.pixelSize/10, portPos[5]) )
      port_4 = UNIT.shapes(l_top).insert( pya.Text('p4', portPos[6]-rrfc.pixelSize/10, portPos[7]) )
      port_5 = UNIT.shapes(l_bottom).insert( pya.Text('p5', portPos[0]+rrfc.pixelSize/10, portPos[1]) ) # Ground Port
      port_6 = UNIT.shapes(l_bottom).insert( pya.Text('p6', portPos[2]+rrfc.pixelSize/10, portPos[3]) ) # Ground Port
      port_7 = UNIT.shapes(l_bottom).insert( pya.Text('p7', portPos[4]-rrfc.pixelSize/10, portPos[5]) ) # Ground Port
      port_8 = UNIT.shapes(l_bottom).insert( pya.Text('p8', portPos[6]-rrfc.pixelSize/10, portPos[7]) ) # Ground Port

    c12, c13, c14, c23, c24, c34 = findConnectivity(pixelMap, rrfc.pixelSize,\
                                   rrfc.ports, rrfc.sides, portPos)
  
    connections = [c12, c13, c14, c23, c24, c34]

    if (rrfc.connect == 0 or rrfc.connect == connections) and rrfc.write == True:
      csvFile = rrfc.outF + ".csv"
      # Export Pixel Map file
      np.savetxt(csvFile, np.flipud(np.transpose(pixelMap)), fmt = '%d', delimiter = ",")
      gdsFile = rrfc.outF + ".gds"
      # Export GDS
      ly.write(gdsFile)
      # View GDS      
      if rrfc.view == True:
        # Load a GDSII file into a new library
        gdsii = gdspy.GdsLibrary(infile=gdsFile)
        gdspy.LayoutViewer(gdsii)
    else:
      csvFile = ''
      gdsFile = ''

  else:
    lib = gdspy.GdsLibrary()

    # Set the database unit to 1 mil
    lib.unit = rrfc.unit

    # Create Cell obj
    cellName = 'RANDOM'
    gdspy.current_library = gdspy.GdsLibrary() # This line of code has to be here to reset the GDS library on every loop
    UNIT = lib.new_cell(cellName)

    port_1, port_2, port_3, port_4, source = randomDesign.addPortsGdspy(rd1,x_total,y_total)
    # Create layer #'s
    l_bottom = {"layer": 10, "datatype": 0}
    l_top = {"layer": 11, "datatype": 0}
    l_sources = {"layer": 5, "datatype": 0}
    if rrfc.sim == 'MEEP':
      if rrfc.ports == 1:
        l_port1 = {"layer": 1, "datatype": 0}
        poly1 = gdspy.Polygon(port_1, **l_port1)
        UNIT.add(poly1)
      if rrfc.ports == 2:
        l_port1 = {"layer": 1, "datatype": 0}
        l_port2 = {"layer": 2, "datatype": 0}
        poly1 = gdspy.Polygon(port_1, **l_port1)
        poly2 = gdspy.Polygon(port_2, **l_port2)
        poly3 = gdspy.Polygon(source, **l_sources)
        UNIT.add(poly1)
        UNIT.add(poly2)
        UNIT.add(poly3)
      if rrfc.ports == 3:
        l_port1 = {"layer": 1, "datatype": 0}
        l_port2 = {"layer": 2, "datatype": 0}
        l_port3 = {"layer": 3, "datatype": 0}
        poly1 = gdspy.Polygon(port_1, **l_port1)
        poly2 = gdspy.Polygon(port_2, **l_port2)
        poly3 = gdspy.Polygon(port_3, **l_port3)
        poly4 = gdspy.Polygon(source, **l_sources)
        UNIT.add(poly1)
        UNIT.add(poly2)
        UNIT.add(poly3)
        UNIT.add(poly4)
      if rrfc.ports == 4:
        l_port1 = {"layer": 1, "datatype": 0}
        l_port2 = {"layer": 2, "datatype": 0}
        l_port3 = {"layer": 3, "datatype": 0}
        l_port4 = {"layer": 4, "datatype": 0}
        poly1 = gdspy.Polygon(port_1, **l_port1)
        poly2 = gdspy.Polygon(port_2, **l_port2)
        poly3 = gdspy.Polygon(port_3, **l_port3)
        poly4 = gdspy.Polygon(port_4, **l_port4)
        poly5 = gdspy.Polygon(source, **l_sources)
        UNIT.add(poly1)
        UNIT.add(poly2)
        UNIT.add(poly3)
        UNIT.add(poly4)
        UNIT.add(poly5)

    pixelMap = randomDesign.addLaunch(rd1,x_total,y_total,pixelMap)
    pixelMap = randomDesign.randomPixRect(rd1, x_dim, y_dim, pixelMap)

    pixMap = (np.transpose(pixelMap))
    rows = np.size(pixMap, 0)
    cols = np.size(pixMap, 1)
    #print(rows, cols)

    caseMap = pixMap

    # When gds is fabricated it naturally creates overlap in the corners because of underetch of the metal. This routine 
    # creates non-overlapping polygons with enough space for manufacture. Space should be 6mil for PCB and typically 2um
    # in thick metal for ICs (A 0 in the pixel map corresponds to no metal and a 1 corresponds to metal
    if rrfc.corner == 'noverlap': 
      for x in range(0,cols):
        for y in range(0,rows):
          if y > 0 and x > 0 and y < rows-1 and x < cols-1: #Do inner portions of map first, edges will be done last
            # First, locate all 1 pixels that are surrounded on all edges by 0 pixels
            # D0D
            # 010
            # D0D
            if caseMap[y,x] == 1 and caseMap[y,x-1] == 0 and caseMap[y,x+1] == 0 and caseMap[y-1,x] == 0 and caseMap[y+1,x] == 0:
              # First, find all pixels that have 0's adjacent and 1's on all diagonals
              if caseMap[y+1,x+1] != 0 and caseMap[y+1,x-1] != 0 and caseMap[y-1,x+1] != 0 and caseMap[y-1,x-1] != 0:
                caseMap[y,x] = 2
              # First, find all pixels that have 0's adjacent and 1's on 3/4 diagonals
              elif caseMap[y+1,x+1] == 0 and caseMap[y+1,x-1] != 0 and caseMap[y-1,x+1] != 0 and caseMap[y-1,x-1] != 0:
                caseMap[y,x] = 3
              elif caseMap[y+1,x+1] != 0 and caseMap[y+1,x-1] == 0 and caseMap[y-1,x+1] != 0 and caseMap[y-1,x-1] != 0:
                caseMap[y,x] = 4
              elif caseMap[y+1,x+1] != 0 and caseMap[y+1,x-1] != 0 and caseMap[y-1,x+1] != 0 and caseMap[y-1,x-1] == 0:
                caseMap[y,x] = 5
              elif caseMap[y+1,x+1] != 0 and caseMap[y+1,x-1] != 0 and caseMap[y-1,x+1] == 0 and caseMap[y-1,x-1] != 0:
                caseMap[y,x] = 6
              # First, find all pixels that have 0's adjacent and 1's on 2/4 diagonals
              elif caseMap[y+1,x+1] == 0 and caseMap[y+1,x-1] != 0 and caseMap[y-1,x+1] != 0 and caseMap[y-1,x-1] == 0:
                caseMap[y,x] = 7
              elif caseMap[y+1,x+1] != 0 and caseMap[y+1,x-1] == 0 and caseMap[y-1,x+1] == 0 and caseMap[y-1,x-1] != 0:
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
            rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize,\
                   y*rrfc.pixelSize)
            UNIT.add(rect)
          if caseMap[y,x] == 2:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 3:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 4:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                      (0, rrfc.pixelSize), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 5:
            points = [(0,0),                                            #  _ 7-sided polygon 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0),                  # / \
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)),                   #|   |
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)),    #|___/
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                      ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 6:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 7:
            points = [(0,0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize),
                      (0,rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 8:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize),
                      (0,(rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 9:
            points = [(0,0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 10:
            points = [(0,0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize)]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 11:
            points = [(0,0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 12:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      (0, rrfc.pixelSize), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 13:
            points = [(0,0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      (0, rrfc.pixelSize)]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 14:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      (0, rrfc.pixelSize), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 15:
            points = [(0,0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize)]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 16:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
      pixelMap[pixelMap > 1] = 1 #Don't know why this variable is changing in the loop above, but this is a temporary fix. 
                                 #only caseMap gets written in the code, but pixelMap changes with caseMap, so this restores
                                 #pixelMap to proper value

    elif rrfc.corner == 'overlap':
      for x in range(0, cols):
        for y in range(0, rows):
          if pixMap[y,x]  == 1:
            rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize,\
                   y*rrfc.pixelSize)
            UNIT.add(rect)

      for x in range(0,cols-1):
        for y in range(0,rows-1):
          if pixMap[y,x] == 1 and pixMap[y+1,x+1] == 1:
            diam = gdspy.Rectangle((-0.707*rrfc.pixelSize/10, -0.707*rrfc.pixelSize/10), \
                                   (0.707*rrfc.pixelSize/10, 0.707*rrfc.pixelSize/10), **l_top).rotate(np.pi/4).\
                                   translate((x+1)*rrfc.pixelSize,(y+1)*rrfc.pixelSize)
            UNIT.add(diam)
      
      for x in range(1,cols):
        for y in range(0,rows-1):
          if pixMap[y,x] == 1 and pixMap[y+1,x-1] == 1:
            diam = gdspy.Rectangle((-0.707*rrfc.pixelSize/10, -0.707*rrfc.pixelSize/10), \
                                   (0.707*rrfc.pixelSize/10, 0.707*rrfc.pixelSize/10), **l_top).rotate(np.pi/4).\
                                   translate((x)*rrfc.pixelSize,(y+1)*rrfc.pixelSize)
            UNIT.add(diam)
    else: # Mo Geometric Modifications
      for x in range(0, cols):
        for y in range(0, rows):
          if pixMap[y,x]  == 1:
            rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize,\
                   y*rrfc.pixelSize)
            UNIT.add(rect)

    # Draw outline
    #pixelMap = np.zeros((int(x_total/rrfc.pixelSize),int(y_total/rrfc.pixelSize)),dtype=int)
    outline = gdspy.Rectangle((0, 0), (x_total, y_total), **l_bottom)
    UNIT.add(outline) 
    c12, c13, c14, c23, c24, c34 = findConnectivity(pixelMap, rrfc.pixelSize,\
                                   rrfc.ports, rrfc.sides, portPos)
  
    connections = [c12, c13, c14, c23, c24, c34]

    if (rrfc.connect == 0 or rrfc.connect == connections) and rrfc.write == True:
      csvFile = rrfc.outF + ".csv"
      # Export Pixel Map file
      np.savetxt(csvFile, np.flipud(np.transpose(pixelMap)), fmt = '%d', delimiter = ",")
      gdsFile = rrfc.outF + ".gds"
      # Export GDS
      lib.write_gds(gdsFile)
    else:
      csvFile = ''
      gdsFile = ''
    # View GDS
    if rrfc.view == True:
      gdspy.LayoutViewer(lib) 

  return portPos, x_total, y_total, csvFile, gdsFile, cellName, launch_l_pixels

def recreateGDS_file(rrfc,sub,imp):

  pixMap = np.flipud(np.loadtxt(rrfc.outF, delimiter=','))
  rows = np.size(pixMap, 0)
  cols = np.size(pixMap, 1)
  if rrfc.sim == 'EMX':

    if rrfc.portPos == '':
      # In EMX sims, ports are defined on metal layers. They are added here, so port position needs to be computed
      width_launch, length_launch = microstrip_calc.synthMicrostrip(sub, imp, rrfc.launchLen)
      width_launch = rrfc.layoutRes*width_launch # Convert from mils to layout resolution. layoutRes defaults to 1 for mils
      length_launch = rrfc.layoutRes*length_launch # Convert from mils to layout resolution. layoutRes defaults to 1 for mils
      if width_launch < rrfc.pixelSize: # Make sure that launch is at least 1 pixel wide
        width_launch = rrfc.pixelSize
      if length_launch < rrfc.pixelSize: # Make sure that launch is at least 1 pixel long
        length_launch = rrfc.pixelSize
      launch_l_pixels = round(length_launch/rrfc.pixelSize) # length of the line to connect to structure in number of pixels
      launch_w_pixels = round(width_launch/rrfc.pixelSize) # width of the line to connect to structure in number of pixels
      rd1 = randomDesign(rrfc.unit, cols-2*launch_l_pixels, rows, rrfc.ports,\
                         rrfc.sides, launch_l_pixels, launch_w_pixels, \
                         rrfc.pixelSize,rrfc.sim,rrfc.sym,rrfc.seed)
      x_total, y_total, portPos = randomDesign.genPortPos(rd1)      
    else:
      portPos = rrfc.portPos
      x_total = cols*rrfc.pixelSize
      y_total = rows*rrfc.pixelSize
    
    ly = pya.Layout()
    # Set the database unit to 1 um. Generally for EMX it is easier to work in um.
    ly.dbu = rrfc.unit*1e6

    # Create Cell obj
    cellName = 'INVDESIGN'
    UNIT = ly.create_cell(cellName)

    # Create layer #'s
    l_top = ly.layer(59, 0) # layer for signal metal, corresponds to OI in FDX process
    l_bottom = ly.layer(15, 0) # Ground Layer corresponds to M1 in FDX process

    # Draw outline
    # Make a mesh shield
    # Assume ground hash is 0.5um strip followed by 1um space, place it on the 
    # bottom metal layer
    for x in range(0, int(y_total/rrfc.pixelSize)):
      rect = UNIT.shapes(l_bottom).insert( pya.Box(0, 0, x_total, 10).moved(0,(x+0.5)*rrfc.pixelSize-5))
    for x in range(0, int(x_total/rrfc.pixelSize)):
      rect = UNIT.shapes(l_bottom).insert( pya.Box(0, 0, 10, y_total).moved((x+0.5)*rrfc.pixelSize-5,0))

    if rrfc.corner == 'overlap':
      for x in range(0, cols):
        for y in range(0, rows):
          if pixMap[y,x]  == 1:
            rect = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize,\
                   y*rrfc.pixelSize)
            UNIT.shapes(l_top).insert(rect)
      
      points = [pya.DPoint(-rrfc.pixelSize/100, 0),
                pya.DPoint(0, rrfc.pixelSize/100),
                pya.DPoint(rrfc.pixelSize/100, 0),
                pya.DPoint(0,-rrfc.pixelSize/100)]
      for x in range(0,cols-1):
        for y in range(0,rows-1):
          if pixMap[y,x] == 1 and pixMap[y+1,x+1] == 1:
            diam = pya.DPolygon(points).moved((x+1)*rrfc.pixelSize/10,(y+1)*rrfc.pixelSize/10)
            UNIT.shapes(l_top).insert(diam)
      
      for x in range(1,cols):
        for y in range(0,rows-1):
          if pixMap[y,x] == 1 and pixMap[y+1,x-1] == 1:
            diam = pya.DPolygon(points).moved((x)*rrfc.pixelSize/10,(y+1)*rrfc.pixelSize/10)
            UNIT.shapes(l_top).insert(diam)

    if rrfc.ports == 1:
      if rrfc.portPos == '':
        port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', portPos[0]+rrfc.pixelSize/10, portPos[1]) )
        port_2 = UNIT.shapes(l_bottom).insert( pya.Text('p2', portPos[0]+rrfc.pixelSize/10, portPos[1]) ) # Ground Port
      else:
        port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', rrfc.scale*portPos[0], rrfc.scale*portPos[1]+rrfc.pixelSize/(2)) )
        port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', rrfc.scale*portPos[0], rrfc.scale*portPos[1]-rrfc.pixelSize/(2)) )
    elif rrfc.ports == 2:
      if rrfc.portPos == '':
        port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', portPos[0]+rrfc.pixelSize/10, portPos[1]) )
        port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', portPos[2]-rrfc.pixelSize/10, portPos[3]) )
        port_3 = UNIT.shapes(l_bottom).insert( pya.Text('p3', portPos[0]+rrfc.pixelSize/10, portPos[1]) ) # Ground Port
        port_4 = UNIT.shapes(l_bottom).insert( pya.Text('p4', portPos[2]-rrfc.pixelSize/10, portPos[3]) ) # Ground Port
      else:
        port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', rrfc.scale*portPos[0], rrfc.scale*portPos[1]+rrfc.pixelSize/(2)) )
        port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', rrfc.scale*portPos[2], rrfc.scale*portPos[3]-rrfc.pixelSize/(2)) )
        port_3 = UNIT.shapes(l_top).insert( pya.Text('p3', rrfc.scale*portPos[0], rrfc.scale*portPos[1]+rrfc.pixelSize/(2)) )
        port_4 = UNIT.shapes(l_top).insert( pya.Text('p4', rrfc.scale*portPos[2], rrfc.scale*portPos[3]-rrfc.pixelSize/(2)) )
    elif rrfc.ports == 3:
      if rrfc.portPos == '':
        port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', portPos[0]+rrfc.pixelSize/10, portPos[1]) )
        port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', portPos[2]+rrfc.pixelSize/10, portPos[3]) )
        port_3 = UNIT.shapes(l_top).insert( pya.Text('p3', portPos[4]-rrfc.pixelSize/10, portPos[5]) )
        port_4 = UNIT.shapes(l_bottom).insert( pya.Text('p4', portPos[0]+rrfc.pixelSize/10, portPos[1]) ) # Ground Port
        port_5 = UNIT.shapes(l_bottom).insert( pya.Text('p5', portPos[2]+rrfc.pixelSize/10, portPos[3]) ) # Ground Port
        port_6 = UNIT.shapes(l_bottom).insert( pya.Text('p6', portPos[4]-rrfc.pixelSize/10, portPos[5]) ) # Ground Port
      else:
        port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', rrfc.scale*portPos[0], rrfc.scale*portPos[1]+rrfc.pixelSize/(2)) )
        port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', rrfc.scale*portPos[2], rrfc.scale*portPos[3]-rrfc.pixelSize/(2)) )
        port_3 = UNIT.shapes(l_top).insert( pya.Text('p3', rrfc.scale*portPos[4], rrfc.scale*portPos[5]+rrfc.pixelSize/(2)) )
        port_4 = UNIT.shapes(l_top).insert( pya.Text('p4', rrfc.scale*portPos[0], rrfc.scale*portPos[1]-rrfc.pixelSize/(2)) )
        port_5 = UNIT.shapes(l_bottom).insert( pya.Text('p5', rrfc.scale*portPos[2], rrfc.scale*portPos[3]+rrfc.pixelSize/(2)) ) # Ground Port
        port_6 = UNIT.shapes(l_bottom).insert( pya.Text('p6', rrfc.scale*portPos[4], rrfc.scale*portPos[5]-rrfc.pixelSize/(2)) ) # Ground Port
    elif rrfc.ports == 4:
      if rrfc.portPos == '':
        port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', portPos[0]+rrfc.pixelSize/10, portPos[1]) )
        port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', portPos[2]+rrfc.pixelSize/10, portPos[3]) )
        port_3 = UNIT.shapes(l_top).insert( pya.Text('p3', portPos[4]-rrfc.pixelSize/10, portPos[5]) )
        port_4 = UNIT.shapes(l_top).insert( pya.Text('p4', portPos[6]-rrfc.pixelSize/10, portPos[7]) )
        port_5 = UNIT.shapes(l_bottom).insert( pya.Text('p5', portPos[0]+rrfc.pixelSize/10, portPos[1]) ) # Ground Port
        port_6 = UNIT.shapes(l_bottom).insert( pya.Text('p6', portPos[2]+rrfc.pixelSize/10, portPos[3]) ) # Ground Port
        port_7 = UNIT.shapes(l_bottom).insert( pya.Text('p7', portPos[4]-rrfc.pixelSize/10, portPos[5]) ) # Ground Port
        port_8 = UNIT.shapes(l_bottom).insert( pya.Text('p8', portPos[6]-rrfc.pixelSize/10, portPos[7]) ) # Ground Port
      else:
        port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', rrfc.scale*portPos[0], rrfc.scale*portPos[1]+rrfc.pixelSize/(2)) )
        port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', rrfc.scale*portPos[2], rrfc.scale*portPos[3]-rrfc.pixelSize/(2)) )
        port_3 = UNIT.shapes(l_top).insert( pya.Text('p3', rrfc.scale*portPos[4], rrfc.scale*portPos[5]+rrfc.pixelSize/(2)) )
        port_4 = UNIT.shapes(l_top).insert( pya.Text('p4', rrfc.scale*portPos[6], rrfc.scale*portPos[7]-rrfc.pixelSize/(2)) )
        port_5 = UNIT.shapes(l_bottom).insert( pya.Text('p5', rrfc.scale*portPos[0], rrfc.scale*portPos[1]+rrfc.pixelSize/(2)) ) # Ground Port
        port_6 = UNIT.shapes(l_bottom).insert( pya.Text('p6', rrfc.scale*portPos[2], rrfc.scale*portPos[3]-rrfc.pixelSize/(2)) ) # Ground Port
        port_7 = UNIT.shapes(l_bottom).insert( pya.Text('p7', rrfc.scale*portPos[4], rrfc.scale*portPos[5]+rrfc.pixelSize/(2)) ) # Ground Port
        port_8 = UNIT.shapes(l_bottom).insert( pya.Text('p8', rrfc.scale*portPos[6], rrfc.scale*portPos[7]-rrfc.pixelSize/(2)) ) # Ground Port

    gdsFile = rrfc.outF.replace('csv', 'gds')
    # Export GDS
    ly.write(gdsFile)
    # View GDS      
    
  else:
    lib = gdspy.GdsLibrary()

    # Set the database unit to 1 mil
    lib.unit = rrfc.unit

    # Create Cell obj
    cellName = 'INVDESIGN'
    gdspy.current_library = gdspy.GdsLibrary() # This line of code has to be here to reset the GDS library on every loop
    UNIT = lib.new_cell(cellName)

    l_bottom = {"layer": 10, "datatype": 0}
    l_top = {"layer": 11, "datatype": 0}
    l_sources = {"layer": 5, "datatype": 0}

    outline = gdspy.Rectangle((0, 0), (cols*rrfc.pixelSize, rows*rrfc.pixelSize), **l_bottom)
    UNIT.add(outline) 

    # When gds is fabricated it naturally creates overlap in the corners because of underetch of the metal. This routine 
    # creates non-overlapping polygons with enough space for manufacture. Space should be 6mil for PCB and typically 2um
    # in thick metal for ICs (A 0 in the pixel map corresponds to no metal and a 1 corresponds to metal
    if rrfc.corner == 'noverlap': 
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
              rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize,\
                     y*rrfc.pixelSize)
              UNIT.add(rect)
            if caseMap[y,x] == 2:
              points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                        (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                        (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                        (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                        (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                        ((rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                        (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                        (0, (rrfc.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 3:
              points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                        (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                        (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                        (rrfc.pixelSize, rrfc.pixelSize), 
                        ((rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                        (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                        (0, (rrfc.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 4:
              points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                        (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                        (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                        (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                        (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                        (0, rrfc.pixelSize), 
                        (0, (rrfc.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 5:
              points = [(0,0),                                            #  _ 7-sided polygon 
                        (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0),                  # / \
                        (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)),                   #|   |
                        (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)),    #|___/
                        (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                        ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                        (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 6:
              points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                        (rrfc.pixelSize,0), 
                        (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                        (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                        ((rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                        (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                        (0, (rrfc.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 7:
              points = [(0,0), 
                        (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                        (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                        (rrfc.pixelSize, rrfc.pixelSize), 
                        ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize),
                        (0,rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 8:
              points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                        (rrfc.pixelSize,0), 
                        (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                        (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                        (0, rrfc.pixelSize),
                        (0,(rrfc.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 9:
              points = [(0,0), 
                        (rrfc.pixelSize,0), 
                        (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                        (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                        ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                        (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 10:
              points = [(0,0), 
                        (rrfc.pixelSize,0), 
                        (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                        (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                        (0, rrfc.pixelSize)]
              poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 11:
              points = [(0,0), 
                        (rrfc.pixelSize,0), 
                        (rrfc.pixelSize, rrfc.pixelSize), 
                        ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                        (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 12:
              points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                        (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                        (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                        (rrfc.pixelSize, rrfc.pixelSize), 
                        (0, rrfc.pixelSize), 
                        (0, (rrfc.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 13:
              points = [(0,0), 
                        (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                        (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                        (rrfc.pixelSize, rrfc.pixelSize), 
                        (0, rrfc.pixelSize)]
              poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 14:
              points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                        (rrfc.pixelSize,0), 
                        (rrfc.pixelSize, rrfc.pixelSize), 
                        (0, rrfc.pixelSize), 
                        (0, (rrfc.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 15:
              points = [(0,0), 
                        (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                        (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                        (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                        (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                        (0, rrfc.pixelSize)]
              poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
              UNIT.add(poly)
            if caseMap[y,x] == 16:
              points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                        (rrfc.pixelSize,0), 
                        (rrfc.pixelSize, rrfc.pixelSize), 
                        ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                        (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                        (0, (rrfc.minPix/2)*np.sqrt(2))]
              poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
              UNIT.add(poly)
      
    elif rrfc.corner == 'overlap':
      for x in range(0, cols):
        for y in range(0, rows):
          if pixMap[y,x]  == 1:
            rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize,\
                   y*rrfc.pixelSize)
            UNIT.add(rect)

      for x in range(0,cols-1):
        for y in range(0,rows-1):
          if pixMap[y,x] == 1 and pixMap[y+1,x+1] == 1:
            diam = gdspy.Rectangle((-0.707*rrfc.pixelSize/10, -0.707*rrfc.pixelSize/10), \
                                   (0.707*rrfc.pixelSize/10, 0.707*rrfc.pixelSize/10), **l_top).rotate(np.pi/4).\
                                   translate((x+1)*rrfc.pixelSize,(y+1)*rrfc.pixelSize)
            UNIT.add(diam)
      
      for x in range(1,cols):
        for y in range(0,rows-1):
          if pixMap[y,x] == 1 and pixMap[y+1,x-1] == 1:
            diam = gdspy.Rectangle((-0.707*rrfc.pixelSize/10, -0.707*rrfc.pixelSize/10), \
                                   (0.707*rrfc.pixelSize/10, 0.707*rrfc.pixelSize/10), **l_top).rotate(np.pi/4).\
                                   translate((x)*rrfc.pixelSize,(y+1)*rrfc.pixelSize)
            UNIT.add(diam)
    else: # Mo Geometric Modifications
      for x in range(0, cols):
        for y in range(0, rows):
          if pixMap[y,x]  == 1:
            rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize,\
                   y*rrfc.pixelSize)
            UNIT.add(rect)

    gdsFile = rrfc.outF.replace('csv', 'gds') 
    ##print(gdsFile) #debug 
    #if rrfc.view == True:
    #  gdspy.LayoutViewer(lib) 
    ## Export GDS
    lib.write_gds(gdsFile)
    
  if rrfc.view == True:
    # Load a GDSII file into a new library
    gdsii = gdspy.GdsLibrary(infile=gdsFile)
    gdspy.LayoutViewer(gdsii)    
  return cellName

def recreateGDS_array(rrfc,pixMap):

  lib = gdspy.GdsLibrary()

  # Set the database unit to 1 mil
  lib.unit = rrfc.unit

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

  outline = gdspy.Rectangle((0, 0), (cols*rrfc.pixelSize, rows*rrfc.pixelSize), **l_bottom)
  UNIT.add(outline) 

  # When gds is fabricated it naturally creates overlap in the corners because of underetch of the metal. This routine 
  # creates non-overlapping polygons with enough space for manufacture. Space should be 6mil for PCB and typically 2um
  # in thick metal for ICs (A 0 in the pixel map corresponds to no metal and a 1 corresponds to metal
  if rrfc.corner == 'noverlap': 
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
            rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize,\
                   y*rrfc.pixelSize)
            UNIT.add(rect)
          if caseMap[y,x] == 2:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 3:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 4:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                      (0, rrfc.pixelSize), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 5:
            points = [(0,0),                                            #  _ 7-sided polygon 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0),                  # / \
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)),                   #|   |
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)),    #|___/
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                      ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 6:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 7:
            points = [(0,0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize),
                      (0,rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 8:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize),
                      (0,(rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 9:
            points = [(0,0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 10:
            points = [(0,0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize)]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 11:
            points = [(0,0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 12:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      (0, rrfc.pixelSize), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 13:
            points = [(0,0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      (0, rrfc.pixelSize)]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 14:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      (0, rrfc.pixelSize), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 15:
            points = [(0,0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize)]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 16:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)

  elif rrfc.corner == 'overlap':
    for x in range(0, cols):
      for y in range(0, rows):
        if pixMap[y,x]  == 1:
          rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize,\
                 y*rrfc.pixelSize)
          UNIT.add(rect)

    for x in range(0,cols-1):
      for y in range(0,rows-1):
        if pixMap[y,x] == 1 and pixMap[y+1,x+1] == 1:
          diam = gdspy.Rectangle((-0.707*rrfc.pixelSize/10, -0.707*rrfc.pixelSize/10), \
                                 (0.707*rrfc.pixelSize/10, 0.707*rrfc.pixelSize/10), **l_top).rotate(np.pi/4).\
                                 translate((x+1)*rrfc.pixelSize,(y+1)*rrfc.pixelSize)
          UNIT.add(diam)
      
    for x in range(1,cols):
      for y in range(0,rows-1):
        if pixMap[y,x] == 1 and pixMap[y+1,x-1] == 1:
          diam = gdspy.Rectangle((-0.707*rrfc.pixelSize/10, -0.707*rrfc.pixelSize/10), \
                                 (0.707*rrfc.pixelSize/10, 0.707*rrfc.pixelSize/10), **l_top).rotate(np.pi/4).\
                                 translate((x)*rrfc.pixelSize,(y+1)*rrfc.pixelSize)
          UNIT.add(diam)
  else: # Mo Geometric Modifications
    for x in range(0, cols):
      for y in range(0, rows):
        if pixMap[y,x]  == 1:
          rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize,\
                 y*rrfc.pixelSize)
          UNIT.add(rect)

  gdsFile = rrfc.outF.replace('csv', 'gds') 
  print(gdsFile)
  gdspy.LayoutViewer(lib) 
  # Export GDS
  lib.write_gds(gdsFile)

def printImage(rrfc):

  lib = gdspy.GdsLibrary()

  # Set the database unit to 1 mil
  lib.unit = rrfc.unit

  # Create Cell obj
  cellName = 'INVDESIGN'
  gdspy.current_library = gdspy.GdsLibrary() # This line of code has to be here to reset the GDS library on every loop
  UNIT = lib.new_cell(cellName)

  l_bottom = {"layer": 10, "datatype": 0}
  l_top = {"layer": 11, "datatype": 0}
  l_sources = {"layer": 5, "datatype": 0}

  pixMap = np.flipud(np.loadtxt(rrfc.outF, delimiter=','))
  print(pixMap)
  rows = np.size(pixMap, 0)
  cols = np.size(pixMap, 1)
  print(rows, cols)

  outline = gdspy.Rectangle((0, 0), (cols*rrfc.pixelSize, rows*rrfc.pixelSize), **l_bottom)
  UNIT.add(outline) 

  # When gds is fabricated it naturally creates overlap in the corners because of underetch of the metal. This routine 
  # creates non-overlapping polygons with enough space for manufacture. Space should be 6mil for PCB and typically 2um
  # in thick metal for ICs (A 0 in the pixel map corresponds to no metal and a 1 corresponds to metal
  if rrfc.corner == 'noverlap': 
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
            rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize,\
                   y*rrfc.pixelSize)
            UNIT.add(rect)
          if caseMap[y,x] == 2:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 3:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 4:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                      (0, rrfc.pixelSize), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 5:
            points = [(0,0),                                            #  _ 7-sided polygon 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0),                  # / \
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)),                   #|   |
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)),    #|___/
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                      ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 6:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2),rrfc.pixelSize),
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 7:
            points = [(0,0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize),
                      (0,rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 8:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize),
                      (0,(rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 9:
            points = [(0,0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 10:
            points = [(0,0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize)]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 11:
            points = [(0,0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 12:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      (0, rrfc.pixelSize), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 13:
            points = [(0,0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      (0, rrfc.pixelSize)]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 14:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      (0, rrfc.pixelSize), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 15:
            points = [(0,0), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize, (rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize)]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)
          if caseMap[y,x] == 16:
            points = [((rrfc.minPix/2)*np.sqrt(2),0), 
                      (rrfc.pixelSize,0), 
                      (rrfc.pixelSize, rrfc.pixelSize), 
                      ((rrfc.minPix/2)*np.sqrt(2), rrfc.pixelSize), 
                      (0, rrfc.pixelSize-(rrfc.minPix/2)*np.sqrt(2)), 
                      (0, (rrfc.minPix/2)*np.sqrt(2))]
            poly = gdspy.Polygon(points, **l_top).translate(x*rrfc.pixelSize,y*rrfc.pixelSize)
            UNIT.add(poly)

  elif rrfc.corner == 'overlap':
    for x in range(0, cols):
      for y in range(0, rows):
        if pixMap[y,x]  == 1:
          rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize,\
                 y*rrfc.pixelSize)
          UNIT.add(rect)

    for x in range(0,cols-1):
      for y in range(0,rows-1):
        if pixMap[y,x] == 1 and pixMap[y+1,x+1] == 1:
          diam = gdspy.Rectangle((-0.707*rrfc.pixelSize/10, -0.707*rrfc.pixelSize/10), \
                                 (0.707*rrfc.pixelSize/10, 0.707*rrfc.pixelSize/10), **l_top).rotate(np.pi/4).\
                                 translate((x+1)*rrfc.pixelSize,(y+1)*rrfc.pixelSize)
          UNIT.add(diam)
      
    for x in range(1,cols):
      for y in range(0,rows-1):
        if pixMap[y,x] == 1 and pixMap[y+1,x-1] == 1:
          diam = gdspy.Rectangle((-0.707*rrfc.pixelSize/10, -0.707*rrfc.pixelSize/10), \
                                 (0.707*rrfc.pixelSize/10, 0.707*rrfc.pixelSize/10), **l_top).rotate(np.pi/4).\
                                 translate((x)*rrfc.pixelSize,(y+1)*rrfc.pixelSize)
          UNIT.add(diam)
  else: # Mo Geometric Modifications
    for x in range(0, cols):
      for y in range(0, rows):
        if pixMap[y,x]  == 1:
          rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize,\
                 y*rrfc.pixelSize)
          UNIT.add(rect)

  svgFile = rrfc.outF.replace('csv', 'svg')  
  UNIT.write_svg(svgFile)
  pngFile = rrfc.outF.replace('csv', 'png')  

  # read svg -> write png
  renderPM.drawToFile(svg2rlg(svgFile), pngFile, fmt='PNG')
  '''
  gdsFile = rrfc.outF.replace('csv', 'gds') 
  #print(gdsFile) #debug 
  if rrfc.view == True:
    gdspy.LayoutViewer(lib) 
  # Export GDS
  lib.write_gds(gdsFile)
  '''
