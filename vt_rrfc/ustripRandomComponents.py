# -*- coding: utf-8 -*-
import csv
import math
import random
import numpy as np
import gdspy
import pya
from vt_rrfc import * 
from vt_rrfc.randomDesigns import *
#from skimage import measure

class rrfc:
  def __init__(rrfc_param, 
               unit,
               ports: int, 
               sides: int,
               corner: str, 
               connect,
               pixelSize: int,
               seed: int, 
               sim: str, 
               view: bool, 
               write: bool, 
               outF: str,
               sym: str):
    """ define the microstrip substrate
    Args:
    unit = grid unit of the layout. Set to 1e-6 for um or 25.4e-6 for mil for example
    ports = number of ports
    sides = the number of sides ports are placed on
    corner = how the corner connection is realize for simulation 
    sym = whether to force symmetry about 'x-axis' or 'y-axis'. Leave empty for no symmetry 
    pixelSize = size of the randomized pixel
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
    rrfc_param.unit = unit
    rrfc_param.ports = ports
    rrfc_param.corner = corner 
    rrfc_param.sides = sides
    rrfc_param.sym = sym
    rrfc_param.pixelSize = pixelSize
    rrfc_param.seed = seed
    rrfc_param.write = write
    rrfc_param.sim = sim
    rrfc_param.view = view
    rrfc_param.outF = outF
    rrfc_param.connect = connect

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
  width_launch, length_launch = microstrip_calc.synthMicrostrip(sub, imp, 30)
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
  pixel_map = np.zeros((int(x_total/rrfc.pixelSize),int(y_total/rrfc.pixelSize)),dtype=int)

  if rrfc.sim == 'EMX':
    ly = pya.Layout()

    # Set the database unit to 1 mil
    ly.dbu = rrfc.unit*1e6

    # Create Cell obj
    cellName = 'RANDOM'
    UNIT = ly.create_cell(cellName)

    # Create layer #'s
    l_top = ly.layer(69, 0) # layer for metal
    l_bottom = ly.layer(6969, 0) # Ground Layer

    # Draw outline
    # pixel_map = np.zeros((int(x_total/rrfc.pixelSize),int(y_total/rrfc.pixelSize)),dtype=int)
    outline = UNIT.shapes(l_bottom).insert( pya.Box(0, 0, x_total, y_total) ) 

    pixel_map = randomDesign.addLaunch(rd1,x_total,y_total,pixel_map)

    pixel_map = randomDesign.randomPixRect(rd1, x_dim, y_dim, pixel_map)

    pixMap = (np.transpose(pixel_map))
    rows = np.size(pixMap, 0)
    cols = np.size(pixMap, 1)
    #print(rows, cols)

    print(pixMap)
    for x in range(0, cols):
      for y in range(0, rows):
        if pixMap[y,x]  == 1:
          rect = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize,\
                 y*rrfc.pixelSize)
          UNIT.shapes(l_top).insert(rect)

    if rrfc.corner == 'overlap':
      points = [pya.DPoint(0, rrfc.pixelSize/10),
                pya.DPoint(rrfc.pixelSize/10, 0),
                pya.DPoint(0, -rrfc.pixelSize/10),
                pya.DPoint(-rrfc.pixelSize/10, 0)]
      for x in range(0,cols-1):
        for y in range(0,rows-1):
          if pixMap[y,x] == 1 and pixMap[y+1,x+1] == 1:
            diam = pya.DPolygon(points).moved((x+1)*rrfc.pixelSize,(y+1)*rrfc.pixelSize)
            UNIT.shapes(l_top).insert(diam)

      for x in range(1,cols):
        for y in range(0,rows-1):
          if pixMap[y,x] == 1 and pixMap[y+1,x-1] == 1:
            diam = pya.DPolygon(points).moved((x)*rrfc.pixelSize,(y+1)*rrfc.pixelSize)
            UNIT.shapes(l_top).insert(diam)

    if rrfc.ports == 2:
      port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', portPos[0], portPos[1]) )
      port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', portPos[2], portPos[3]) )
    elif rrfc.ports == 3:
      port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', portPos[0], portPos[1]) )
      port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', portPos[2], portPos[3]) )
      port_3 = UNIT.shapes(l_top).insert( pya.Text('p3', portPos[4], portPos[5]) )
    elif rrfc.ports == 4:
      port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', portPos[0], portPos[1]) )
      port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', portPos[2], portPos[3]) )
      port_3 = UNIT.shapes(l_top).insert( pya.Text('p3', portPos[4], portPos[5]) )
      port_4 = UNIT.shapes(l_top).insert( pya.Text('p4', portPos[6], portPos[7]) )

    c12, c13, c14, c23, c24, c34 = findConnectivity(pixel_map, rrfc.pixelSize,\
                                   rrfc.ports, rrfc.sides, portPos)
  
    connections = [c12, c13, c14, c23, c24, c34]

    if rrfc.connect == connections and rrfc.write == True:
      csvFile = rrfc.outF + ".csv"
      # Export Pixel Map file
      np.savetxt(csvFile, np.flipud(np.transpose(pixel_map)), fmt = '%d', delimiter = ",")
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

    pixel_map = randomDesign.addLaunch(rd1,x_total,y_total,pixel_map)

    pixel_map = randomDesign.randomPixRect(rd1, x_dim, y_dim, pixel_map)

    pixMap = (np.transpose(pixel_map))
    rows = np.size(pixMap, 0)
    cols = np.size(pixMap, 1)
    #print(rows, cols)

    for x in range(0, cols):
      for y in range(0, rows):
        if pixMap[y,x]  == 1:
          rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize,\
                 y*rrfc.pixelSize)
          UNIT.add(rect)

    if rrfc.corner == 'overlap':
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
      

    #joined = gdspy.boolean(UNIT,None,'or')
    #joined.fillet(0.25)
    #UNIT.add(joined)

    # Draw outline
    #pixel_map = np.zeros((int(x_total/rrfc.pixelSize),int(y_total/rrfc.pixelSize)),dtype=int)
    outline = gdspy.Rectangle((0, 0), (x_total, y_total), **l_bottom)
    UNIT.add(outline) 

    c12, c13, c14, c23, c24, c34 = findConnectivity(pixel_map, rrfc.pixelSize,\
                                   rrfc.ports, rrfc.sides, portPos)
  
    connections = [c12, c13, c14, c23, c24, c34]

    if rrfc.connect == connections and rrfc.write == True:
      csvFile = rrfc.outF + ".csv"
      # Export Pixel Map file
      np.savetxt(csvFile, np.flipud(np.transpose(pixel_map)), fmt = '%d', delimiter = ",")
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

def recreateGDS_file(rrfc):

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

  #print(rrfc.outF) #debug
  pixMap = np.flipud(np.loadtxt(rrfc.outF, delimiter=','))
  rows = np.size(pixMap, 0)
  cols = np.size(pixMap, 1)
  #print(rows, cols) #debug

  outline = gdspy.Rectangle((0, 0), (cols*rrfc.pixelSize, rows*rrfc.pixelSize), **l_bottom)
  UNIT.add(outline) 

  for x in range(0, cols):
    for y in range(0, rows):
      if pixMap[y,x]  == 1:
        rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize,\
               y*rrfc.pixelSize)
        UNIT.add(rect)

    # These sections are added to mimic overlap that is left in fabrication at corner joints between pixels due 
    # to underetch or finite rotational tool dimensions on a CNC milling machine.
    if rrfc.corner == 'overlap':
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

  gdsFile = rrfc.outF.replace('csv', 'gds') 
  #print(gdsFile) #debug
  if rrfc.view == True:
    gdspy.LayoutViewer(lib) 
  # Export GDS
  lib.write_gds(gdsFile)
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

  for x in range(0, cols):
    for y in range(0, rows):
      if pixMap[y,x]  == 1:
        rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize,\
               y*rrfc.pixelSize)
        UNIT.add(rect)

    # These sections are added to mimic overlap that is left in fabrication at corner joints between pixels due 
    # to underetch or finite rotational tool dimensions on a CNC milling machine.
    if rrfc.corner == 'overlap':
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

  gdsFile = rrfc.outF.replace('csv', 'gds') 
  print(gdsFile)
  gdspy.LayoutViewer(lib) 
  # Export GDS
  lib.write_gds(gdsFile)
