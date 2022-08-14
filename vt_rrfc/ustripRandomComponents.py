import csv
import math
import random
import numpy as np
import gdspy
import pya
from vt_rrfc import microstrip as ustrip 
from skimage import measure

class rrfc:
  def __init__(rrfc_param, 
               unit,
               ports: int, 
               sides: int, 
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
  random.seed(rrfc.seed)

  # Pixel Size
  width_launch, length_launch = ustrip.microstrip_calc.synthMicrostrip(sub, imp, 30)
  launch_l_pixels = round(length_launch/rrfc.pixelSize) # length of the line to connect to structure in number of pixels
  launch_w_pixels = round(width_launch/rrfc.pixelSize) # width of the line to connect to structure in number of pixels
  length_launch = launch_l_pixels*rrfc.pixelSize
  width_launch = launch_w_pixels*rrfc.pixelSize

  # Set horizontal and vertical pixel limits
  y_pixels = math.ceil(y_dim/rrfc.pixelSize) 
  x_pixels = math.ceil(x_dim/rrfc.pixelSize)
  if y_pixels % 2 == 0:
    y_pixels += 1
  if x_pixels % 2 == 0:
    x_pixels += 1

  if rrfc.ports == 2:
    if rrfc.sides == 2:
      x_total = x_pixels*rrfc.pixelSize
      y_total = y_pixels*rrfc.pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, y_total/2, x_total, y_total/2, 0, 0, 0, 0] 
    else: 
      print('For a 2-port network, the number of sides must be equal to 2.')
      quit()
  elif rrfc.ports == 3:
    if rrfc.sides == 2:
      x_total = x_pixels*rrfc.pixelSize
      y_total = y_pixels*rrfc.pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, (math.ceil(int(y_pixels/4))+0.5)*rrfc.pixelSize, 0, (math.ceil(int(3*y_pixels/4))+0.5)*rrfc.pixelSize, \
                 x_total, (math.ceil(int(y_pixels/2))+0.5)*rrfc.pixelSize, 0, 0]
    elif rrfc.sides == 3:
      x_total = x_pixels*rrfc.pixelSize
      y_total = y_pixels*rrfc.pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, (math.ceil(int(y_pixels/2))+0.5)*rrfc.pixelSize, x_total, (math.ceil(int(y_pixels/2))+0.5)*rrfc.pixelSize, \
                 (int(x_total/rrfc.pixelSize)/2)*rrfc.pixelSize, y_total, 0, 0]
    else:
      print('For a 3-port network, the number of sides must be equal to either 2 or 3.')
      quit()
  elif rrfc.ports == 4:
    if rrfc.sides == 2:
      x_total = x_pixels*rrfc.pixelSize
      y_total = y_pixels*rrfc.pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, (math.ceil(int(y_pixels/4))+0.5)*rrfc.pixelSize, 0, (math.ceil(int(3*y_pixels/4))+0.5)*rrfc.pixelSize, \
                 x_total, (math.ceil(int(y_pixels/4))+0.5)*rrfc.pixelSize, x_total, (math.ceil(int(3*y_pixels/4))\
                 +0.5)*rrfc.pixelSize]
    elif rrfc.sides == 4:
      x_total = x_pixels*rrfc.pixelSize
      y_total = y_pixels*rrfc.pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, (math.ceil(int((y_total/rrfc.pixelSize)/2))+0.5)*rrfc.pixelSize, \
                 (int(x_total/rrfc.pixelSize)/2)*rrfc.pixelSize, y_total, \
                 x_total, (math.ceil(int((y_total/rrfc.pixelSize)/2))+0.5)*rrfc.pixelSize, \
                 (int(x_total/rrfc.pixelSize)/2)*rrfc.pixelSize, 0]
    else:
      print('For a 4-port network, the number of sides must be equal to either 2 or 4.')
      quit()

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
    pixel_map = np.zeros((int(x_total/rrfc.pixelSize),int(y_total/rrfc.pixelSize)),dtype=int)
    outline = UNIT.shapes(l_bottom).insert( pya.Box(0, 0, x_total, y_total) ) 

    if rrfc.ports == 2:
      # Add launches: assume a rectangle with port 1 = West, 2 = East
      for x in range(launch_l_pixels):
        for y in range(launch_w_pixels):
          launch_1 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize, \
                     (y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
          pixel_map[x,y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)] = 1
          launch_2 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x_total - \
                     length_launch + x*rrfc.pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                     math.floor(launch_w_pixels/2))*rrfc.pixelSize)
          pixel_map[int(x_pixels) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                    math.floor(launch_w_pixels/2)] = 1
          UNIT.shapes(l_top).insert(launch_1)
          UNIT.shapes(l_top).insert(launch_2)
        y += 1
      x += 1
    elif rrfc.ports ==3:
      if rrfc.sides == 2:
        # Add launches: assume a rectangle with port 1 = Southwest, 2 = Northwest, 3 = East
        for x in range(launch_l_pixels):
          for y in range(launch_w_pixels):
            launch_1 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize, \
                       (y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[x,y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
            launch_2 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize, \
                       (y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[x,y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
            launch_3 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x_total - \
                       length_launch + x*rrfc.pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                       math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[int(x_total/rrfc.pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                      math.floor(launch_w_pixels/2)] = 1
            UNIT.shapes(l_top).insert(launch_1)
            UNIT.shapes(l_top).insert(launch_2)
            UNIT.shapes(l_top).insert(launch_3)
          y += 1
        x += 1
      elif rrfc.sides == 3:
        # Add launches: assume a rectangle with port 1 = west, 2 = east, 3 = north
        for x in range(launch_l_pixels):
          for y in range(launch_w_pixels):
            launch_1 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize, \
                       (y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[x,y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)] = 1
            launch_2 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x_total - \
                       length_launch + x*rrfc.pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                       math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[int(x_total/rrfc.pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                      math.floor(launch_w_pixels/2)] = 1
            launch_3 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(\
                       (y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2)) * \
                       rrfc.pixelSize, y_total - length_launch + x*rrfc.pixelSize)
            pixel_map[y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2), \
                      int(y_total/rrfc.pixelSize)-launch_l_pixels + x] = 1
            UNIT.shapes(l_top).insert(launch_1)
            UNIT.shapes(l_top).insert(launch_2)
            UNIT.shapes(l_top).insert(launch_3)
          y += 1
        x += 1
    elif rrfc.ports == 4:
      if rrfc.sides == 2:
        # Add launches: assume a rectangle with port 1 = southwest, 2 = northwest, 3 = southeast, 4 = northeast
        for x in range(launch_l_pixels):
          for y in range(launch_w_pixels):
            launch_1 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize, \
                       (y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[x,y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
            launch_2 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize, \
                       (y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[x,y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
            launch_3 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x_total - \
                       length_launch + x*rrfc.pixelSize, (y+math.ceil(int(y_pixels/4)) - \
                       math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[int(x_total/rrfc.pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/4)) - \
                      math.floor(launch_w_pixels/2)] = 1
            launch_4 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x_total - \
                       length_launch + x*rrfc.pixelSize, (y+math.ceil(int(3*y_pixels/4)) - \
                       math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[int(x_total/rrfc.pixelSize) - launch_l_pixels + x, y+math.ceil(int(3*y_pixels/4)) - \
                      math.floor(launch_w_pixels/2)] = 1
            UNIT.shapes(l_top).insert(launch_1)
            UNIT.shapes(l_top).insert(launch_2)
            UNIT.shapes(l_top).insert(launch_3)
            UNIT.shapes(l_top).insert(launch_4)
          y += 1
        x += 1
      elif rrfc.sides == 4:
        # Add launches: assume a rectangle with port 1 = west, 2 = north, 3 = east, 4 = south
        for x in range(launch_l_pixels):
          for y in range(launch_w_pixels):
            launch_1 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize, \
                       (y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)+launch_l_pixels)*rrfc.pixelSize)
            pixel_map[x,y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)+launch_l_pixels] = 1
            launch_2 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(\
                       (y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2)) * \
                       rrfc.pixelSize, y_total - length_launch + x*rrfc.pixelSize)
            pixel_map[y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2), \
                      int(y_total/rrfc.pixelSize)-launch_l_pixels + x] = 1
            launch_3 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x_total - \
                       length_launch + x*rrfc.pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                       math.floor(launch_w_pixels/2)+launch_l_pixels)*rrfc.pixelSize)
            pixel_map[int(x_total/rrfc.pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                      math.floor(launch_w_pixels/2) + launch_l_pixels] = 1
            launch_4 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(\
                       (y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2)) * \
                       rrfc.pixelSize, x*rrfc.pixelSize)
            pixel_map[y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2), x] = 1
            UNIT.shapes(l_top).insert(launch_1)
            UNIT.shapes(l_top).insert(launch_2)
            UNIT.shapes(l_top).insert(launch_3)
            UNIT.shapes(l_top).insert(launch_4)
          y += 1
        x += 1

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

    # Design Random Structure
    # pixel_map = np.zeros((x_pixels,y_pixels),dtype=int)
    for x in range(x_pixels-2*launch_l_pixels):
      for y in range(y_pixels):
        bit = random.random()
        if bit >= 0.5:
        # creates a rectangle for the "on" pixel
          if rrfc.ports == 2:
            if rrfc.sym == 'x-axis':
              if y >= y_pixels/2 - 0:
                break
              else:
                rect1 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize \
                        +length_launch, y*rrfc.pixelSize)
                rect2 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize \
                        +length_launch, y_dim - (y+1)*rrfc.pixelSize)
                UNIT.shapes(l_top).insert(rect1)
                UNIT.shapes(l_top).insert(rect2)
                pixel_map[x+launch_l_pixels,y] = 1;
                pixel_map[x+launch_l_pixels,y_pixels - (y+1)] = 1;
            elif rrfc.sym == 'y-axis':
              if x >= (x_pixels-2*launch_l_pixels)/2 - 0:
                break
              else:
                rect1 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize \
                        +length_launch, y*rrfc.pixelSize)
                rect2 = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x_dim \
                        -x*rrfc.pixelSize-length_launch, y*rrfc.pixelSize)
                UNIT.shapes(l_top).insert(rect1)
                UNIT.shapes(l_top).insert(rect2)
            else: #Asymmetric
              rect = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize \
                     + length_launch, y*rrfc.pixelSize)
              UNIT.shapes(l_top).insert(rect)
              pixel_map[x+launch_l_pixels,y] = 1;
          elif rrfc.ports == 3:
            if rrfc.sides == 2:
              rect = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize \
                     + length_launch, y*rrfc.pixelSize)
              UNIT.shapes(l_top).insert(rect)
              pixel_map[x+launch_l_pixels,y] = 1;
            elif rrfc.sides == 3:
              rect = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize \
                     + length_launch, y*rrfc.pixelSize)
              UNIT.shapes(l_top).insert(rect)
              pixel_map[x+launch_l_pixels,y] = 1;
          elif rrfc.ports == 4:
            if rrfc.sides == 2:
              rect = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize \
                     + length_launch, y*rrfc.pixelSize)
              UNIT.shapes(l_top).insert(rect)
              pixel_map[x+launch_l_pixels,y] = 1;
            elif rrfc.sides == 4:
              rect = pya.Box(0, 0, rrfc.pixelSize, rrfc.pixelSize).moved(x*rrfc.pixelSize \
                     + length_launch, y*rrfc.pixelSize + length_launch)
              UNIT.shapes(l_top).insert(rect)
              pixel_map[x+launch_l_pixels,y+launch_l_pixels] = 1;
          # pixel_map[x+launch_l_pixels,y] = 1;
        y += 1  
      x += 1

    c12, c13, c14, c23, c24, c34 = ustrip.findConnectivity(pixel_map, rrfc.pixelSize, rrfc.ports, portPos)
  
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

    # Create layer #'s
    if rrfc.ports == 1:
      l_port1 = {"layer": 1, "datatype": 0}
    if rrfc.ports == 2:
      l_port1 = {"layer": 1, "datatype": 0}
      l_port2 = {"layer": 2, "datatype": 0}
    if rrfc.ports == 3:
      l_port1 = {"layer": 1, "datatype": 0}
      l_port2 = {"layer": 2, "datatype": 0}
      l_port3 = {"layer": 3, "datatype": 0}
    if rrfc.ports == 4:
      l_port1 = {"layer": 1, "datatype": 0}
      l_port2 = {"layer": 2, "datatype": 0}
      l_port3 = {"layer": 3, "datatype": 0}
      l_port4 = {"layer": 4, "datatype": 0}
    l_bottom = {"layer": 10, "datatype": 0}
    l_top = {"layer": 11, "datatype": 0}
    l_sources = {"layer": 5, "datatype": 0}

    # Draw outline
    pixel_map = np.zeros((int(x_total/rrfc.pixelSize),int(y_total/rrfc.pixelSize)),dtype=int)
    outline = gdspy.Rectangle((0, 0), (x_total, y_total), **l_bottom)
    UNIT.add(outline) 
  
    if rrfc.ports == 2:
      # Add launches: assume a rectangle with port 1 = West, 2 = East
      for x in range(launch_l_pixels):
        for y in range(launch_w_pixels):
          launch_1 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize, \
                     (y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
          pixel_map[x,y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)] = 1
          launch_2 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x_total - \
                     length_launch + x*rrfc.pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                     math.floor(launch_w_pixels/2))*rrfc.pixelSize)
          pixel_map[int(x_pixels) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                    math.floor(launch_w_pixels/2)] = 1
          UNIT.add(launch_1)
          UNIT.add(launch_2)
        y += 1
      x += 1
    elif rrfc.ports ==3:
      if rrfc.sides == 2:
        # Add launches: assume a rectangle with port 1 = Southwest, 2 = Northwest, 3 = East
        for x in range(launch_l_pixels):
          for y in range(launch_w_pixels):
            launch_1 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize, \
                       (y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[x,y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
            launch_2 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize, \
                       (y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[x,y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
            launch_3 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x_total - \
                       length_launch + x*rrfc.pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                       math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[int(x_total/rrfc.pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                      math.floor(launch_w_pixels/2)] = 1
            UNIT.add(launch_1)
            UNIT.add(launch_2)
            UNIT.add(launch_3)
          y += 1
        x += 1
      elif rrfc.sides == 3:
        # Add launches: assume a rectangle with port 1 = west, 2 = east, 3 = north
        for x in range(launch_l_pixels):
          for y in range(launch_w_pixels):
            launch_1 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize, \
                       (y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[x,y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)] = 1
            launch_2 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x_total - \
                       length_launch + x*rrfc.pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                       math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[int(x_total/rrfc.pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                      math.floor(launch_w_pixels/2)] = 1
            launch_3 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(\
                       (y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2)) * \
                       rrfc.pixelSize, y_total - length_launch + x*rrfc.pixelSize)
            pixel_map[y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2), \
                      int(y_total/rrfc.pixelSize)-launch_l_pixels + x] = 1
            UNIT.add(launch_1)
            UNIT.add(launch_2)
            UNIT.add(launch_3)
          y += 1
        x += 1
    elif rrfc.ports == 4:
      if rrfc.sides == 2:
        # Add launches: assume a rectangle with port 1 = southwest, 2 = northwest, 3 = southeast, 4 = northeast
        for x in range(launch_l_pixels):
          for y in range(launch_w_pixels):
            launch_1 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize, \
                       (y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[x,y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
            launch_2 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize, \
                       (y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[x,y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
            launch_3 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x_total - \
                       length_launch + x*rrfc.pixelSize, (y+math.ceil(int(y_pixels/4)) - \
                       math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[int(x_total/rrfc.pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/4)) - \
                      math.floor(launch_w_pixels/2)] = 1
            launch_4 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x_total - \
                       length_launch + x*rrfc.pixelSize, (y+math.ceil(int(3*y_pixels/4)) - \
                       math.floor(launch_w_pixels/2))*rrfc.pixelSize)
            pixel_map[int(x_total/rrfc.pixelSize) - launch_l_pixels + x, y+math.ceil(int(3*y_pixels/4)) - \
                      math.floor(launch_w_pixels/2)] = 1
            UNIT.add(launch_1)
            UNIT.add(launch_2)
            UNIT.add(launch_3)
            UNIT.add(launch_4)
          y += 1
        x += 1
      elif rrfc.sides == 4:
        # Add launches: assume a rectangle with port 1 = west, 2 = north, 3 = east, 4 = south
        for x in range(launch_l_pixels):
          for y in range(launch_w_pixels):
            launch_1 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize, \
                       (y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)+launch_l_pixels)*rrfc.pixelSize)
            pixel_map[x,y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)+launch_l_pixels] = 1
            launch_2 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(\
                       (y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2)) * \
                       rrfc.pixelSize, y_total - length_launch + x*rrfc.pixelSize)
            pixel_map[y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2), \
                      int(y_total/rrfc.pixelSize)-launch_l_pixels + x] = 1
            launch_3 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x_total - \
                       length_launch + x*rrfc.pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                       math.floor(launch_w_pixels/2)+launch_l_pixels)*rrfc.pixelSize)
            pixel_map[int(x_total/rrfc.pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                      math.floor(launch_w_pixels/2) + launch_l_pixels] = 1
            launch_4 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(\
                       (y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2)) * \
                       rrfc.pixelSize, x*rrfc.pixelSize)
            pixel_map[y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2), x] = 1
            UNIT.add(launch_1)
            UNIT.add(launch_2)
            UNIT.add(launch_3)
            UNIT.add(launch_4)
          y += 1
        x += 1

    # Add ports and sources to the gds
    if rrfc.ports == 2:
      if rrfc.sim == 'MEEP':
        port_1 = [(length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_1 = gdspy.Polygon(port_1, **l_port1)
        UNIT.add(poly_1)
        port_2 = [(x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_2 = gdspy.Polygon(port_2, **l_port2)
        UNIT.add(poly_2)
        source = [(1, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (1, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_3 = gdspy.Polygon(source, **l_sources)
        UNIT.add(poly_3)
      else:
        t = 1 #print('If you are using ADS, you do not need to add ports. They are added with AEL')
    elif rrfc.ports == 3:
      if rrfc.sides == 2:
        if rrfc.sim == 'MEEP':
          port_1 = [(length_launch*0.5, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                    (length_launch*0.5, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
          poly_1 = gdspy.Polygon(port_1, **l_port1)
          UNIT.add(poly_1)
          port_2 = [(length_launch*0.5, (math.ceil(int(3*y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                    (length_launch*0.5, (math.ceil(int(3*y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
          poly_2 = gdspy.Polygon(port_2, **l_port2)
          UNIT.add(poly_2)
          port_3 = [(x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                    (x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
          poly_3 = gdspy.Polygon(port_3, **l_port3)
          UNIT.add(poly_3)
          source = [(1, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                    (1, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
          poly_4 = gdspy.Polygon(source, **l_sources)
          UNIT.add(poly_4)
        else:
          t = 1 #print('If you are using ADS, you do not need to add ports. They are added with AEL')
      elif rrfc.sides == 3:
        if rrfc.sim == 'MEEP':
          port_1 = [(length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                    (length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
          poly_1 = gdspy.Polygon(port_1, **l_port1)
          UNIT.add(poly_1)
          port_2 = [(x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                    (x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
          poly_2 = gdspy.Polygon(port_2, **l_port2)
          UNIT.add(poly_2)
          port_3 = [((int(x_total/rrfc.pixelSize)/2-1.05*(launch_w_pixels/2))*rrfc.pixelSize, y_total - length_launch*0.5), \
                    ((int(x_total/rrfc.pixelSize)/2+1.05*(launch_w_pixels/2))*rrfc.pixelSize, y_total - length_launch*0.5)]
          poly_3 = gdspy.Polygon(port_3, **l_port3)
          UNIT.add(poly_3)
          source = [(1, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                    (1, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
          poly_4 = gdspy.Polygon(source, **l_sources)
          UNIT.add(poly_4)
        else:
          t = 1 #print('If you are using ADS, you do not need to add ports. They are added with AEL')
    elif rrfc.ports == 4:
      if rrfc.sides == 2:
        if rrfc.sim == 'MEEP':
          port_1 = [(length_launch*0.5, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                    (length_launch*0.5, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
          poly_1 = gdspy.Polygon(port_1, **l_port1)
          UNIT.add(poly_1)
          port_2 = [(length_launch*0.5, (math.ceil(int(3*y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                    (length_launch*0.5, (math.ceil(int(3*y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
          poly_2 = gdspy.Polygon(port_2, **l_port2)
          UNIT.add(poly_2)
          port_3 = [(x_total - length_launch*0.5, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                    (x_total - length_launch*0.5, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
          poly_3 = gdspy.Polygon(port_3, **l_port3)
          UNIT.add(poly_3)
          port_4 = [(x_total - length_launch*0.5, (math.ceil(int(3*y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                    (x_total - length_launch*0.5, (math.ceil(int(3*y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
          poly_4 = gdspy.Polygon(port_4, **l_port4)
          UNIT.add(poly_4)
          source = [(1, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                    (1, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
          poly_5 = gdspy.Polygon(source, **l_sources)
          UNIT.add(poly_5)
        else:
          t = 1 #print('If you are using ADS, you do not need to add ports. They are added with AEL')
      elif rrfc.sides == 4:
        if rrfc.sim == 'MEEP':
          port_1 = [(length_launch*0.5, (math.ceil(int((y_total/rrfc.pixelSize)/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                    (length_launch*0.5, (math.ceil(int((y_total/rrfc.pixelSize)/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
          poly_1 = gdspy.Polygon(port_1, **l_port1)
          UNIT.add(poly_1)
          port_2 = [((int(x_total/rrfc.pixelSize)/2-1.05*(launch_w_pixels/2))*rrfc.pixelSize, y_total - length_launch*0.5), \
                    ((int(x_total/rrfc.pixelSize)/2+1.05*(launch_w_pixels/2))*rrfc.pixelSize, y_total - length_launch*0.5)]
          poly_2 = gdspy.Polygon(port_2, **l_port2)
          UNIT.add(poly_2)
          port_3 = [(x_total - length_launch*0.5, (math.ceil(int((y_total/rrfc.pixelSize)/2))-\
                     1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                    (x_total - length_launch*0.5, (math.ceil(int((y_total/rrfc.pixelSize)/2))+\
                     1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
          poly_3 = gdspy.Polygon(port_3, **l_port3)
          UNIT.add(poly_3)
          port_4 = [((int(x_total/rrfc.pixelSize)/2-1.05*(launch_w_pixels/2))*rrfc.pixelSize, length_launch*0.5), \
                    ((int(x_total/rrfc.pixelSize)/2+1.05*(launch_w_pixels/2))*rrfc.pixelSize, length_launch*0.5)]
          poly_4 = gdspy.Polygon(port_4, **l_port4)
          UNIT.add(poly_4)
          source = [(1, (math.ceil(int((y_total/rrfc.pixelSize)/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                    (1, (math.ceil(int((y_total/rrfc.pixelSize)/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
          poly_5 = gdspy.Polygon(source, **l_sources)
          UNIT.add(poly_5)
        else:
          t = 1 #print('If you are using ADS, you do not need to add ports. They are added with AEL')
  
    # Design Random Structure
    # pixel_map = np.zeros((x_pixels,y_pixels),dtype=int)
    for x in range(x_pixels-2*launch_l_pixels):
      for y in range(y_pixels):
        bit = random.random()
        if bit >= 0.5:
        # creates a rectangle for the "on" pixel
          if rrfc.ports == 2:
            if rrfc.sym == 'x-axis':
              if y >= y_pixels/2 - 0:
                break
              rect1 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize \
                      +length_launch, y*rrfc.pixelSize)
              rect2 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize \
                      +length_launch, y_dim - (y+1)*rrfc.pixelSize)
              UNIT.add(rect1)
              UNIT.add(rect2)
              pixel_map[x+launch_l_pixels,y] = 1;
              pixel_map[x+launch_l_pixels,y_pixels - (y+1)] = 1;
            if rrfc.sym == 'y-axis':
              if x >= (x_pixels-2*launch_l_pixels)/2 - 0:
                break
              rect1 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize \
                      +length_launch, y*rrfc.pixelSize)
              rect2 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x_dim \
                      -x*rrfc.pixelSize-length_launch, y*rrfc.pixelSize)
              UNIT.add(rect1)
              UNIT.add(rect2)
            else: #Asymmetric
              rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize \
                     + length_launch, y*rrfc.pixelSize)
              UNIT.add(rect)
              pixel_map[x+launch_l_pixels,y] = 1;
          elif rrfc.ports == 3:
            if rrfc.sides == 2:
              rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize \
                     + length_launch, y*rrfc.pixelSize)
              UNIT.add(rect)
              pixel_map[x+launch_l_pixels,y] = 1;
            elif rrfc.sides == 3:
              rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize \
                     + length_launch, y*rrfc.pixelSize)
              UNIT.add(rect)
              pixel_map[x+launch_l_pixels,y] = 1;
          elif rrfc.ports == 4:
            if rrfc.sides == 2:
              rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize \
                     + length_launch, y*rrfc.pixelSize)
              UNIT.add(rect)
              pixel_map[x+launch_l_pixels,y] = 1;
            elif rrfc.sides == 4:
              rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize \
                     + length_launch, y*rrfc.pixelSize + length_launch)
              UNIT.add(rect)
              pixel_map[x+launch_l_pixels,y+launch_l_pixels] = 1;
          # pixel_map[x+launch_l_pixels,y] = 1;
        y += 1  
      x += 1

    c12, c13, c14, c23, c24, c34 = ustrip.findConnectivity(pixel_map, rrfc.pixelSize, rrfc.ports, portPos)
  
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

  return portPos, x_total, y_total, csvFile, gdsFile, cellName

def randomGDS(sub,
              rrfc,
              x_el_units: int,
              y_el_units: int,
              imp: float,
              eLength: float):
  """
  This script will create a rectangular pixel grid that is x_el_units*eLength wide by y_el_units*eLength tall 
  according to the following rules:
  eLength = Electrical Length in Degrees
  x_el_units = Number of electrical length units wide (east-to-west)
  y_el_units = Number of electrical length units tall (north-to-south)
  pixelSize = Size of the pixel in library units (code defaults to mils, but can be modified by changing lib.unit below)
  imp = characteristic impedance of the launch paths from the ports
  
  The output of the script is a pixel map, which is a CSV file with ones and zeros in the positions of the pixels in the grid, 
  and a GDS file that can be imported into an EM simulator. It will also output an array with port positions that can be used 
  for placing ports in an ADS simulation. 
  """
  random.seed(rrfc.seed)
  lib = gdspy.GdsLibrary()

  # Set the database unit to 1 mil
  lib.unit = rrfc.unit

  # Create Cell obj
  cellName = 'RANDOM'
  gdspy.current_library = gdspy.GdsLibrary() # This line of code has to be here to reset the GDS library on every loop
  UNIT = lib.new_cell(cellName)

  # Create layer #'s
  if rrfc.ports == 1:
    l_port1 = {"layer": 1, "datatype": 0}
  if rrfc.ports == 2:
    l_port1 = {"layer": 1, "datatype": 0}
    l_port2 = {"layer": 2, "datatype": 0}
  if rrfc.ports == 3:
    l_port1 = {"layer": 1, "datatype": 0}
    l_port2 = {"layer": 2, "datatype": 0}
    l_port3 = {"layer": 3, "datatype": 0}
  if rrfc.ports == 4:
    l_port1 = {"layer": 1, "datatype": 0}
    l_port2 = {"layer": 2, "datatype": 0}
    l_port3 = {"layer": 3, "datatype": 0}
    l_port4 = {"layer": 4, "datatype": 0}
  l_bottom = {"layer": 10, "datatype": 0}
  l_top = {"layer": 11, "datatype": 0}
  l_sources = {"layer": 5, "datatype": 0}

  # Metal dimensions
  #width_50 = 37.6 # 50 Ohm line width at 2.4 GHz
  #length_90 = 687.5 # Quarter-wavelength line at 2.4 GHz
  width_50, elec_length = ustrip.microstrip_calc.synthMicrostrip(sub, imp, eLength)
#  launch_l_pixels = 10 # length of the line to connect to structure in number of pixels

  # Pixel Size
  width_launch, length_launch = ustrip.microstrip_calc.synthMicrostrip(sub, imp, 30)
  launch_l_pixels = round(length_launch/rrfc.pixelSize) # length of the line to connect to structure in number of pixels
  launch_w_pixels = round(width_launch/rrfc.pixelSize) # width of the line to connect to structure in number of pixels
  length_launch = launch_l_pixels*rrfc.pixelSize
  width_launch = launch_w_pixels*rrfc.pixelSize

  # Set horizontal and vertical pixel limits
  y_units = y_el_units # integer units of electrical length
  x_units = x_el_units # integer units of electrical length
  y_dim = math.floor(y_units*elec_length)
  x_dim = math.floor(x_units*elec_length)
  y_pixels = math.ceil(y_dim/rrfc.pixelSize)
  x_pixels = math.ceil(x_dim/rrfc.pixelSize)
  if y_pixels % 2 == 0:
    y_pixels += 1
  if x_pixels % 2 == 0:
    x_pixels += 1

  if rrfc.ports == 2:
    if rrfc.sides == 2:
      x_total = (2*launch_l_pixels + x_pixels)*rrfc.pixelSize
      y_total = y_pixels*rrfc.pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, y_total/2, x_total, y_total/2, 0, 0, 0, 0] 
    else: 
      print('For a 2-port network, the number of sides must be equal to 2.')
      quit()
  elif rrfc.ports == 3:
    if rrfc.sides == 2:
      x_total = (2*launch_l_pixels + x_pixels)*rrfc.pixelSize
      y_total = y_pixels*rrfc.pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, (math.ceil(int(y_pixels/4))+0.5)*rrfc.pixelSize, 0, (math.ceil(int(3*y_pixels/4))+0.5)*rrfc.pixelSize, \
                 x_total, (math.ceil(int(y_pixels/2))+0.5)*rrfc.pixelSize, 0, 0]
    elif rrfc.sides == 3:
      x_total = (2*launch_l_pixels + x_pixels)*rrfc.pixelSize
      y_total = (launch_l_pixels + y_pixels)*rrfc.pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, (math.ceil(int(y_pixels/2))+0.5)*rrfc.pixelSize, x_total, (math.ceil(int(y_pixels/2))+0.5)*rrfc.pixelSize, \
                 (int(x_total/rrfc.pixelSize)/2)*rrfc.pixelSize, y_total, 0, 0]
    else:
      print('For a 3-port network, the number of sides must be equal to either 2 or 3.')
      quit()
  elif rrfc.ports == 4:
    if rrfc.sides == 2:
      x_total = (2*launch_l_pixels + x_pixels)*rrfc.pixelSize
      y_total = y_pixels*rrfc.pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, (math.ceil(int(y_pixels/4))+0.5)*rrfc.pixelSize, 0, (math.ceil(int(3*y_pixels/4))+0.5)*rrfc.pixelSize, \
                 x_total, (math.ceil(int(y_pixels/4))+0.5)*rrfc.pixelSize, x_total, (math.ceil(int(3*y_pixels/4))\
                 +0.5)*rrfc.pixelSize]
    elif rrfc.sides == 4:
      x_total = (2*launch_l_pixels + x_pixels)*rrfc.pixelSize
      y_total = (2*launch_l_pixels + y_pixels)*rrfc.pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, (math.ceil(int((y_total/rrfc.pixelSize)/2))+0.5)*rrfc.pixelSize, \
                 (int(x_total/rrfc.pixelSize)/2)*rrfc.pixelSize, y_total, \
                 x_total, (math.ceil(int((y_total/rrfc.pixelSize)/2))+0.5)*rrfc.pixelSize, \
                 (int(x_total/rrfc.pixelSize)/2)*rrfc.pixelSize, 0]
    else:
      print('For a 4-port network, the number of sides must be equal to either 2 or 4.')
      quit()

  # Draw outline
  pixel_map = np.zeros((int(x_total/rrfc.pixelSize),int(y_total/rrfc.pixelSize)),dtype=int)
  outline = gdspy.Rectangle((0, 0), (x_total, y_total), **l_bottom)
  UNIT.add(outline) 
  
  if rrfc.ports == 2:
    # Add launches: assume a rectangle with port 1 = West, 2 = East
    for x in range(launch_l_pixels):
      for y in range(launch_w_pixels):
        launch_1 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize, \
                   (y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
        pixel_map[x,y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)] = 1
        launch_2 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x_total - \
                   length_launch + x*rrfc.pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                   math.floor(launch_w_pixels/2))*rrfc.pixelSize)
        pixel_map[int(x_pixels) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                  math.floor(launch_w_pixels/2)] = 1
        UNIT.add(launch_1)
        UNIT.add(launch_2)
      y += 1
    x += 1
  elif rrfc.ports ==3:
    if rrfc.sides == 2:
      # Add launches: assume a rectangle with port 1 = Southwest, 2 = Northwest, 3 = East
      for x in range(launch_l_pixels):
        for y in range(launch_w_pixels):
          launch_1 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize, \
                     (y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
          pixel_map[x,y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
          launch_2 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize, \
                     (y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
          pixel_map[x,y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
          launch_3 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x_total - \
                     length_launch + x*rrfc.pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                     math.floor(launch_w_pixels/2))*rrfc.pixelSize)
          pixel_map[int(x_total/rrfc.pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                    math.floor(launch_w_pixels/2)] = 1
          UNIT.add(launch_1)
          UNIT.add(launch_2)
          UNIT.add(launch_3)
        y += 1
      x += 1
    elif rrfc.sides == 3:
      # Add launches: assume a rectangle with port 1 = west, 2 = east, 3 = north
      for x in range(launch_l_pixels):
        for y in range(launch_w_pixels):
          launch_1 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize, \
                     (y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
          pixel_map[x,y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)] = 1
          launch_2 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x_total - \
                     length_launch + x*rrfc.pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                     math.floor(launch_w_pixels/2))*rrfc.pixelSize)
          pixel_map[int(x_total/rrfc.pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                    math.floor(launch_w_pixels/2)] = 1
          launch_3 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(\
                     (y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2)) * \
                     rrfc.pixelSize, y_total - length_launch + x*rrfc.pixelSize)
          pixel_map[y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2), \
                    int(y_total/rrfc.pixelSize)-launch_l_pixels + x] = 1
          UNIT.add(launch_1)
          UNIT.add(launch_2)
          UNIT.add(launch_3)
        y += 1
      x += 1
  elif rrfc.ports == 4:
    if rrfc.sides == 2:
      # Add launches: assume a rectangle with port 1 = southwest, 2 = northwest, 3 = southeast, 4 = northeast
      for x in range(launch_l_pixels):
        for y in range(launch_w_pixels):
          launch_1 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize, \
                     (y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
          pixel_map[x,y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
          launch_2 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize, \
                     (y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2))*rrfc.pixelSize)
          pixel_map[x,y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
          launch_3 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x_total - \
                     length_launch + x*rrfc.pixelSize, (y+math.ceil(int(y_pixels/4)) - \
                     math.floor(launch_w_pixels/2))*rrfc.pixelSize)
          pixel_map[int(x_total/rrfc.pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/4)) - \
                    math.floor(launch_w_pixels/2)] = 1
          launch_4 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x_total - \
                     length_launch + x*rrfc.pixelSize, (y+math.ceil(int(3*y_pixels/4)) - \
                     math.floor(launch_w_pixels/2))*rrfc.pixelSize)
          pixel_map[int(x_total/rrfc.pixelSize) - launch_l_pixels + x, y+math.ceil(int(3*y_pixels/4)) - \
                    math.floor(launch_w_pixels/2)] = 1
          UNIT.add(launch_1)
          UNIT.add(launch_2)
          UNIT.add(launch_3)
          UNIT.add(launch_4)
        y += 1
      x += 1
    elif rrfc.sides == 4:
      # Add launches: assume a rectangle with port 1 = west, 2 = north, 3 = east, 4 = south
      for x in range(launch_l_pixels):
        for y in range(launch_w_pixels):
          launch_1 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize, \
                     (y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)+launch_l_pixels)*rrfc.pixelSize)
          pixel_map[x,y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)+launch_l_pixels] = 1
          launch_2 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(\
                     (y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2)) * \
                     rrfc.pixelSize, y_total - length_launch + x*rrfc.pixelSize)
          pixel_map[y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2), \
                    int(y_total/rrfc.pixelSize)-launch_l_pixels + x] = 1
          launch_3 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x_total - \
                     length_launch + x*rrfc.pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                     math.floor(launch_w_pixels/2)+launch_l_pixels)*rrfc.pixelSize)
          pixel_map[int(x_total/rrfc.pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                    math.floor(launch_w_pixels/2) + launch_l_pixels] = 1
          launch_4 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(\
                     (y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2)) * \
                     rrfc.pixelSize, x*rrfc.pixelSize)
          pixel_map[y+math.ceil(int(x_total / rrfc.pixelSize/2))-math.floor(launch_w_pixels/2), x] = 1
          UNIT.add(launch_1)
          UNIT.add(launch_2)
          UNIT.add(launch_3)
          UNIT.add(launch_4)
        y += 1
      x += 1

  # Add ports and sources to the gds
  if rrfc.ports == 2:
    if rrfc.sim == 'ADS':
      port_1 = [(length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                (length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
      poly_1 = gdspy.Polygon(port_1, **l_port1)
      UNIT.add(poly_1)
      port_2 = [(x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                (x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
      poly_2 = gdspy.Polygon(port_2, **l_port2)
      UNIT.add(poly_2)
      source = [(1, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                (1, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
      poly_3 = gdspy.Polygon(source, **l_sources)
      UNIT.add(poly_3)
    elif rrfc.sim == 'EMX':
      port_1 = gdspy.Label("p1", (portPos[0], portPos[1]), "w", layer=11)
      UNIT.add(port_1)
      port_2 = gdspy.Label("p2", (portPos[2], portPos[3]), "e", layer=11)
      UNIT.add(port_2)
    else:
      print('You must choose an available simulator')
      quit()
  elif rrfc.ports == 3:
    if rrfc.sides == 2:
      if rrfc.sim == 'ADS':
        port_1 = [(length_launch*0.5, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (length_launch*0.5, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_1 = gdspy.Polygon(port_1, **l_port1)
        UNIT.add(poly_1)
        port_2 = [(length_launch*0.5, (math.ceil(int(3*y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (length_launch*0.5, (math.ceil(int(3*y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_2 = gdspy.Polygon(port_2, **l_port2)
        UNIT.add(poly_2)
        port_3 = [(x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_3 = gdspy.Polygon(port_3, **l_port3)
        UNIT.add(poly_3)
        source = [(1, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (1, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_4 = gdspy.Polygon(source, **l_sources)
        UNIT.add(poly_4)
      elif rrfc.sim == 'EMX':
        port_1 = gdspy.Label("p1", (portPos[0], portPos[1]), "w", layer=11)
        UNIT.add(port_1)
        port_2 = gdspy.Label("p2", (portPos[2], portPos[3]), "w", layer=11)
        UNIT.add(port_2)
        port_3 = gdspy.Label("p3", (portPos[4], portPos[5]), "e", layer=11)
        UNIT.add(port_3)
      else:
        print('You must choose an available simulator')
        quit()
    elif rrfc.sides == 3:
      if rrfc.sim == 'ADS':
        port_1 = [(length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_1 = gdspy.Polygon(port_1, **l_port1)
        UNIT.add(poly_1)
        port_2 = [(x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_2 = gdspy.Polygon(port_2, **l_port2)
        UNIT.add(poly_2)
        port_3 = [((int(x_total/rrfc.pixelSize)/2-1.05*(launch_w_pixels/2))*rrfc.pixelSize, y_total - length_launch*0.5), \
                  ((int(x_total/rrfc.pixelSize)/2+1.05*(launch_w_pixels/2))*rrfc.pixelSize, y_total - length_launch*0.5)]
        poly_3 = gdspy.Polygon(port_3, **l_port3)
        UNIT.add(poly_3)
        source = [(1, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (1, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_4 = gdspy.Polygon(source, **l_sources)
        UNIT.add(poly_4)
      elif rrfc.sim == 'EMX':
        port_1 = gdspy.Label("p1", (portPos[0], portPos[1]), "w", layer=11)
        UNIT.add(port_1)
        port_2 = gdspy.Label("p2", (portPos[2], portPos[3]), "e", layer=11)
        UNIT.add(port_2)
        port_3 = gdspy.Label("p3", (portPos[4], portPos[5]), "n", layer=11)
        UNIT.add(port_3)
      else:
        print('You must choose an available simulator')
        quit()
  elif rrfc.ports == 4:
    if rrfc.sides == 2:
      if rrfc.sim == 'ADS':
        port_1 = [(length_launch*0.5, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (length_launch*0.5, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_1 = gdspy.Polygon(port_1, **l_port1)
        UNIT.add(poly_1)
        port_2 = [(length_launch*0.5, (math.ceil(int(3*y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (length_launch*0.5, (math.ceil(int(3*y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_2 = gdspy.Polygon(port_2, **l_port2)
        UNIT.add(poly_2)
        port_3 = [(x_total - length_launch*0.5, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (x_total - length_launch*0.5, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_3 = gdspy.Polygon(port_3, **l_port3)
        UNIT.add(poly_3)
        port_4 = [(x_total - length_launch*0.5, (math.ceil(int(3*y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (x_total - length_launch*0.5, (math.ceil(int(3*y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_4 = gdspy.Polygon(port_4, **l_port4)
        UNIT.add(poly_4)
        source = [(1, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (1, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_5 = gdspy.Polygon(source, **l_sources)
        UNIT.add(poly_5)
      elif rrfc.sim == 'EMX':
        port_1 = gdspy.Label("p1", (portPos[0], portPos[1]), "w", layer=11)
        UNIT.add(port_1)
        port_2 = gdspy.Label("p2", (portPos[2], portPos[3]), "w", layer=11)
        UNIT.add(port_2)
        port_3 = gdspy.Label("p3", (portPos[4], portPos[5]), "e", layer=11)
        UNIT.add(port_3)
        port_4 = gdspy.Label("p4", (portPos[6], portPos[7]), "e", layer=11)
        UNIT.add(port_4)
      else:
        print('You must choose an available simulator')
        quit()
    elif rrfc.sides == 4:
      if rrfc.sim == 'ADS':
        port_1 = [(length_launch*0.5, (math.ceil(int((y_total/rrfc.pixelSize)/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (length_launch*0.5, (math.ceil(int((y_total/rrfc.pixelSize)/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_1 = gdspy.Polygon(port_1, **l_port1)
        UNIT.add(poly_1)
        port_2 = [((int(x_total/rrfc.pixelSize)/2-1.05*(launch_w_pixels/2))*rrfc.pixelSize, y_total - length_launch*0.5), \
                  ((int(x_total/rrfc.pixelSize)/2+1.05*(launch_w_pixels/2))*rrfc.pixelSize, y_total - length_launch*0.5)]
        poly_2 = gdspy.Polygon(port_2, **l_port2)
        UNIT.add(poly_2)
        port_3 = [(x_total - length_launch*0.5, (math.ceil(int((y_total/rrfc.pixelSize)/2))-\
                   1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (x_total - length_launch*0.5, (math.ceil(int((y_total/rrfc.pixelSize)/2))+\
                   1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_3 = gdspy.Polygon(port_3, **l_port3)
        UNIT.add(poly_3)
        port_4 = [((int(x_total/rrfc.pixelSize)/2-1.05*(launch_w_pixels/2))*rrfc.pixelSize, length_launch*0.5), \
                  ((int(x_total/rrfc.pixelSize)/2+1.05*(launch_w_pixels/2))*rrfc.pixelSize, length_launch*0.5)]
        poly_4 = gdspy.Polygon(port_4, **l_port4)
        UNIT.add(poly_4)
        source = [(1, (math.ceil(int((y_total/rrfc.pixelSize)/2))-1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize), \
                  (1, (math.ceil(int((y_total/rrfc.pixelSize)/2))+1.05*(launch_w_pixels/2)+0.5)*rrfc.pixelSize)]
        poly_5 = gdspy.Polygon(source, **l_sources)
        UNIT.add(poly_5)
      elif rrfc.sim == 'EMX':
        port_1 = gdspy.Label("p1", (portPos[0], portPos[1]), "w", layer=11)
        UNIT.add(port_1)
        port_2 = gdspy.Label("p2", (portPos[2], portPos[3]), "n", layer=11)
        UNIT.add(port_2)
        port_3 = gdspy.Label("p3", (portPos[4], portPos[5]), "e", layer=11)
        UNIT.add(port_3)
        port_4 = gdspy.Label("p4", (portPos[6], portPos[7]), "s", layer=11)
        UNIT.add(port_4)
      else:
        print('You must choose an available simulator')
        quit()

  # Design Random Structure
  # pixel_map = np.zeros((x_pixels,y_pixels),dtype=int)
  for x in range(x_pixels):
    for y in range(y_pixels):
      bit = random.random()
      if bit >= 0.5:
      # creates a rectangle for the "on" pixel
        if rrfc.ports == 2:
          if rrfc.sym == 'x-axis':
            if y >= y_pixels/2 - 0:
              break
            rect1 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize \
                    +length_launch, y*rrfc.pixelSize)
            rect2 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize \
                    +length_launch, y_dim - (y+1)*rrfc.pixelSize)
            UNIT.add(rect1)
            UNIT.add(rect2)
            pixel_map[x+launch_l_pixels,y] = 1;
            pixel_map[x+launch_l_pixels,y_pixels - (y+1)] = 1;
          if rrfc.sym == 'y-axis':
            if x >= (x_pixels-2*launch_l_pixels)/2 - 0:
              break
            rect1 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize \
                    +length_launch, y*rrfc.pixelSize)
            rect2 = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x_dim \
                    -x*rrfc.pixelSize-length_launch, y*rrfc.pixelSize)
            UNIT.add(rect1)
            UNIT.add(rect2)
          else: #Asymmetric
            rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize \
                   + length_launch, y*rrfc.pixelSize)
            UNIT.add(rect)
            pixel_map[x+launch_l_pixels,y] = 1;
        elif rrfc.ports == 3:
          if rrfc.sides == 2:
            rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize \
                   + length_launch, y*rrfc.pixelSize)
            pixel_map[x+launch_l_pixels,y] = 1;
          elif rrfc.sides == 3:
            rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize \
                   + length_launch, y*rrfc.pixelSize)
            pixel_map[x+launch_l_pixels,y] = 1;
        elif rrfc.ports == 4:
          if rrfc.sides == 2:
            rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize \
                   + length_launch, y*rrfc.pixelSize)
            pixel_map[x+launch_l_pixels,y] = 1;
          elif rrfc.sides == 4:
            rect = gdspy.Rectangle((0, 0), (rrfc.pixelSize, rrfc.pixelSize), **l_top).translate(x*rrfc.pixelSize \
                   + length_launch, y*rrfc.pixelSize + length_launch)
            pixel_map[x+launch_l_pixels,y+launch_l_pixels] = 1;
        UNIT.add(rect)
        # pixel_map[x+launch_l_pixels,y] = 1;
      y += 1  
    x += 1

  c12, c13, c14, c23, c24, c34 = ustrip.findConnectivity(pixel_map, rrfc.pixelSize, ports, portPos)
  
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


  return portPos, x_total, y_total, csvFile, gdsFile, cellName

def recreateGDS(rrfc):

  lib = gdspy.GdsLibrary()

  # Set the database unit to 1 mil
  lib.unit = rrfc.unit

  # Create Cell obj
  cellName = 'DESIGN'
  gdspy.current_library = gdspy.GdsLibrary() # This line of code has to be here to reset the GDS library on every loop
  UNIT = lib.new_cell(cellName)

  l_bottom = {"layer": 10, "datatype": 0}
  l_top = {"layer": 11, "datatype": 0}
  l_sources = {"layer": 5, "datatype": 0}

  pixMap = np.flipud(np.loadtxt(rrfc.outF, delimiter=','))
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

  gdsFile = rrfc.outF.replace('csv', 'gds') 
  print(gdsFile)
  gdspy.LayoutViewer(lib) 
  # Export GDS
  lib.write_gds(gdsFile)
