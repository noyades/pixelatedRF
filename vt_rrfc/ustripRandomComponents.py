import csv
import math
import random
import numpy as np
import gdspy
from vt_rrfc import microstrip as ustrip 
from skimage import measure

def randomGDS_dim(
              sub,
              filtOrder: int,
              ports: int,
              sides: int,
              x_dim: int,
              y_dim: int,
              pixelSize: int,
              imp: float,
              eLength: float,
              seedNum: int,
              connect, 
              sim: str,
              view: bool,
              write: bool,
              outF: str):
  """
  This script will create a rectangular pixel grid that is x_el_units*eLength wide by y_el_units*eLength tall 
  according to the following rules:
  eLength = Electrical Length in Degrees
  x_el_units = Number of electrical length units wide (east-to-west)
  y_el_units = Number of electrical length units tall (north-to-south)
  pixelSize = Size of the pixel in library units (code defaults to mils, but can be modified by changing lib.unit below)
  freq = frequency at which the wavelength is defined
  imp = characteristic impedance of the launch paths from the ports
  thick = thickness of the conductor in library units (code defaults to mils, but can be modified by changing lib.units below)
  height = thickness of the substrate in library units (code defaults to mils, but can be modified by changing lib.unit below) 
           Currently this code assumes that the design is a two layer design with conductors above and below the substrate
  relPerm = relative permittivity of the substrate material
  seedNum = fixed seed number for the random number generator. This allows the same design to be reproduced each time.
  sim = The type of simulator the GDS that is produced is meant to support. There are currently two options, EMX and 
        everything else. The primary reason is the way that ports are defined in EMX, where they must be defined by a
        label in the conductor layer.
  view = Do you want to open up a view window after each gds is created for inspection. This is primarily for de-bugging as new
         features are added. Value must be boolean (either True or False)
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
  
  The output of the script is a pixel map, which is a CSV file with ones and zeros in the positions of the pixels in the grid, 
  and a GDS file that can be imported into an EM simulator. It will also output an array with port positions that can be used 
  for placing ports in an ADS simulation. 
  """
  random.seed(seedNum)
  lib = gdspy.GdsLibrary()

  # Set the database unit to 1 mil
  lib.unit = 25.4e-6

  # Create Cell obj
  cellName = 'RANDOM'
  gdspy.current_library = gdspy.GdsLibrary() # This line of code has to be here to reset the GDS library on every loop
  UNIT = lib.new_cell(cellName)

  # Create layer #'s
  if ports == 1:
    l_port1 = {"layer": 1, "datatype": 0}
  if ports == 2:
    l_port1 = {"layer": 1, "datatype": 0}
    l_port2 = {"layer": 2, "datatype": 0}
  if ports == 3:
    l_port1 = {"layer": 1, "datatype": 0}
    l_port2 = {"layer": 2, "datatype": 0}
    l_port3 = {"layer": 3, "datatype": 0}
  if ports == 4:
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
  launch_l_pixels = round(length_launch/pixelSize) # length of the line to connect to structure in number of pixels
  launch_w_pixels = round(width_launch/pixelSize) # width of the line to connect to structure in number of pixels
  length_launch = launch_l_pixels*pixelSize
  width_launch = launch_w_pixels*pixelSize

  # Set horizontal and vertical pixel limits
  y_pixels = math.ceil(y_dim/pixelSize) 
  x_pixels = math.ceil(x_dim/pixelSize)
  if y_pixels % 2 == 0:
    y_pixels += 1
  if x_pixels % 2 == 0:
    x_pixels += 1

  if ports == 2:
    if sides == 2:
      x_total = x_pixels*pixelSize
      y_total = y_pixels*pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, y_total/2, x_total, y_total/2, 0, 0, 0, 0] 
    else: 
      print('For a 2-port network, the number of sides must be equal to 2.')
      quit()
  elif ports == 3:
    if sides == 2:
      x_total = x_pixels*pixelSize
      y_total = y_pixels*pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, (math.ceil(int(y_pixels/4))+0.5)*pixelSize, 0, (math.ceil(int(3*y_pixels/4))+0.5)*pixelSize, \
                 x_total, (math.ceil(int(y_pixels/2))+0.5)*pixelSize, 0, 0]
    elif sides == 3:
      x_total = x_pixels*pixelSize
      y_total = y_pixels*pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, (math.ceil(int(y_pixels/2))+0.5)*pixelSize, x_total, (math.ceil(int(y_pixels/2))+0.5)*pixelSize, \
                 (int(x_total/pixelSize)/2)*pixelSize, y_total, 0, 0]
    else:
      print('For a 3-port network, the number of sides must be equal to either 2 or 3.')
      quit()
  elif ports == 4:
    if sides == 2:
      x_total = x_pixels*pixelSize
      y_total = y_pixels*pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, (math.ceil(int(y_pixels/4))+0.5)*pixelSize, 0, (math.ceil(int(3*y_pixels/4))+0.5)*pixelSize, \
                 x_total, (math.ceil(int(y_pixels/4))+0.5)*pixelSize, x_total, (math.ceil(int(3*y_pixels/4))\
                 +0.5)*pixelSize]
    elif sides == 4:
      x_total = x_pixels*pixelSize
      y_total = y_pixels*pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, (math.ceil(int((y_total/pixelSize)/2))+0.5)*pixelSize, \
                 (int(x_total/pixelSize)/2)*pixelSize, y_total, \
                 x_total, (math.ceil(int((y_total/pixelSize)/2))+0.5)*pixelSize, \
                 (int(x_total/pixelSize)/2)*pixelSize, 0]
    else:
      print('For a 4-port network, the number of sides must be equal to either 2 or 4.')
      quit()

  # Draw outline
  pixel_map = np.zeros((int(x_total/pixelSize),int(y_total/pixelSize)),dtype=int)
  outline = gdspy.Rectangle((0, 0), (x_total, y_total), **l_bottom)
  UNIT.add(outline) 
  
  if ports == 2:
    # Add launches: assume a rectangle with port 1 = West, 2 = East
    for x in range(launch_l_pixels):
      for y in range(launch_w_pixels):
        launch_1 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize, \
                   (y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2))*pixelSize)
        pixel_map[x,y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)] = 1
        launch_2 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x_total - \
                   length_launch + x*pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                   math.floor(launch_w_pixels/2))*pixelSize)
        pixel_map[int(x_pixels) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                  math.floor(launch_w_pixels/2)] = 1
        UNIT.add(launch_1)
        UNIT.add(launch_2)
      y += 1
    x += 1
  elif ports ==3:
    if sides == 2:
      # Add launches: assume a rectangle with port 1 = Southwest, 2 = Northwest, 3 = East
      for x in range(launch_l_pixels):
        for y in range(launch_w_pixels):
          launch_1 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize, \
                     (y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[x,y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
          launch_2 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize, \
                     (y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[x,y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
          launch_3 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x_total - \
                     length_launch + x*pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                     math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[int(x_total/pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                    math.floor(launch_w_pixels/2)] = 1
          UNIT.add(launch_1)
          UNIT.add(launch_2)
          UNIT.add(launch_3)
        y += 1
      x += 1
    elif sides == 3:
      # Add launches: assume a rectangle with port 1 = west, 2 = east, 3 = north
      for x in range(launch_l_pixels):
        for y in range(launch_w_pixels):
          launch_1 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize, \
                     (y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[x,y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)] = 1
          launch_2 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x_total - \
                     length_launch + x*pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                     math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[int(x_total/pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                    math.floor(launch_w_pixels/2)] = 1
          launch_3 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(\
                     (y+math.ceil(int(x_total / pixelSize/2))-math.floor(launch_w_pixels/2)) * \
                     pixelSize, y_total - length_launch + x*pixelSize)
          pixel_map[y+math.ceil(int(x_total / pixelSize/2))-math.floor(launch_w_pixels/2), \
                    int(y_total/pixelSize)-launch_l_pixels + x] = 1
          UNIT.add(launch_1)
          UNIT.add(launch_2)
          UNIT.add(launch_3)
        y += 1
      x += 1
  elif ports == 4:
    if sides == 2:
      # Add launches: assume a rectangle with port 1 = southwest, 2 = northwest, 3 = southeast, 4 = northeast
      for x in range(launch_l_pixels):
        for y in range(launch_w_pixels):
          launch_1 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize, \
                     (y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[x,y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
          launch_2 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize, \
                     (y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[x,y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
          launch_3 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x_total - \
                     length_launch + x*pixelSize, (y+math.ceil(int(y_pixels/4)) - \
                     math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[int(x_total/pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/4)) - \
                    math.floor(launch_w_pixels/2)] = 1
          launch_4 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x_total - \
                     length_launch + x*pixelSize, (y+math.ceil(int(3*y_pixels/4)) - \
                     math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[int(x_total/pixelSize) - launch_l_pixels + x, y+math.ceil(int(3*y_pixels/4)) - \
                    math.floor(launch_w_pixels/2)] = 1
          UNIT.add(launch_1)
          UNIT.add(launch_2)
          UNIT.add(launch_3)
          UNIT.add(launch_4)
        y += 1
      x += 1
    elif sides == 4:
      # Add launches: assume a rectangle with port 1 = west, 2 = north, 3 = east, 4 = south
      for x in range(launch_l_pixels):
        for y in range(launch_w_pixels):
          launch_1 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize, \
                     (y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)+launch_l_pixels)*pixelSize)
          pixel_map[x,y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)+launch_l_pixels] = 1
          launch_2 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(\
                     (y+math.ceil(int(x_total / pixelSize/2))-math.floor(launch_w_pixels/2)) * \
                     pixelSize, y_total - length_launch + x*pixelSize)
          pixel_map[y+math.ceil(int(x_total / pixelSize/2))-math.floor(launch_w_pixels/2), \
                    int(y_total/pixelSize)-launch_l_pixels + x] = 1
          launch_3 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x_total - \
                     length_launch + x*pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                     math.floor(launch_w_pixels/2)+launch_l_pixels)*pixelSize)
          pixel_map[int(x_total/pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                    math.floor(launch_w_pixels/2) + launch_l_pixels] = 1
          launch_4 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(\
                     (y+math.ceil(int(x_total / pixelSize/2))-math.floor(launch_w_pixels/2)) * \
                     pixelSize, x*pixelSize)
          pixel_map[y+math.ceil(int(x_total / pixelSize/2))-math.floor(launch_w_pixels/2), x] = 1
          UNIT.add(launch_1)
          UNIT.add(launch_2)
          UNIT.add(launch_3)
          UNIT.add(launch_4)
        y += 1
      x += 1

  # Add ports and sources to the gds
  if ports == 2:
    if sim == 'ADS':
      port_1 = [(length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                (length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
      poly_1 = gdspy.Polygon(port_1, **l_port1)
      UNIT.add(poly_1)
      port_2 = [(x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                (x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
      poly_2 = gdspy.Polygon(port_2, **l_port2)
      UNIT.add(poly_2)
      source = [(1, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                (1, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
      poly_3 = gdspy.Polygon(source, **l_sources)
      UNIT.add(poly_3)
    elif sim == 'EMX':
      port_1 = gdspy.Label("p1", (portPos[0], portPos[1]), "w", layer=11)
      UNIT.add(port_1)
      port_2 = gdspy.Label("p2", (portPos[2], portPos[3]), "e", layer=11)
      UNIT.add(port_2)
    else:
      print('You must choose an available simulator')
      quit()
  elif ports == 3:
    if sides == 2:
      if sim == 'ADS':
        port_1 = [(length_launch*0.5, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (length_launch*0.5, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_1 = gdspy.Polygon(port_1, **l_port1)
        UNIT.add(poly_1)
        port_2 = [(length_launch*0.5, (math.ceil(int(3*y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (length_launch*0.5, (math.ceil(int(3*y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_2 = gdspy.Polygon(port_2, **l_port2)
        UNIT.add(poly_2)
        port_3 = [(x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_3 = gdspy.Polygon(port_3, **l_port3)
        UNIT.add(poly_3)
        source = [(1, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (1, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_4 = gdspy.Polygon(source, **l_sources)
        UNIT.add(poly_4)
      elif sim == 'EMX':
        port_1 = gdspy.Label("p1", (portPos[0], portPos[1]), "w", layer=11)
        UNIT.add(port_1)
        port_2 = gdspy.Label("p2", (portPos[2], portPos[3]), "w", layer=11)
        UNIT.add(port_2)
        port_3 = gdspy.Label("p3", (portPos[4], portPos[5]), "e", layer=11)
        UNIT.add(port_3)
      else:
        print('You must choose an available simulator')
        quit()
    elif sides == 3:
      if sim == 'ADS':
        port_1 = [(length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_1 = gdspy.Polygon(port_1, **l_port1)
        UNIT.add(poly_1)
        port_2 = [(x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_2 = gdspy.Polygon(port_2, **l_port2)
        UNIT.add(poly_2)
        port_3 = [((int(x_total/pixelSize)/2-1.05*(launch_w_pixels/2))*pixelSize, y_total - length_launch*0.5), \
                  ((int(x_total/pixelSize)/2+1.05*(launch_w_pixels/2))*pixelSize, y_total - length_launch*0.5)]
        poly_3 = gdspy.Polygon(port_3, **l_port3)
        UNIT.add(poly_3)
        source = [(1, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (1, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_4 = gdspy.Polygon(source, **l_sources)
        UNIT.add(poly_4)
      elif sim == 'EMX':
        port_1 = gdspy.Label("p1", (portPos[0], portPos[1]), "w", layer=11)
        UNIT.add(port_1)
        port_2 = gdspy.Label("p2", (portPos[2], portPos[3]), "e", layer=11)
        UNIT.add(port_2)
        port_3 = gdspy.Label("p3", (portPos[4], portPos[5]), "n", layer=11)
        UNIT.add(port_3)
      else:
        print('You must choose an available simulator')
        quit()
  elif ports == 4:
    if sides == 2:
      if sim == 'ADS':
        port_1 = [(length_launch*0.5, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (length_launch*0.5, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_1 = gdspy.Polygon(port_1, **l_port1)
        UNIT.add(poly_1)
        port_2 = [(length_launch*0.5, (math.ceil(int(3*y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (length_launch*0.5, (math.ceil(int(3*y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_2 = gdspy.Polygon(port_2, **l_port2)
        UNIT.add(poly_2)
        port_3 = [(x_total - length_launch*0.5, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (x_total - length_launch*0.5, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_3 = gdspy.Polygon(port_3, **l_port3)
        UNIT.add(poly_3)
        port_4 = [(x_total - length_launch*0.5, (math.ceil(int(3*y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (x_total - length_launch*0.5, (math.ceil(int(3*y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_4 = gdspy.Polygon(port_4, **l_port4)
        UNIT.add(poly_4)
        source = [(1, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (1, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_5 = gdspy.Polygon(source, **l_sources)
        UNIT.add(poly_5)
      elif sim == 'EMX':
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
    elif sides == 4:
      if sim == 'ADS':
        port_1 = [(length_launch*0.5, (math.ceil(int((y_total/pixelSize)/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (length_launch*0.5, (math.ceil(int((y_total/pixelSize)/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_1 = gdspy.Polygon(port_1, **l_port1)
        UNIT.add(poly_1)
        port_2 = [((int(x_total/pixelSize)/2-1.05*(launch_w_pixels/2))*pixelSize, y_total - length_launch*0.5), \
                  ((int(x_total/pixelSize)/2+1.05*(launch_w_pixels/2))*pixelSize, y_total - length_launch*0.5)]
        poly_2 = gdspy.Polygon(port_2, **l_port2)
        UNIT.add(poly_2)
        port_3 = [(x_total - length_launch*0.5, (math.ceil(int((y_total/pixelSize)/2))-\
                   1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (x_total - length_launch*0.5, (math.ceil(int((y_total/pixelSize)/2))+\
                   1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_3 = gdspy.Polygon(port_3, **l_port3)
        UNIT.add(poly_3)
        port_4 = [((int(x_total/pixelSize)/2-1.05*(launch_w_pixels/2))*pixelSize, length_launch*0.5), \
                  ((int(x_total/pixelSize)/2+1.05*(launch_w_pixels/2))*pixelSize, length_launch*0.5)]
        poly_4 = gdspy.Polygon(port_4, **l_port4)
        UNIT.add(poly_4)
        source = [(1, (math.ceil(int((y_total/pixelSize)/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (1, (math.ceil(int((y_total/pixelSize)/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_5 = gdspy.Polygon(source, **l_sources)
        UNIT.add(poly_5)
      elif sim == 'EMX':
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
  for x in range(x_pixels-2*launch_l_pixels):
    for y in range(y_pixels):
      bit = random.random()
      if bit >= 0.5:
      # creates a rectangle for the "on" pixel
        if ports == 2:
          rect = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize \
                 + length_launch, y*pixelSize)
          pixel_map[x+launch_l_pixels,y] = 1;
        elif ports == 3:
          if sides == 2:
            rect = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize \
                   + length_launch, y*pixelSize)
            pixel_map[x+launch_l_pixels,y] = 1;
          elif sides == 3:
            rect = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize \
                   + length_launch, y*pixelSize)
            pixel_map[x+launch_l_pixels,y] = 1;
        elif ports == 4:
          if sides == 2:
            rect = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize \
                   + length_launch, y*pixelSize)
            pixel_map[x+launch_l_pixels,y] = 1;
          elif sides == 4:
            rect = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize \
                   + length_launch, y*pixelSize + length_launch)
            pixel_map[x+launch_l_pixels,y+launch_l_pixels] = 1;
        UNIT.add(rect)
        # pixel_map[x+launch_l_pixels,y] = 1;
      y += 1  
    x += 1

  c12, c13, c14, c23, c24, c34 = ustrip.findConnectivity(pixel_map, pixelSize, ports, portPos)
  
  connections = [c12, c13, c14, c23, c24, c34]

  if connect == connections and write == True:
    csvFile = outF + ".csv"
    # Export Pixel Map file
    np.savetxt(csvFile, np.flipud(np.transpose(pixel_map)), fmt = '%d', delimiter = ",")
    gdsFile = outF + ".gds"
    # Export GDS
    lib.write_gds(gdsFile)
  else:
    csvFile = ''
    gdsFile = ''
  # View GDS
  if view == True:
    gdspy.LayoutViewer(lib) 

  return portPos, x_total, y_total, csvFile, gdsFile, cellName

def randomGDS(sub,
              ports: int,
              sides: int,
              x_el_units: int,
              y_el_units: int,
              pixelSize: int,
              imp: float,
              eLength: float,
              seedNum: int,
              connect, 
              sim: str,
              view: bool,
              write: bool,
              outF: str):
  """
  This script will create a rectangular pixel grid that is x_el_units*eLength wide by y_el_units*eLength tall 
  according to the following rules:
  eLength = Electrical Length in Degrees
  x_el_units = Number of electrical length units wide (east-to-west)
  y_el_units = Number of electrical length units tall (north-to-south)
  pixelSize = Size of the pixel in library units (code defaults to mils, but can be modified by changing lib.unit below)
  freq = frequency at which the wavelength is defined
  imp = characteristic impedance of the launch paths from the ports
  thick = thickness of the conductor in library units (code defaults to mils, but can be modified by changing lib.units below)
  height = thickness of the substrate in library units (code defaults to mils, but can be modified by changing lib.unit below) 
           Currently this code assumes that the design is a two layer design with conductors above and below the substrate
  relPerm = relative permittivity of the substrate material
  seedNum = fixed seed number for the random number generator. This allows the same design to be reproduced each time.
  sim = The type of simulator the GDS that is produced is meant to support. There are currently two options, EMX and 
        everything else. The primary reason is the way that ports are defined in EMX, where they must be defined by a
        label in the conductor layer.
  view = Do you want to open up a view window after each gds is created for inspection. This is primarily for de-bugging as new
         features are added. Value must be boolean (either True or False)
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
  
  The output of the script is a pixel map, which is a CSV file with ones and zeros in the positions of the pixels in the grid, 
  and a GDS file that can be imported into an EM simulator. It will also output an array with port positions that can be used 
  for placing ports in an ADS simulation. 
  """
  random.seed(seedNum)
  lib = gdspy.GdsLibrary()

  # Set the database unit to 1 mil
  lib.unit = 25.4e-6

  # Create Cell obj
  cellName = 'RANDOM'
  gdspy.current_library = gdspy.GdsLibrary() # This line of code has to be here to reset the GDS library on every loop
  UNIT = lib.new_cell(cellName)

  # Create layer #'s
  if ports == 1:
    l_port1 = {"layer": 1, "datatype": 0}
  if ports == 2:
    l_port1 = {"layer": 1, "datatype": 0}
    l_port2 = {"layer": 2, "datatype": 0}
  if ports == 3:
    l_port1 = {"layer": 1, "datatype": 0}
    l_port2 = {"layer": 2, "datatype": 0}
    l_port3 = {"layer": 3, "datatype": 0}
  if ports == 4:
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
  launch_l_pixels = round(length_launch/pixelSize) # length of the line to connect to structure in number of pixels
  launch_w_pixels = round(width_launch/pixelSize) # width of the line to connect to structure in number of pixels
  length_launch = launch_l_pixels*pixelSize
  width_launch = launch_w_pixels*pixelSize

  # Set horizontal and vertical pixel limits
  y_units = y_el_units # integer units of electrical length
  x_units = x_el_units # integer units of electrical length
  y_dim = math.floor(y_units*elec_length)
  x_dim = math.floor(x_units*elec_length)
  y_pixels = math.ceil(y_dim/pixelSize)
  x_pixels = math.ceil(x_dim/pixelSize)
  if y_pixels % 2 == 0:
    y_pixels += 1
  if x_pixels % 2 == 0:
    x_pixels += 1

  if ports == 2:
    if sides == 2:
      x_total = (2*launch_l_pixels + x_pixels)*pixelSize
      y_total = y_pixels*pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, y_total/2, x_total, y_total/2, 0, 0, 0, 0] 
    else: 
      print('For a 2-port network, the number of sides must be equal to 2.')
      quit()
  elif ports == 3:
    if sides == 2:
      x_total = (2*launch_l_pixels + x_pixels)*pixelSize
      y_total = y_pixels*pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, (math.ceil(int(y_pixels/4))+0.5)*pixelSize, 0, (math.ceil(int(3*y_pixels/4))+0.5)*pixelSize, \
                 x_total, (math.ceil(int(y_pixels/2))+0.5)*pixelSize, 0, 0]
    elif sides == 3:
      x_total = (2*launch_l_pixels + x_pixels)*pixelSize
      y_total = (launch_l_pixels + y_pixels)*pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, (math.ceil(int(y_pixels/2))+0.5)*pixelSize, x_total, (math.ceil(int(y_pixels/2))+0.5)*pixelSize, \
                 (int(x_total/pixelSize)/2)*pixelSize, y_total, 0, 0]
    else:
      print('For a 3-port network, the number of sides must be equal to either 2 or 3.')
      quit()
  elif ports == 4:
    if sides == 2:
      x_total = (2*launch_l_pixels + x_pixels)*pixelSize
      y_total = y_pixels*pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, (math.ceil(int(y_pixels/4))+0.5)*pixelSize, 0, (math.ceil(int(3*y_pixels/4))+0.5)*pixelSize, \
                 x_total, (math.ceil(int(y_pixels/4))+0.5)*pixelSize, x_total, (math.ceil(int(3*y_pixels/4))\
                 +0.5)*pixelSize]
    elif sides == 4:
      x_total = (2*launch_l_pixels + x_pixels)*pixelSize
      y_total = (2*launch_l_pixels + y_pixels)*pixelSize
      #Define the position of the ports for BEM simulators. This is a one column array with:
      #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
      portPos = [0, (math.ceil(int((y_total/pixelSize)/2))+0.5)*pixelSize, \
                 (int(x_total/pixelSize)/2)*pixelSize, y_total, \
                 x_total, (math.ceil(int((y_total/pixelSize)/2))+0.5)*pixelSize, \
                 (int(x_total/pixelSize)/2)*pixelSize, 0]
    else:
      print('For a 4-port network, the number of sides must be equal to either 2 or 4.')
      quit()

  # Draw outline
  pixel_map = np.zeros((int(x_total/pixelSize),int(y_total/pixelSize)),dtype=int)
  outline = gdspy.Rectangle((0, 0), (x_total, y_total), **l_bottom)
  UNIT.add(outline) 
  
  if ports == 2:
    # Add launches: assume a rectangle with port 1 = West, 2 = East
    for x in range(launch_l_pixels):
      for y in range(launch_w_pixels):
        launch_1 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize, \
                   (y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2))*pixelSize)
        pixel_map[x,y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)] = 1
        launch_2 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x_total - \
                   length_launch + x*pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                   math.floor(launch_w_pixels/2))*pixelSize)
        pixel_map[int(x_pixels) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                  math.floor(launch_w_pixels/2)] = 1
        UNIT.add(launch_1)
        UNIT.add(launch_2)
      y += 1
    x += 1
  elif ports ==3:
    if sides == 2:
      # Add launches: assume a rectangle with port 1 = Southwest, 2 = Northwest, 3 = East
      for x in range(launch_l_pixels):
        for y in range(launch_w_pixels):
          launch_1 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize, \
                     (y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[x,y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
          launch_2 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize, \
                     (y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[x,y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
          launch_3 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x_total - \
                     length_launch + x*pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                     math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[int(x_total/pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                    math.floor(launch_w_pixels/2)] = 1
          UNIT.add(launch_1)
          UNIT.add(launch_2)
          UNIT.add(launch_3)
        y += 1
      x += 1
    elif sides == 3:
      # Add launches: assume a rectangle with port 1 = west, 2 = east, 3 = north
      for x in range(launch_l_pixels):
        for y in range(launch_w_pixels):
          launch_1 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize, \
                     (y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[x,y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)] = 1
          launch_2 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x_total - \
                     length_launch + x*pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                     math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[int(x_total/pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                    math.floor(launch_w_pixels/2)] = 1
          launch_3 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(\
                     (y+math.ceil(int(x_total / pixelSize/2))-math.floor(launch_w_pixels/2)) * \
                     pixelSize, y_total - length_launch + x*pixelSize)
          pixel_map[y+math.ceil(int(x_total / pixelSize/2))-math.floor(launch_w_pixels/2), \
                    int(y_total/pixelSize)-launch_l_pixels + x] = 1
          UNIT.add(launch_1)
          UNIT.add(launch_2)
          UNIT.add(launch_3)
        y += 1
      x += 1
  elif ports == 4:
    if sides == 2:
      # Add launches: assume a rectangle with port 1 = southwest, 2 = northwest, 3 = southeast, 4 = northeast
      for x in range(launch_l_pixels):
        for y in range(launch_w_pixels):
          launch_1 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize, \
                     (y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[x,y+math.ceil(int(y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
          launch_2 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize, \
                     (y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[x,y+math.ceil(int(3*y_pixels/4))-math.floor(launch_w_pixels/2)] = 1
          launch_3 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x_total - \
                     length_launch + x*pixelSize, (y+math.ceil(int(y_pixels/4)) - \
                     math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[int(x_total/pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/4)) - \
                    math.floor(launch_w_pixels/2)] = 1
          launch_4 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x_total - \
                     length_launch + x*pixelSize, (y+math.ceil(int(3*y_pixels/4)) - \
                     math.floor(launch_w_pixels/2))*pixelSize)
          pixel_map[int(x_total/pixelSize) - launch_l_pixels + x, y+math.ceil(int(3*y_pixels/4)) - \
                    math.floor(launch_w_pixels/2)] = 1
          UNIT.add(launch_1)
          UNIT.add(launch_2)
          UNIT.add(launch_3)
          UNIT.add(launch_4)
        y += 1
      x += 1
    elif sides == 4:
      # Add launches: assume a rectangle with port 1 = west, 2 = north, 3 = east, 4 = south
      for x in range(launch_l_pixels):
        for y in range(launch_w_pixels):
          launch_1 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize, \
                     (y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)+launch_l_pixels)*pixelSize)
          pixel_map[x,y+math.ceil(int(y_pixels/2))-math.floor(launch_w_pixels/2)+launch_l_pixels] = 1
          launch_2 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(\
                     (y+math.ceil(int(x_total / pixelSize/2))-math.floor(launch_w_pixels/2)) * \
                     pixelSize, y_total - length_launch + x*pixelSize)
          pixel_map[y+math.ceil(int(x_total / pixelSize/2))-math.floor(launch_w_pixels/2), \
                    int(y_total/pixelSize)-launch_l_pixels + x] = 1
          launch_3 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x_total - \
                     length_launch + x*pixelSize, (y+math.ceil(int(y_pixels/2)) - \
                     math.floor(launch_w_pixels/2)+launch_l_pixels)*pixelSize)
          pixel_map[int(x_total/pixelSize) - launch_l_pixels + x, y+math.ceil(int(y_pixels/2)) - \
                    math.floor(launch_w_pixels/2) + launch_l_pixels] = 1
          launch_4 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(\
                     (y+math.ceil(int(x_total / pixelSize/2))-math.floor(launch_w_pixels/2)) * \
                     pixelSize, x*pixelSize)
          pixel_map[y+math.ceil(int(x_total / pixelSize/2))-math.floor(launch_w_pixels/2), x] = 1
          UNIT.add(launch_1)
          UNIT.add(launch_2)
          UNIT.add(launch_3)
          UNIT.add(launch_4)
        y += 1
      x += 1

  # Add ports and sources to the gds
  if ports == 2:
    if sim == 'ADS':
      port_1 = [(length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                (length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
      poly_1 = gdspy.Polygon(port_1, **l_port1)
      UNIT.add(poly_1)
      port_2 = [(x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                (x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
      poly_2 = gdspy.Polygon(port_2, **l_port2)
      UNIT.add(poly_2)
      source = [(1, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                (1, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
      poly_3 = gdspy.Polygon(source, **l_sources)
      UNIT.add(poly_3)
    elif sim == 'EMX':
      port_1 = gdspy.Label("p1", (portPos[0], portPos[1]), "w", layer=11)
      UNIT.add(port_1)
      port_2 = gdspy.Label("p2", (portPos[2], portPos[3]), "e", layer=11)
      UNIT.add(port_2)
    else:
      print('You must choose an available simulator')
      quit()
  elif ports == 3:
    if sides == 2:
      if sim == 'ADS':
        port_1 = [(length_launch*0.5, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (length_launch*0.5, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_1 = gdspy.Polygon(port_1, **l_port1)
        UNIT.add(poly_1)
        port_2 = [(length_launch*0.5, (math.ceil(int(3*y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (length_launch*0.5, (math.ceil(int(3*y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_2 = gdspy.Polygon(port_2, **l_port2)
        UNIT.add(poly_2)
        port_3 = [(x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_3 = gdspy.Polygon(port_3, **l_port3)
        UNIT.add(poly_3)
        source = [(1, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (1, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_4 = gdspy.Polygon(source, **l_sources)
        UNIT.add(poly_4)
      elif sim == 'EMX':
        port_1 = gdspy.Label("p1", (portPos[0], portPos[1]), "w", layer=11)
        UNIT.add(port_1)
        port_2 = gdspy.Label("p2", (portPos[2], portPos[3]), "w", layer=11)
        UNIT.add(port_2)
        port_3 = gdspy.Label("p3", (portPos[4], portPos[5]), "e", layer=11)
        UNIT.add(port_3)
      else:
        print('You must choose an available simulator')
        quit()
    elif sides == 3:
      if sim == 'ADS':
        port_1 = [(length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_1 = gdspy.Polygon(port_1, **l_port1)
        UNIT.add(poly_1)
        port_2 = [(x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (x_total - length_launch*0.5, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_2 = gdspy.Polygon(port_2, **l_port2)
        UNIT.add(poly_2)
        port_3 = [((int(x_total/pixelSize)/2-1.05*(launch_w_pixels/2))*pixelSize, y_total - length_launch*0.5), \
                  ((int(x_total/pixelSize)/2+1.05*(launch_w_pixels/2))*pixelSize, y_total - length_launch*0.5)]
        poly_3 = gdspy.Polygon(port_3, **l_port3)
        UNIT.add(poly_3)
        source = [(1, (math.ceil(int(y_pixels/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (1, (math.ceil(int(y_pixels/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_4 = gdspy.Polygon(source, **l_sources)
        UNIT.add(poly_4)
      elif sim == 'EMX':
        port_1 = gdspy.Label("p1", (portPos[0], portPos[1]), "w", layer=11)
        UNIT.add(port_1)
        port_2 = gdspy.Label("p2", (portPos[2], portPos[3]), "e", layer=11)
        UNIT.add(port_2)
        port_3 = gdspy.Label("p3", (portPos[4], portPos[5]), "n", layer=11)
        UNIT.add(port_3)
      else:
        print('You must choose an available simulator')
        quit()
  elif ports == 4:
    if sides == 2:
      if sim == 'ADS':
        port_1 = [(length_launch*0.5, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (length_launch*0.5, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_1 = gdspy.Polygon(port_1, **l_port1)
        UNIT.add(poly_1)
        port_2 = [(length_launch*0.5, (math.ceil(int(3*y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (length_launch*0.5, (math.ceil(int(3*y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_2 = gdspy.Polygon(port_2, **l_port2)
        UNIT.add(poly_2)
        port_3 = [(x_total - length_launch*0.5, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (x_total - length_launch*0.5, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_3 = gdspy.Polygon(port_3, **l_port3)
        UNIT.add(poly_3)
        port_4 = [(x_total - length_launch*0.5, (math.ceil(int(3*y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (x_total - length_launch*0.5, (math.ceil(int(3*y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_4 = gdspy.Polygon(port_4, **l_port4)
        UNIT.add(poly_4)
        source = [(1, (math.ceil(int(y_pixels/4))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (1, (math.ceil(int(y_pixels/4))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_5 = gdspy.Polygon(source, **l_sources)
        UNIT.add(poly_5)
      elif sim == 'EMX':
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
    elif sides == 4:
      if sim == 'ADS':
        port_1 = [(length_launch*0.5, (math.ceil(int((y_total/pixelSize)/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (length_launch*0.5, (math.ceil(int((y_total/pixelSize)/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_1 = gdspy.Polygon(port_1, **l_port1)
        UNIT.add(poly_1)
        port_2 = [((int(x_total/pixelSize)/2-1.05*(launch_w_pixels/2))*pixelSize, y_total - length_launch*0.5), \
                  ((int(x_total/pixelSize)/2+1.05*(launch_w_pixels/2))*pixelSize, y_total - length_launch*0.5)]
        poly_2 = gdspy.Polygon(port_2, **l_port2)
        UNIT.add(poly_2)
        port_3 = [(x_total - length_launch*0.5, (math.ceil(int((y_total/pixelSize)/2))-\
                   1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (x_total - length_launch*0.5, (math.ceil(int((y_total/pixelSize)/2))+\
                   1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_3 = gdspy.Polygon(port_3, **l_port3)
        UNIT.add(poly_3)
        port_4 = [((int(x_total/pixelSize)/2-1.05*(launch_w_pixels/2))*pixelSize, length_launch*0.5), \
                  ((int(x_total/pixelSize)/2+1.05*(launch_w_pixels/2))*pixelSize, length_launch*0.5)]
        poly_4 = gdspy.Polygon(port_4, **l_port4)
        UNIT.add(poly_4)
        source = [(1, (math.ceil(int((y_total/pixelSize)/2))-1.05*(launch_w_pixels/2)+0.5)*pixelSize), \
                  (1, (math.ceil(int((y_total/pixelSize)/2))+1.05*(launch_w_pixels/2)+0.5)*pixelSize)]
        poly_5 = gdspy.Polygon(source, **l_sources)
        UNIT.add(poly_5)
      elif sim == 'EMX':
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
        if ports == 2:
          rect = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize \
                 + length_launch, y*pixelSize)
          pixel_map[x+launch_l_pixels,y] = 1;
        elif ports == 3:
          if sides == 2:
            rect = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize \
                   + length_launch, y*pixelSize)
            pixel_map[x+launch_l_pixels,y] = 1;
          elif sides == 3:
            rect = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize \
                   + length_launch, y*pixelSize)
            pixel_map[x+launch_l_pixels,y] = 1;
        elif ports == 4:
          if sides == 2:
            rect = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize \
                   + length_launch, y*pixelSize)
            pixel_map[x+launch_l_pixels,y] = 1;
          elif sides == 4:
            rect = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize \
                   + length_launch, y*pixelSize + length_launch)
            pixel_map[x+launch_l_pixels,y+launch_l_pixels] = 1;
        UNIT.add(rect)
        # pixel_map[x+launch_l_pixels,y] = 1;
      y += 1  
    x += 1

  c12, c13, c14, c23, c24, c34 = ustrip.findConnectivity(pixel_map, pixelSize, ports, portPos)
  
  connections = [c12, c13, c14, c23, c24, c34]

  if connect == connections and write == True:
    csvFile = outF + ".csv"
    # Export Pixel Map file
    np.savetxt(csvFile, np.flipud(np.transpose(pixel_map)), fmt = '%d', delimiter = ",")
    gdsFile = outF + ".gds"
    # Export GDS
    lib.write_gds(gdsFile)
  else:
    csvFile = ''
    gdsFile = ''
  # View GDS
  if view == True:
    gdspy.LayoutViewer(lib) 


  return portPos, x_total, y_total, csvFile, gdsFile, cellName

def recreateGDS(pixelSize: int,
                inF: str):

  lib = gdspy.GdsLibrary()

  # Set the database unit to 1 mil
  lib.unit = 25.4e-6

  # Create Cell obj
  cellName = 'DESIGN'
  gdspy.current_library = gdspy.GdsLibrary() # This line of code has to be here to reset the GDS library on every loop
  UNIT = lib.new_cell(cellName)

  l_bottom = {"layer": 10, "datatype": 0}
  l_top = {"layer": 11, "datatype": 0}
  l_sources = {"layer": 5, "datatype": 0}

  pixMap = np.flipud(np.loadtxt(inF, delimiter=','))
  rows = np.size(pixMap, 0)
  cols = np.size(pixMap, 1)
  print(rows, cols)

  outline = gdspy.Rectangle((0, 0), (cols*pixelSize, rows*pixelSize), **l_bottom)
  UNIT.add(outline) 

  for x in range(0, cols):
    for y in range(0, rows):
      if pixMap[y,x]  == 1:
        rect = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize,\
               y*pixelSize)
        UNIT.add(rect)

  gdsFile = inF.replace('csv', 'gds') 
  # Export GDS
  lib.write_gds(gdsFile)
  gdspy.LayoutViewer(lib) 
