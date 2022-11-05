import csv
import math
import random
import numpy as np
import gdspy
import pya
from vt_rrfc import * 
from skimage import measure

class rfc:
  def __init__(rfc_param, 
               unit,
               pixelSize: int,
               sim: str, 
               view: bool, 
               write: bool, 
               outF: str):
    """ define the microstrip substrate
    Args:
    unit = grid unit of the layout. Set to 1e-6 for um or 25.4e-6 for mil for example
    pixelSize = size of the randomized pixel
    seed = random seed number
    sim = boolean flag on what simulator to assume (this is mainly for port placement and setup)
    view = boolean flag on whether to view gds (do not use in sim mode)
    write = boolean flag on whether to write files
    outF = Output file string
    """    
    rfc_param.unit = unit
    rfc_param.pixelSize = pixelSize
    rfc_param.write = write
    rfc_param.sim = sim
    rfc_param.view = view
    rfc_param.outF = outF

def uStripSteppedImpFilterGDS(
              sub,
              rfc,
              filtType: str,
              order: int,
              w_h: int,
              w_l: int,
              Zo_h: float,
              Zo_l: float,
              imp: float):

  # Setup filter coefficients
  if order == 1:
    if filtType == 'butter':
      coeffs = [2.0000]
    elif filtType == 'chebyp5':
      coeffs = [0.6986]
    elif filtType == 'cheby3':
      coeffs = [1.9953]
    elif filtType == 'bessel':
      coeffs = [2.0000]
    else:
      print('The only valid filter types are butter, chebyp5, cheby3 or bessel')
  elif order == 2:
    if filtType == 'butter':
      coeffs = [1.4142, 1.4142]
    elif filtType == 'chebyp5':
      coeffs = [1.4029, 0.7071]
    elif filtType == 'cheby3':
      coeffs = [3.1013, 0.5339]
    elif filtType == 'bessel':
      coeffs = [1.5774, 0.4226]
    else:
      print('The only valid filter types are butter, chebyp5, cheby3 or bessel')
  elif order == 3:
    if filtType == 'butter':
      coeffs = [1.0000, 2.0000, 1.0000]
    elif filtType == 'chebyp5':
      coeffs = [1.5963, 1.0967, 1.5963]
    elif filtType == 'cheby3':
      coeffs = [3.3487, 0.7117, 3.3487]
    elif filtType == 'bessel':
      coeffs = [1.2550, 0.5528, 0.1922]
    else:
      print('The only valid filter types are butter, chebyp5, cheby3 or bessel')
  elif order == 4:
    if filtType == 'butter':
      coeffs = [0.7654, 1.8478, 1.8478, 0.7654]
    elif filtType == 'chebyp5':
      coeffs = [1.6703, 1.1926, 2.3661, 0.8419]
    elif filtType == 'cheby3':
      coeffs = [3.4389, 0.7483, 4.3471, 0.5920]
    elif filtType == 'bessel':
      coeffs = [1.0598, 0.5116, 0.3181, 0.1104]
    else:
      print('The only valid filter types are butter, chebyp5, cheby3 or bessel')
  elif order == 5:
    if filtType == 'butter':
      coeffs = [0.6180, 1.6180, 2.0000, 1.6180, 0.6180]
    elif filtType == 'chebyp5':
      coeffs = [1.7058, 1.2296, 2.5408, 1.2296, 1.7058]
    elif filtType == 'cheby3':
      coeffs = [3.4817, 0.7618, 4.5381, 0.7618, 3.4817]
    elif filtType == 'bessel':
      coeffs = [0.9303, 0.4577, 0.3312, 0.2090, 0.0718]
    else:
      print('The only valid filter types are butter, chebyp5, cheby3 or bessel')
  elif order == 6:
    if filtType == 'butter':
      coeffs = [0.5176, 1.4142, 1.9318, 1.9318, 1.4142, 0.5176]
    elif filtType == 'chebyp5':
      coeffs = [1.7254, 1.2479, 2.6064, 1.3137, 2.4758, 0.8696]
    elif filtType == 'cheby3':
      coeffs = [3.5045, 0.7685, 4.6061, 0.7929, 4.4641, 0.6033]
    elif filtType == 'bessel':
      coeffs = [0.8377, 0.4116, 0.3158, 0.2364, 0.1480, 0.0505]
    else:
      print('The only valid filter types are butter, chebyp5, cheby3 or bessel')
  elif order == 7:
    if filtType == 'butter':
      coeffs = [0.4450, 1.2470, 1.8019, 2.0000, 1.8019, 1.2470, 0.4450]
    elif filtType == 'chebyp5':
      coeffs = [1.7372, 1.2583, 2.6381, 1.3444, 2.6381, 1.2583, 1.7372]
    elif filtType == 'cheby3':
      coeffs = [3.5182, 0.7723, 4.6386, 0.8039, 4.6386, 0.7723, 3.5182]
    elif filtType == 'bessel':
      coeffs = [0.7677, 0.3744, 0.2944, 0.2378, 0.1778, 0.1104, 0.0375]
    else:
      print('The only valid filter types are butter, chebyp5, cheby3 or bessel')
  elif order == 8:
    if filtType == 'butter':
      coeffs = [0.3902, 1.1111, 1.6629, 1.9615, 1.9615, 1.6629, 1.1111, 0.3902]
    elif filtType == 'chebyp5':
      coeffs = [1.7451, 1.2647, 2.6564, 1.3590, 2.6964, 1.3389, 2.5093, 0.8796]
    elif filtType == 'cheby3':
      coeffs = [3.5277, 0.7745, 4.6575, 0.8089, 4.6990, 0.8018, 4.4990, 0.6073]
    elif filtType == 'bessel':
      coeffs = [0.7125, 0.3446, 0.2735, 0.2297, 0.1867, 0.1387, 0.0855, 0.0289]
    else:
      print('The only valid filter types are butter, chebyp5, cheby3 or bessel')
  elif order == 9:
    if filtType == 'butter':
      coeffs = [0.3473, 1.0000, 1.5321, 1.8794, 2.0000, 1.8794, 1.5321, 1.0000, 0.3473]
    elif filtType == 'chebyp5':
      coeffs = [1.7504, 1.2690, 2.6678, 1.3673, 2.7239, 1.3673, 2.6678, 1.2690, 1.7504]
    elif filtType == 'cheby3':
      coeffs = [3.5340, 0.7760, 4.6692, 0.8118, 4.7272, 0.8118, 4.6692, 0.7760, 3.5340]
    elif filtType == 'bessel':
      coeffs = [0.6678, 0.3203, 0.2547, 0.2184, 0.1859, 0.1506, 0.1111, 0.0682, 0.0230]
    else:
      print('The only valid filter types are butter, chebyp5, cheby3 or bessel')
  elif order == 10:
    if filtType == 'butter':
      coeffs = [0.3129, 0.9080, 1.4142, 1.7820, 1.9754, 1.9754, 1.7820, 1.4142, 0.9080, 0.3129]
    elif filtType == 'chebyp5':
      coeffs = [1.7543, 1.2721, 2.6754, 1.3725, 2.7392, 1.3806, 2.7231, 1.3485, 2.5239, 0.8842]
    elif filtType == 'cheby3':
      coeffs = [3.5384, 0.7771, 4.6768, 0.8136, 4.7425, 0.8164, 4.7260, 0.8051, 4.5142, 0.6091]
    elif filtType == 'bessel':
      coeffs = [0.6305, 0.3002, 0.2384, 0.2066, 0.1808, 0.1539, 0.1240, 0.0911, 0.0557, 0.0187]
    else:
      print('The only valid filter types are butter, chebyp5, cheby3 or bessel')
  else:
    print('The only valid filter orders are 1-10')

  # design filter sections
  x = 0
  filtElecLength = np.zeros((len(coeffs),1),dtype=float)
  for coeff in coeffs:
    if x % 2 == 0:
      filtElecLength[x] = (180/math.pi) * coeff * Zo_l/imp
    else:
      filtElecLength[x] = (180/math.pi) * coeff * imp/Zo_h
    x += 1

  # design the dimensions of the filter sections
  x = 0
  widthSec = np.zeros((len(coeffs),1),dtype=float)
  lengthSec = np.zeros((len(coeffs),1),dtype=float)
  x_pixels = np.zeros((len(coeffs),1),dtype=float)
  y_pixels = np.zeros((len(coeffs),1),dtype=float)
  for filt in filtElecLength:
    if x % 2 == 0:
      widthSec[x], lengthSec[x] = microstrip_calc.synthMicrostrip(sub, Zo_l, filt)
    else:
      widthSec[x], lengthSec[x] = microstrip_calc.synthMicrostrip(sub, Zo_h, filt)
    y_pixels[x] = np.around(widthSec[x]/rfc.pixelSize)
    x_pixels[x] = np.around(lengthSec[x]/rfc.pixelSize)
    x += 1

  # Pixel Size
  width_launch, length_launch = microstrip_calc.synthMicrostrip(sub, imp, 30)
  launch_l_pixels = round(length_launch/rfc.pixelSize) # length of the line to connect to structure in number of pixels
  launch_w_pixels = round(width_launch/rfc.pixelSize) # length of the line to connect to structure in number of pixels
  length_launch = launch_l_pixels*rfc.pixelSize

  # Set horizontal and vertical pixel limits
  # The conditional below pads the structure, but prevents a fair geometric comparison hence, it is being deprecated
  """
  if (2*max(y_pixels) + 1) < launch_w_pixels:
    y_dim = 3*launch_w_pixels
  else:
    y_dim = (2*max(y_pixels) + 1)
  """
  y_dim = max(y_pixels)
  x_dim = (2*launch_l_pixels + sum(x_pixels))
  x_total = int(x_dim*rfc.pixelSize)
  y_total = int(y_dim*rfc.pixelSize)
  portPos = [0, y_total/2, x_total, y_total/2, 0, 0, 0, 0]

  if rfc.sim == 'EMX':
    ly = pya.Layout()

    # Set the database unit to 1 mil
    ly.dbu = rfc.unit*1e6

    # Create Cell obj
    cellName = 'uStripFilter'
    UNIT = ly.create_cell(cellName)

    # Create layer #'s
    l_top = ly.layer(69, 0) # layer for metal
    l_bottom = ly.layer(6969, 0) # Ground Layer

    # Draw outline
    pixel_map = np.zeros((int(x_dim),int(y_dim)),dtype=int)
    outline = UNIT.shapes(l_bottom).insert( pya.Box(0, 0, x_total, y_total) ) 

    # Add port labels
    port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', 0, y_total/2) )
    port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', x_total, y_total/2) )

    # Add launches: assume a rectangle with port 1 = west, 2 = east
    for x in range(launch_l_pixels):
      for y in range(launch_w_pixels):
        launch_1 = pya.Box(0, 0, rfc.pixelSize, rfc.pixelSize).moved(x*rfc.pixelSize, \
                 (y+math.ceil(int(y_dim/2))-math.floor(launch_w_pixels/2))*rfc.pixelSize)
        pixel_map[x,y+math.ceil(int(y_dim/2))-math.floor(launch_w_pixels/2)] = 1
        launch_2 = pya.Box(0, 0, rfc.pixelSize, rfc.pixelSize).moved(x_total - \
                   length_launch + x*rfc.pixelSize, (y+math.ceil(int(y_dim/2)) - \
                   math.floor(launch_w_pixels/2))*rfc.pixelSize)
        pixel_map[int(x_dim) - launch_l_pixels + x, y+math.ceil(int(y_dim/2)) - \
                  math.floor(launch_w_pixels/2)] = 1
        UNIT.shapes(l_top).insert(launch_1)
        UNIT.shapes(l_top).insert(launch_2)
      y += 1
    x += 1

    # Add the filter elements
    for z in range(len(coeffs)):
      for x in range(int(x_pixels[z])):
        for y in range(int(y_pixels[z])):
          rect = pya.Box(0, 0, rfc.pixelSize, rfc.pixelSize).moved(length_launch + \
                 x*rfc.pixelSize + int(sum(x_pixels[0:z])*rfc.pixelSize), (y+math.ceil(int(y_dim/2)) - \
                 math.floor(int(y_pixels[z])/2))*rfc.pixelSize)
          pixel_map[launch_l_pixels + x + int(sum(x_pixels[0:z])),y+math.ceil(int(y_dim/2)) - \
                 math.floor(int(y_pixels[z])/2)] = 1
          UNIT.shapes(l_top).insert(rect)
        y += 1
      x += 1
    z += 1

    """
    if rfc.corner == 'overlap':
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
    """

    if rfc.write == True:
      csvFile = rfc.outF + ".csv"
      # Export Pixel Map file
      np.savetxt(csvFile, np.transpose(pixel_map), fmt = '%d', delimiter = ",")
      gdsFile = rfc.outF + ".gds"
      # Export GDS
      ly.write(gdsFile)
    else:
      csvFile = ''
      gdsFile = ''
    if rfc.view == True:
      # Load a GDSII file into a new library
      gdsii = gdspy.GdsLibrary(infile=gdsFile)
      gdspy.LayoutViewer(gdsii)

  else:
    lib = gdspy.GdsLibrary()

    # Set the database unit to 1 mil
    lib.unit = rfc.unit

    # Create Cell obj
    cellName = 'uStripFilter'
    gdspy.current_library = gdspy.GdsLibrary() # This line of code has to be here to reset the GDS library on every loop
    UNIT = lib.new_cell(cellName)
  
    # Create layer #'s
    l_port1 = {"layer": 1, "datatype": 0}
    l_port2 = {"layer": 2, "datatype": 0}
    l_bottom = {"layer": 10, "datatype": 0}
    l_top = {"layer": 11, "datatype": 0}
    l_sources = {"layer": 5, "datatype": 0}
  
    # Draw outline
    pixel_map = np.zeros((int(x_dim),int(y_dim)),dtype=int)
    outline = gdspy.Rectangle((0, 0), (x_total, y_total), **l_bottom)
    UNIT.add(outline) 

    # Add ports and sources to the gds
    if rfc.sim == 'MEEP':
      port_1 = [(length_launch*0.5, y_total/2 + 1.05*width_launch/2), \
                (length_launch*0.5, y_total/2 - 1.05*width_launch/2)]
      poly_1 = gdspy.Polygon(port_1, **l_port1)
      UNIT.add(poly_1)
      port_2 = [(x_total - length_launch*0.5, y_total/2 + 1.05*width_launch/2), \
                (x_total - length_launch*0.5, y_total/2 - 1.05*width_launch/2)]
      poly_2 = gdspy.Polygon(port_2, **l_port2)
      UNIT.add(poly_2)
      source = [(1, y_total/2 + 1.05*width_launch/2), (1, y_total/2 - 1.05*width_launch/2)]
      poly_3 = gdspy.Polygon(source, **l_sources)
      UNIT.add(poly_3)
    else:
      t = 1 #print('If you are using ADS, you do not need to add ports. They are added with AEL')


    # Add launches: assume a rectangle with port 1 = west, 2 = east
    for x in range(launch_l_pixels):
      for y in range(launch_w_pixels):
        launch_1 = gdspy.Rectangle((0, 0), (rfc.pixelSize, rfc.pixelSize), **l_top).translate(x*rfc.pixelSize, \
                   (y+math.ceil(int(y_dim/2))-math.floor(launch_w_pixels/2))*rfc.pixelSize)
        pixel_map[x,y+math.ceil(int(y_dim/2))-math.floor(launch_w_pixels/2)] = 1
        launch_2 = gdspy.Rectangle((0, 0), (rfc.pixelSize, rfc.pixelSize), **l_top).translate(x_total - \
                   length_launch + x*rfc.pixelSize, (y+math.ceil(int(y_dim/2)) - \
                   math.floor(launch_w_pixels/2))*rfc.pixelSize)
        pixel_map[int(x_dim) - launch_l_pixels + x, y+math.ceil(int(y_dim/2)) - \
                  math.floor(launch_w_pixels/2)] = 1
        UNIT.add(launch_1)
        UNIT.add(launch_2)
      y += 1
    x += 1

    # Add the filter elements
    for z in range(len(coeffs)):
      for x in range(int(x_pixels[z])):
        for y in range(int(y_pixels[z])):
          rect = gdspy.Rectangle((0, 0), (rfc.pixelSize, rfc.pixelSize), **l_top).translate(length_launch + \
                 x*rfc.pixelSize + int(sum(x_pixels[0:z])*rfc.pixelSize), (y+math.ceil(int(y_dim/2)) - \
                 math.floor(int(y_pixels[z])/2))*rfc.pixelSize)
          pixel_map[launch_l_pixels + x + int(sum(x_pixels[0:z])),y+math.ceil(int(y_dim/2)) - \
                 math.floor(int(y_pixels[z])/2)] = 1
          UNIT.add(rect)
        y += 1
      x += 1
    z += 1

    """
    if rfc.corner == 'overlap':
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
    """

    if rfc.write == True:
      csvFile = rfc.outF + ".csv"
      # Export Pixel Map file
      np.savetxt(csvFile, np.transpose(pixel_map), fmt = '%d', delimiter = ",")
      gdsFile = rfc.outF + ".gds"
      # Export GDS
      lib.write_gds(gdsFile)
    else:
      csvFile = ''
      gdsFile = ''
    if rfc.view == True:
      gdspy.LayoutViewer(lib) 
  
  return portPos, x_total, y_total, csvFile, gdsFile, cellName, launch_l_pixels

def uStripTeeJunctionGDS(
              sub,
              rfc,
              x_pixels: float,
              y_pixels: float,
              imp: float):

  # Pixel Size
  width_launch, length_launch = microstrip_calc.synthMicrostrip(sub, imp, 30)
  launch_l_pixels = round(length_launch/rfc.pixelSize) # length of the line to connect to structure in number of pixels
  launch_w_pixels = round(width_launch/rfc.pixelSize) # length of the line to connect to structure in number of pixels
  length_launch = launch_l_pixels*rfc.pixelSize

  # Set horizontal and vertical pixel limits
  # The conditional below pads the structure, but prevents a fair geometric comparison hence, it is being deprecated
  y_dim = y_pixels + launch_l_pixels
  x_dim = (2*launch_l_pixels + x_pixels)
  x_total = int(x_dim*rfc.pixelSize)
  y_total = int(y_dim*rfc.pixelSize)
  portPos = [0, (y_total-length_launch)/2, x_total, (y_total-length_launch)/2, x_total/2, y_total, 0, 0]

  if rfc.sim == 'EMX':
    ly = pya.Layout()

    # Set the database unit to 1 mil
    ly.dbu = rfc.unit*1e6

    # Create Cell obj
    cellName = 'uStripTee'
    UNIT = ly.create_cell(cellName)

    # Create layer #'s
    l_top = ly.layer(69, 0) # layer for metal
    l_bottom = ly.layer(6969, 0) # Ground Layer

    # Draw outline
    pixel_map = np.zeros((int(x_dim),int(y_dim)),dtype=int)
    outline = UNIT.shapes(l_bottom).insert( pya.Box(0, 0, x_total, y_total) ) 

    # Add port labels
    port_1 = UNIT.shapes(l_top).insert( pya.Text('p1', 0, (y_total-length_launch)/2) )
    port_2 = UNIT.shapes(l_top).insert( pya.Text('p2', x_total, (y_total-length_launch)/2) )
    port_3 = UNIT.shapes(l_top).insert( pya.Text('p3', x_total/2, y_total) )

    # Add the filter elements
    widthSec, lengthSec = microstrip_calc.synthMicrostrip(sub, imp, 90)
    tee_w_pixels = round(widthSec/rfc.pixelSize)
    tee_l_pixels = int(lengthSec/rfc.pixelSize)
    for x in range(tee_l_pixels):
      for y in range(tee_w_pixels):
        rect1 = pya.Box(0, 0, rfc.pixelSize, rfc.pixelSize).moved(length_launch + \
               x*rfc.pixelSize, (y+math.ceil(int((y_dim-launch_l_pixels)/2))-math.floor(launch_w_pixels/2))*rfc.pixelSize)
        pixel_map[launch_l_pixels + x,(y+math.ceil(int((y_dim-launch_l_pixels)/2))-math.floor(launch_w_pixels/2))] = 1
        UNIT.shapes(l_top).insert(rect1)
      y += 1
    x += 1

    for x in range(tee_w_pixels):
      for y in range(int((y_pixels - tee_w_pixels)/2)):
        rect2 = pya.Box(0, 0, rfc.pixelSize, rfc.pixelSize).moved((x+math.ceil(int(x_dim/2)) - \
                   math.floor(launch_w_pixels/2))*rfc.pixelSize, (y+math.ceil(int((y_dim-launch_l_pixels)/2)) + \
                   math.floor(launch_w_pixels/2)+1)*rfc.pixelSize)
        pixel_map[(x+math.ceil(int(x_dim/2)) - math.floor(launch_w_pixels/2)),(y+math.ceil(int((y_dim-launch_l_pixels)/2)) + \
                   math.floor(launch_w_pixels/2)+1)] = 1
        UNIT.shapes(l_top).insert(rect2)
      y += 1
    x += 1

    # Add launches: assume a rectangle with port 1 = west, 2 = east, 3 = north
    for x in range(launch_l_pixels):
      for y in range(launch_w_pixels):
        launch_1 = pya.Box(0, 0, rfc.pixelSize, rfc.pixelSize).moved(x*rfc.pixelSize, \
                 (y+math.ceil(int((y_dim-launch_l_pixels)/2))-math.floor(launch_w_pixels/2))*rfc.pixelSize)
        pixel_map[x,y+math.ceil(int((y_dim-launch_l_pixels)/2))-math.floor(launch_w_pixels/2)] = 1
        launch_2 = pya.Box(0, 0, rfc.pixelSize, rfc.pixelSize).moved(x_total - \
                   length_launch + x*rfc.pixelSize, (y+math.ceil(int((y_dim-launch_l_pixels)/2)) - \
                   math.floor(launch_w_pixels/2))*rfc.pixelSize)
        pixel_map[int(x_dim) - launch_l_pixels + x, y+math.ceil(int((y_dim-launch_l_pixels)/2)) - \
                  math.floor(launch_w_pixels/2)] = 1
        launch_3 = pya.Box(0, 0, rfc.pixelSize, rfc.pixelSize).moved((y+math.ceil(int(x_dim/2)) - \
                   math.floor(launch_w_pixels/2))*rfc.pixelSize, y_total - \
                   length_launch + x*rfc.pixelSize)
        pixel_map[(y+math.ceil(int(x_dim/2)) - \
                   math.floor(launch_w_pixels/2)), y_dim - launch_l_pixels + x] = 1
        UNIT.shapes(l_top).insert(launch_1)
        UNIT.shapes(l_top).insert(launch_2)
        UNIT.shapes(l_top).insert(launch_3)
      y += 1
    x += 1

    print(rfc.write)
    if rfc.write == True:
      csvFile = rfc.outF + ".csv"

      # Export Pixel Map file
      np.savetxt(csvFile, np.flipud(np.transpose(pixel_map)), fmt = '%d', delimiter = ",")
      gdsFile = rfc.outF + ".gds"
      # Export GDS
      ly.write(gdsFile)
    else:
      csvFile = ''
      gdsFile = ''
    if rfc.view == True:
      # Load a GDSII file into a new library
      gdsii = gdspy.GdsLibrary(infile=gdsFile)
      gdspy.LayoutViewer(gdsii)

  else:
    lib = gdspy.GdsLibrary()

    # Set the database unit to 1 mil
    lib.unit = rfc.unit

    # Create Cell obj
    cellName = 'uStripTee'
    gdspy.current_library = gdspy.GdsLibrary() # This line of code has to be here to reset the GDS library on every loop
    UNIT = lib.new_cell(cellName)
  
    # Create layer #'s
    l_port1 = {"layer": 1, "datatype": 0}
    l_port2 = {"layer": 2, "datatype": 0}
    l_port3 = {"layer": 3, "datatype": 0}
    l_bottom = {"layer": 10, "datatype": 0}
    l_top = {"layer": 11, "datatype": 0}
    l_sources = {"layer": 5, "datatype": 0}
  
    # Draw outline
    pixel_map = np.zeros((int(x_dim),int(y_dim)),dtype=int)
    outline = gdspy.Rectangle((0, 0), (x_total, y_total), **l_bottom)
    UNIT.add(outline)

    # Add ports and sources to the gds
    if rfc.sim == 'MEEP':
      port_1 = [(length_launch*0.5, (y_total - length_launch)/2 + 1.05*width_launch/2), \
                (length_launch*0.5, (y_total - length_launch)/2 - 1.05*width_launch/2)]
      poly_1 = gdspy.Polygon(port_1, **l_port1)
      UNIT.add(poly_1)
      port_2 = [(x_total - length_launch*0.5, (y_total - length_launch)/2 + 1.05*width_launch/2), \
                (x_total - length_launch*0.5, (y_total - length_launch)/2 - 1.05*width_launch/2)]
      poly_2 = gdspy.Polygon(port_2, **l_port2)
      UNIT.add(poly_2)
      port_3 = [(x_total/2 + 1.05*width_launch/2, y_total - length_launch*0.5), \
                (x_total/2 - 1.05*width_launch/2, y_total - length_launch*0.5)]
      poly_3 = gdspy.Polygon(port_3, **l_port3)
      UNIT.add(poly_3)
      source = [(1, (y_total - length_launch)/2 + 1.05*width_launch/2), \
                (1, (y_total - length_launch)/2 - 1.05*width_launch/2)]
      poly_4 = gdspy.Polygon(source, **l_sources)
      UNIT.add(poly_4)
    else:
      t = 1 #print('If you are using ADS, you do not need to add ports. They are added with AEL')

    # Add the filter elements
    widthSec, lengthSec = microstrip_calc.synthMicrostrip(sub, imp, 90)
    tee_w_pixels = round(widthSec/rfc.pixelSize)
    tee_l_pixels = int(lengthSec/rfc.pixelSize)
    for x in range(tee_l_pixels):
      for y in range(tee_w_pixels):
        rect = gdspy.Rectangle((0, 0), (rfc.pixelSize, rfc.pixelSize), **l_top).translate(length_launch + \
               x*rfc.pixelSize, (y+math.ceil(int((y_dim-launch_l_pixels)/2))-math.floor(launch_w_pixels/2))*rfc.pixelSize)
        pixel_map[launch_l_pixels + x,(y+math.ceil(int((y_dim-launch_l_pixels)/2))-math.floor(launch_w_pixels/2))] = 1
        UNIT.add(rect)
      y += 1
    x += 1

    for x in range(tee_w_pixels):
      for y in range(int((y_pixels - tee_w_pixels)/2)):
        rect = gdspy.Rectangle((0, 0), (rfc.pixelSize, rfc.pixelSize), **l_top).translate((x+math.ceil(int(x_dim/2)) - \
                   math.floor(launch_w_pixels/2))*rfc.pixelSize, (y+math.ceil(int((y_dim-launch_l_pixels)/2)) + \
                   math.floor(launch_w_pixels/2)+1)*rfc.pixelSize)
        pixel_map[(x+math.ceil(int(x_dim/2)) - math.floor(launch_w_pixels/2)),(y+math.ceil(int((y_dim-launch_l_pixels)/2)) + \
                   math.floor(launch_w_pixels/2)+1)] = 1
        UNIT.add(rect)
      y += 1
    x += 1

    # Add launches: assume a rectangle with port 1 = west, 2 = east
    for x in range(launch_l_pixels):
      for y in range(launch_w_pixels):
        launch_1 = gdspy.Rectangle((0, 0), (rfc.pixelSize, rfc.pixelSize), **l_top).translate(x*rfc.pixelSize, \
                 (y+math.ceil(int((y_dim-launch_l_pixels)/2))-math.floor(launch_w_pixels/2))*rfc.pixelSize)
        pixel_map[x,y+math.ceil(int((y_dim-launch_l_pixels)/2))-math.floor(launch_w_pixels/2)] = 1
        launch_2 = gdspy.Rectangle((0, 0), (rfc.pixelSize, rfc.pixelSize), **l_top).translate(x_total - \
                   length_launch + x*rfc.pixelSize, (y+math.ceil(int((y_dim-launch_l_pixels)/2)) - \
                   math.floor(launch_w_pixels/2))*rfc.pixelSize)
        pixel_map[int(x_dim) - launch_l_pixels + x, y+math.ceil(int((y_dim-launch_l_pixels)/2)) - \
                  math.floor(launch_w_pixels/2)] = 1
        launch_3 = gdspy.Rectangle((0, 0), (rfc.pixelSize, rfc.pixelSize), **l_top).translate((y+math.ceil(int(x_dim/2)) - \
                   math.floor(launch_w_pixels/2))*rfc.pixelSize, y_total - \
                   length_launch + x*rfc.pixelSize)
        pixel_map[(y+math.ceil(int(x_dim/2)) - \
                   math.floor(launch_w_pixels/2)), y_dim - launch_l_pixels + x] = 1
        UNIT.add(launch_1)
        UNIT.add(launch_2)
        UNIT.add(launch_3)
      y += 1
    x += 1

    if rfc.write == True:
      csvFile = rfc.outF + ".csv"
      # Export Pixel Map file
      np.savetxt(csvFile, np.flipud(np.transpose(pixel_map)), fmt = '%d', delimiter = ",")
      gdsFile = rfc.outF + ".gds"
      # Export GDS
      lib.write_gds(gdsFile)
    else:
      csvFile = ''
      gdsFile = ''
    if rfc.view == True:
      gdspy.LayoutViewer(lib) 

  return portPos, x_total, y_total, csvFile, gdsFile, cellName, launch_l_pixels
