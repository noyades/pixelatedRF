import csv
import math
import random
import numpy as np
import gdspy
from vt_rrfc import microstrip as ustrip 
from skimage import measure

def uStripSteppedImpFilterGDS(
              sub,
              filtType: str,
              order: int,
              w_h: int,
              w_l: int,
              Zo_h: float,
              Zo_l: float,
              pixelSize: int,
              imp: float,
              sim: str,
              view: bool,
              write: bool,
              outF: str):

  lib = gdspy.GdsLibrary()

  # Set the database unit to 1 mil
  lib.unit = 25.4e-6

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
  
  # Setup filter coefficients
  if order == 1:
    if filtType == 'butter':
      coeffs = [2.0000]
  if order == 2:
    if filtType == 'butter':
      coeffs = [1.4142, 1.4142]
  elif order == 3:
    if filtType == 'butter':
      coeffs = [1.0000, 2.0000, 1.0000]
  elif order == 4:
    if filtType == 'butter':
      coeffs = [0.7654, 1.8478, 1.8478, 0.7654]
  elif order == 5:
    if filtType == 'butter':
      coeffs = [0.6180, 1.6180, 2.0000, 1.6180, 0.6180]
  elif order == 6:
    if filtType == 'butter':
      coeffs = [0.5176, 1.4142, 1.9318, 1.9318, 1.4142, 0.5176]
  elif order == 7:
    if filtType == 'butter':
      coeffs = [0.4450, 1.2470, 1.8019, 2.0000, 1.8019, 1.2470, 0.4450]
  elif order == 8:
    if filtType == 'butter':
      coeffs = [0.3902, 1.1111, 1.6629, 1.9615, 1.9615, 1.6629, 1.1111, 0.3902]
  elif order == 9:
    if filtType == 'butter':
      coeffs = [0.3473, 1.0000, 1.5321, 1.8794, 2.0000, 1.8794, 1.5321, 1.0000, 0.3473]
  elif order == 10:
    if filtType == 'butter':
      coeffs = [0.3129, 0.9080, 1.4142, 1.7820, 1.9754, 1.9754, 1.7820, 1.4142, 0.9080, 0.3129]

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
      widthSec[x], lengthSec[x] = ustrip.microstrip_calc.synthMicrostrip(sub, Zo_l, filt)
    else:
      widthSec[x], lengthSec[x] = ustrip.microstrip_calc.synthMicrostrip(sub, Zo_h, filt)
    y_pixels[x] = np.around(widthSec[x]/pixelSize)
    x_pixels[x] = np.around(lengthSec[x]/pixelSize)
    x += 1

  # Pixel Size
  width_launch, length_launch = ustrip.microstrip_calc.synthMicrostrip(sub, imp, 30)
  launch_l_pixels = round(length_launch/pixelSize) # length of the line to connect to structure in number of pixels
  launch_w_pixels = round(width_launch/pixelSize) # length of the line to connect to structure in number of pixels
  length_launch = launch_l_pixels*pixelSize

  # Set horizontal and vertical pixel limits
  if (2*max(y_pixels) + 1) < launch_w_pixels:
    y_dim = 3*launch_w_pixels
  else:
    y_dim = (2*max(y_pixels) + 1)
  x_dim = (2*launch_l_pixels + sum(x_pixels))
  x_total = int(x_dim*pixelSize)
  y_total = int(y_dim*pixelSize)
  portPos = [0, y_total/2, x_total, y_total/2, 0, 0, 0, 0]
  
  # Draw outline
  pixel_map = np.zeros((int(x_dim),int(y_dim)),dtype=int)
  outline = gdspy.Rectangle((0, 0), (x_total, y_total), **l_bottom)
  UNIT.add(outline) 

  # Add ports and sources to the gds
  if sim == 'ADS':
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
  elif sim == 'EMX':
    port_1 = gdspy.Label("p1", (0, y_total/2), "w", layer=11)
    UNIT.add(port_1)
    port_2 = gdspy.Label("p2", (x_total, y_total/2), "e", layer=11)
    UNIT.add(port_2)
  else:
    print('You must choose an available simulator')
    quit()

  # Add launches: assume a rectangle with port 1 = west, 2 = east
  for x in range(launch_l_pixels):
    for y in range(launch_w_pixels):
      launch_1 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x*pixelSize, \
               (y+math.ceil(int(y_dim/2))-math.floor(launch_w_pixels/2))*pixelSize)
      pixel_map[x,y+math.ceil(int(y_dim/2))-math.floor(launch_w_pixels/2)] = 1
      launch_2 = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(x_total - \
                 length_launch + x*pixelSize, (y+math.ceil(int(y_dim/2)) - \
                 math.floor(launch_w_pixels/2))*pixelSize)
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
        rect = gdspy.Rectangle((0, 0), (pixelSize, pixelSize), **l_top).translate(length_launch + \
               x*pixelSize + int(sum(x_pixels[0:z])*pixelSize), (y+math.ceil(int(y_dim/2)) - \
               math.floor(int(y_pixels[z])/2))*pixelSize)
        pixel_map[launch_l_pixels + x + int(sum(x_pixels[0:z])),y+math.ceil(int(y_dim/2)) - \
               math.floor(int(y_pixels[z])/2)] = 1
        UNIT.add(rect)
      y += 1
    x += 1
  z += 1

  if write == True:
    csvFile = outF + ".csv"
    # Export Pixel Map file
    np.savetxt(csvFile, np.transpose(pixel_map), fmt = '%d', delimiter = ",")
    gdsFile = outF + ".gds"
    # Export GDS
    lib.write_gds(gdsFile)
  else:
    csvFile = ''
    gdsFile = ''
  if view == True:
    gdspy.LayoutViewer(lib) 
  
  return portPos, x_total, y_total, csvFile, gdsFile, cellName
