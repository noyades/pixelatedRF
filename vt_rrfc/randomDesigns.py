import math
import random
import numpy as np
import gdspy
import pya

class randomDesign:
  def __init__(self, 
               unit,
               x_pixels: int,
               y_pixels: int,
               ports: int, 
               sides: int, 
               launch_pixels_l: int,
               launch_pixels_w: int,
               pixelSize: int,
               sim: str,
               sym: str,
               seed: int):
               #connect,
               #seed: int, 
               #sim: str, 
               #view: bool, 
               #write: bool, 
               #outF: str,
               #sym: str):
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
    self.unit = unit
    self.x_pixels = x_pixels
    self.y_pixels = y_pixels
    self.ports = ports
    self.sides = sides
    self.launch_pixels_l = launch_pixels_l
    self.launch_pixels_w = launch_pixels_w
    self.sym = sym
    self.pixelSize = pixelSize
    #self.seed = seed
    #self.write = write
    self.sim = sim
    #self.view = view
    #self.outF = outF
    #self.connect = connect
    random.seed(seed)

  def genPortPos(self):
    if self.ports == 2:
      if self.sides == 2:
        x_total = (self.x_pixels+2*self.launch_pixels_l)*self.pixelSize
        y_total = self.y_pixels*self.pixelSize
        #Define the position of the ports for BEM simulators. This is a one column array with:
        #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
        portPos = [0, y_total/2, x_total, y_total/2, 0, 0, 0, 0] 
      else: 
        print('For a 2-port network, the number of sides must be equal to 2.')
        quit()
    elif self.ports == 3:
      if self.sides == 2:
        x_total = (self.x_pixels+2*self.launch_pixels_l)*self.pixelSize
        y_total = self.y_pixels*self.pixelSize
        #Define the position of the ports for BEM simulators. This is a one column array with:
        #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
        portPos = [0, (math.ceil(int(self.y_pixels/4))+0.5)*self.pixelSize,\
                   0, (math.ceil(int(3*self.y_pixels/4))+0.5)*self.pixelSize, \
                   x_total, (math.ceil(int(self.y_pixels/2))+0.5)*self.pixelSize, 0, 0]
      elif self.sides == 3:
        x_total = (self.x_pixels+2*self.launch_pixels_l)*self.pixelSize
        y_total = (self.y_pixels+self.launch_pixels_l)*self.pixelSize
        #Define the position of the ports for BEM simulators. This is a one column array with:
        #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
        portPos = [0, (math.ceil(int(self.y_pixels/2))+0.5)*self.pixelSize,\
                   x_total, (math.ceil(int(self.y_pixels/2))+0.5)*self.pixelSize, \
                   (int(x_total/self.pixelSize)/2)*self.pixelSize, y_total, 0, 0]
      else:
        print('For a 3-port network, the number of sides must be equal to either 2 or 3.')
        quit()
    elif self.ports == 4:
      if self.sides == 2:
        x_total = (self.x_pixels+2*self.launch_pixels_l)*self.pixelSize
        y_total = self.y_pixels*self.pixelSize
        #Define the position of the ports for BEM simulators. This is a one column array with:
        #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
        portPos = [0, (math.ceil(int(self.y_pixels/4))+0.5)*self.pixelSize,\
                   0, (math.ceil(int(3*self.y_pixels/4))+0.5)*self.pixelSize,\
                   x_total, (math.ceil(int(self.y_pixels/4))+0.5)*self.pixelSize,\
                   x_total, (math.ceil(int(3*self.y_pixels/4))+0.5)*self.pixelSize]
      elif self.sides == 4:
        x_total = (self.x_pixels+2*self.launch_pixels_l)*self.pixelSize
        y_total = (self.y_pixels+2*self.launch_pixels_l)*self.pixelSize
        #Define the position of the ports for BEM simulators. This is a one column array with:
        #[port1x, port1y, port2x, port2y, port3x, port3y, port4x, port4y]
        portPos = [0, (math.ceil(int((y_total/self.pixelSize)/2))+0.5)*self.pixelSize, \
                   (int(x_total/self.pixelSize)/2)*self.pixelSize, y_total, \
                   x_total, (math.ceil(int((y_total/self.pixelSize)/2))+0.5)*self.pixelSize, \
                   (int(x_total/self.pixelSize)/2)*self.pixelSize, 0]
      else:
        print('For a 4-port network, the number of sides must be equal to either 2 or 4.')
        quit()

    return x_total, y_total, portPos

  def addLaunch(self,x_total,y_total,pixel_map):
    if self.ports == 2:
      # Add launches: assume a rectangle with port 1 = West, 2 = East
      for x in range(self.launch_pixels_l):
        for y in range(self.launch_pixels_w):
          pixel_map[x,y+math.ceil(int(self.y_pixels/2))-math.floor(self.launch_pixels_w/2)] = 1
          pixel_map[int(x_total/self.pixelSize) - self.launch_pixels_l + x, y+math.ceil(int(self.y_pixels/2)) - \
                    math.floor(self.launch_pixels_w/2)] = 1
        y += 1
      x += 1
    elif self.ports ==3:
      if self.sides == 2:
        # Add launches: assume a rectangle with port 1 = Southwest, 2 = Northwest, 3 = East
        for x in range(self.launch_pixels_l):
          for y in range(self.launch_pixels_w):
            pixel_map[x,y+math.ceil(int(self.y_pixels/4))-math.floor(self.launch_pixels_w/2)] = 1
            pixel_map[x,y+math.ceil(int(3*self.y_pixels/4))-math.floor(self.launch_pixels_w/2)] = 1
            pixel_map[int(x_total/self.pixelSize) - self.launch_pixels_l + x, y+math.ceil(int(self.y_pixels/2)) - \
                      math.floor(self.launch_pixels_w/2)] = 1
          y += 1
        x += 1
      elif self.sides == 3:
        # Add launches: assume a rectangle with port 1 = west, 2 = east, 3 = north
        for x in range(self.launch_pixels_l):
          for y in range(self.launch_pixels_w):
            pixel_map[x,int(y+math.ceil(int(self.y_pixels/2))-math.floor(self.launch_pixels_w/2))] = 1
            pixel_map[int(x_total/self.pixelSize) - self.launch_pixels_l + x, int(y+math.ceil(int(self.y_pixels/2)) - \
                      math.floor(self.launch_pixels_w/2))] = 1
            pixel_map[y+math.ceil(int(x_total / self.pixelSize/2))-math.floor(self.launch_pixels_w/2), \
                      int(y_total/self.pixelSize)-self.launch_pixels_l + x] = 1
          y += 1
        x += 1
    elif self.ports == 4:
      if self.sides == 2:
        # Add launches: assume a rectangle with port 1 = southwest, 2 = northwest, 3 = southeast, 4 = northeast
        for x in range(self.launch_pixels_l):
          for y in range(self.launch_pixels_w):
            pixel_map[x,y+math.ceil(int(self.y_pixels/4))-math.floor(self.launch_pixels_w/2)] = 1
            pixel_map[x,y+math.ceil(int(3*self.y_pixels/4))-math.floor(self.launch_pixels_w/2)] = 1
            pixel_map[int(x_total/self.pixelSize) - self.launch_pixels_l + x, y+math.ceil(int(self.y_pixels/4)) - \
                      math.floor(self.launch_pixels_w/2)] = 1
            pixel_map[int(x_total/self.pixelSize) - self.launch_pixels_l + x, y+math.ceil(int(3*self.y_pixels/4)) - \
                      math.floor(self.launch_pixels_w/2)] = 1
          y += 1
        x += 1
      elif self.sides == 4:
        # Add launches: assume a rectangle with port 1 = west, 2 = north, 3 = east, 4 = south
        for x in range(self.launch_pixels_l):
          for y in range(self.launch_pixels_w):
            pixel_map[x,y+math.ceil(int(self.y_pixels/2))-math.floor(self.launch_pixels_w/2)+self.launch_pixels_l] = 1
            pixel_map[y+math.ceil(int(x_total / self.pixelSize/2))-math.floor(self.launch_pixels_w/2), \
                      int(y_total/self.pixelSize)-self.launch_pixels_l + x] = 1
            pixel_map[int(x_total/self.pixelSize) - self.launch_pixels_l + x, y+math.ceil(int(self.y_pixels/2)) - \
                      math.floor(self.launch_pixels_w/2) + self.launch_pixels_l] = 1
            pixel_map[y+math.ceil(int(x_total / self.pixelSize/2))-math.floor(self.launch_pixels_w/2), x] = 1
          y += 1
        x += 1
      
    return pixel_map

  def addPortsGdspy(self,x_total,y_total):
    # Add ports and sources to the gds
    port_1 = []
    port_2 = []
    port_3 = []
    port_4 = []
    port_5 = []
    source = []
    if self.ports == 2:
      port_1 = [(self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(self.y_pixels/2))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                (self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(self.y_pixels/2))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
        #poly_1 = gdspy.Polygon(port_1, **self.l_port1)

      port_2 = [(x_total - self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(self.y_pixels/2))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                (x_total - self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(self.y_pixels/2))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
        #poly_2 = gdspy.Polygon(port_2, **self.l_port2)

      source = [(1, (math.ceil(int(self.y_pixels/2))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                (1, (math.ceil(int(self.y_pixels/2))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
        #poly_3 = gdspy.Polygon(source, **self.l_sources)
      
    elif self.ports == 3:
      if self.sides == 2:
        port_1 = [(self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(self.y_pixels/4))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                  (self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(self.y_pixels/4))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
          #poly_1 = gdspy.Polygon(port_1, **self.l_port1)

        port_2 = [(self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(3*self.y_pixels/4))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                  (self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(3*self.y_pixels/4))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
          #poly_2 = gdspy.Polygon(port_2, **self.l_port2)

        port_3 = [(x_total - self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(self.y_pixels/2))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                  (x_total - self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(self.y_pixels/2))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
          #poly_3 = gdspy.Polygon(port_3, **self.l_port3)

        source = [(1, (math.ceil(int(self.y_pixels/4))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                  (1, (math.ceil(int(self.y_pixels/4))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
          #poly_4 = gdspy.Polygon(source, **self.l_sources)
      elif self.sides == 3:
        port_1 = [(self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(self.y_pixels/2))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                  (self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(self.y_pixels/2))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
          #poly_1 = gdspy.Polygon(port_1, **self.l_port1)

        port_2 = [(x_total - self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(self.y_pixels/2))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                  (x_total - self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(self.y_pixels/2))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
          #poly_2 = gdspy.Polygon(port_2, **self.l_port2)

        port_3 = [((int(x_total/self.pixelSize)/2-1.05*(self.launch_pixels_w/2))*self.pixelSize, y_total - self.launch_pixels_l*self.pixelSize*0.5), \
                  ((int(x_total/self.pixelSize)/2+1.05*(self.launch_pixels_w/2))*self.pixelSize, y_total - self.launch_pixels_l*self.pixelSize*0.5)]
          #poly_3 = gdspy.Polygon(port_3, **self.l_port3)

        source = [(1, (math.ceil(int(self.y_pixels/2))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                  (1, (math.ceil(int(self.y_pixels/2))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
          #poly_4 = gdspy.Polygon(source, **self.l_sources)

    elif self.ports == 4:
      if self.sides == 2:
        port_1 = [(self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(self.y_pixels/4))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                  (self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(self.y_pixels/4))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
          #poly_1 = gdspy.Polygon(port_1, **self.l_port1)

        port_2 = [(self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(3*self.y_pixels/4))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                  (self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(3*self.y_pixels/4))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
          #poly_2 = gdspy.Polygon(port_2, **self.l_port2)

        port_3 = [(x_total - self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(self.y_pixels/4))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                  (x_total - self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(self.y_pixels/4))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
          #poly_3 = gdspy.Polygon(port_3, **self.l_port3)

        port_4 = [(x_total - self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(3*self.y_pixels/4))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                  (x_total - self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int(3*self.y_pixels/4))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
          #poly_4 = gdspy.Polygon(port_4, **self.l_port4)

        source = [(1, (math.ceil(int(self.y_pixels/4))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                  (1, (math.ceil(int(self.y_pixels/4))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
          #poly_5 = gdspy.Polygon(source, **self.l_sources)

      elif self.sides == 4:
        port_1 = [(self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int((y_total/self.pixelSize)/2))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                  (self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int((y_total/self.pixelSize)/2))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
          #poly_1 = gdspy.Polygon(port_1, **self.l_port1)

        port_2 = [((int(x_total/self.pixelSize)/2-1.05*(self.launch_pixels_w/2))*self.pixelSize, y_total - self.launch_pixels_l*self.pixelSize*0.5), \
                  ((int(x_total/self.pixelSize)/2+1.05*(self.launch_pixels_w/2))*self.pixelSize, y_total - self.launch_pixels_l*self.pixelSize*0.5)]
          #poly_2 = gdspy.Polygon(port_2, **self.l_port2)

        port_3 = [(x_total - self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int((y_total/self.pixelSize)/2))-\
                   1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                  (x_total - self.launch_pixels_l*self.pixelSize*0.5, (math.ceil(int((y_total/self.pixelSize)/2))+\
                   1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
          #poly_3 = gdspy.Polygon(port_3, **self.l_port3)

        port_4 = [((int(x_total/self.pixelSize)/2-1.05*(self.launch_pixels_w/2))*self.pixelSize, self.launch_pixels_l*self.pixelSize*0.5), \
                  ((int(x_total/self.pixelSize)/2+1.05*(self.launch_pixels_w/2))*self.pixelSize, self.launch_pixels_l*self.pixelSize*0.5)]
          #poly_4 = gdspy.Polygon(port_4, **self.l_port4)

        source = [(1, (math.ceil(int((y_total/self.pixelSize)/2))-1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize), \
                  (1, (math.ceil(int((y_total/self.pixelSize)/2))+1.05*(self.launch_pixels_w/2)+0.5)*self.pixelSize)]
          #poly_5 = gdspy.Polygon(source, **self.l_sources)
    
    return port_1, port_2, port_3, port_4, source

  def randomPixRect(self, x_dim, y_dim, pixel_map):
    # Design Random Structure
    # pixel_map = np.zeros((self.x_pixels,self.y_pixels),dtype=int)
    for x in range(self.x_pixels):
      for y in range(self.y_pixels):
        bit = random.random()
        if bit >= 0.5:
        # creates a rectangle for the "on" pixel
          if self.ports == 2:
            if self.sym == 'x-axis':
              if y >= self.y_pixels/2 - 0:
                break
              pixel_map[x+self.launch_pixels_l,y] = 1;
              pixel_map[x+self.launch_pixels_l,self.y_pixels - (y+1)] = 1;
            elif self.sym == 'y-axis':
              if x >= (self.x_pixels)/2 - 0:
                break
              pixel_map[x+self.launch_pixels_l,y] = 1;
              pixel_map[self.x_pixels + self.launch_pixels_l - (x+1),y] = 1;
            elif self.sym == 'xy-axis':
              if x>= (self.x_pixels)/2 - 0 or y >= self.y_pixels/2 - 0:
                break
              pixel_map[x+self.launch_pixels_l,y] = 1;
              pixel_map[self.x_pixels + self.launch_pixels_l - (x+1),y] = 1;
              pixel_map[x+self.launch_pixels_l,self.y_pixels - (y+1)] = 1;
              pixel_map[self.x_pixels + self.launch_pixels_l - (x+1),self.y_pixels - (y+1)] = 1;
            else: #Asymmetric
              pixel_map[x+self.launch_pixels_l,y] = 1;
          elif self.ports == 3:
            if self.sym == 'x-axis':
              if y >= self.y_pixels/2 - 0:
                break
              pixel_map[x+self.launch_pixels_l,y] = 1;
              pixel_map[x+self.launch_pixels_l,self.y_pixels - (y+1)] = 1;
            elif self.sym == 'y-axis':
              if x >= (self.x_pixels)/2 - 0:
                break
              if self.sides == 2:
                pixel_map[x+self.launch_pixels_l,y] = 1;
                pixel_map[self.x_pixels + self.launch_pixels_l - (x+1),y] = 1; 
              elif self.sides == 3:
                pixel_map[x+self.launch_pixels_l,y] = 1;
                pixel_map[self.x_pixels + self.launch_pixels_l - (x+1),y] = 1; 
            elif self.sym == 'xy-axis':
              if x>= (self.x_pixels)/2 - 0 or y >= self.y_pixels/2 - 0:
                break
              pixel_map[x+self.launch_pixels_l,y] = 1;
              pixel_map[self.x_pixels + self.launch_pixels_l - (x+1),y] = 1;
              pixel_map[x+self.launch_pixels_l,self.y_pixels - (y+1)] = 1;
              pixel_map[self.x_pixels + self.launch_pixels_l - (x+1),self.y_pixels - (y+1)] = 1;
            else: #Asymmetric
              if self.sides == 2:
                pixel_map[x+self.launch_pixels_l,y] = 1;
              elif self.sides == 3:
                pixel_map[x+self.launch_pixels_l,y] = 1;
          elif self.ports == 4:
            if self.sides == 2:
              pixel_map[x+self.launch_pixels_l,y] = 1;
            elif self.sides == 4:
              pixel_map[x+self.launch_pixels_l,y+self.launch_pixels_l] = 1;
          # pixel_map[x+self.launch_pixels_l,y] = 1;
        y += 1  
      x += 1
    return pixel_map
