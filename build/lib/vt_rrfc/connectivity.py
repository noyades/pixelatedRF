from skimage import measure 

class Connectivity:
  def find_connectivity(arr, pixelSize, ports, sides, portPos, connectivity=1):
    """
    Parameters:
    arr (ndarray): The input array to be labeled.
    pixelSize (int): The size of each pixel.
    ports (int): The number of ports.
    sides (int): The number of sides.
    portPos (list): The positions of the ports.

    Returns:
    tuple: A tuple containing boolean values indicating the connectivity between ports.
    """

    labeled = measure.label(arr, background=False, connectivity=connectivity)

    if ports < 2:
      c12 = False
      c13 = False
      c14 = False
      c23 = False
      c24 = False
      c34 = False
    if ports == 2:
      c13 = False
      c14 = False
      c23 = False
      c24 = False
      c34 = False
      if labeled[int(portPos[0] / pixelSize), int(portPos[1] / pixelSize)] == \
         labeled[int(portPos[2] / pixelSize) - 1, int(portPos[3] / pixelSize)]:
         c12 = True
      else:
         c12 = False
    if ports == 3:
      if sides == 2:
        c14 = False
        c24 = False
        c34 = False
        if labeled[int(portPos[0] / pixelSize), int(portPos[1] / pixelSize)] == \
           labeled[int(portPos[2] / pixelSize), int(portPos[3] / pixelSize)]:
           c12 = True
        else:
           c12 = False
        if labeled[int(portPos[0] / pixelSize), int(portPos[1] / pixelSize)] == \
           labeled[int(portPos[4] / pixelSize) - 1, int(portPos[5] / pixelSize)]:
           c13 = True
        else:
           c13 = False
        if labeled[int(portPos[2] / pixelSize), int(portPos[1] / pixelSize)] == \
           labeled[int(portPos[4] / pixelSize) - 1, int(portPos[5] / pixelSize)]:
           c23 = True
        else:
           c23 = False
      if sides == 3:
        c14 = False
        c24 = False
        c34 = False
        if labeled[int(portPos[0] / pixelSize), int(portPos[1] / pixelSize)] == \
           labeled[int(portPos[2] / pixelSize) - 1, int(portPos[3] / pixelSize)]:
           c12 = True
        else:
           c12 = False
        if labeled[int(portPos[0] / pixelSize), int(portPos[1] / pixelSize)] == \
           labeled[int(portPos[4] / pixelSize), int(portPos[5] / pixelSize) - 1]:
           c13 = True
        else:
           c13 = False
        if labeled[int(portPos[2] / pixelSize) - 1, int(portPos[1] / pixelSize)] == \
           labeled[int(portPos[4] / pixelSize), int(portPos[5] / pixelSize)-1]:
           c23 = True
        else:
          c23 = False
    if ports == 4:
      if sides == 2:
        if labeled[int(portPos[0] / pixelSize), int(portPos[1] / pixelSize)] == \
           labeled[int(portPos[2] / pixelSize), int(portPos[3] / pixelSize)]:
           c12 = True
        else:
           c12 = False
        if labeled[int(portPos[0] / pixelSize), int(portPos[1] / pixelSize)] == \
           labeled[int(portPos[4] / pixelSize) - 1, int(portPos[5] / pixelSize)]:
           c13 = True
        else:
           c13 = False
        if labeled[int(portPos[0] / pixelSize), int(portPos[1] / pixelSize)] == \
           labeled[int(portPos[6] / pixelSize) - 1, int(portPos[7] / pixelSize)]:
           c14 = True
        else:
          c14 = False
        if labeled[int(portPos[2] / pixelSize), int(portPos[3] / pixelSize)] == \
           labeled[int(portPos[4] / pixelSize) - 1, int(portPos[5] / pixelSize)]:
           c23 = True
        else:
           c23 = False
        if labeled[int(portPos[2] / pixelSize), int(portPos[3] / pixelSize)] == \
           labeled[int(portPos[6] / pixelSize) - 1, int(portPos[7] / pixelSize)]:
           c24 = True
        else:
           c24 = False
        if labeled[int(portPos[4] / pixelSize) - 1, int(portPos[5] / pixelSize)] == \
           labeled[int(portPos[6] / pixelSize) - 1, int(portPos[7] / pixelSize)]:
           c34 = True
        else:
           c34 = False
      if sides == 4:
        if labeled[int(portPos[0] / pixelSize), int(portPos[1] / pixelSize)] == \
           labeled[int(portPos[2] / pixelSize), int(portPos[3] / pixelSize) - 1]:
           c12 = True
        else:
          c12 = False
        if labeled[int(portPos[0] / pixelSize), int(portPos[1] / pixelSize)] == \
           labeled[int(portPos[4] / pixelSize) - 1, int(portPos[5] / pixelSize)]:
           c13 = True
        else:
           c13 = False
        if labeled[int(portPos[0] / pixelSize), int(portPos[1] / pixelSize)] == \
           labeled[int(portPos[6] / pixelSize), int(portPos[7] / pixelSize) - 1]:
           c14 = True
        else:
           c14 = False
        if labeled[int(portPos[2] / pixelSize), int(portPos[3] / pixelSize) - 1] == \
           labeled[int(portPos[4] / pixelSize) - 1, int(portPos[5] / pixelSize)]:
           c23 = True
        else:
           c23 = False
        if labeled[int(portPos[2] / pixelSize), int(portPos[3] / pixelSize) - 1] == \
           labeled[int(portPos[6] / pixelSize), int(portPos[7] / pixelSize) - 1]:
           c24 = True
        else:
           c24 = False
        if labeled[int(portPos[4] / pixelSize) - 1, int(portPos[5] / pixelSize)] == \
           labeled[int(portPos[6] / pixelSize), int(portPos[7] / pixelSize) - 1]:
           c34 = True
        else:
           c34 = False
    return c12, c13, c14, c23, c24, c34
