import numpy as np
import os
import time, sys
from vt_rrfc import microstrip as ustrip
from vt_rrfc import ustripComponents as usc
from vt_rrfc import ustripRandomComponents as usrc
from vt_rrfc import adsAelGeneration as ael

fc = 5.6e9
bw = 120e6
pixelSize = 18 # the size of the randomized pixel in mils. Typically contrained by a PCB manufacturer.

pathName = '/home/jswalling/pythonWork/rrfc/data/designMaps/' # Base path for file creation
inF = pathName + 'arlaBPF_pixelSize' + str(pixelSize) + '_fc=' + \
       str(np.round(fc*1e-9,1)) + 'GHz_BW=' + str(int(bw*1e-6)) + 'MHz.csv'

#os.chdir(pathName)

usrc.recreateGDS(pixelSize, inF)


