import numpy as np
import os
import time, sys
from vt_rrfc import *

fc  = 6e9
fc2 = 2e9

pixelSize = 18 # the size of the randomized pixel in mils. Typically contrained by a PCB manufacturer.

pathName = '/home/jswalling/pythonWork/rrfc/data/designMaps/' # Base path for file creation
inF = pathName + 'arlaBPF_pixelSize' + str(pixelSize) + 'notch_fc=' + \
       str(np.round(fc2*1e-9,1)) + 'GHz_pass_fc=' + str(np.round(fc*1e-9,1)) + 'MHz_v3.csv'

rrfc1 = rfc(pixelSize=pixelSize,outF=inF,unit=25.4e-6,sim=0,view=0,write=0)
#os.chdir(pathName)

recreateGDS(rrfc1)


