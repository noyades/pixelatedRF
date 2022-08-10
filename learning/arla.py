# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import skrf as rf
import tensorflow as tf
import sys
from ch02 import ch02

from tensorflow.keras import datasets, layers, models

# Path to data files 
spPath = 'H:/My Drive/Research/randomRFComponents/data/spFiles/'
pixPath = 'H:/My Drive/Research/randomRFComponents/data/pixelMaps/'
desPath = 'H:/My Drive/Research/randomRFComponents/data/designMaps/'
pixelSize = 18
# This code loads the pixel maps, that map the geometry of the randomized 
# components. Each binary number in the map represents whether a pixel is 
# populated with metal or not. The pixel size is in the file name in units of 
# mils, so the map can be used to reconstruct a gds file; gds files are also 
# saved and available in H:/My Drive/Research/randomRFComponents/data/gdsFiles
csvFile = pixPath + 'randomGDSSteppedImpFilter_order=3_pixelSize=' + \
          str(pixelSize) + '_sim=2.csv'
initPix = np.loadtxt(csvFile, delimiter=',')
rows = np.size(initPix,0)
cols = np.size(initPix,1)
n = 1357
pixMap = np.zeros((n,rows,cols),dtype=int)
for x in range(0,n):
    csvFile = pixPath + 'randomGDSSteppedImpFilter_order=3_pixelSize=' + \
              str(pixelSize) + '_sim='\
        + str(x) + '.csv'
    pixMap[x,:,:] = np.loadtxt(csvFile, delimiter=',')
  
# All S-parameter data is collected with the same frequency spacing (0 - 10GHz)
# in 50 MHz increments. This loads a frequency variable and pre-determines the 
# size of the data arrays
initSP = rf.Network(spPath + \
      'randomGDSSteppedImpFilter_Type=butter_order=3_pixelSize=' + \
                str(pixelSize) + '_sim=1.s2p')
spFreq = initSP.frequency._f
m = np.size(spFreq)    
spFreq = spFreq.reshape(1,m)
initSP 

# Load the S-parameter data for the randomized component simulations
absSp11 = np.zeros((n,m),dtype=float)
phaSp11 = np.zeros((n,m),dtype=float)
absSp12 = np.zeros((n,m),dtype=float)
phaSp12 = np.zeros((n,m),dtype=float)
absSp21 = np.zeros((n,m),dtype=float)
phaSp21 = np.zeros((n,m),dtype=float)
absSp22 = np.zeros((n,m),dtype=float)
phaSp22 = np.zeros((n,m),dtype=float)
rl_in = np.zeros((n,m),dtype=float)
rl_out = np.zeros((n,m),dtype=float)
il = np.zeros((n,m),dtype=float)
for x in range(0,n):
    se_ntwk = rf.Network(spPath + \
          'randomGDSSteppedImpFilter_Type=butter_order=3_pixelSize=' + \
                    str(pixelSize) + '_sim=' + \
          str(x) + '.s2p')
    spt = se_ntwk.s11._s
    absSp11[x,:] = abs(spt.reshape(1,m))
    phaSp11[x,:] = np.angle(spt.reshape(1,m))
    spt = se_ntwk.s12._s
    absSp12[x,:] = abs(spt.reshape(1,m))
    phaSp12[x,:] = np.angle(spt.reshape(1,m))
    spt = se_ntwk.s21._s
    absSp21[x,:] = abs(spt.reshape(1,m))
    phaSp21[x,:] = np.angle(spt.reshape(1,m))
    spt = se_ntwk.s22._s
    absSp22[x,:] = abs(spt.reshape(1,m))
    phaSp22[x,:] = np.angle(spt.reshape(1,m))

    rl_in[x,:] = 20*np.log10(absSp11[x,:])
    rl_out[x,:] = 20*np.log10(absSp22[x,:])
    il[x,:] = 20*np.log10(absSp21[x,:])
del spt

# In the next part of the algorithm, a reinforcement parameter will be defined.
# The parameter will take into account the desire for good input/output match,
# and desire for strong transmission in particular band or bands. This will be
# used to calculated a cumulative reward matrix that will drive the evaluation
# of future structures; hence this portion of the design is to use a set of 
# random simulation data to train the reward function. The algorithm is 
# presented in:
# Sourangsu Banerji, Apratim Majumder, Alexander Hamrick, Rajesh Menon, 
# Berardi Sensale-Rodriguez,
# Machine learning enables design of on-chip integrated silicon T-junctions 
# with footprint of 1.2 μm×1.2μm,
# Nano Communication Networks,
# Volume 25,
# 2020,
# 100312,
# ISSN 1878-7789,
# https://doi.org/10.1016/j.nancom.2020.100312.

# For this example, we will assume that there is one objective, to design a 
# passband filter with a chosen center frequency and bandwidth. 
fc = 5.6e9 
bw = 120e6
fs = fc - bw/2
fe = fc + bw/2
spIntFreq = np.arange(0,10.005e9,5e6) # frequency vector to interpolate data to
# ideally will ensure that fs and fe land on grid
spIntFreq = spIntFreq.tolist()

fs_index = spIntFreq.index(fs)
fe_index = spIntFreq.index(fe)
    

from scipy.interpolate import interp1d
spIntFreq = np.arange(0,10.005e9,5e6) # set a frequency range to interpolate to
il_new = np.zeros((n,np.size(spIntFreq)),dtype=float)
rw = np.zeros((n),dtype=float)
for x in range(0,n):
    f = interp1d(spFreq[0,:],il[x,:])
    il_new = f(spIntFreq)
    il_max = -1.5 # Target insertion loss for the filter
    rw1 = (1/(fe_index-fs_index))*np.sum(2**(10**((il_new[fs_index:fe_index]-il_max)/20)))
    f = interp1d(spFreq[0,:],rl_in[x,:])
    rl_in_new = f(spIntFreq)
    rl_in_min = -10 # Target input return loss for the filter
    rw2 = (1/(fe_index-fs_index))*np.sum(2**(10**((rl_in_min-rl_in_new[fs_index:fe_index])/20)))
    f = interp1d(spFreq[0,:],rl_out[x,:])
    rl_out_new = f(spIntFreq)
    rl_out_min = -10 # Target output return loss for the filter
    rw3 = (1/(fe_index-fs_index))*np.sum(2**(10**((rl_out_min-rl_out_new[fs_index:fe_index])/20)))    
    rw[x] = rw1 + rw2 + rw3 # Reward Function, assumes equal weight for input, 
    # output match and insertion loss 

# Reshape the pixel maps for each of the training designs so that the map is in
# a single row for each training set
rePixMap = pixMap.reshape((n,rows*cols))
crm = (np.dot(rw,rePixMap) - np.mean(np.dot(rw,rePixMap)))/np.std(np.dot(rw,rePixMap))

# Implement the ARLA algorithm
invDesPixMap = np.zeros((1,np.size(rePixMap,1)),dtype=int)
for x in range(0,np.size(rePixMap,1)):
    if crm[x] > 0:
        invDesPixMap[0,x] = 1

outF = desPath + 'arlaBPF_pixelSize' + str(pixelSize) + '_fc=' + \
       str(np.round(fc*1e-9,1)) + 'GHz_BW=' + str(int(bw*1e-6)) + 'MHz'
invDesPixMap = invDesPixMap.reshape(rows,cols)
csvFile = outF + ".csv"
# Export Pixel Map file
np.savetxt(csvFile, invDesPixMap, fmt = '%d', delimiter = ",")
