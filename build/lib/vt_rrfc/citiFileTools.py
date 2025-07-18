import re
import os
from datetime import datetime
import numpy as np


def readCiti(file_path, verbose=False):
    """
    Load the contents of a CITI file into an NPort
    
    :returns: NPort holding data contained in the CITI file
    :rtype: :class:`nport.NPort`
    
    """
    file_path = os.path.abspath(file_path)
    citifile = CITIFile(file_path)
    assert citifile.params[0][0][0].lower() == "freq"
    freqs = citifile.data[0][0]
    ports = np.sqrt(len(citifile.params[0]) - 1)
    ports = int(ports)
    
    #re_param = re.compile(r"^S\[(\d+),(\d+)\]$")
    indices = []
    for param in citifile.params[0][1:(ports**2)+1]:
        name = param[0]
        m = re.match(r"^S\[(\d+),(\d+)\]$", name)
        port1 = int(m.group(1))
        port2 = int(m.group(2))
        indices.append((port1, port2))

    S = np.zeros((len(freqs),ports**2),dtype=complex)
    for index in range(len(freqs)):
      k = 0
      for k in range(0,ports**2):
        S[index,k] = citifile.data[0][k+1][index]
    return freqs, S

# Collection of object classes for reading calibration lab data file types
#
# Author:           J. Wayde Allen
# Creation Date:    2001-05-22
# Revised:          2001-05-23 JWA
#                   2010-01-28 Brecht Machiels
#                               * made parsing more robust
#                               * changed indentation from 3 to 4 spaces
#
#  The software was developed and is owned by ITS/NTIA, an agency
#  of the Federal Government.  Pursuant to title 15 United States
#  Code Section 105, works of Federal employees are not subject to
#  copyright protection in the United States.  This software is
#  provided by ITS as a service and is expressly provided "AS IS".
#  NEITHER ITS NOR NTIA MAKES ANY WARRANTY OF ANY KIND, EXPRESS,
#  IMPLIED OR STATUTORY, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
#  WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
#  NON-INFRINGEMENT AND DATA ACCURACY.  ITS/NTIA does warrant or
#  make any representations regarding the use of the software or
#  the results thereof, including but not limited to the
#  correctness, accuracy, reliability or usefulness of the
#  software.
#
#  This software is free software;  you can use, copy, modify and
#  redistribute it upon your acceptance of these terms and
#  conditions and upon your express agreement to provide
#  appropriate acknowledgements of ITS/NTIA's ownership of and
#  development of this software by keeping this exact text present
#  in any copied or derivative works.

import string, sys

class CITIFile:
    def __init__(self, filename):
        self.filename = filename

        # The following are the main data structures
        self.packages = {}        
        self.constants = []
        self.params = []
        self.data = []
        self.instrmnt = []

        # Open the citifile 
        myfile = open(self.filename, 'r')

        # Define some special control and book keeping variables
        packagecounter = -1 # Index to the number of Citifile packages
        packagenames = []   # List of the package names

        while 1:
            line = myfile.readline()

            if not line:
                break
            
            linetxt = str.strip(line)
            line = str.split(linetxt)

            #This line starts a new Citifile data package
            #update the package counter and create blank indices
            if len(line) > 0:

                if line[0] == 'CITIFILE':
                    packagecounter = packagecounter + 1
                    packagenames.append("")  #Create a blank name entry

                    self.constants.append([])
                    self.params.append([])
                    self.data.append([])
                    self.instrmnt.append([])

                    indata = 'NO'       #Not reading data
                    invarlist = 'NO'    #Not reading independant variable data
                    datacount = 0       #Index to package data blocks

                #Skip device-specific variables
                if line[0][0] == '#':
                    continue

                #Should be one name per package
                elif line[0] == 'NAME':  
                    packagenames[packagecounter] = line[1]

                elif line[0] == 'CONSTANT':
                    self.constants[packagecounter].append((line[1],line[2]))

                elif line[0] == 'VAR':
                    self.params[packagecounter].append((line[1],line[2],line[3]))

                elif line[0] == 'SEG_LIST_BEGIN':
                    invarlist = 'SEG'
                    self.data[packagecounter].append([])

                elif line[0] == 'SEG' and invarlist == 'SEG':

                    #Decode the start, stop and number of points entries
                    start = float(line[1])
                    stop  = float(line[2])
                    numpoints  = int(line[3])

                    #Compute the actual data values from this information
                    #and put it in the data block
                    step = (stop - start) / (numpoints - 1)
                    next = start
                    count = 0
                    while next <= stop:
                        count = count + 1
                        self.data[packagecounter][datacount].append(next)
                        next = next + step

                elif line[0] == 'SEG_LIST_END':
                    invarlist = 'NO'
                    #We've filled this data bin so point to the next one
                    datacount = datacount + 1

                elif line[0] == 'VAR_LIST_BEGIN':
                    invarlist = 'VARLIST'
                    self.data[packagecounter].append([])

                elif line[0] != 'VAR_LIST_END' and invarlist == 'VARLIST':
                    datum = float(line[0])
                    self.data[packagecounter][datacount].append(datum)

                elif line[0] == 'VAR_LIST_END':
                    invarlist = 'NO'
                    datacount = datacount + 1

                elif line[0] == 'DATA':
                    self.params[packagecounter].append((line[1],line[2]))

                elif line[0] == 'BEGIN': 
                    indata = 'YES'
                    self.data[packagecounter].append([])

                elif line[0] != 'END' and indata == 'YES':

                    if self.params[packagecounter][datacount][1] == 'RI':
                        real,imag = str.split(linetxt,',')
                        value = complex(float(real),float(imag))

                    elif self.params[packagecounter][datacount][1] == 'MAG':
                        value = float(line[0])

                    self.data[packagecounter][datacount].append(value)

                elif line[0] == 'END': 
                    indata = 'NO'
                    datacount = datacount + 1

                else:
                    #Anything else must be instrument specific so make these
                    #lines available for parsing by the user
                    self.instrmnt[packagecounter].append(line)

        #We've read and sorted all of these data
        #Create dictionary of package index and names
        for values in range(0,packagecounter+1):
            self.packages[values] = packagenames[values]
