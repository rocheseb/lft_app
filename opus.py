# -*- coding: utf-8 -*-

from __future__ import print_function # allows the use of Python 3.x print function in python 2.x code so that print('a','b') prints 'a b' and not ('a','b')

"""
This code contain a class that can be used to read Bruker OPUS files.

Acknowledgment:

Part of this code was developed by the Infrared Observations team at BIRA-IASB who received funding from the AGACC-2 project  
(Advanced exploitation of Ground-based measurements for Atmospheric Chemistry and Climate applications-II), 
CONTRAT N° SD/CS/07A of the Belgian Science for Sustainable Development research programme ), 
the Belgian federal support to ICOS (ministerial decree FR/35/IC1), 
the EU FP7 ICOS_Inwire project and BIRA-IASB own resources.
"""

import struct

import os.path

import numpy as np

import datetime as dt

class Opus:

    """Retrieve data from a Bruker OPUS file

    Arguments:

    opus_file_path -- the path to a Bruker OPUS file

    """

    def __init__(self, opus_file_path):

        if not os.path.isfile(opus_file_path):

            self.file = ''

        elif os.path.getsize(opus_file_path) < 16:  # (size in bytes) check if file is not completely empty

            self.file = ''

        else:

            self.file = str(opus_file_path)

        self.param = []

        self.xdata = []

        self.ydata = []

        self.subID = []

        self.time  = []

    def get_data(self, request=''):

        """Read parameters and/or data from the Bruker Opus file

        Arguments:

        request -- a string, detailing the data that has to be read from the file, multiple parameters can be combined

                    "p" - get parameters only, no data

                    "r" - return only reduced a subset of the data in the file (e.g. for a quick visual quality-check)

        """

        if len(self.file) == 0:  # check if the file entry is not empty

            return

        wantsdata = True

        reduced = False

        # see what the user has requested

        if 'p' in request:

            wantsdata = False

        if 'r' in request:

            reduced = True

        file_size = os.path.getsize(self.file)

        with open(self.file, 'rb') as fid:  # open the file

            # The first 4 bytes contain the magic number of the OPUS file

            magic = -16905718  # this magic number identifies the file as being an OPUS file

            cbin = fid.read(4)

            cdat = struct.unpack('i', cbin)

            if not (cdat[0] == magic):

                print("Magic number mismatch")

            #print(cdat,cdat[0],magic)

            fid.read(8)  # skip 8 bytes of version info

            # The next 3 integers contain the structure of the directory block

            cbin = fid.read(3 * 4)

            cdat = struct.unpack('iii', cbin)

            pdb = cdat[0]  # start of the directory block (bytes)

            mndb = cdat[1]  # max number of entries in directory block

            ndb = cdat[2]  # actual number of entries in directory block

            if ndb > mndb:
                print("Invalid file structure: ndb > mndb")
                fid.close()
                return

            elif ndb < 2:
                print("The file appears to contain no data")
                fid.close()
                return

            # Now we read each entry in the directory block

            opusdir = []

            file_size_from_dir = 0  # size of the file in bytes, based on the directory structure

            for bid in range(ndb):

                fid.seek(pdb + bid * 12)  # each entry has a length of 12 bytes (BBHII i.e. 11244)

                cbin = fid.read(12)  # Entry contains (type, subtype, ???, length(units of int32), pointer(in bytes))

                if len(cbin) != 12:
                    print("Invalid directory entry in file: {file}".format(file=self.file))
                    return

                cdat = struct.unpack('BBHII', cbin)  # each entry refers to a block further down in the binary file

                opusdir.append(cdat)

                cur_size = (cdat[3] * 4 + cdat[4])  # byte pointer + block size (in units in int32)

                if cur_size > file_size_from_dir:

                    file_size_from_dir = cur_size

            if file_size_from_dir > file_size:  # check if the file is big enough to hold all the promised data
                print("File size is too small: {file}".format(file=self.file))
                return

            parameters = []

            pairs = {}  # will store the parameter, value pairs

            for bid in range(1, ndb):  # we skip the first block since this is the directory block which we have just parsed

                if opusdir[bid][0] in [160, 48, 32, 96, 23, 64]:

                    if opusdir[bid][0] in [23]:  # we treat data parameter blocks differently

                        datapairs = {'subID': opusdir[bid][1]}

                    # move to the beginning of a parameter block

                    fid.seek(opusdir[bid][4])

                    while fid.tell() < opusdir[bid][4] + opusdir[bid][3]*4:

                        # First we identify the type of the current entry in the block

                        cbin = fid.read(8)

                        if len(cbin) != 8:
                            print("Invalid parameter entry, skipping remaider of block")
                            break

                        cdat = struct.unpack('ccccHH', cbin)  # H = unsighned short (2 bytes)

                        pname = ''.join(cdat[0:3])  # parameter name

                        if pname.lower() == "end":  # skip to the next block when we reach then END entry

                            break

                        ptype = cdat[4]  # parameter type

                        plen = cdat[5]  # half-length in bytes of the entry i.e. plen * 2 is the length in bytes

                        if ptype in [2, 3, 4]:  # String of Enum

                            cbin = fid.read(2*plen)

                            cdat = struct.unpack('c'*2*plen, cbin)

                            pval = ''.join(cdat).split('\x00')[0]  # String ends at first ASCII char 0

                        elif ptype == 0 and plen == 2:  # Int32

                            cbin = fid.read(4)

                            cdat = struct.unpack('i', cbin)

                            pval = cdat[0]

                        elif ptype == 1 and plen == 4:  # Real64

                            cbin = fid.read(8)

                            cdat = struct.unpack('d', cbin)

                            pval = cdat[0]

                        elif plen == 0:

                            pval = 0

                        elif ptype > 16:
                            print("Invalid parameter type, skipping remainder of block")
                            break

                        else:  # Unknown type

                            fid.read(2*plen)  # skip the entry
                            print("Unknown parameter type, it will be ignored")
                            continue

                        if opusdir[bid][0] in [23]:  # we treat data parameter blocks differently

                            datapairs[pname] = pval

                        else:

                            pairs[pname] = pval

                    if opusdir[bid][0] in [23]:

                        parameters.append(datapairs)  # Add each individual data parameters block to the output

            parameters.insert(0, pairs)  # Add all the 'non data' parameter blocks together to the output

            # find the measurement times from the data parameter blocks

            try: [self.time.append(dt.datetime.strptime((p['DAT']+' '+p['TIM'])[:19],'%d/%m/%Y %H:%M:%S')) for p in parameters[1:]]

            except ValueError: [self.time.append(dt.datetime.strptime((p['DAT']+' '+p['TIM'])[:19],'%Y/%m/%d %H:%M:%S')) for p in parameters[1:]]

            # store the parameters in the class

            self.param = parameters

            txdata = []

            tydata = []

            tsubID = []

            if wantsdata:  # If the user has requested data, we get it

                for bid in range(1, ndb):  # skip the first block: this is the directory block which we have just parsed

                    if opusdir[bid][0] in [7]:  # This time we focus only on the data blocks

                        # find the corresponding data parameter block based on the subID

                        dind = 0

                        for sis in range(1, len(parameters)):

                            if parameters[sis]['subID'] == opusdir[bid][1]:

                                dind = sis

                                break

                        if dind == 0:

                            fid.close()

                            return

                        fid.seek(opusdir[bid][4])

                        # read the binary data from file, and convert it to floating point data

                        cdat = np.fromfile(fid, dtype="f", count=parameters[dind]['NPT'])

                        # create a linearly spaced x-axis based on the parameters in the data block

                        xtemp = np.linspace(start=parameters[dind]['FXV'],

                                            stop=parameters[dind]['LXV'],

                                            num=parameters[dind]['NPT'])

                        if reduced:

                            if parameters[dind]['NPT'] > 1e5:

                                step = int(parameters[dind]['NPT']/1e5 + 1)

                                idxs = np.arange(0, parameters[dind]['NPT'], step)

                                cdat = [cdat[int(j)] for j in idxs]

                                xtemp = [xtemp[int(j)] for j in idxs]

                        tydata.append(cdat)

                        txdata.append(xtemp)

                        tsubID.append(parameters[dind]['subID'])  # changed sis into dind

                self.xdata = txdata

                self.ydata = tydata

                self.subID = tsubID