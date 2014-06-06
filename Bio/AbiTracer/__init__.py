# Copyright (C) 2014 David Bulger
# This is modified version of seqtrace-0.9.0 and Biopython packages
# Both copyrights for those packages can be found below

# Copyright 2011 by Wibowo Arindrarto (w.arindrarto@gmail.com)
# Revisions copyright 2011 by Peter Cock.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

# Copyright (C) 2014 Brian J. Stucky
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Bio.SeqIO parser for the ABI format.

ABI is the format used by Applied Biosystem's sequencing machines to store
sequencing results.

For more details on the format specification, visit:
http://www.appliedbiosystem.com/support/software_community/ABIF_File_Format.pdf

"""
__docformat__ = "epytext en"

# Import
import datetime
import struct
from struct import unpack, pack
import os.path
from datetime import datetime
from os.path import basename

# Initialization of Key Values
tracesamps = {}
fname = ''
max_traceval = -1
comments = {}

def AbiIterator(input_file_name, trim=False):
    #Iterator for the Abi file format

    abi = open(input_file_name, mode = 'rb')

    # check if input file is a valid Abi file # Read ABI magic number
    abi.seek(0)
    abinum = abi.read(4)
    #print abinum

    # check the major version number
    try:
        version = struct.unpack('>H', abi.read(2))[0]
    except struct.error:
        raise ABIError('The ABI file header is invalid.  The file appears to be damaged.')
    #print version
    if (version / 100) != 1:
        raise ABIVersionError(version / 100, version % 100)
        
    # skip the next 10 bytes
    abi.read(10)
        
    # get the file index information
    try:
        index_entry_len = struct.unpack('>h', abi.read(2))[0]
        num_index_entries = struct.unpack('>i', abi.read(4))[0]
        total_index_size = struct.unpack('>i', abi.read(4))[0]
        index_offset = struct.unpack ('>i', abi.read(4))[0]
    except struct.error:
        raise ABIError('The ABI file header is invalid.  The file appears to be damaged.')
        
    #index_entry_len, num_index_entries, total_index_size, index_offset

    # read the ABI index block
    abi.seek(index_offset, 0)
    abiindex = list()

    for cnt in range(num_index_entries):
        try:
            abiindex.append(dict(did=0, idv=0, dformat=0, fsize=0, dcnt=0, dlen=0, offset=0))
            abiindex[cnt]['did'] = abi.read(4)
            abiindex[cnt]['idv'] = struct.unpack('>I', abi.read(4))[0]
            abiindex[cnt]['dformat'] = struct.unpack('>H', abi.read(2))[0]
            abiindex[cnt]['fsize'] = struct.unpack('>H', abi.read(2))[0]
            abiindex[cnt]['dcnt'] = struct.unpack('>I', abi.read(4))[0]
            abiindex[cnt]['dlen'] = struct.unpack('>I', abi.read(4))[0]
            abiindex[cnt]['offset'] = struct.unpack('>I', abi.read(4))[0]
            # skip 4 bytes (the unused "data handle" field)
            abi.read(4)
            #print abiindex
        except struct.error:
            raise ABIIndexError(cnt, num_index_entries)
    #printABIIndex(abiindex,'FWO_')
    

    # Read Base Calls from file
    row = getIndexEntry(abiindex,'PBAS', 1)
    if row is None:
        raise ABIError('No base call data were found in the ABI file.  The file might be damaged.')
    basecalls = readString(abi,row).upper()
    #print basecalls
    
    # Read Base Call Confidence Scores (Quality Scores) from file
    row = getIndexEntry(abiindex,'PCON', 1)
    if row is None:
        raise ABIError('No confidence score data were found in the ABI file.  SeqTrace requires confidence scores for all base calls.')
    bcconf = read1ByteInts(abi,row)
    #print bcconf

    # Read Trace Data
    base_order = getBaseDataOrder(abiindex)
    maxval = 0    
    # This is the ID for the first 'DATA' index entry that points to the processed
    # trace data.  The man page for the Staden program convert_trace states
    # that IDs 9-12 contain the processed data; IDs 1-4 contain the raw data.
    # The ABIF documentation from ABI also suggests that IDs 1-8 will always contain
    # raw data, and 9-12 will contain the processed data.  Is this always correct?
    start_id = 9
    for cnt in range(0, 4):
        row = getIndexEntry(abiindex,'DATA', start_id + cnt)
        if row == None:
            raise ABIError('Could not find trace data index entries for all bases.  The file might be damaged.')
        # read the trace data from the file
        lst = read2ByteInts(abi,row)
        tmpmax = max(lst)
        if tmpmax > maxval:
            maxval = tmpmax
        tracesamps[base_order[cnt]] = lst
        #print base_order
        #print "tmpmax: " + str(tmpmax)
    max_traceval = maxval
    #print "max_traceval: " + str(max_traceval)

    # read the base call locations from the file
    row = getIndexEntry(abiindex,'PLOC', 1)
    if row is None:
        raise ABIError('No base location data were found in the ABI file.  The file might be damaged.')
    basepos = read2ByteInts(abi,row)
    #print basepos

    # read the comments from the file
    #readComments(abi,abiindex)

   
    # use the file name as SeqRecord.name if available
    try:
        file_name = basename(abi.name).replace('.ab1', '')
    except:
        file_name = ""

    record = output_data(seq=basecalls,
                       name=file_name,
                       phd=bcconf,
                       DATA1=tracesamps[base_order[0]], 
                       DATA2=tracesamps[base_order[1]], 
                       DATA3=tracesamps[base_order[2]], 
                       DATA4=tracesamps[base_order[3]], 
                       FWO = base_order,
                       POS = basepos)                                 
    yield record

class output_data(object):
    def __init__(self, seq, name,
                 phd, 
                 DATA1,
                 DATA2,
                 DATA3,
                 DATA4,
                 FWO,
                 POS):
        self.seq = seq
        self.name = name
        self.phd = phd
        self.DATA1 = DATA1
        self.DATA2 = DATA2
        self.DATA3 = DATA3
        self.DATA4 = DATA4
        self.FWO = FWO
        self.POS = POS

class TraceFileError(Exception):
    pass

class ABIError(TraceFileError):
    pass

class ABIVersionError(ABIError):
    def __init__(self, ver_major, ver_minor):
        self.ver_major = ver_major
        self.ver_minor = ver_minor

    def __str__(self):
        return 'This file uses version ' + str(self.ver_major) + '.' + str(self.ver_minor) + ' of the ABI format.  This software only supports version 1.x of the format.'

class ABIIndexError(ABIError):
    def __init__(self, indexnum, indextotal):
        self.indexnum = indexnum
        self.indextotal = indextotal

    def __str__(self):
        return 'Error reading ABI file index entry ' + str(self.indexnum) + ' of ' + str(self.indextotal) + ' expected entries.  The file might be damaged.'

class ABIDataError(ABIError):
    def __init__(self, expectedlen, actuallen):
        self.expectedlen = expectedlen
        self.actuallen = actuallen

    def __str__(self):
        return 'Error reading ABI file data.  Expected ' + str(self.expectedlen) + ' bytes but only got ' + str(self.actuallen) + ' bytes.  The file appears to be damaged.'

def printABIIndex(abiindex, data_id):
    for entry in abiindex:
        if entry['did'] == data_id:
            print 'entry ID:', entry['did']
            print 'idv:', entry['idv']
            print 'data format:', entry['dformat']
            print 'format size:', entry['fsize']
            print 'data count:', entry['dcnt']
            print 'total data length:', entry['dlen']
            print 'data offset:', entry['offset']

def getIndexEntry(abiindex,data_id, number):
        for row in abiindex:
            if (row['did'] == data_id) and (row['idv'] == number):
                return row
        return None
    
def getIndexEntriesById(abiindex,data_id):
        entries = list()
        for row in abiindex:
            if row['did'] == data_id:
                entries.append(row)
        return entries

# Attempts to get a bunch of information about the sequencing run from the ABI file.  As much as possible,
# the keys used for individual comment values correspond with the keys used for the same values by the
# Staden software package.  However, this method also retrieves some comments that are not read by the
# Staden package.  To avoid confusion, these additional comment values are not given 4-letter keys.
def readComments(abi,abiindex):
    # get the sample name
    entry = getIndexEntry(abiindex,'SMPL', 1)
    if entry:
        comments['NAME'] = readString(abi,entry)
    # get the run name
    entry = getIndexEntry(abiindex,'RunN', 1)
    if entry:
        comments['Run name'] = readString(abi,entry)

    # get the lane number
    entry = getIndexEntry(abiindex,'LANE', 1)
    if entry:
        comments['LANE'] = str(read2ByteInts(abi,entry)[0])

    # get the signal strengths for each dye
    entry = getIndexEntry(abiindex,'S/N%', 1)
    if entry:
        stvals = read2ByteInts(abi,entry)

        # use the "filter wheel order" to determine the base/value pairings
        order = getBaseDataOrder(abiindex)
        sigst = {}
        for cnt in range(0, len(order)):
            sigst[order[cnt]] = stvals[cnt]

        comments['SIGN'] = 'A={0},C={1},G={2},T={3}'.format(sigst['A'], sigst['C'], sigst['G'], sigst['T'])

    # get the average peak spacing
    entry = getIndexEntry(abiindex,'SPAC', 1)
    if entry:
        spacing = read4ByteFloats(abi,entry)[0]
        # if spacing is invalid, estimate it ourselves (the Staden code [seqIOABI.c] indicates this is a possibility)
        if spacing < 0:
            spacing = float(basepos[-1] - basepos[0]) / (len(basepos) - 1)
        comments['SPAC'] = '{0:.2f}'.format(spacing)

    # get the run dates and times
    d_entries = getIndexEntriesById(abiindex,'RUND')
    t_entries = getIndexEntriesById(abiindex,'RUNT')
    if (len(d_entries) > 1) and (len(t_entries) > 1):
        sdate = readDateTime(abi,getIndexEntry(abiindex,'RUND', 1), getIndexEntry(abiindex,'RUNT', 1))
        edate = readDateTime(abi,getIndexEntry(abiindex,'RUND', 2), getIndexEntry(abiindex,'RUNT', 2))
        #print sdate, edate
        comments['RUND'] = sdate.strftime('%Y%m%d.%H%M%S') + ' - ' + edate.strftime('%Y%m%d.%H%M%S')
        comments['DATE'] = sdate.strftime('%a %d %b %H:%M:%S %Y') + ' to ' + edate.strftime('%a %d %b %H:%M:%S %Y')

    # get the data collection dates and times
    if (len(d_entries) == 4) and (len(t_entries) == 4):
        sdate = readDateTime(abi,getIndexEntry(abiindex,'RUND', 3), getIndexEntry(abiindex,'RUNT', 3))
        edate = readDateTime(abi,getIndexEntry(abiindex,'RUND', 4), getIndexEntry(abiindex,'RUNT', 4))
        #print sdate, edate
        comments['Data coll. dates/times'] = sdate.strftime('%a %d %b %H:%M:%S %Y') + ' to ' + edate.strftime('%a %d %b %H:%M:%S %Y')

    # get the dye set/primer (mobility) file
    entry = getIndexEntry(abiindex,'PDMF', 1)
    if entry:
        comments['DYEP'] = readString(abi,entry)

    # get the sequencing machine name and serial number
    entry = getIndexEntry(abiindex,'MCHN', 1)
    if entry:
        comments['MACH'] = readString(abi,entry)

    # get the sequencing machine model
    entry = getIndexEntry(abiindex,'MODL', 1)
    if entry:
        comments['MODL'] = readString(abi,entry)

    # get the basecaller name
    entry = getIndexEntry(abiindex,'SPAC', 2)
    if entry:
        comments['BCAL'] = readString(abi,entry)

    # get the data collection software version
    entry = getIndexEntry(abiindex,'SVER', 1)
    if entry:
        comments['VER1'] = readString(abi,entry)

    # get the basecaller version
    entry = getIndexEntry(abiindex,'SVER', 2)
    if entry:
        comments['VER2'] = readString(abi,entry)

    # get the plate size
    entry = getIndexEntry(abiindex,'PSZE', 1)
    if entry:
        comments['Plate size'] = str(read4ByteInts(abi,entry)[0])

    # get the gel name
    # This is included here because it is read by the Staden package, but it does not appear to be
    # included in the modern ABIF documentation.
    entry = getIndexEntry(abiindex,'GELN', 1)
    if entry:
        comments['GELN'] = readString(abi,entry)

    # get the instrument (matrix) file
    # This is included here because it is read by the Staden package, but it does not appear to be
    # included in the modern ABIF documentation.
    entry = getIndexEntry(abiindex,'MTXF', 1)
    if entry:
        comments['MTXF'] = readString(abi,entry)

    # 'APrX' points to a long XML string with detailed information about the analysis protocol used
    #entry = self.getIndexEntry('APrX', 1)
    #if entry:
    #    self.readString(entry)

def readDateTime(abi, dateindexrow, timeindexrow):
    # date format:
    #   bits 31-16: year
    #   bits 15-8: month
    #   bits 7-0: day of month
    # time format:
    #   bits 31-24: hour
    #   bits 23-16: minutes
    #   bits 15-8: seconds
    datenum = read4ByteInts(abi,dateindexrow)[0]
    timenum = read4ByteInts(abi,timeindexrow)[0]
    dateobj = datetime(year=(datenum >> 16), month=((datenum >> 8) & 0xff), day=(datenum & 0xff),
            hour=(timenum >> 24), minute=((timenum >> 16) & 0xff), second=((timenum >> 8) & 0xff))

    return dateobj

def readString(abi, indexrow):
    if indexrow['fsize'] != 1:
        raise ABIError('Index entry contains an invalid format size for string data.')
    if indexrow['dformat'] not in (2, 18, 19):
        raise ABIError('Index entry contains an invalid data type for character data.')
   
    if indexrow['dlen'] <= 4:
        # The actual data are stored in the offset field of the index entry.  Because the offset
        # was read as an unsigned, big-endian integer, the bytes should be in the correct order for
        # the following bit shift operations.
        lst = list()
        for cnt in range(0, indexrow['dcnt']):
            val = (indexrow['offset'] >> ((3 - cnt) * 8)) & 0xff
            lst.append(chr(val))
    
        strval = ''.join(lst)
    else:
        # get the data from the file
        abi.seek(indexrow['offset'], 0)
        strval = abi.read(indexrow['dcnt'])
    
    if indexrow['dlen'] != len(strval):
        raise ABIDataError(indexrow['dlen'], len(strval))

    # If this is a Pascal-style string (format 18), then remove the first character (which specifies
    # the string length).  If this is a C-style string (format 19), then remove the trailing
    # null character.
    if indexrow['dformat'] == 18:
        strval = strval[1:]
    elif indexrow['dformat'] == 19:
        strval = strval[:-1]

    return strval
    
def read1ByteInts(abi, indexrow):
    if indexrow['fsize'] != 1:
        raise ABIError('Index entry contains an invalid format size for 1-byte integers.')
    
    # see if the data format is signed or unsigned
    if indexrow['dformat'] == 1:
        formatstr = 'B'
    elif indexrow['dformat'] == 2:
        formatstr = 'b'
    else:
        raise ABIError('Index entry contains an invalid data type ID for 1-byte integers.')

    lst = list()
    
    if indexrow['dlen'] <= 4:
        # The actual data are stored in the offset field of the index entry.  Because the offset
        # was read as an unsigned, big-endian integer, the bytes should be in the correct order for
        # the following bit shift operations.
        # First, repack the integer to deal with the possibility of signed integers (shift operations
        # would only return positive values).
        data = pack('>I', indexrow['offset'])
        for cnt in range(0, indexrow['dcnt']):
            val = unpack(formatstr, data[cnt:cnt+1])[0]
            lst.append(val)
    else:
        # get the data from the file
        abi.seek(indexrow['offset'], 0)
        for cnt in range(0, indexrow['dcnt']):
            lst.append(unpack(formatstr, abi.read(1))[0])
    
    if indexrow['dlen'] != len(lst):
        raise ABIDataError(indexrow['dlen'], len(lst))

    return lst
   
def read2ByteInts(abi, indexrow):
    if indexrow['fsize'] != 2:
        raise ABIError('Index entry contains an invalid format size for 2-byte integers.')

    # see if the data format is signed or unsigned
    if indexrow['dformat'] == 3:
        formatstr = '>H'
    elif indexrow['dformat'] == 4:
        formatstr = '>h'
    else:
        raise ABIError('Index entry contains an invalid data type ID for 2-byte integers.')
   
    lst = list()
   
    if indexrow['dlen'] <= 4:
        # The actual data are stored in the offset field of the index entry.  Because the offset
        # was read as an unsigned, big-endian integer, the bytes should be in the correct order for
        # the following operations.
        # First, repack the integer to deal with the possibility of signed integers (shift operations
        # would only return positive values).
        data = pack('>I', indexrow['offset'])
        for cnt in range(0, indexrow['dcnt']):
            val = unpack(formatstr, data[cnt*2:cnt*2+2])[0]
            lst.append(val)
    else:
        # get the data from the file
        abi.seek(indexrow['offset'], 0)
        for cnt in range(0, indexrow['dcnt']):
            lst.append(unpack(formatstr, abi.read(2))[0])
    
    if indexrow['dlen'] != (len(lst) * 2):
        raise ABIDataError(indexrow['dlen'], (len(lst) * 2))
    
    return lst
    
def read4ByteInts(abi, indexrow):
    if indexrow['fsize'] != 4:
        raise ABIError('Index entry contains an invalid format size for 4-byte integers.')
    if indexrow['dformat'] not in (5, 10, 11):
        raise ABIError('Index entry contains an invalid data type ID for 4-byte integers.')
    
    lst = list()
    
    if indexrow['dlen'] == 4:
        # The actual data are stored in the offset field of the index entry.  In the case of 4-byte
        # ints, the offset value is the data value.  It must be repacked, though, to reinterpret it
        # as a signed integer.
        data = pack('>I', indexrow['offset'])
        val = unpack('>i', data)[0]
        lst.append(val)
    else:
        # get the data from the file
        abi.seek(indexrow['offset'], 0)
        for cnt in range(0, indexrow['dcnt']):
            lst.append(unpack('>i', abi.read(4))[0])
 
    if indexrow['dlen'] != (len(lst) * 4):
        raise ABIDataError(indexrow['dlen'], (len(lst) * 4))
    
    return lst

def read4ByteFloats(abi, indexrow):
    if indexrow['fsize'] != 4:
        raise ABIError('Index entry contains an invalid format size for 4-byte floating point numbers.')
    if indexrow['dformat'] != 7:
        raise ABIError('Index entry contains an invalid data type ID for 4-byte floating point numbers.')
    
    lst = list()
   
    if indexrow['dlen'] <= 4:
        # The actual data are stored in the offset field of the index entry.
        data = pack('>I', indexrow['offset'])
        lst.append(unpack('>f', data)[0])
    else:
        # get the data from the file
        abi.seek(indexrow['offset'], 0)
        for cnt in range(0, indexrow['dcnt']):
            lst.append(unpack('>f', abi.read(4))[0])
   
    if indexrow['dlen'] != (len(lst) * 4):
        raise ABIDataError(indexrow['dlen'], (len(lst) * 4))
    
    return lst

def getBaseDataOrder(abiindex):
    # retrieve the "filter wheel order" row from the file index
    rows = getIndexEntriesById(abiindex,'FWO_')
    
    if len(rows) > 1:
        raise ABIError('Found multiple filter wheel order index entries in ABI file.')
    if rows[0]['dlen'] != 4:
        raise ABIError('Incorrect data length for filter wheel order index entry.')
    
    # the data length is only 4 bytes, so the actual data is stored in the offset
    val = rows[0]['offset']
    
    base_order = list()
    
    base_order.append(chr((val >> 24) & 0xff))
    base_order.append(chr((val >> 16) & 0xff))
    base_order.append(chr((val >> 8) & 0xff))
    base_order.append(chr(val & 0xff))
    
    return base_order

if __name__ == '__main__':
    pass
