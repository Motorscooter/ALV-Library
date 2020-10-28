# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 19:52:50 2018

@author: scott.downard
ALVTestLib
This library is designed to help users create scripts for analyzing test bay data
12/03/2018
ImpaxRead(self,direct)
ImpaxRead() is a class that creates an object based on a directory that points to 
            an IMPAX file.
    ImpaxRead.file()
            Output: String
            String related to input directory
    ImpaxRead.filetype()
            Output: String
            String related to file type
            File Types: Force, Acceleration, Displacement, Pressure, Primary/Secondary Voltage
                        Primary/Secondary Current
    ImpaxRead.Datapoints()
            Output: Integer
            Number of data points inside the input directory related to the vertical
            datapoint (filetype)
    ImpaxRead.Sample()
            Output: Integer
            Sample rate that was used for testing
    ImpaxRead.BuildSheet()
            Output: String
            Sample ID for file that was input
    ImpaxRead.vUnits()
            Output: String
            Units of vertical data
    ImpaxRead.hUnits()
            Output: String
            Units of horizontal data (majority of the time it is msec)
    ImpaxRead.vData()
            Output: Tuple
            An array of data points calculated from .verOff object, and verScale object
    ImpaxRead.hData()
            Output: Tuple
            An array of data points calculated from horOff object, and horScale object
    ImpaxRead.vOff()
            Output: Double
            The vertical offset of the data from 0
    ImpaxRead.hOff()
            Output: Double
            The horizontal offset of the data from 0
    ImpaxRead.vScale()
            Output: Double
            The vertical scale for the data
    ImpaxRead.hScale()
            Output: Double
            The horizontal scale for the data
    ImpaxRead.vdataPoints(verData,verScale,verOff)
            Output: Numpy Array
            calculates the converted data into correct data
    ImpaxRead.hdataPoints(hData,hscale,hoff)
            Output: Numpy Array
            calculates the horizontal data
DiademRead(self,direct) is a class that creates an object with an directory input that points
            a Diadem file
TankRead(self,direct) is a class that creates an object with a directory input that points
           a tank data file
"""
import os
from scipy.signal import filtfilt
import struct
import numpy as np
import pandas as pd


class ImpaxRead:

    def __init__(self, direct):
        self.file = direct
        filename, file_ext = os.path.splitext(direct)
        file_ext = file_ext.strip('.') 
        try:
            file_ext = int(file_ext)
        except ValueError:
            try:
                file_ext == '.aift'
            except ValueError:
                print('File is not an IMPAX file please use another file')
        if "CU0" in filename:
            if "O01" in filename:    
                self.file_type = "Primary Current"
            elif "O02" in filename:
                self.file_type = "Secondary Current"
        elif "VO0" in filename:
            if "O01" in filename:   
                self.file_type = "Primary Voltage"
            elif "O02" in filename:
                self.file_type = "Secondary Voltage"
        elif "ACX" in filename:
            self.file_type = "Acceleration"
        elif "DSX" in filename:
            self.file_type = "Displacement"
        elif "PR0" in filename:
            self.file_type = "Pressure"
        elif "Bag" in filename:
            self.file_type = "Pressure"
        elif "FOX" in filename:
            self.file_type = "X Force"
        elif "FOY" in filename:
            self.file_type = "Y Force"
        elif "FOZ" in filename:
            self.file_type = "Z Force"
        elif "CMB" in filename:
            if "CMB1" in filename:
                self.file_type = "Combustion 01"
            elif "CMB2" in filename:
                self.file_type = "Combustion 02"
        elif "VLT" in filename:
            if "VLT1" in filename:
                self.file_type = "Primary Voltage"
            elif "VLT2" in filename:
                self.file_type = "Secondary Voltage"
        elif "CUR" in filename:
            if "CUR1" in filename:
                self.file_type = "Primary Current"
            elif "CUR2" in filename:
                self.file_type = "Secondary Current"
        elif "TNK" in filename:
            if "TNK1" in filename:
                self.file_type = "Tank 01"
            elif "TNK2" in filename:
                self.file_type = "Tank 02"
        else:
            self.file_type = ""
        self.header = {"Filename" : filename}
        with open(direct, encoding="latin1") as fup:
            read_flag = True
            data_pts_flag = True
            while read_flag:
                line = fup.readline()
                if "=" in line:
                    line = line.strip()
                    split_line = line.split("=")
                    self.header[split_line[0]] = self.parse_header_data(split_line[1].strip())
                    if 'BuildSheet' in line:
                        read_flag = True
            fup.close()
            self.BS = str(self.header['BuildSheet'])
            self.sampleID = str(self.header['Sample'])
            self.vScale = self.header['VERTSCALE']
            self.vOffset = self.header['VERTOFFSET']
            self.vUnits = self.header['VERTUNITS']
            self.hOffset = self.header['HORZOFFSET']
            self.hScale = self.header['HORZSCALE']
            self.hUPS = self.header['HUNITPERSEC']
            self.numDataPts = int(self.header['RECLEN'])
            self.data_type = self.header['DATATYPE']
            self.test_group = self.header['TestGroup']
            if 'CLOCK_RATE' in self.header:
                self.sampleRate = self.header['CLOCK_RATE']
            elif 'DDR_RATE' in self.header:
                self.sampleRate = self.header['DDR_RATE']
            fup.close()

    def vert_data_points(self, sae_filter=0):
        with open(self.file, 'rb') as fup:
            data_type_dict = dict(DOUBLE=['d', 8], SHORT=['h', 2], LONG=['l', 4])
            size = fup.read()
            vert_data = struct.unpack(data_type_dict[self.data_type][0]*self.numDataPts,
                                      size[len(size)-(self.numDataPts*data_type_dict[self.data_type][1]):])
            vert_data = np.asarray(vert_data)
            vert_data.astype(float)
            vert_data = vert_data * self.vScale + self.vOffset
        fup.close()
        if sae_filter < 0:
            sae_filter = 0
            print('SAE filter cannot be negative filter was set to 0 (RAW). Please enter number greater than 0')           
        if sae_filter == 0:
            return vert_data
        else:
            vert_data = filter_processing(vert_data, sae_filter, self.sampleRate)
        return vert_data
        
    def h_data_points(self):
        time = np.arange(self.hOffset, ((self.numDataPts*self.hScale)+self.hOffset), self.hScale)
        return time

    def parse_header_data(self, dataString): # will turn a string into a float if possible.
        try:
            data = float(dataString.strip())
            return data
        except ValueError:
            return dataString

class DiademRead:
    """
    1, OS
    2, DiaDem Version
    105, Time
    110, Test Date
    Local Channel Header
    200, Channel Name
    201, Description
    202, Units
    210, Channel Type
    211, Source File
    213, Channel
    214, Data Type
    220, Number of Values for Data
    221, Pointer to first value in data
    240, Start Value offset
    241, Scale factor for Offset
    250, Min Value
    251, Max Value
    252, No Values
    253, Monotonous
    260, Display data type
    261, Visibility
    262, Display Order
    263, Content Width
    265, Property Width """


    def __init__(self, direct):
        filename, file_ext = os.path.splitext(direct)
        location, file = os.path.split(direct)
        self.file = file
        self.loca = location
        file_ext = file_ext.strip('.')
     #Enter error handling code
     #Read .DAT file to build dictionary with needed data from above
        with open(direct,'r') as fup:
            self.fileName = filename
            self.channel = {}
            channel_header = False
            for line in fup:
                if channel_header:
                    linetxt = line.split(',')
                    if linetxt[0] == '200':
                        self.channel[linetxt[1].rstrip()] = {}
                        channelname = linetxt[1].rstrip()
                    elif linetxt[0] == '201':
                        self.channel[channelname]['Description'] = linetxt[1].rstrip()
                    elif linetxt[0] == '202':
                        self.channel[channelname]['Units'] = linetxt[1].rstrip()
                    elif linetxt[0] == '211':
                        self.channel[channelname]['R32'] = (self.loca + '\\' + linetxt[1]).rstrip()
                    elif linetxt[0] == '214':
                        self.channel[channelname]['Data Type'] = linetxt[1].rstrip()
                    elif linetxt[0] == '220':
                        self.channel[channelname]['NumDatPts'] = int(linetxt[1].rstrip())
                    elif linetxt[0] == '221':
                        self.channel[channelname]['DataLoc'] = int(linetxt[1].rstrip())
                    elif linetxt[0] == '240':
                        self.channel[channelname]['OffSet'] = float(linetxt[1].rstrip())
                    elif linetxt[0] == '241':
                        self.channel[channelname]['Scale'] = float(linetxt[1].rstrip())
                    elif linetxt[0] == '250':
                        self.channel[channelname]['Min'] = float(linetxt[1].rstrip())
                    elif linetxt[0] == '251':
                        self.channel[channelname]['Max'] = float(linetxt[1].rstrip())

                if 'BEGINCHANNELHEADER' in line:
                    channel_header = True
                elif 'ENDCHANNGELHEADER' in line:
                    channel_header = False
                    channelname = None
            fup.close()
        data_type_dict = dict(REAL32=['f', 4], SHORT=['h', 2], LONG=['l', 4])
        key_file = list(self.channel.keys())[-1]
        self.data = pd.DataFrame()
        with open(self.channel[key_file]['R32'], 'rb') as fup:
            size = fup.read()
            data_array = struct.unpack('f' * int(len(size) / 4), size)
            for key in self.channel.keys():
                raw_data = data_array[self.channel[key]['DataLoc']-1:self.channel[key]['NumDatPts']+self.channel[key]['DataLoc']-1]
                raw_array = np.asarray(raw_data)
                float_array = raw_array.astype(float)
                calc_data = float_array * self.channel[key]['Scale'] + self.channel[key]['OffSet']
                self.data[key] = pd.Series(calc_data)
                fup.close()

def filter_processing(data, cfc, sample_rate):
    # calculate sample rate
    # Filter RAW data using J211 SAE Filtering
    sample_rate = int(sample_rate)
    t = 1/sample_rate
    wd = 2 * np.pi * cfc * 1.25*5.0/3.0
    x = wd * t/2
    wa = np.sin(x)/np.cos(x)
    a0 = wa**2.0/(1.0+np.sqrt(2.0)*wa+(wa**2.0))
    a1 = 2.0*a0
    a2 = a0
    b0 = 1
    b1 = (-2.0*((wa**2.0)-1.0))/(1.0+np.sqrt(2.0)*wa+(wa**2.0))
    b2 = (-1.0+np.sqrt(2.0)*wa-(wa**2.0))/(1.0+np.sqrt(2.0)*wa+(wa**2.0))
    b = [a0, a1, a2]
    a = [b0, -1*b1, -1*b2]
#    CFC = 5/3*CFC
#    wn = CFC/sample_rate * 2
#    b,a = butter(2,wn,'low')
    y = filtfilt(b, a, data)
    return y

# def linear_analysis(impactor_distance, dist_data, accel_data, time):
#     bottomOut = False
#     max_displacement = max(dist_data)
#     if impactor_distance < max_displacement:
#         bottomOut = True
#     if bottomOut:
#
#     max_accel = min(accel_data)


def filter_analysis(data, sample_rate, cfc_list = None ):
    filter_analysis_dict = {}
    if cfc_list is None:
        cfc_list = [0,60,180,1000,5000]
    for i in cfc_list:
        if i == 0:
            filter_analysis_dict['RAW'] = filter_processing(data, i, sample_rate)
        else:
            filter_analysis_dict[str(i)] = filter_processing(data, i, sample_rate)
    return filter_analysis_dict

