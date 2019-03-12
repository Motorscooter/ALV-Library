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
            
        with open(direct, 'rb') as fup:
            read_flag = True
            data_pts_flag = True
            while read_flag:
                line = fup.readline()
                if bytes('VERTSCALE', 'latin-1') in line:
                    self.vScale = float(line[10:-1].rstrip())
                if bytes('VERTOFFSET', 'latin-1') in line:
                    self.vOffset = float(line[11:-1].rstrip())
                if bytes('VERTUNITS', 'latin-1') in line:
                    self.vUnits = line[10:-1].rstrip()
                if bytes('HORZSCALE', 'latin-1') in line:
                    self.hScale = float(line[10:-1].rstrip())
                if bytes('HORZOFFSET', 'latin-1') in line:
                    self.hOffset = float(line[11:-1].rstrip())
                if bytes('HORZUNITS', 'latin-1') in line:
                    self.hUnits = line[10:-1].rstrip()
                if bytes('HUNITPERSEC', 'latin-1') in line:
                    self.hUPS = float(line[12:-1].rstrip())
                if bytes('RECLEN', 'latin-1') in line and data_pts_flag:
                    data_pts_flag = False
                    self.numDataPts = int(line[7:-1].rstrip())
                if bytes('CLOCK_RATE', 'latin-1') in line:
                    self.sampleRate = float(line[11:-1].rstrip())
                if bytes('DDR_RATE', 'latin-1') in line:
                    self.sampleRate = float(line[9:-1].rstrip())
                if bytes('InflatorID', 'latin-1') in line:
                    self.sampleID = (line[11:-1].rstrip()).decode()
                if bytes('DATATYPE', 'latin-1') in line:
                    self.data_type = (line[9:-1].rstrip()).decode('latin-1')
                if bytes('BuildSheet', 'latin-1') in line:
                    self.BS = line[11:-1].rstrip()
                    read_flag = False
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


class DiademRead:
    x = 0
 

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
