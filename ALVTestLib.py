
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
track_image(directory) convert a tri file into a dictionary for data manipulation through python
track_to_excel(dictionary, title, save_directory) takes in the dictionary from track_image method, and outputs a user friendly excel file with the data from each tool used in track image.
"""
import os
from scipy.signal import filtfilt
from scipy import integrate
import struct
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
import xlsxwriter

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
                        read_flag = False
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


    def __init__(self, direct,channel_list = []):
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
                        channelname = linetxt[1].rstrip()
                        if channel_list:
                            if channelname not in channel_list:
                                continue
                        self.channel[channelname] = {}
                        
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

class TrackImage:
    """Track image class will create a class that internally has all the data from the TRI files.

    Attributes:

    Sequence Dictionary: Level 1 Key Filename
                                Level 2 Frames
                                        Frame Duration
                                        Spaital Resolution
                                        Frame Rate
                                        Time Offset
    Trajectory Dictionary: Level 1 Name/Title (Whatever the title was for the trajectorgraphy)
                                Level 2 Pandas Dataframe Index is Frames columns (x_world, y_world, x_image, y_image)
    Airbag 2D Dictionary: Level 1 Name/Title
                                Level 2 Tracking
                                    Level 3 (x (list), y (list)) Tuple

                                Level 2 Parameters
                                    Level 3 Parameter Name
                                        Level 4 Pandas Data Frame (x world, y world, x image, y image)

    Methods:
        Velocity Tracking Inputs(Sequence, edge, tracking_direction)
        Threshhold(sequence, edge, track_direction, threshhold)
    """
    def __init__(self, direct):
        filename, file_ext = os.path.splitext(direct)
        location, file = os.path.split(direct)
        self.file = file
        self.loca = location
        file_ext = file_ext.strip('.')
        xml_data = open(direct, 'r').read()
        root = ET.XML(xml_data)
        self.sequence = {}
        self.trajectory = {}
        self.AB2D = {}

        for child in root.iter('Object'):
            if child.get('type') == 'Sequence':
                file_path = child.find('./Content/Filepath')
                filename, file_ext = os.path.splitext(file_path)
                self.sequence[filename] = {}
                self.sequence[filename]['Frames'] = child.find('./Content/FramesCount')
                self.sequence[filename]['Frame Duration'] = child.find('./Content/FrameDuration')
                self.sequence[filename]['Spatial Resolution'] = child.find('./Content/SpatialResolution')
                self.sequence[filename]['Frame Rate'] = int(child.find('./Content/FrameRate'))
                self.sequence[filename]['Time Offset'] = child.find('./Content/TimeOffset')
            elif child.get('type') == 'Trajectory2D':
                title = child.get('name')
                temp_dict = {}
                self.trajectory[title] = {}
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
    #CFC = 5/3*CFC
    #wn = CFC/sample_rate * 2
    #b,a = butter(2,wn,'low')
    y = filtfilt(b, a, data)
    return y


def linear_calc(data):
    linear_dict = {}
    return linear_dict


def filter_analysis(data, sample_rate, cfc_list=None):
    filter_analysis_dict = {}
    if cfc_list is None:
        cfc_list = [0, 60, 180, 1000, 5000]
    for i in cfc_list:
        if i == 0:
            filter_analysis_dict['RAW'] = filter_processing(data, i, sample_rate)
        else:
            filter_analysis_dict[str(i)] = filter_processing(data, i, sample_rate)
    return filter_analysis_dict


def nij_calc(dummy_type, m_y, f_x, f_z):

    """Purpose of this method is to calculate the Neck Injury Criterion per FMVSS208
       The inputs for this are the dummy type. Right now 6 Year Old and 3 Year Old are the only dummys available
       3 Year Old = 0
       6 Year Old = 1
       5th In position = 2
       50th In position = 3
       Moment Y Channel
       Force in X direction
       Force in Z Direction
       Then using the methods defined by FMVSS208 the Nij is calculated"""
    # Setting up critical forces dictionary. Depending if the force is negative or positive depends on what the
    fcz = {0: {'Tension': 2120, 'Compression': -2120}, 1: {'Tension': 2800, 'Compression': -2800},
           2: {'Tension': 4287, 'Compression': -3880}, 3: {'Tension': 6806, 'Compression': -6160}}
    mcy = {0: {'Flexion': 68, 'Extension': -27}, 1: {'Flexion': 93, 'Extension': -37},
           2: {'Flexion': 155, 'Extension': -67}, 3: {'Flexion': 310, 'Extension': -135}}
    d = {0: 0, 1: 0.01778, 2: 0.01778, 3: 0.01778}
    # Data needs to be filtered at CFC 600
    # m_y_filter = filter_processing(M_y, 600, sample_rate)
    # f_z_filter = filter_processing(F_z, 600, sample_rate)
    # f_x_filter = filter_processing(F_x, 600, sample_rate)
    Nij = []
    for i in range(0, len(m_y)):
        moc_y = m_y[i] - (d[dummy_type] * f_x[i])
        if f_z[i] < 0:
            a = f_z[i] / fcz[dummy_type]['Compression']
        else:
            a = f_z[i] / fcz[dummy_type]['Tension']
        if m_y[i] < 0:
            b = moc_y/mcy[dummy_type]['Extension']
        else:
            b = moc_y/mcy[dummy_type]['Flexion']
        Nij.append(a+b)
    return Nij


def hic_calc(a_x, a_y, a_z, time, delta_hic):
    # Purpose of this method is to calculate HIC (Head Injury Criterion). This is done using integration over specific
    # curves. Defined by delta_hic (unlimited, 15ms, or 30ms)
    # time is used to calculate the sample rate
    # Calculate sample rate
    sample = int(sample_rate(time))
    # filter acceleration
    a_x_filter = filter_processing(a_x, 1000, sample)
    a_y_filter = filter_processing(a_y, 1000, sample)
    a_z_filter = filter_processing(a_z, 1000, sample)
    # Convert acceleration to m/s^2
    # a_x_filter = a_x_filter * 9.81
    # a_y_filter = a_y_filter * 9.81
    # a_z_filter = a_z_filter * 9.81
    # Calculate resultant of acceleration
    a_r = np.sqrt((np.power(a_x_filter, 2) + np.power(a_y_filter, 2) + np.power(a_z_filter, 2)))
    # Null data
    null_idx = np.where(time >= 0)
    null_idx = null_idx[0][0]
    idx2 = np.where(time >= 0.2)[0][0]
    a_r_null = a_r[null_idx:idx2]
    time_null = time[null_idx:idx2]
    # Calculate the number of data points in delta_hic convert delta_hic to seconds.
    delta_hic = delta_hic/1000
    time_len = int(delta_hic/(time[1]-time[0]))
    hic_list = []
    for i in time[:-time_len]:
        t1 = np.where(time == i)[0][0]
        t2 = t1 + time_len
        time_array = time_null[t1:t2]
        accel_array = a_r_null[t1:t2]
        a = time[t2] - time[t1]
        b = 1 / a
        hic_a = (np.power((b * np.cumtrapz(accel_array, time_array)), 2.5)*a)
        hic_list.append(hic_a)
        del time_array
        del accel_array
    hic = max(hic_list)
    return hic, a_r


def linear_calc(data):
    linear_dict = {}

    return linear_dict


def track_image(directory):
    # Method is for taking track image TRI files and outputting the Trajectory2D and Airbag2D into a dictionary for
    # easier data manipulation with Python.
    xml_data = open(directory, 'r').read()
    root = ET.XML(xml_data)
    frames = {}
    frames_2d = {}
    track_data_x = {}
    track_data_y = {}
    ab2_d_x = {}
    ab2_d_y = {}
    parameters = {}
    for child in root.iter('Object'):
        # Find info for current file Frames per Second, Frames, time step, offset time
        if child.get('type') == 'Sequence':
            time_off = child.find('./Content/TimeOffset')
            frame_rate = child.find('./Content/FrameRate')
            frame_count = child.find('./Content/FramesCount')
            time_step = child.find('./Content/FrameDuration')
        # Find Trajectory2D data inside xml file
        if child.get('type') == 'Trajectory2D':
            track_data_x[child.get('name')] = []
            track_data_y[child.get('name')] = []
            frames[child.get('name')] = []
            for subchild in child.iter('Content'):
                for data in subchild.iter('xworlds'):
                    track_data_x[child.get('name')] = [float(i) for i in data.text.split(' ')]
                for data in subchild.iter('yworlds'):
                    track_data_y[child.get('name')] = [float(i) for i in data.text.split(' ')]
                for data in subchild.iter('images'):
                    frames[child.get('name')] = [int(i) for i in data.text.split(' ')]
        # Find Airbag2D data
        if child.get('type') == 'Airbag2D':
            parameters[child.get('name')] = {}
            str_data_x = []
            str_data_y = []
            l_edge = {}
            r_edge = {}
            u_edge = {}
            d_edge = {}
            cg = {}
            ab2_d_x[child.get('name')] = []
            ab2_d_y[child.get('name')] = []
            frames_2d[child.get('name')] = []
            l_edge['X'] = [float(i) for i in child.find('.//Parameters//LeftEdges//xworlds').text.split(' ')]
            r_edge['X'] = [float(i) for i in child.find('.//Parameters//RightEdges//xworlds').text.split(' ')]
            u_edge['X'] = [float(i) for i in child.find('.//Parameters//UpperEdges//xworlds').text.split(' ')]
            d_edge['X'] = [float(i) for i in child.find('.//Parameters//LowerEdges//xworlds').text.split(' ')]
            l_edge['Y'] = [float(i) for i in child.find('.//Parameters//LeftEdges//yworlds').text.split(' ')]
            r_edge['Y'] = [float(i) for i in child.find('.//Parameters//RightEdges//yworlds').text.split(' ')]
            u_edge['Y'] = [float(i) for i in child.find('.//Parameters//UpperEdges//yworlds').text.split(' ')]
            d_edge['Y'] = [float(i) for i in child.find('.//Parameters//LowerEdges//yworlds').text.split(' ')]
            cg['Y'] = [float(i) for i in child.find('.//Parameters//Centers//yworlds').text.split(' ')]
            cg['X'] = [float(i) for i in child.find('.//Parameters//Centers//xworlds').text.split(' ')]
            parameters[child.get('name')] = {'Left Edge': l_edge, 'Right Edge': r_edge, 'Upper Edge': u_edge, 'Lower Edge': d_edge, 'CG': cg}

            for subchild in child.find('.//Curves'):
                str_data_x.append(subchild.find('.//xworlds').text.split(' '))
                str_data_y.append(subchild.find('.//yworlds').text.split(' '))
                frames_2d[child.get('name')].append(int(subchild.get('frame')))
            for i in range(0, len(str_data_x)):
                ab2_d_x[child.get('name')].append([float(j) for j in str_data_x[i]])
                ab2_d_y[child.get('name')].append([float(k) for k in str_data_y[i]])
    # Calculate time array for file
    offset = float(time_off.text)/1000
    delta = float(time_step.text)/1000
    end_time = int(frame_count.text)/int(float(frame_rate.text))*1000
    time = np.arange(offset, end_time, delta)
    track_dict = {"Trajectory X": track_data_x, "Trajectory Y": track_data_y, "Trajectory Frame": frames,
                  "AB2D X": ab2_d_x, "AB2D Y": ab2_d_y, "AB2D Frame": frames_2d,
                  "Frame rate": int(float(frame_rate.text)), 'Time Offset': float(time_off.text)/1000,
                  "Number of Frames": int(frame_count.text), "Time_Step": float(time_step.text)/1000,
                  "Time (msec)": time, "Parameters": parameters }
    return track_dict


def track_to_excel(dictionary, title, save_direct):
    save_direct = save_direct + '\\'+title + '.xlsx'
    workbook = xlsxwriter.Workbook(save_direct)
    cell_format = workbook.add_format({'bold': 1, 'align': 'center', 'valign': 'vcenter','border': 1})
    merge_format = workbook.add_format({'bold': 1, 'align': 'center', 'valign': 'vcenter','border': 1,
                                        'rotation': 90})
    for key in dictionary.keys():
        if 'Trajectory X' in key and dictionary[key]:
            frame_start = []
            frame_end = []
            for frame_key in dictionary['Trajectory Frame'].keys():
                frame_start.append(dictionary['Trajectory Frame'][frame_key][0])
                frame_end.append(dictionary['Trajectory Frame'][frame_key][-1])
            frame_list = list(range(1, dictionary['Number of Frames'] + 1))
            start_idx = frame_list.index(min(frame_start))
            end_idx = frame_list.index(max(frame_end))
            traj_work = workbook.add_worksheet('2D Trajectory')
            traj_work.merge_range('A1:A2', 'Frames', cell_format)
            traj_work.merge_range('B1:B2', 'Time', cell_format)
            row = 2
            col = 0
            for i in frame_list[start_idx:end_idx+1]:
                time_idx = frame_list.index(i)
                traj_work.write(row, col, i)
                traj_work.write(row, col+1, dictionary['Time (msec)'][time_idx])
                row += 1
            row = 0
            col = 2
            for sub_key in dictionary[key].keys():
                traj_work.merge_range(row, col, row, col + 1, sub_key, cell_format)
                traj_work.write(1, col, 'X', cell_format)
                traj_work.write(1, col+1, 'Y', cell_format)
                count = 0
                for i in dictionary[key][sub_key]:
                    data_row = dictionary['Trajectory Frame'][sub_key][count] - start_idx + 1
                    traj_work.write(data_row, col, i)
                    traj_work.write(data_row, col+1, dictionary['Trajectory Y'][sub_key][count])
                    count += 1
                col += 2

        if 'AB2D X' in key and dictionary[key]:
            row = 0
            col = 0
            frame_list = list(range(1,dictionary['Number of Frames']+1))
            data_points_length = {}
            ab2D_work = workbook.add_worksheet('Airbag 2D')
            row_position = 0
            for sub_key in dictionary[key].keys():
                #find max length of list data to merge rows
                row = row_position
                col = 0
                ab2D_work.write(row,col,'Frame',cell_format)
                row += 1
                ab2D_work.write(row,col,'Time',cell_format)
                row += 1
                array_length = [len(i) for i in dictionary[key][sub_key]]
                merge_row = max(array_length)
                ab2D_work.merge_range(row,0,row + merge_row,0, sub_key, merge_format)
                col += 1
                count = 0
                list_count = 0
                for i in dictionary[key][sub_key]: 
                    row = row_position  
                    ab2D_work.merge_range(row,col,row, col + 1,dictionary['AB2D Frame'][sub_key][count],cell_format)
                    row += 1
                    time = dictionary['Time (msec)'][dictionary['AB2D Frame'][sub_key][count]]
                    ab2D_work.merge_range(row,col,row,col+1, time ,cell_format)
                    row += 1
                    ab2D_work.write(row,col,'X',cell_format)
                    ab2D_work.write(row,col+1,'Y',cell_format)
                    row += 1
                    matrix_count = 0
                    for j in i:
                        ab2D_work.write(row,col,j)
                        ab2D_work.write(row,col+1,dictionary['AB2D Y'][sub_key][list_count][matrix_count])
                        matrix_count += 1
                        row += 1
                    row = row
                    matrix_count = 0
                    list_count += 1
                    col += 2
                    count += 1
                row_position = row_position + merge_row + 4
        if 'Parameters' in key:
            param = workbook.add_worksheet('Parameters')
            col = 0
            for sub_key in dictionary[key].keys():
                row = 0
                param.merge_range(row,col,row,col+11,sub_key,cell_format)
                row += 1
                param.merge_range(row, col, row + 1, col, 'Frame', cell_format)
                param.merge_range(row, col+1, row + 1, col + 1, 'Time (msec)', cell_format)
                row += 2
                for i in dictionary['AB2D Frame'][sub_key]:
                    param.write(row,col,i)
                    param.write(row,col+1, dictionary['Time (msec)'][i - 1])
                    row += 1
                col += 2
                row = 1
                for type_key in dictionary[key][sub_key].keys():
                    param.merge_range(row, col, row, col + 1, type_key, cell_format)
                    row += 1
                    param.write(row,col,'X',cell_format)
                    param.write(row,col+1,'Y',cell_format)
                    row += 1
                    count = 0
                    for i in dictionary[key][sub_key][type_key]['X']:
                        param.write(row,col,i)
                        param.write(row,col+1,dictionary[key][sub_key][type_key]['Y'][count])
                        row += 1
                        count += 1
                    col += 2
                    row = 1

    workbook.close()


def sample_rate(time):

    """Sample Rate Calculator

    Args:
        time (Pandas Series): Need to input a series of the time values from the calculated data. data should be converted to msec if its not in those units. 

    Returns:
        [float]: [calculated sample rate by 1/(max time value - min time value) then  multiplying it by the number of samples]
    """
    sam_rate = 1/(time.max()/1000 - time.min()/1000) * (len(time)-1)
    return sam_rate


def contact_time(accel_list, time_list, negative_accel = True):
    """Calculate contact time

    Args:
        accel_list (list of acceleration data): [list or series of acceleration]
        time_list (float units (msec)): [time list/pd.series to find contact time]
        negative_accel (bool, optional): [If acceleration is sloping positive this value should be false]. Defaults to True.

    Returns:
        [float, integer]: [returns the time at which the impactor makes contact with the PAB based off of ATS201]
    """
    samp_rate = sample_rate(time_list)
    #calculate number of samples per msec. For finding rate of change greater than 1ms
    sample = sample_rate/1000
    sample = int(sample_rate + (sample_rate/2))
    count = 0
    for i in accel_list:
        if negative_accel == True and i <= -1:
            accel_rate = (accel_list[count + sample] - i)/(time_list[count + sample]/1000 - time_list[count]/1000)
            if accel_rate >= 100.00:
                contact_time = time_list[count]
                time_index = count
                break
        elif negative_accel == False and i >= 1:
            accel_rate = (accel_list[count + sample] - i)/(time_list[count + sample]/1000 - time_list[count]/1000)
            if accel_rate >= 100.00:
                contact_time = time_list[count]
                time_index = count
                break            
    return contact_time, time_index


def calc_force(accel, impactor_mass):
    """Calculate force

    Args:
        accel (list or series): Acceleration array read from diadem or Impax file
        impactor_mass (float): Mass of the impactor. This is needed to calculate the force

    Returns:
        : Pandas series with calculated force in (N)
    """
    force = []
    for i in accel:
        force.append(i * 9.81 * impactor_mass)
    force_df = pd.Series(force)
    return force_df


def velocity(accel, time):
    """calculate the velocity of the test using acceleration and time and integrating the curve. Make sure the the acceleration is in m/s^2 and the time is in msec.

    Args:
        accel (array): [acceleration array in m/s^2]
        time (array): [time array in msec]
    """
    vel = integrate.cumtrapz(accel, time/1000)
    return vel


def kinetic_energy(velocity, impactor_mass):
    """Calculate kinetic energy of the impactor

    Args:
        velocity ([Array]): [velocity, either calculated with post processed files, or calculated using velocity method in this library]
        impactor_mass ([float]): [mass of the impactor in kg]
    """
    ke = []
    for i in velocity:
        ke.append(0.5 * impactor_mass * i**2)
    return ke


def null_data(time, data, null_at = 0.0):
    """Null data at a given time point (default is 0msec)

    Returns:
        data structures: [Pandas data frame that outputs the nulled data]
    """
    idx = time[time == null_at].index[0]
    trun_time = time.truncate(before = idx)
    null_data = data - data[idx]
    null_data = null_data.truncate(before = idx)
    return null_data


def ride_down(acceleration, time, negative_curve = True):
    """Calculates the ride down value for the acceleration curve. The ride down is the slope of the acceleration from contact to max acceleration

    Args:
        acceleration (array): [acceleration array units are not taken into account for this. ]
        time ([array]): [time units are not taken into account]
        negative_curve (BOOLEAN) = [if value is True (default) the acceleration curve is negative so the max value will actually the min value.]

    Returns:
        [Ride down]: [Takes max acceleration, and accleration at contact time and outputs the slope of that curve. ]
    """
    if negative_curve == True:
        con_time, con_idx = contact_time(acceleration, time)
        peak_accel = acceleration.min()
        peak_time = acceleration[acceleration == peak_accel].index[0]
        contact_accel = acceleration[con_idx]
        ride_dwn_value = (peak_accel - contact_accel)/(peak_time - con_time)
    elif negative_curve == False:
        con_time, con_idx = contact_time(acceleration, time)
        peak_accel = acceleration.max()
        peak_time = acceleration[acceleration == peak_accel].index[0]
        contact_accel = acceleration[con_idx]
        ride_dwn_value = (peak_accel - contact_accel)/(peak_time - con_time)
    
    return ride_dwn_value


