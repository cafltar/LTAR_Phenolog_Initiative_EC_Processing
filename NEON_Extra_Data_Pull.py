# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 15:08:15 2018

@author: Eric Russell
Dept. Civil and Environmental Engineering
Washington State University
eric.s.russell@wsu.edu


"""

import pandas as pd
import glob
import numpy as np

def indx_fill(df):   
    df.index = pd.to_datetime(df.index)
        # Sort index in case it came in out of order, a possibility depending on filenames and naming scheme
    df = df.sort_index()
        # Remove any duplicate times, can occur if files from mixed sources and have overlapping endpoints
    df = df[~df.index.duplicated(keep='first')]
        # Fill in missing times due to tower being down and pad dataframe to midnight of the first and last day
    idx = pd.date_range(df.index[0].floor('D'),df.index[len(df.index)-1].ceil('D'),freq = '30min')
    df = df.reindex(idx, fill_value=np.NaN)
    return df
# Function to read in files assuming a consistent header and column structure
def Fast_Read(filenames):
    if len(filenames) == 0:
        print('No Files in directory, check the path name.')
        return  # 'exit' function and return error
    else:
        #Initialize dataframe used within function
        Final = [];Final = pd.DataFrame(Final)
        for k in range (0,len(filenames)):
            #Read in data and concat to one dataframe; no processing until data all read in
            df = pd.read_csv(filenames[k],index_col = 'startDateTime',header= 0,low_memory=False)
            Final = pd.concat([Final,df])
        # Convert time index
        Out = indx_fill(Final)
    return Out # Return dataframe to main function.    

#Columns from the different *.csv files that are needed from the different file sets
Cols = pd.read_csv('F:\\LTAR\\LTAR_National_Projects\\PhenologyInitiative\\NEON_Met_Data_Cols.csv',header=0)
#Read-in the precipitation data
fnames =glob.glob('F:\\LTAR\\LTAR_National_Projects\\PhenologyInitiative\\Jornada\\NEON\\NEON_precipitation\\2017\\*.csv')
Precip = Fast_Read(fnames)
#Read-in the PAR data
fnames = glob.glob('F:\\LTAR\\LTAR_National_Projects\\PhenologyInitiative\\Jornada\\NEON\\NEON_par\\2017\\*.csv')
par = Fast_Read(fnames)
#Read-in the net radiation data
fnames =glob.glob('F:\\LTAR\\LTAR_National_Projects\\PhenologyInitiative\\Jornada\\NEON\\NEON_rad-net\\2017\\*.csv')
rad_net = Fast_Read(fnames)
#Read-in the relative humidity data
fnames =glob.glob('F:\\LTAR\\LTAR_National_Projects\\PhenologyInitiative\\Jornada\\NEON\\NEON_rel-humidity\\2017\\*.csv')
RH = Fast_Read(fnames)
#Create index for the dataframe from the relative humidity files
data = []; data= pd.DataFrame(data, RH.index)
#Concat the different files into one dataframe will all column headers
data = pd.concat([Precip,par,rad_net,RH],axis = 1, sort = True)
#Drop columns not wanted based on the initial column list read-in
cls = data.columns 
s = cls.isin(Cols['Keep'])
data = data.drop(data[cls[~s]],axis = 1)
#Output data to desired location
data.to_csv('F:\\LTAR\\LTAR_National_Projects\\PhenologyInitiative\\Jornada\\NEON\\Non_Flux_Data_Test.csv')