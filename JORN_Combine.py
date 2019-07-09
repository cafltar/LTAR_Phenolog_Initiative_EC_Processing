# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 14:14:57 2019

@author: eric.russell
Laboratory for Atmospheric Research
Dept. Civil and Environ. Engineering
Washington State University
eric.s.russell@wsu.edu
"""

import pandas as pd
import datetime

#Read in the flux data from the EC dataset
flux = pd.read_csv(r'C:\Users\Eric\Desktop\LTAR\LTAR_National_Projects\PhenologyInitiative\Konza\KONA_Neon_EC_Data.csv', header = 0, index_col = 'TimeStamp')
flux.index = pd.to_datetime(flux.index)
#Read in the meteorology data from the non-flux dataset
extra = pd.read_csv(r'C:\Users\Eric\Desktop\LTAR\LTAR_National_Projects\PhenologyInitiative\Konza\KONA_Non_Flux_Data_Test.csv', header = 0, index_col = 'TimeStamp')
extra.index = pd.to_datetime(extra.index)
#Concat the data into one dataframe from the two sources
f = pd.concat([flux,extra], axis = 1)
f = f.drop(columns = ['endDateTime']) # Drop erroneous column
#Output full combined dataset
f.to_csv(r'C:\Users\Eric\Desktop\LTAR\LTAR_National_Projects\PhenologyInitiative\Konza\KONA_NEON_Extract_2018.csv')
#Format timestamp into the AmeriFlux format to match other datasets used
f['TIMESTAMP_END'] = f.index.shift(1, '30T')
f['TIMESTAMP_START'] = f.index.shift(0, '30T')    
f= f.drop(f.index[0])
f['TIMESTAMP_START']= f.TIMESTAMP_START.map(lambda x: datetime.datetime.strftime(x, '%Y%m%d%H%M'))
f['TIMESTAMP_END']= f.TIMESTAMP_END.map(lambda x: datetime.datetime.strftime(x, '%Y%m%d%H%M'))
#Read-in columns to keep from separate list
Cols = pd.read_csv(r'C:\Users\Eric\Desktop\LTAR\LTAR_National_Projects\PhenologyInitiative\Jornada\NEON\Cols_Save_NEON_Jorn.csv',header = 0)
fcols = f.columns
#Rename columns
for k in range (0,len(Cols)):
        if Cols['Old'][k] in fcols:
            qn = Cols['Old'][k] == fcols
            f= f.rename(columns={fcols[qn][0]:Cols['Keep'][k]})
s = Cols.isin(f.columns)
#Drop columns not in the needed data list
f = f.drop(f.columns[~s['AF']],axis = 1)
#Output formatted dataset
f = f['01-01-2018':'12-31-2018']
f.to_csv(r'C:\Users\Eric\Desktop\LTAR\LTAR_National_Projects\PhenologyInitiative\Konza\EC_LTAR_KONA_neon_20180101_20181231.csv')
