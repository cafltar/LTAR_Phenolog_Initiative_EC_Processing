# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 14:32:20 2019

@author: Eric Russell
Laboratory for Atmospheric Research
Dept. Civil and Environmental Engineering
Washington State University
eric.s.russell@wsu.edu

"""

import pandas as pd
import calendar
import datetime
import os
os.chdir(r'C:\Users\Eric\Documents\GitHub\LTAR_Phenology_Data_Processing')       
import LTAR_Pheno_QC_Functions as LLT
#Function to convert DoY into more usable format; not quite ISO8601
def JulianDate_to_MMDDYYY(y,jd):
    month = 1
    while jd - calendar.monthrange(y,month)[1] > 0 and month <= 12:
        jd = jd - calendar.monthrange(y,month)[1]
        month = month + 1
    return month,jd,y
#Read in data; UCB data was from the AmeriFlux output option from EddyPro and is inconsistent with the needed format
datae = pd.read_csv(r'C:\Users\Eric\Desktop\LTAR\LTAR_National_Projects\PhenologyInitiative\EC Data\Processed\Unprocessed\EC_LTAR_UCB_hawbeckereddy_201700101_20171231_redone.csv',header=0, skiprows = [1])
Y = datae['YEAR'];D = datae['DOY'];H = datae['HRMIN']
#Conversion of time into more usable format from Doy and decimal formats
idx = []
for i in range(0,len(Y)):
    month, jd, y = JulianDate_to_MMDDYYY(Y[i],D[i])
    dt = str(y)+'-'+str(month)+'-'+str(jd)+' '+str(H[i])
    idx.append(dt)
#Create new index and fill in missing 30-minute time periods
datae.index = idx
datae.index = pd.to_datetime(datae.index)
datae = LLT.indx_fill(datae, '30min')
#Shift time to match correct start and end times and generate correct form
datae['TIMESTAMP_END'] = datae.index.shift(0, '30T')
datae['TIMESTAMP_START'] = datae.index.shift(-1, '30T')    
datae['TIMESTAMP_START']= datae.TIMESTAMP_START.map(lambda x: datetime.datetime.strftime(x, '%Y%m%d%H%M'))
datae['TIMESTAMP_END']= datae.TIMESTAMP_END.map(lambda x: datetime.datetime.strftime(x, '%Y%m%d%H%M'))
cols = datae.columns.tolist()
#Move TimeStamps to the first two columns of the dataset
cols.insert(0,cols.pop(cols.index('TIMESTAMP_START')))
cols.insert(1,cols.pop(cols.index('TIMESTAMP_END')))
datae = datae.reindex(columns = cols)   
#Output data following naming convention
datae.to_csv(r'C:\Users\Eric\Desktop\LTAR\LTAR_National_Projects\PhenologyInitiative\EC Data\Processed\Unprocessed\EC_LTAR_UCB_hawbeckereddy_201700101_20171231_redoneTime.csv', index = False)