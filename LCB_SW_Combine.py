# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 15:40:25 2019

@author: Eric
"""

import pandas as pd

df = pd.read_csv(r'C:\Users\Eric\Desktop\LTAR\LTAR_National_Projects\PhenologyInitiative\EC Data\Processed\Unprocessed\MissingSW\LTAR_EC_LCB_ope3_20170101_20171231.csv',header = 0)
SW = pd.read_csv(r'C:\Users\Eric\Desktop\LTAR\LTAR_National_Projects\PhenologyInitiative\EC Data\Processed\Unprocessed\MissingSW\OPE3_IncidentSolar.csv',header = 0, index_col = 'Time Stamp')
SW.index = pd.to_datetime(SW.index)
dt = []
for k in range (0,len(df)):
    Y =  str(int(df['TIMESTAMP_END'][k]))[0:4]
    M =  str(int(df['TIMESTAMP_END'][k]))[4:6]  
    D =  str(int(df['TIMESTAMP_END'][k]))[6:8]
    hh = str(int(df['TIMESTAMP_END'][k]))[8:10] 
    mm = str(int(df['TIMESTAMP_END'][k]))[10:12]
    dt.append(Y+'-'+M+'-'+D+' '+hh+':'+mm)
dt = pd.DataFrame(dt);df.index = dt[0]
df.index=pd.to_datetime(df.index) # Time-based index
#%
df['SW_IN'] = SW['Incident Solar']

df.to_csv(r'C:\Users\Eric\Desktop\LTAR\LTAR_National_Projects\PhenologyInitiative\EC Data\Processed\Unprocessed\EC_LTAR_LCB_ope3_20170101_20171001.csv',index= False)