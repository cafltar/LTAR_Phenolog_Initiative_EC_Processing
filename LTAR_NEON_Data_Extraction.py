# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 09:30:07 2018

@author: Eric
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 14:51:25 2018

@author: Eric Russell
Dept. Civil and Environmental Engineering
Washington State University
eric.s.russell@wsu.edu

This script pulls out selected columns (Cols) from the NEON dataset into a pandas dataframe with predetermined column headers.
The index used is for the beginning time of the averaging period but could be changed as needed.
Data can be output or manipulated as needed once in the pandas dataframe for the users usage. 
The file names and variable names are free to be changed by the user as they see fit; this is my naming convention
This script pulls individual lists of data from the *.h5 files and not whole subsections of particular variables...
a script like that could be created or this script could be modified as needed to pull whole objects but thought most likely...
use would be to pull individual variables so wrote it this way.
"""
import pandas as pd
import h5py
import glob

#Read-in data column names needed from the HDF5 file from separate 
Cols = pd.read_csv('C:\\Users\\Eric\\Desktop\\LTAR\\LTAR_National_Projects\\PhenologyInitiative\\NEON_Data_IDs_CLBJ.csv',header=0)
#Collect all the *.h5 files within the directory listed; can add an extra line to get only specific files based on files names if needed.
fnames =glob.glob('C:\\Users\\Eric\\Desktop\\LTAR\\LTAR_National_Projects\\PhenologyInitiative\\LBJ\\*.h5')
#Initialize blank dataframes prior to drawing in the 
Final = []; Final = pd.DataFrame(Final, columns=[Cols['Var_Name']])
Final_m = []; Final_m = pd.DataFrame(Final)
#%%
#Loop to read in the *h5 files grabbed in the glob.glob statement and over the different columns to be extracted
for K in range(0,len(fnames)):
    f = h5py.File(fnames[K], 'r')
    for k in range(0,len(Cols)):
        if k ==0: # Grabs the Begin-time from the first time around the loop for timestamp inclusion
            tme = serc_refl = pd.DataFrame(f[Cols['NEON'][k]][0:]['timeBgn'])
        serc_refl = pd.DataFrame(f[Cols['NEON'][k]][0:][Cols['Col_Num'][k]], columns =[Cols['Var_Name'][k]])
        Final_m[Cols['Var_Name'][k]] = serc_refl[Cols['Var_Name'][k]]
    #Decode time from full ISO8601 to YYYY-MM-DD hh:mm:ss
    tme = tme.astype(object)
    tme[0] = tme[0].str.decode("utf-8")
    Final_m.index = pd.to_datetime(tme[0])
    #Final concat before the final output of the dataset
    Final = pd.concat([Final,Final_m])
    Final_m = []; Final_m = pd.DataFrame(Final_m, columns=[Cols['Var_Name']])
# Output data
Final.to_csv(r'C:\Users\Eric\Desktop\CLBJ_Neon_EC_Data.csv')