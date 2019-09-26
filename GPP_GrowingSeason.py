# -*- coding: utf-8 -*-
"""
Created on Wed May 29 10:12:24 2019

@author: Eric
"""


import pandas as pd
import matplotlib.pyplot as plt
import glob as glob
import numpy as np
import os
os.chdir(r'C:\Users\Eric\Desktop\PyScripts\Phenology_LTAR')
import Phenology_Readins as PR

Tair, H, LE, Fc, Reco, GPP, ET, GDD_10_F, GDD_5_F = PR.Phen_ECData_Read() # Add the file paths as inputs to the functions once get things more settled
GCC_M, RCC_M, BCC_M, Names = PR.Phen_CC_Readin()

NME = pd.read_excel(r'C:\Users\Eric\Desktop\LTAR\LTAR_National_Projects\PhenologyInitiative\Camera_Data\EC_PhenoCam_Site_Names.xlsx', header = 0)
SES = pd.read_csv(r'C:\Users\Eric\Desktop\LTAR\LTAR_National_Projects\PhenologyInitiative\Camera_Data\dt_ltar_pheno_dates_all_sites_year.csv',header = 0)
SES = SES[SES['site_name'].isin(NME['PhenoCam'])]
SES = SES[SES['year'] > 2016]

Site_SES = SES.groupby('site_name')
h_col = ['Site','Date','SE','Year']

#%%
datesum = []
for k in range (5,6):
    if NME['PhenoCam'][k] != -9999:
        a = Site_SES.get_group(NME['PhenoCam'][k]).values
        a = pd.DataFrame(a, columns = h_col)
        fig = plt.figure(k+1, figsize=(9,6))
        ax = plt.subplot(211)
        plt.plot(GPP[NME['EC Tower'][k]],'g')
        plt.plot(Fc[NME['EC Tower'][k]],'k')
        yu = round(np.max(GPP[NME['EC Tower'][k]])); yu = yu +(10-yu%10)
        yl = round(np.min(GPP[NME['EC Tower'][k]])); yl = yl - (10-yl%10)
        for ses in range (0, len(a)):
            plt.vlines(a['Date'][ses],yl,yu)
            if ses%2 == 0:
                ax.fill_between(GPP[NME['EC Tower'][k]].index,yl,yu, where = (GPP[NME['EC Tower'][k]].index>=a['Date'][ses])& (GPP[NME['EC Tower'][k]].index<=a['Date'][ses+1]), facecolor='grey',alpha = 0.25, interpolate=True)
                datesum.append([a['Date'][ses], a['Date'][ses+1], np.sum(GPP[NME['EC Tower'][k]][a['Date'][ses]:a['Date'][ses+1]]),NME['EC Tower'][k],np.sum(GCC_M[NME['PhenoCam'][k]][a['Date'][ses]:a['Date'][ses+1]]),NME['EC Tower'][k]])
        if ~np.isnan(GPP[NME['EC Tower'][k]][-1]):
            PR.format_plot(ax, 12,12,3,5, yl,yu,'2017-01-01','2019-01-01')
        else:
            PR.format_plot(ax, 12,12,3,5, yl,yu,'2017-01-01','2018-01-01')
        plt.ylabel('[gC m$^{-2}$ day$^{-1}$]',fontsize = 12)
        plt.legend(['GEP','NEE'],loc = 2,fontsize = 12)
        ax = plt.twinx()
        plt.plot(GCC_M[NME['PhenoCam'][k]],'m.', alpha = 0.65)
        if ~np.isnan(GPP[NME['EC Tower'][k]][-1]):
            PR.format_plot(ax, 12,12,3,5, 0.3,0.5,'2017-01-01','2019-01-01')
        else:
            PR.format_plot(ax, 12,12,3,5, 0.3,0.5,'2017-01-01','2018-01-01')
        plt.ylabel('GCC 90%',fontsize = 12);plt.legend(['Camera'],loc =1,fontsize = 12)
        
        ax = plt.subplot(212)
        plt.plot(ET[NME['EC Tower'][k]],'b')
        yu = round(np.max(ET[NME['EC Tower'][k]]))
        yu = yu +(10-yu%10)
        yl = round(np.min(ET[NME['EC Tower'][k]]))
        for ses in range (0, len(a)):
            plt.vlines(a['Date'][ses],yl,yu)
            if ses%2 == 0:
                ax.fill_between(ET[NME['EC Tower'][k]].index,yl,yu, where = (ET[NME['EC Tower'][k]].index>=a['Date'][ses])& (ET[NME['EC Tower'][k]].index<=a['Date'][ses+1]), facecolor='grey',alpha = 0.25, interpolate=True)
        if ~np.isnan(ET[NME['EC Tower'][k]][-1]):
            PR.format_plot(ax, 12,12,3,5, yl,yu,'2017-01-01','2019-01-01')
        else:
            PR.format_plot(ax, 12,12,3,5, yl,yu,'2017-01-01','2018-01-01')
        plt.ylabel('[mm H$_2$O day$^{-1}$]',fontsize = 12)
        plt.legend(['ET'],loc = 2,fontsize = 12)
        ax = plt.twinx()
        plt.plot(GCC_M[NME['PhenoCam'][k]],'m.', alpha = 0.65)
        if ~np.isnan(ET[NME['EC Tower'][k]][-1]):
            PR.format_plot(ax, 12,12,3,5, 0.3,0.5,'2017-01-01','2019-01-01')
        else:
            PR.format_plot(ax, 12,12,3,5, 0.3,0.5,'2017-01-01','2018-01-01')
        plt.ylabel('GCC 90%',fontsize = 12);plt.legend(['Camera'],loc =1,fontsize = 12)
        st = fig.suptitle(NME['EC Tower'][k],fontsize = 16)
        plt.tight_layout()
        st.set_y(0.975)
        fig.subplots_adjust(top=0.90)
#        file = r'C:\Users\Eric\Desktop\LTAR\LTAR_National_Projects\PhenologyInitiative\Camera_Data\Figures'
#        plt.savefig(file+ '\\' + str(NME['EC Tower'][k])+'_20172018_TS_StartEndSeason_Plots.png',dpi = 300)
#        plt.close()
        
#%%
import re
ds = pd.DataFrame(datesum)
qn = ds[4]>0
df = ds[qn]
a_s = []
for k in df.index:
    aa = re.findall(r'[A-Z]',df[3][k])
    a_s.append(''.join(aa))

df[6] = a_s

g = df.groupby(df[6])
plt.figure(1000, figsize=(8,6))
ax = plt.subplot(111)
for k in g.groups.keys():
    q = g.get_group(k)
    plt.plot(q[2],q[4],'.', label = k, ms = 20)
plt.ylabel('iGCC', fontsize = 16)
plt.xlabel('MGS GEP', fontsize = 16)
PR.format_plot(ax, 16,16,5,4,0,120,0,2500)
plt.yticks(fontsize = 16);plt.xticks(fontsize = 16)
plt.legend(fontsize = 14, ncol = 3)