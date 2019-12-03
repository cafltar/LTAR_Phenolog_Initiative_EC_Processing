# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 10:58:29 2018

@author: Eric Russell
Laboratory for Atmospheric Research
Dept. Civil and Environmental Engineering
Washington State University
eric.s.russell@wsu.edu

"""
import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from matplotlib.ticker import AutoMinorLocator

def format_plot(ax,yf,xf,xminor,yminor,yl,yu,xl,xu):
    #subplot has to have ax as the axis handle
    plt.yticks(fontsize = yf);plt.xticks(fontsize = xf)
    minor_locator = AutoMinorLocator(xminor)
    ax.xaxis.set_minor_locator(minor_locator)
    minor_locator = AutoMinorLocator(yminor)
    ax.yaxis.set_minor_locator(minor_locator)
    ax.tick_params(axis='both',direction='in',length=12.5,width=2)
    ax.tick_params(axis='both',which = 'minor',direction='in',length=5)
    plt.ylim([yl,yu])
    plt.xlim([xl,xu])  
    return

Path = 'C:\\Users\\Eric\\Desktop\\LTAR\\LTAR_National_Projects\\PhenologyInitiative\\EC Data\\Final_Data\\Daily\\'
files = glob.glob(Path+'*Daily.csv')
Fc = []; Fc = pd.DataFrame(Fc)
Reco = []; Reco = pd.DataFrame(Reco)
GPP = []; GPP = pd.DataFrame(GPP)
ET = []; ET = pd.DataFrame(ET)
Tair = []; Tair = pd.DataFrame(Tair)
LE = []; LE = pd.DataFrame(LE)
H = []; H = pd.DataFrame(H)
year, Site = [],[]
GDD_10_F, GDD_5_F, GDD_0_F = [], [],[]
ds = '2017-01-01'
de = '2018-12-31'
for K in range(0, len(files)):
    datae = pd.read_csv(files[K],header=0,index_col='TimeStamp',skiprows = [1])
    datae.index = pd.to_datetime(datae.index)
    datae = datae[ds:de]
    Site.append(files[K][113:-43])
    year.append(files[K][-33:-29])
    Fc = pd.concat([Fc,datae['Fc_Gapfilled_L3_Sum_Carbon']],axis = 1)
    Reco = pd.concat([Reco,datae['Reco_L3_Sum_Carbon']],axis = 1)
    GPP = pd.concat([GPP,datae['GPP_L3_Sum_Carbon']],axis = 1)
    ET = pd.concat([ET,datae['ET_Gapfilled_L3_Sum']],axis=1)
#    H= pd.concat([H,datae['H_Gapfilled_L3_Sum']],axis=1)
    LE= pd.concat([LE,datae['LE_Gapfilled_L3_Sum']],axis=1)
    Tair = pd.concat([Tair,datae['Tair_Gapfilled_L3']],axis=1)
    qn = datae['Tair_Gapfilled_L3_Max'] >50
    for m in range (0,len(qn)):
        if qn[m]==True:
            datae['Tair_Gapfilled_L3_Max'][m] = (datae['Tair_Gapfilled_L3_Max'][m+1]+datae['Tair_Gapfilled_L3_Max'][m-1])/2
    # Insert a step to check that Tmax isn't some obscene value like at Cook; not sure why that is the case.
    G_Tmin = pd.DataFrame(datae['Tair_Gapfilled_L3_Min']); G_Tmin.index = pd.to_datetime(G_Tmin.index)
    G_Tmin['Tair_Gapfilled_L3_Min'][G_Tmin['Tair_Gapfilled_L3_Min']<10] = 10
    GDD_10 = ((datae['Tair_Gapfilled_L3_Max']+G_Tmin['Tair_Gapfilled_L3_Min'])/2)-10
    G_Tmin['Tair_Gapfilled_L3_Min'][G_Tmin['Tair_Gapfilled_L3_Min']<5] = 5
    GDD_5 = ((datae['Tair_Gapfilled_L3_Max']+G_Tmin['Tair_Gapfilled_L3_Min'])/2)-5
    G_Tmin['Tair_Gapfilled_L3_Min'][G_Tmin['Tair_Gapfilled_L3_Min']<0] = 0
    GDD_0 = ((datae['Tair_Gapfilled_L3_Max']+G_Tmin['Tair_Gapfilled_L3_Min'])/2)-0
    GDD_0[GDD_5<0] = 0; GDD_0_F.append(GDD_0)
    GDD_5[GDD_5<0] = 0; GDD_5_F.append(GDD_5)
    GDD_10[GDD_10<0] = 0; GDD_10_F.append(GDD_10)
    
cls = pd.DataFrame(Site)
year = pd.DataFrame(year, columns = ['year'])
Fc.index=pd.to_datetime(Fc.index)
Reco.index=pd.to_datetime(Reco.index)
GPP.index=pd.to_datetime(GPP.index)
ET.index=pd.to_datetime(ET.index)

GDD_10_F = pd.DataFrame(GDD_10_F)
GDD_10_F = GDD_10_F.transpose()
GDD_10_F.columns = Site

GDD_5_F = pd.DataFrame(GDD_5_F)
GDD_5_F = GDD_5_F.transpose()
GDD_5_F.columns = Site

GDD_0_F = pd.DataFrame(GDD_0_F)
GDD_0_F = GDD_0_F.transpose()
GDD_0_F.columns = Site

Tair.columns = Site
#H.columns = Site
LE.columns = Site

Fc.columns = cls[0]
Reco.columns = cls[0]
GPP.columns = cls[0]
ET.columns = cls[0]
Fc_W = Fc.resample('W').sum()
Reco_W = Reco.resample('W').sum()
GPP_W = GPP.resample('W').sum()
ET_W = ET.resample('W').sum()
GPP_Y = GPP.resample('Y').sum()
Reco_Y = Reco.resample('Y').sum()
ET_Y = ET.resample('Y').sum()


GPP_Y = GPP_Y.groupby(level=0, axis=1).sum()
ET_Y = ET_Y.groupby(level=0, axis=1).sum()
Reco_Y = Reco_Y.groupby(level=0, axis=1).sum()
GPP_W = GPP_W.groupby(level=0, axis=1).sum()
Reco_W = Reco_W.groupby(level=0, axis=1).sum()
ET_W = ET_W.groupby(level=0, axis=1).sum()



a = []
for K in range(0,len(cls)):
    a.append(cls[0][K].split('_'))
a = pd.DataFrame(a, columns = ['ABV','Tower','Nope'])
ass = a.ABV.unique()
#%%
GPP_Site = []; GPP_Site = pd.DataFrame()
ET_Site = []; ET_Site = pd.DataFrame()
for k in range(0,len(ass)):
    tempp= GPP_W.filter(like=ass[k],axis=1)
    ss = len(tempp.transpose())
    if ass[k] != 'GACP' and ass[k] !='GB':
        tempp = np.sum(tempp, axis = 1)/ss
    elif ass[k] == 'GB': 
        front = np.sum(tempp['2017-01-01':'2017-12-31'], axis = 1); front.index = pd.to_datetime(front.index)
        back = np.sum(tempp['2018-01-01':'2018-12-31'], axis = 1)/ss; back.index = pd.to_datetime(back.index)
        tempp = pd.concat([front,back])
    elif ass[k] == 'GACP':
        front = np.sum(tempp['2017-01-01':'2017-12-31'], axis = 1)/ss; front.index = pd.to_datetime(front.index)
        back = np.sum(tempp['2018-01-01':'2018-12-31'], axis = 1); back.index = pd.to_datetime(back.index)
        tempp = pd.concat([front,back])
    GPP_Site = pd.concat([GPP_Site,tempp], axis = 1)
GPP_Site.columns = ass
GPP_Site[GPP_Site == 0] = np.NaN

for k in range(0,len(ass)):
    tempp= ET_W.filter(like=ass[k],axis=1)
    ss = len(tempp.transpose())
    if ass[k] != 'GACP' and ass[k] !='GB':
        tempp = np.sum(tempp, axis = 1)/ss
    elif ass[k] == 'GB': 
        front = np.sum(tempp['2017-01-01':'2017-12-31'], axis = 1); front.index = pd.to_datetime(front.index)
        back = np.sum(tempp['2018-01-01':'2018-12-31'], axis = 1)/ss; back.index = pd.to_datetime(back.index)
        tempp = pd.concat([front,back])
    elif ass[k] == 'GACP':
        front = np.sum(tempp['2017-01-01':'2017-12-31'], axis = 1)/ss; front.index = pd.to_datetime(front.index)
        back = np.sum(tempp['2018-01-01':'2018-12-31'], axis = 1); back.index = pd.to_datetime(back.index)
        tempp = pd.concat([front,back])
    ET_Site = pd.concat([ET_Site,tempp], axis = 1)
ET_Site.columns = ass
ET_Site[ET_Site == 0] = np.NaN

#%%
colors = pl.cm.jet(np.linspace(0,1,len(cls)))

plt.figure(1,figsize=(11,8))
ax = plt.subplot(211)
for K in range (0,len(cls)):
    plt.plot(Fc_W[cls[0][K]],linewidth = 1.5,color = colors[K],alpha = 0.75)
plt.ylabel('NEE \n [gC week$^{-1}$ m$^{-2}$]',fontsize=14)
plt.legend(cls[0],ncol=3,loc= 'lower left',bbox_to_anchor=(0,1.02,1,0.2),mode = 'expand')
format_plot(ax,16,16,1,5,-150,50,ds,de)
plt.xticks(rotation = 5)
ax = plt.subplot(212)
for K in range (0,len(cls)):
    plt.plot(ET_W[cls[0][K]],linewidth = 1.5,color = colors[K],alpha = 0.75)
plt.ylabel('ET \n [mm-H$_2$O week$^{-1}$ m$^{-2}$]',fontsize=14)
format_plot(ax,16,16,1,5,-10,40,ds,de)
plt.xticks(rotation = 5)
plt.tight_layout(rect=[0,0,1,0.925])
#plt.savefig('C:\\Users\\Eric\\Desktop\\LTAR\\LTAR_National_Projects\\PhenologyInitiative\\EC Data\\Final_Data\\Figures\\FC_ET_Weekly_Pheno_LTAR.png',dpi=600)

plt.figure(2,figsize=(11,8))
ax = plt.subplot(211)
for K in range (0,len(cls)):
    plt.plot(GPP_W[cls[0][K]],linewidth = 1.5,color = colors[K], alpha = 0.75)
plt.legend(cls[0],ncol=3,loc= 'lower left',bbox_to_anchor=(0,1.02,1,0.2),mode = 'expand')
plt.ylabel('GPP \n [gC week$^{-1}$ m$^{-2}$]',fontsize=14)
format_plot(ax,16,16,2,5,-50,200,ds,de)
plt.xticks(rotation = 5)
ax = plt.subplot(212)
for K in range (0, len(cls)):
    plt.plot(Reco_W[cls[0][K]],linewidth = 1.5,color = colors[K],alpha = 0.75)
plt.ylabel('Reco \n [gC week$^{-1}$ m$^{-2}$]',fontsize=14)
format_plot(ax,16,16,2,5,-25,100,ds,de)
plt.xticks(rotation =5)
plt.tight_layout(rect=[0,0,1,0.925])
#plt.savefig('C:\\Users\\Eric\\Desktop\\LTAR\\LTAR_National_Projects\\PhenologyInitiative\\EC Data\\Final_Data\\Figures\\GPP_Reco_Weekly_Pheno_LTAR.png',dpi=600)

#%% 
a = []
for K in range(0,len(cls)):
    a.append(cls[0][K].split('_'))
a = pd.DataFrame(a)
year.index = a[0]
G_S = np.sum(GPP_W[ds:de]); G_S.index = a[0]; G_S = pd.DataFrame(G_S, columns = ['G_S'])
G_S = pd.concat([G_S, year],axis = 1); G_S.sort_values(by=['year'])
R_S = np.sum(Reco_W[ds:de]); R_S.index = a[0];R_S = pd.DataFrame(R_S, columns = ['R_S'])

df = pd.concat([G_S,R_S],axis = 1)
df['site'] = df.index
colors = {'ABS':'green','CAF':'brown','CMR':'midnightblue', 'CPE':'black','GAC':'palegreen','GB_':'magenta','JORN':'teal','KONA':'lime',
          'NP':'orange','PRHPA':'blue','SEG':'yellowgreen','TG':'tan','UCB':'cyan', 'UMRB':'yellow','WGEW':'grey' }
grp = df.groupby('site')
mrk = {'2017':'x','2018':'.'}
xy_line = (0,2500)
plt.figure(3,figsize=(11,8))
ax = plt.subplot(111)
for key, group in grp:
    grp_yr = grp.get_group(key).groupby('year')
    for key_yr, group in grp_yr:
        if key_yr=='2018':
            group.plot(ax=ax,kind='scatter',x='G_S',y='R_S',marker = mrk[key_yr], color = colors[key], s= 500, alpha = 0.75)    
            if key == 'JORN' or key== 'KONA':
                group.plot(ax=ax,kind='scatter',x='G_S',y='R_S',label = key,marker = mrk[key_yr], color = colors[key], s= 500, alpha = 0.75)    
        else: group.plot(ax=ax,kind='scatter',x='G_S',y='R_S',marker = mrk[key_yr],label = key, color = colors[key],linewidth = 3, s= 500, alpha = 0.75)    
#plt.ylim([0,2500]); plt.xlim([0,2500])
#%%
GPP_Y = GPP_Y.transpose();Reco_Y = Reco_Y.transpose()
GPP_Y[abs(GPP_Y) <=2] = np.NaN; Reco_Y[abs(Reco_Y) <=2] = np.NaN
Site= pd.DataFrame(Site)
Site.columns = ['Site'] ; xs = Site.Site.unique()
#%%
colors = {'ABS':'green','CAF':'brown','CMR':'midnightblue', 'CPE':'black','GAC':'tomato','GB_':'magenta','JOR':'teal','KON':'lime',
          'NP_':'orange','PRH':'blue','SEG':'yellowgreen','TG_':'tan','UCB':'cyan', 'UMR':'yellow','WGE':'grey' }

import re
xy_line = (0,2500)
plt.figure(3,figsize=(11,8))
ax = plt.subplot(111)
for k in range(0, len(xs)):
    r = re.findall('([A-Z])', GPP_Y.index[k])
    plt.plot(GPP_Y['2017-12-31'][xs[k]],Reco_Y['2017-12-31'][xs[k]],'x', color = colors[str(xs[k])[0:3]],ms = 17,mew = 3, alpha = 0.7)
    plt.plot(GPP_Y['2018-12-31'][xs[k]],Reco_Y['2018-12-31'][xs[k]],'.', color = colors[str(xs[k])[0:3]],ms = 25, alpha = 0.7,
             label= ''.join(r) if ''.join(r) not in plt.gca().get_legend_handles_labels()[1] else '')
plt.plot(xy_line,xy_line,'r--', linewidth = 3)
plt.legend(fontsize=16, ncol = 4)
plt.ylabel('Total Reco [gC m$^{-2}$]', fontsize = 18); plt.yticks(np.linspace(0,2500,11),fontsize = 16)
plt.xlabel('Total GEP [gC m$^{-2}$]', fontsize = 18); plt.xticks(np.linspace(0,2500,11),fontsize = 16)
ax.annotate('More Uptake',xy=(1000,1000),xytext=(1500,500), arrowprops=dict(arrowstyle='<-'),ha = 'center',va='center',fontsize = 24)
ax.annotate('More Loss',xy=(1000,1000),xytext=(500,1500), arrowprops=dict(arrowstyle='<-'),ha = 'center',va='center',fontsize = 24)
plt.tight_layout()
#plt.savefig('C:\\Users\\Eric\\Desktop\\LTAR\\LTAR_National_Projects\\PhenologyInitiative\\EC Data\\Final_Data\\Figures\\Annual_Reco_GPP_Comparison_PI_2017_2018.png',dpi=600)
#%%
colors = {'ABS-UF':'green','CAF':'brown','CMRB':'midnightblue', 'CPER':'black','GACP':'tomato','GB':'magenta','JORN':'teal','KBS':'purple','KONA':'lime',
          'NP':'orange','PRHPA':'blue','SEG':'yellowgreen','TG':'tan','UCB':'cyan', 'UMRB':'yellow','WGEW':'grey' }
plt.figure(4, figsize=(11,8))

for k in range (0,len(colors)):
    ax= plt.subplot(211)
    plt.plot(GPP_Site[ass[k]],linewidth = 2, color = colors[ass[k]],label = ass[k])
    ax = plt.subplot(212)
    plt.plot(ET_Site[ass[k]],linewidth = 2,color = colors[ass[k]])
    plt.ylabel('Normalized ET\n [mm H$_2$O week$^{-1}$]',fontsize = 14)
ax= plt.subplot(211)
format_plot(ax, 14,14,3,5,-10,160,'2017-01-01','2019-01-01')
plt.legend(ncol = 5, fontsize = 12)
plt.ylabel('Normalized GEP \n[gC m$^{-2}$ week$^{-1}$]',fontsize = 14)
ax= plt.subplot(212)
format_plot(ax, 14,14,3,5,-5,40,'2017-01-01','2019-01-01')
plt.tight_layout()
#plt.savefig('C:\\Users\\Eric\\Desktop\\LTAR\\LTAR_National_Projects\\PhenologyInitiative\\EC Data\\Final_Data\\Figures\\Pheno_EC_Weekly_GEP_ET_TowerNormal_20172018.png',dpi=400)

#%%
a = []
for K in range(0,len(cls)):
    a.append(cls[0][K].split('_'))
a = pd.DataFrame(a)
year.index = a[0]
G_S = np.sum(GPP_W); G_S.index = a[0]; G_S = pd.DataFrame(G_S, columns = ['G_S'])
G_S = pd.concat([G_S, year],axis = 1); G_S.sort_values(by=['year'])
R_S = np.sum(GDD_10_F); R_S.index = a[0];R_S = pd.DataFrame(R_S, columns = ['R_S'])
df = pd.concat([G_S,R_S],axis = 1)
df['site'] = df.index
colors = {'ABS-UF':'green','CAF':'brown','CMRB':'midnightblue', 'CPER':'black','GB':'magenta','JORN':'teal',
          'NP':'orange','PRHPA':'blue','UCB':'cyan', 'UMRB':'yellow','WGEW':'grey'}
mrk = {'2017':'x','2018':'.'}
grp = df.groupby('site')
plt.figure(4,figsize=(11,8))
ax = plt.subplot(311)
for key, group in grp:
    grp_yr = grp.get_group(key).groupby('year')
    for key_yr, group in grp_yr:
        if key_yr=='2018':
            group.plot(ax=ax,kind='scatter',x='R_S',y='G_S',marker = mrk[key_yr], color = colors[key], s= 500, alpha = 0.75)    
            if key == 'JORN': group.plot(ax=ax,kind='scatter',x='R_S',y='G_S',label = key,marker = mrk[key_yr], color = colors[key], s= 500, alpha = 0.75)    
        else: group.plot(ax=ax,kind='scatter',x='R_S',y='G_S',marker = mrk[key_yr],label = key, color = colors[key],linewidth = 3, s= 500, alpha = 0.75)    
plt.ylim([0,3000]); plt.xlim([0,6000])
plt.xlabel('Total GDD [10 degC]', fontsize = 18); plt.xticks(np.linspace(0,6000,11),fontsize = 16)
plt.ylabel('Total GPP [gC m$^{-2}$]', fontsize = 18); plt.yticks(np.linspace(0,3000,11),fontsize = 16)
plt.legend(fontsize=14, ncol = 3, loc = 1)

for K in range(0,len(cls)):
    a.append(cls[0][K].split('_'))
a = pd.DataFrame(a)
G_S = np.sum(GPP_W); G_S.index = a[0]; G_S = pd.DataFrame(G_S, columns = ['G_S'])
G_S = pd.concat([G_S, year],axis = 1); G_S.sort_values(by=['year'])
R_S = np.sum(GDD_5_F); R_S.index = a[0];R_S = pd.DataFrame(R_S, columns = ['R_S'])
df = pd.concat([G_S,R_S],axis = 1)
df['site'] = df.index
colors = {'ABS-UF':'green','CAF':'brown','CMRB':'midnightblue', 'CPER':'black','GB':'magenta','JORN':'teal',
          'NP':'orange','PRHPA':'blue','UCB':'cyan', 'UMRB':'yellow','WGEW':'grey'}
grp = df.groupby('site')
ax = plt.subplot(312)
for key, group in grp:
    grp_yr = grp.get_group(key).groupby('year')
    for key_yr, group in grp_yr:
        if key_yr=='2018':
            group.plot(ax=ax,kind='scatter',x='R_S',y='G_S',marker = mrk[key_yr], color = colors[key], s= 500, alpha = 0.75)    
            if key == 'JORN': group.plot(ax=ax,kind='scatter',x='R_S',y='G_S',marker = mrk[key_yr], color = colors[key], s= 500, alpha = 0.75)    
        else: group.plot(ax=ax,kind='scatter',x='R_S',y='G_S',marker = mrk[key_yr], color = colors[key],linewidth = 3, s= 500, alpha = 0.75)    
plt.ylim([0,3000]); plt.xlim([0,7500])
plt.xlabel('Total GDD [5 degC]', fontsize = 18); plt.xticks(np.linspace(0,7500,11),fontsize = 16)
plt.ylabel('Total GPP [gC m$^{-2}$]', fontsize = 18); plt.yticks(np.linspace(0,3000,11),fontsize = 16)

for K in range(0,len(cls)):
    a.append(cls[0][K].split('_'))
a = pd.DataFrame(a)
G_S = np.sum(GPP_W); G_S.index = a[0]; G_S = pd.DataFrame(G_S, columns = ['G_S'])
G_S = pd.concat([G_S, year],axis = 1); G_S.sort_values(by=['year'])
R_S = np.sum(GDD_0_F); R_S.index = a[0];R_S = pd.DataFrame(R_S, columns = ['R_S'])
df = pd.concat([G_S,R_S],axis = 1)
df['site'] = df.index
colors = {'ABS-UF':'green','CAF':'brown','CMRB':'midnightblue', 'CPER':'black','GB':'magenta','JORN':'teal',
          'NP':'orange','PRHPA':'blue','UCB':'cyan', 'UMRB':'yellow','WGEW':'grey'}
grp = df.groupby('site')
ax = plt.subplot(313)
for key, group in grp:
    grp_yr = grp.get_group(key).groupby('year')
    for key_yr, group in grp_yr:
        if key_yr=='2018':
            group.plot(ax=ax,kind='scatter',x='R_S',y='G_S',marker = mrk[key_yr], color = colors[key], s= 500, alpha = 0.75)    
            if key == 'JORN': group.plot(ax=ax,kind='scatter',x='R_S',y='G_S',marker = mrk[key_yr], color = colors[key], s= 500, alpha = 0.75)    
        else: group.plot(ax=ax,kind='scatter',x='R_S',y='G_S',marker = mrk[key_yr], color = colors[key],linewidth = 3, s= 500, alpha = 0.75)    
#plt.ylim([0,3000]); plt.xlim([0,7500])
plt.xlabel('Total GDD [5 degC]', fontsize = 18); plt.xticks(np.linspace(0,7500,11),fontsize = 16)
plt.ylabel('Total GPP [gC m$^{-2}$]', fontsize = 18); plt.yticks(np.linspace(0,3000,11),fontsize = 16)
plt.tight_layout()
#plt.savefig('C:\\Users\\Eric\\Desktop\\LTAR\\LTAR_National_Projects\\PhenologyInitiative\\EC Data\\Final_Data\\Figures\\Annual_GDD_GPP_Comparison_PI_2018_2017.png',dpi=600)
#%%
#Path = 'C:\\Users\\Eric\\Desktop\\LTAR\\LTAR_National_Projects\\PhenologyInitiative\\EC Data\\Final_Data\\'
#files = glob.glob(Path+'*30Mins.csv')
#Fc_UC = []; Fc_UC = pd.DataFrame(Fc_UC)
#H_UC = []; H_UC = pd.DataFrame(H_UC)
#ET_UC = []; ET_UC = pd.DataFrame(ET_UC)
#Site = []
#for K in range(0, len(files)):
#    datae = pd.read_csv(files[K],header=0,index_col='TimeStamp')
#    # Need an Fc, GPP, Reco, H, LE, VPD, Tair, variable
#    Fc_UC = pd.concat([Fc_UC,datae['Fc_Gapfilled_SD']],axis = 1, sort = 'True')
#    ET_UC = pd.concat([ET_UC,datae['ET_Gapfilled_SD']],axis = 1, sort = 'True')
#    H_UC = pd.concat([H_UC,datae['H_GapFilled_SD']],axis = 1, sort = 'True')
#    Site.append(files[K][107:-44])
#
#Fc_UC.index=pd.to_datetime(Fc_UC.index)
#ET_UC.index=pd.to_datetime(ET_UC.index)
#Fc_Sum_UC = ((np.sum(Fc_UC**2))**0.5)*12*10**(-6)*60*30 
#ET_Sum_UC = ((np.sum(ET_UC**2))**0.5)
#
#Tots = []; Tots = pd.DataFrame(Tots, index = Site)
#Tots['NEE'] = np.sum(Fc_W)
#Tots['GPP'] = np.sum(GPP_W*-1)
#Tots['Reco'] = np.sum(Reco_W)
#Tots['ET'] = np.sum(ET_W)
#Tots['GDD_5'] = np.sum(GDD_5_F)
#Tots['GDD_10'] = np.sum(GDD_10_F)
#Tots['H'] = np.sum(H*0.0018)
#Tots['Tair'] = np.mean(Tair)
#Tots['Tair_min'] = np.min(Tair)
#Tots['Tair_max'] = np.max(Tair)
#
#BR = H/LE
#qn = (BR<50)&(BR>0)
#Tots['BR'] = np.mean(BR[qn])


