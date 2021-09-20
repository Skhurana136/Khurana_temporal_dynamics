# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 15:25:45 2021

@author: khurana
"""

# Native libraries
import os

# Third party libraries for data retrieval and processing
import h5py
import numpy as np
import pandas as pd

# Visualisation libraries
import matplotlib.pyplot as plt
import seaborn as sns

# In house libraries
import data_reader.data_processing as proc

#Define directories
sourcedatadir = "E:/Saturated_flow/EGUGoldschmidtdataset6"
hdf5directory = "Y:/Home/khurana/4. Publications/Restructuring"

# Define commonly used scenarios
Regimes = ["Slow", "Equal", "Fast"]
reglist = ["Slow", "Medium", "Fast"]
gvarnames = ["DO", "DOC", "Nitrate", "Ammonium"]#, "TOC", "Nitrogen"]
imposedtimeseries = ["1","2","5"]
Trial = proc.masterscenarios("Saturated").keys()

# Define plotting related parameters
mymarklist = ["^", "o", "s", "d"]
colorlist = ["indianred", "g", "steelblue"]
PeDapalette = {0: "blue", 1: "orange", 2: "g", 3: "indianred"}
colormaps = [plt.cm.Blues, plt.cm.Oranges, plt.cm.Greens, plt.cm.Reds]

# Generate Dat from Da and Pe
data  = pd.read_csv(os.path.join(hdf5directory,"Paper1","Figurecodes", "Da_29012021_95pcloss.csv"), sep=",")
data["logDa"] = np.log10(data.Da63)

#Consider time series in terms of Dat.
#Sort all the values in Dat
finaldata = data[data['Chem'].isin (gvarnames)]

finaldata.loc[finaldata["logDa"] > 0.5, "PeDamark"] = 3
finaldata.loc[(finaldata["logDa"] > 0) & (finaldata["logDa"] < 0.5), "PeDamark"] = 2
finaldata.loc[(finaldata["logDa"] > -1) & (finaldata["logDa"] < 0), "PeDamark"] = 1
finaldata.loc[finaldata["logDa"] < -1, "PeDamark"] = 0

labels = {0 : 'log$_{10}$Da < -1',
          1 : '-1 < log$_{10}$Da < 0',
          2 : '0 <log$_{10}$Da < 0.5',
          3 : 'log$_{10}$Da > 0.5'}

subfinal = finaldata[["Trial", "Regime", "Chem", "Time_series","PeDamark"]]
subfinal['key'] = subfinal.Trial + subfinal.Regime + subfinal.Chem + subfinal.Time_series.astype(str)
colorcriteria = subfinal[["PeDamark", "key"]].to_dict('records')
colorcriteria = dict(zip(subfinal['key'], subfinal['PeDamark']))

# Combine velocities and normalised concentrations in a single dataframe
basedata = np.load(os.path.join(sourcedatadir, "SlowAR_0/NS-AH/NS-AH_df.npy"))
basevelocity = np.mean(basedata[2, -1, :, :])
hr = h5py.File(os.path.join(hdf5directory, "Paper2","Figurecodes","Temporal_analysis_full_data.h5"), mode = 'r')
alldata = pd.DataFrame(columns = ['Vel', 'Conc', 'Dat'])
for danum in [0,1,2,3]:
    n_danum = np.empty(1)
    vel_all3 = np.empty(1)
    for t in imposedtimeseries:
        n_danumt = []
        directory = os.path.join(sourcedatadir, "SlowAR_" + t)
        for j in ["37"]:
            data = np.load(os.path.join(directory,"NS-A"+j,"NS-A"+j+"_df.npy"))
            velocity = np.mean(data[2, :, :, :], axis = (-1,-2))/basevelocity
        for Reg in Regimes:
            if Reg == "Equal":
                r = "Medium"
            else:
                r = Reg    
            for g in gvarnames:
                for j in Trial:
                    if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                        pass
                    elif int(colorcriteria[j+r+g+"0"]) == danum:
                        n = hr.get(t + "/" + Reg + "/" + j + "/" + g).value
                        n_danumt.append(n)
                    else:
                        pass
        datapv = np.asarray(n_danumt).transpose()
        vel_t = np.asarray(list(velocity[1:])*np.shape(datapv)[1])
        dataflat = datapv.flatten('F')
        n_danum = np.append(n_danum,dataflat)
        vel_all3 = np.append(vel_all3,vel_t)
    df = pd.DataFrame(vel_all3,columns=['Vel'])
    df["Conc"]= n_danum
    df["Dat"] = danum
    alldata = pd.concat([alldata, df], axis=0)
    
#for danum in [0]:#,1,2,3]:
#    d=alldata[alldata.Dat==danum]
#    print(d.shape)

# Visualise distributions
#fig, ax = plt.subplots(2,2,sharex=True,figsize=(6,6))
#for danum in [0,1,2,3]:
#    axe = ax.flat[danum]
#    data = alldata[alldata.Dat==danum]
#    sns.jointplot(data.Vel, data.Conc, color = PeDapalette[danum], ax= axe)
#g = sns.FacetGrid(alldata, col = 'Dat', col_wrap=2,  hue = 'Dat')
#g = (g.map(sns.jointplot, "Vel", "Conc"))

sns.jointplot("Vel", "Conc", data=alldata, hue="Dat")
#plt.yscale("log")

kindtry = 'hist'
JG0 = sns.jointplot("Vel", "Conc", data=alldata[alldata.Dat==0], kind=kindtry,space=0, color="steelblue")
JG1 = sns.jointplot("Vel", "Conc", data=alldata[alldata.Dat==1], kind=kindtry,space=0, color="orange")
JG2 = sns.jointplot("Vel", "Conc", data=alldata[alldata.Dat==2], kind=kindtry,space=0, color="g")
JG3 = sns.jointplot("Vel", "Conc", data=alldata[alldata.Dat==3], kind=kindtry,space=0, color="indianred")

import matplotlib.gridspec as gridspec
import SeabornFig2Grid as sfg

fig = plt.figure(figsize=(8,8))
gs = gridspec.GridSpec(2, 2)

mg0 = sfg.SeabornFig2Grid(JG0, fig, gs[0])
mg1 = sfg.SeabornFig2Grid(JG1, fig, gs[1])
mg2 = sfg.SeabornFig2Grid(JG2, fig, gs[3])
mg3 = sfg.SeabornFig2Grid(JG3, fig, gs[2])
gs.tight_layout(fig)
#gs.update(top=0.7)

plt.show()

#subplots migration
f = plt.figure()
for J in [JG0,JG1, JG2, JG3]:
    for A in J.fig.axes:
        f._axstack.add(f._make_key(A), A)

#subplots size adjustment
f.axes[0].set_position([0.025, 0.025, 0.2,  0.2])
f.axes[1].set_position([0.025, 0.225, 0.2,  0.025])
f.axes[2].set_position([0.225, 0.025, 0.025, 0.2])
f.axes[3].set_position([0.275, 0.025, 0.2,  0.2])
f.axes[4].set_position([0.275, 0.225, 0.2,  0.025])
f.axes[5].set_position([0.475, 0.025, 0.025, 0.2])
f.axes[6].set_position([0.025, 0.025, 0.2,  0.2])
f.axes[7].set_position([0.025, 0.225, 0.2,  0.025])
f.axes[8].set_position([0.225, 0.025, 0.025, 0.2])
f.axes[9].set_position([0.275, 0.025, 0.2,  0.1])
f.axes[10].set_position([0.275, 0.225, 0.2,  0.025])
f.axes[11].set_position([0.475, 0.025, 0.025, 0.2])