# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 14:02:56 2020

@author: khurana
"""
#This script analyses the effect of velocity changes at the end of long dry spells or long wet spells

import numpy as np
import pandas as pd
import data_reader.data_processing as proc
import matplotlib.pyplot as plt
import analyses.saturated_transient as sta
import h5py

#set up basic constants 
Regimes = ["Slow", "Equal", "Fast"]
velocities = [0.00038, 0.0038, 0.038]
scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run

# Scenarios to investigate:
Trial = list(scdict.keys())
#Trial = ["H", "37", "38", "39", "40", "41", "42", "43", "44", "45"]
reginvest = Regimes
imposedtimeseries = ["1", "2", "5"]

vardict = proc.speciesdict("Saturated")
gvarnames = list(t for t in vardict.keys() if vardict[t]["Location"]=="Mobile") + ["Nitrogen", "TOC"]

#Identify indices where velocity is higher or lower than average
row = []
for Reg in reginvest:
    #basedata = np.load("E:/Saturated_flow/EGUGoldschmidtdataset6/" + Reg + "AR_0/NS-AH/NS-AH_df.npy")
    #basevelocity = np.mean(basedata[2, -1, :, :])
    for t in imposedtimeseries:
        directory = "E:/Saturated_flow/EGUGoldschmidtdataset6/" + Reg + "AR_" + t + "/"
        print (Reg, t)
        for j in [37, 45, 79]:
            basedata = np.load("E:/Saturated_flow/EGUGoldschmidtdataset6/" + Reg + "AR_0/NS-A" + str(j) + "/NS-A" + str(j) + "_df.npy")
            basevelocity = np.mean(basedata[2, -1, :, :])
            data = np.load(directory + "NS-A"+str(j)+"/NS-A"+str(j)+"_df.npy")
            velocity = list(np.mean(data[2,:,:,:], axis = (-1,-2)))
            slowind = [velocity.index(v) for v in velocity if v > basevelocity]
            fastind = [velocity.index(v) for v in velocity if v < basevelocity]            
            row.append([Reg, t, j, fastind, slowind])

#Velocities are calculated in negative (downward y direction) so the comparison sign is switched to consider this.

flowdata = pd.DataFrame.from_records (row, columns = ["Regime", "Time_series", "Trial", "Faster", "Slower"])

#indices of interest
def gensubtimeseries (series, minlength):
    subseries = []
    for idx in range(len(series)-1):
        diff = series[idx+1] - series[idx]
        if diff == 1:
            subseries.append(series[idx])
        if diff > 1:
            subseries.append(series[idx])
            #print(len(subseries))
            if len(subseries)>minlength:
                break
            else:
                subseries = []
    return subseries

#Dry periods
datasize_threshold = 50
dry_all_time_series = []
for t in imposedtimeseries:
    subset = flowdata[flowdata["Time_series"]==t].reset_index()
    data = list(subset["Slower"][0])
    subdata = gensubtimeseries(data, datasize_threshold)
    dry_all_time_series.append(subdata)
df = pd.DataFrame.from_records(dry_all_time_series)
df.to_pickle("Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/dry_all_time_series")

#Wet periods
wet_all_time_series = []
for t in imposedtimeseries:
    subset = flowdata[flowdata["Time_series"]==t].reset_index()
    data = list(subset["Faster"][0])
    subdata = gensubtimeseries(data, datasize_threshold)
    wet_all_time_series.append(subdata)
df = pd.DataFrame.from_records(wet_all_time_series)
df.to_pickle("Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/wet_all_time_series")        

#Visualized velocities in all domains at these indices
palette = {0: "grey", 0.1: "orange", 1: "indianred", 5: "g", 10: "steelblue"}
style = {1: "solid", 2: "dashed", 5:"dotted", 10: "dashdot"}
fig, ax = plt.subplots(3,3,figsize = (8,8), sharey = 'row', sharex = 'col')
#plt.suptitle("Regime: " + Reg + "& Time series: " + t, fontsize = 18)
for Reg in reginvest:
    basedata = np.load("E:/Saturated_flow/EGUGoldschmidtdataset6/" + Reg + "AR_0/NS-AH/NS-AH_df.npy")
    basevelocity = np.mean(basedata[2, -1, :, :])
    for t in imposedtimeseries:
        directory = "E:/Saturated_flow/EGUGoldschmidtdataset6/" + Reg + "AR_" + t + "/"
        for j in Trial:
                a = ax.flat[reginvest.index(Reg)*len(imposedtimeseries) + imposedtimeseries.index(t)]
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    data = np.load(directory + "NS-A"+j+"/NS-A"+j+"_df.npy")
                    velocity = np.mean(data[2, wet_all_time_series[imposedtimeseries.index(t)][0]:wet_all_time_series[imposedtimeseries.index(t)][-1]+10, :, :], axis = (-1,-2))/basevelocity
                    a.plot(np.abs(velocity[1:]), color=palette[scdict[j]["Het"]], linestyle = style[scdict[j]["Anis"]], label = scdict[j]["Het"])
picname = "Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/wet_periods_velocities.png"
plt.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.1) 

#Records concentration at outlet for these periods and save in a dataframe

h5file = h5py.File("Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/Temporal_analysis_ratios.h5", mode = 'w')
for Reg in reginvest:
    basedata = np.load("E:/Saturated_flow/EGUGoldschmidtdataset6/" + Reg + "AR_0/NS-AH/NS-AH_df.npy")
    basevelocity = np.mean(basedata[2, -1, :, :])
    for t in imposedtimeseries:
        directory = "E:/Saturated_flow/EGUGoldschmidtdataset6/" + Reg + "AR_" + t + "/"
        print (Reg, t)
        for j in Trial:
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    basedata = np.load("E:/Saturated_flow/EGUGoldschmidtdataset6/" + Reg + "AR_0/NS-A"+j+"/NS-A"+j+"_df.npy")
                    baseconcs, flow, heads = sta.calcconcmasstimenew (basedata,0,-1,0,-1, 51, gvarnames, "Saturated")
                    data = np.load(directory + "NS-A"+j+"/NS-A"+j+"_df.npy")
                    concs, flow, heads = sta.calcconcmasstimenew (data,0,-1,0,-1, 51, gvarnames, "Saturated")
                    drysubconcs = concs[dry_all_time_series[imposedtimeseries.index(t)][0]:dry_all_time_series[imposedtimeseries.index(t)][-1]+10, -1, :]/baseconcs[-1, -1,:]
                    wetsubconcs = concs[wet_all_time_series[imposedtimeseries.index(t)][0]:wet_all_time_series[imposedtimeseries.index(t)][-1]+10, -1, :]/baseconcs[-1, -1,:]
                    for g in gvarnames:
                        h5file.create_dataset("dry_periods/"+ t + "/" + Reg + "/" + j + "/" + g, data=drysubconcs[:, gvarnames.index(g)])
                        h5file.create_dataset("wet_periods/"+ t + "/" + Reg + "/" + j + "/" + g, data=wetsubconcs[:, gvarnames.index(g)])                
h5file.close()

toplot = ["DOC", "TOC", "DO", "Nitrogen"]
hf = h5py.File("Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/Temporal_analysis.h5", mode = 'r')
print(hf.keys())
for Reg in reginvest:
    for t in imposedtimeseries:
        fig, ax = plt.subplots(2,2,figsize = (8,8), sharex = True)
        plt.suptitle("Regime: " + Reg + "& Time series: " + t, fontsize = 18)
        for g in toplot:
            a = ax.flat[toplot.index(g)]
            a.set_title(g, fontsize = 15)
            #num = 0
            for j in Trial:
                col = scdict[j]['Het']
                sty = scdict[j]["Anis"]
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    n = hf.get("dry_periods/"+ t + "/" + Reg + "/" + j + "/" + g)
                    #x_axis = list(range(len(n)))
                    a.plot(n[1:], color=palette[col], linestyle = style[sty], label = col)
        #a.legend(fontsize = 12)
        picname = "Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/dry_periods_" + Reg + "_" + t + ".png"
        plt.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)
hf.close()

#There is an anomaly in DO concentrations. Find it.

hf = h5py.File("Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/Temporal_analysis_ratios.h5", mode = 'r')
print(hf.keys())
row = []
for Reg in reginvest:
    if Reg == "Equal":
        r = "Medium"
    else:
        r = Reg
    for t in imposedtimeseries:
        for g in ["DO", "DOC"]:
            a = ax.flat[toplot.index(g)]
            a.set_title(g, fontsize = 15)
            #num = 0
            for j in Trial:
                col = scdict[j]['Het']
                sty = scdict[j]["Anis"]
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    n = hf.get("dry_periods/"+ t + "/" + Reg + "/" + j + "/" + g)
                    #if any(n[:]) > 2:
                    row.append([Reg, t, j, g, max(n[:])])
df = pd.DataFrame.from_records(row)

#Consider time series in terms of Dat.
#Sort all the values in Dat
#Then take average at each time point
#Then take rolling mean for the averaged time series
data  = pd.read_csv("Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/mass_flux_sensitivity_generalized.csv", sep="\t")
gvarnames = ["DOC", "DO", "Nitrogen", "TOC"]
finaldata = data[data['Chem'].isin (gvarnames)]
mymarklist = ["^", "o", "s", "d"]
reglist = ["Slow", "Medium", "Fast"]
colorlist = ["indianred", "g", "steelblue"]

finaldata.loc[finaldata["PeDa"] > 40, "PeDamark"] = 3
finaldata.loc[(finaldata["PeDa"] > 15) & (finaldata["PeDa"] < 40), "PeDamark"] = 2
finaldata.loc[(finaldata["PeDa"] > 1) & (finaldata["PeDa"] < 15), "PeDamark"] = 1
finaldata.loc[finaldata["PeDa"] < 1, "PeDamark"] = 0
labels = {0 : r'$Da_t < 1$',
          1 : r'$1 < Da_t < 15$',
          2 : r'$15 < Da_t < 40$',
          3 : r'$Da_t > 40$'}

subfinal = finaldata[["Trial", "Regime", "Chem", "Time_series","PeDamark"]]
subfinal['key'] = subfinal.Trial + subfinal.Regime + subfinal.Chem + subfinal.Time_series.astype(str)
colorcriteria = subfinal[["PeDamark", "key"]].to_dict('records')
colorcriteria = dict(zip(subfinal['key'], subfinal['PeDamark']))

PeDapalette = {0: "grey", 1: "orange", 2: "g", 3: "indianred"}

hr = h5py.File("Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/Temporal_analysis_ratios.h5", mode = 'r')
hw = h5py.File("Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/Temporal_analysis_ratios_Dat.h5", mode = 'w')
print(hr.keys())
for timdesc in ["wet", "dry"]:
    for t in imposedtimeseries:
        n0=[]
        n1=[]
        n2=[]
        n3=[]
        for Reg in reginvest:
            if Reg == "Equal":
                r = "Medium"
            else:
                r = Reg    
            for g in gvarnames:
                for j in Trial:
                    if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                        pass
                    else:                   
                        n = hr.get(timdesc+"_periods/"+ t + "/" + Reg + "/" + j + "/" + g).value
                        if int(colorcriteria[j+r+g+t]) == 0:
                            n0.append(n)
                        elif int(colorcriteria[j+r+g+t]) == 1: 
                            n1.append(n)
                        elif int(colorcriteria[j+r+g+t]) == 2:
                            n2.append(n)
                        elif int(colorcriteria[j+r+g+t]) == 3:
                            n3.append(n)
        for Dat, k in zip([0,1,2,3],[n0, n1, n2, n3]):
            df = pd.DataFrame.from_records(k)
            hw.create_dataset(timdesc + "_periods/"+ t + "/Dat" + str(Dat) + "/", data=df.mean())
hr.close()
hw.close()

#Calculate rolling average for each Dat in each time series
periods = ["wet", "dry"]
hr = h5py.File("Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/Temporal_analysis_ratios_Dat.h5", mode = 'r')
basedata = np.load("E:/Saturated_flow/EGUGoldschmidtdataset6/SlowAR_0/NS-AH/NS-AH_df.npy")
basevelocity = np.mean(basedata[2, -1, :, :])
fig, ax = plt.subplots(2,3, figsize = (10,8), sharey = True)
plt.suptitle("Select time windows")
for timdesc in periods:
    if timdesc == "wet":
        scdata = wet_all_time_series
    else:
        scdata = dry_all_time_series
    for t in imposedtimeseries:
        a = ax.flat[periods.index(timdesc)*len(imposedtimeseries) + imposedtimeseries.index(t)]
        a.set_title(t)
        for Dat in [0,1,2,3]:
            n = hr.get(timdesc + "_periods/"+ t + "/Dat" + str(Dat) + "/").value
            print(np.shape(n))
            a.plot(n[1:], label = Dat)
            #a.plot(np.convolve(n[1:], np.ones((3,))/3,mode = 'valid'), label = Dat)
        directory = "E:/Saturated_flow/EGUGoldschmidtdataset6/SlowAR_" + t + "/"
        for j in ["37"]:
            a = ax.flat[periods.index(timdesc)*len(imposedtimeseries) + imposedtimeseries.index(t)]
            data = np.load(directory + "NS-A"+j+"/NS-A"+j+"_df.npy")
            velocity = np.mean(data[2, scdata[imposedtimeseries.index(t)][0]:scdata[imposedtimeseries.index(t)][-1]+10, :, :], axis = (-1,-2))/basevelocity
            a.plot(np.abs(velocity[1:]), color= "gray")
a.legend()
hr.close()

#Do the above steps for the entire time series
directory1 = "C:/Users/khurana/Documents/Holiday_homework/"
h5file = h5py.File(directory1 + "Paper2/Figurecodes/Temporal_analysis_full.h5", mode = 'w')
for Reg in reginvest:
    basedata = np.load("E:/Saturated_flow/EGUGoldschmidtdataset6/" + Reg + "AR_0/NS-AH/NS-AH_df.npy")
    basevelocity = np.mean(basedata[2, -1, :, :])
    for t in imposedtimeseries:
        directory = "E:/Saturated_flow/EGUGoldschmidtdataset6/" + Reg + "AR_" + t + "/"
        print (Reg, t)
        for j in Trial:
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    basedata = np.load("E:/Saturated_flow/EGUGoldschmidtdataset6/" + Reg + "AR_0/NS-A"+j+"/NS-A"+j+"_df.npy")
                    baseconcs, flow, heads = sta.calcconcmasstimenew (basedata,0,-1,0,-1, 51, gvarnames, "Saturated")
                    data = np.load(directory + "NS-A"+j+"/NS-A"+j+"_df.npy")
                    concs, flow, heads = sta.calcconcmasstimenew (data,0,-1,0,-1, 51, gvarnames, "Saturated")
                    subconcs = concs[1:, -1, :]/baseconcs[-1, -1,:]
                    for g in gvarnames:
                        h5file.create_dataset(t + "/" + Reg + "/" + j + "/" + g, data=subconcs[:, gvarnames.index(g)])
h5file.close()

#Consider time series in terms of Dat.
#Sort all the values in Dat
#Then take average at each time point
#Then take rolling mean for the averaged time series
data  = pd.read_csv(directory1 + "Paper1/Figurecodes/mass_flux_sensitivity_generalized.csv", sep="\t")
gvarnames = ["DOC", "DO", "Nitrogen", "TOC"]
finaldata = data[data['Chem'].isin (gvarnames)]
mymarklist = ["^", "o", "s", "d"]
reglist = ["Slow", "Medium", "Fast"]
colorlist = ["indianred", "g", "steelblue"]

finaldata.loc[finaldata["PeDa"] > 40, "PeDamark"] = 3
finaldata.loc[(finaldata["PeDa"] > 15) & (finaldata["PeDa"] < 40), "PeDamark"] = 2
finaldata.loc[(finaldata["PeDa"] > 1) & (finaldata["PeDa"] < 15), "PeDamark"] = 1
finaldata.loc[finaldata["PeDa"] < 1, "PeDamark"] = 0
labels = {0 : r'$Da_t < 1$',
          1 : r'$1 < Da_t < 15$',
          2 : r'$15 < Da_t < 40$',
          3 : r'$Da_t > 40$'}

subfinal = finaldata[["Trial", "Regime", "Chem", "Time_series","PeDamark"]]
subfinal['key'] = subfinal.Trial + subfinal.Regime + subfinal.Chem + subfinal.Time_series.astype(str)
colorcriteria = subfinal[["PeDamark", "key"]].to_dict('records')
colorcriteria = dict(zip(subfinal['key'], subfinal['PeDamark']))

PeDapalette = {0: "grey", 1: "orange", 2: "g", 3: "indianred"}

hr = h5py.File(directory1 + "Paper2/Figurecodes/Temporal_analysis_full.h5", mode = 'r')
hw = h5py.File(directory1 + "Paper2/Figurecodes/Temporal_analysis_full_max_Dat.h5", mode = 'w')
print(hr.keys())

for t in imposedtimeseries:
    n0=[]
    n1=[]
    n2=[]
    n3=[]
    for Reg in reginvest:
        if Reg == "Equal":
            r = "Medium"
        else:
            r = Reg    
        for g in gvarnames:
            for j in Trial:
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:                   
                    n = hr.get(t + "/" + Reg + "/" + j + "/" + g).value
                    if int(colorcriteria[j+r+g+t]) == 0:
                        n0.append(n)
                    elif int(colorcriteria[j+r+g+t]) == 1: 
                        n1.append(n)
                    elif int(colorcriteria[j+r+g+t]) == 2:
                        n2.append(n)
                    elif int(colorcriteria[j+r+g+t]) == 3:
                        n3.append(n)
    for Dat, k in zip([0,1,2,3],[n0, n1, n2, n3]):
        df = pd.DataFrame.from_records(k)
        hw.create_dataset(t + "/Dat" + str(Dat) + "/", data=df.mean())
hr.close()
hw.close()

#Calculate rolling average for each Dat in each time series
hr = h5py.File(directory1 + "Paper2/Figurecodes/Temporal_analysis_full_max_Dat.h5", mode = 'r')
basedata = np.load("E:/Saturated_flow/EGUGoldschmidtdataset6/SlowAR_0/NS-AH/NS-AH_df.npy")
basevelocity = np.mean(basedata[2, -1, :, :])
fig, ax = plt.subplots(3, 1, figsize = (10,8), sharey = True)
for t in imposedtimeseries:
    a = ax.flat[imposedtimeseries.index(t)]
    a.set_title("Time series: T"+ t, fontsize = 16)
    for Dat in [0,1,2,3]:
        n = hr.get(t + "/Dat" + str(Dat) + "/").value
        print(np.shape(n))
        a.plot(n[1:], label = labels[Dat])
        #a.plot(np.convolve(n[1:], np.ones((3,))/3,mode = 'valid'), label = Dat)
    directory = "E:/Saturated_flow/EGUGoldschmidtdataset6/SlowAR_" + t + "/"
    for j in ["37"]:
        a = ax.flat[imposedtimeseries.index(t)]
        data = np.load(directory + "NS-A"+j+"/NS-A"+j+"_df.npy")
        velocity = np.mean(data[2, :, :, :], axis = (-1,-2))/basevelocity
        a.plot(np.abs(velocity[1:]), color= "gray")
a.legend()
hr.close()
picname = directory1 + "Paper2/Figurecodes/Temporal_analysis_full_data.png"
plt.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)