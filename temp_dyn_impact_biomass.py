# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 21:18:49 2020

@author: khurana
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 12:26:03 2019

@author: khurana
"""
import os
import numpy as np
import h5py
import pandas as pd
import data_reader.data_processing as proc
import analyses.transient as ta
import matplotlib.pyplot as plt

sourcedatadir = "E:/Saturated_flow/EGUGoldschmidtdataset6"
hdf5directory = "Y:/Home/khurana/4. Publications/Restructuring"

#set up basic constants 
Regimes = ["Slow", "Equal", "Fast"]
domainodes = {"Original": {'ynodes' : 51},
              "Big" : {'ynodes' : 126},
              "Double" : {'ynodes' : 101},
              "Half" : {'ynodes' : 26}}
scdict = proc.masterscenarios("Saturated") #master dictionary of all spatially heterogeneous scenarios that were run

# Scenarios to investigate:
Trial = list(scdict.keys())
reginvest = Regimes
domaininvest = list(domainodes.keys())[:1]
imposedtimeseries = ["1","2","5"]

vardict = proc.speciesdict("Saturated")
States = ["Active", "Inactive"]
gvarnames = list(t for t in vardict.keys() if vardict[t]["State"] in (States))
regcolors = {"Slow":"indianred", "Equal":"g", "Fast":"steelblue"}
regmapcolors = {"Slow":plt.cm.Reds, "Equal":plt.cm.Greens, "Fast":plt.cm.Blues}

h5file = h5py.File(os.path.join(hdf5directory, "Paper2", "Figurecodes","Temporal_analysis_biomass.h5"), mode = 'w')
for Reg in Regimes:
    basedata = np.load(os.path.join(sourcedatadir, Reg + "AR_0/NS-AH/NS-AH_df.npy"))
    basevelocity = np.mean(basedata[2, -1, :, :])
    for t in imposedtimeseries:
        directory = os.path.join(sourcedatadir,Reg + "AR_" + t)
        print (Reg, t)
        for j in Trial:
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    basefile = os.path.join(Reg + "AR_0","NS-A"+j,"NS-A"+j+"_df.npy")
                    basedata = np.load(os.path.join(sourcedatadir, basefile))
                    baseconcs = ta.biomass_time (basedata,0,-1,0,-1, gvarnames, "Saturated")
                    data = np.load(os.path.join(directory, "NS-A"+j,"NS-A"+j+"_df.npy"))
                    concs = ta.biomass_time (data,0,-1,0,-1, gvarnames, "Saturated")
                    subconcs = concs[1:, :]/baseconcs[-1, :]
                    for g in gvarnames:
                        dataset_name = t + "/" + Reg + "/" + j + "/" + g
                        print(dataset_name)
                        h5file.create_dataset(dataset_name, data=subconcs[:, gvarnames.index(g)])
h5file.close()

basedata = np.load(os.path.join(sourcedatadir, "SlowAR_0/NS-AH/NS-AH_df.npy"))
basevelocity = np.mean(basedata[2, -1, :, :])
nosulfate = list(g for g in gvarnames if "sulphate" not in g)
species = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
position = ["Immobile", "Mobile"]
state = ["Active", "Inactive"]

hr = h5py.File(os.path.join(hdf5directory, "Paper2","Figurecodes","Temporal_analysis_biomass.h5"), mode = 'r')
for Reg in Regimes:
    if Reg == "Equal":
        r = "Medium"
    else:
        r = Reg 
    fig,ax = plt.subplots(nrows = 4, ncols = 3, figsize = (10,10), sharex = True, sharey = True)
    for g in nosulfate:
        axeidx = nosulfate.index(g)
        axe = ax.flat[axeidx]
        axe.axhline(y=1.2, color = "black", linestyle = "dashed")
        axe.axhline(y=0.80, color = "black", linestyle = "dashed")
        axe.axvline(x=1.2, color = "black", linestyle = "dashed")
        axe.axvline(x=0.80, color = "black", linestyle = "dashed")
        for t in imposedtimeseries:
            directory = os.path.join(sourcedatadir, "SlowAR_" + t)
            for j in ["37"]:
                data = np.load(os.path.join(directory,"NS-A"+j,"NS-A"+j+"_df.npy"))
                velocity = np.mean(data[2, :, :, :], axis = (-1,-2))/basevelocity
            n_reg = []
            for j in Trial:
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    n = hr.get(t + "/" + Reg + "/" + j + "/" + g).value
                    n_reg.append(n)
            datapv = np.asarray(n_reg).transpose()
            for col in list(range(np.shape(datapv)[1])):
                columndata = datapv[:,col]
                colorindex = np.argwhere((columndata<=0.8) | (columndata>=1.2))
                grayindex = np.argwhere((columndata>0.8) & (columndata<1.2))
                axe.scatter(velocity[colorindex],columndata[colorindex], color = regcolors[Reg], alpha = 0.5, s = 0.5, edgecolor = None)
                axe.scatter(velocity[grayindex],columndata[grayindex], color = "gainsboro", alpha = 0.2,s = 0.5, edgecolor = "white")
        axe.tick_params(labelsize = 14)
        plt.yscale("log")
        if axeidx>=(len(nosulfate)-3):
            axe.set_xlabel("Velocity (-)", fontsize = 14)
        plt.xticks ((0,0.5,1,1.5,2),(0,0.5,1,1.5,2),fontsize = 14)
    for a in [0,1,2,3]:
        ax[a,0].set_ylabel(state[a%2], fontsize = 14)
    for a in [0,1]:
        ix = [a*2,a*2 + 1]
        for i in ix:
            ax[i,-1].annotate(position[a],
                              xy=(0, 0.5),
                              xytext=(180, 0),
                              xycoords="axes fraction",
                              textcoords="offset points",
                              size="large",
                              ha="left",
                              va="center",
                              rotation="vertical",
                              fontsize=14)
    plt.annotate("Normalized biomass",
        xy=(-0.1, 2.5),
        xytext=(-450, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        rotation="vertical",
        fontsize=16)
    plt.yticks ((0.02,0.1,1,10),(0.02,0.1,1,10),fontsize = 14)
    plt.savefig(os.path.join(hdf5directory,"Paper2","Figurecodes",r+"_Temp_variation_biomass.png"),
    dpi=300,
    bbox_inches="tight",
    pad=0.1)

fig,ax = plt.subplots(nrows = 4, ncols = 3, figsize = (10,10), sharex = True, sharey = True)
for g in nosulfate:
    axeidx = nosulfate.index(g)
    axe = ax.flat[axeidx]
    if axeidx<3:
        axe.set_title(species[axeidx], fontsize = 14)
    axe.axhline(y=1.2, color = "black", linestyle = "dashed")
    axe.axhline(y=0.80, color = "black", linestyle = "dashed")
    axe.axvline(x=1.2, color = "black", linestyle = "dashed")
    axe.axvline(x=0.80, color = "black", linestyle = "dashed")
    for t in imposedtimeseries:
        directory = os.path.join(sourcedatadir, "SlowAR_" + t)
        for j in ["37"]:
            data = np.load(os.path.join(directory,"NS-A"+j,"NS-A"+j+"_df.npy"))
            velocity = np.mean(data[2, :, :, :], axis = (-1,-2))/basevelocity
        for Reg in Regimes:
            if Reg == "Equal":
                r = "Medium"
            else:
                r = Reg 
            n_reg = []
            for j in Trial:
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    n = hr.get(t + "/" + Reg + "/" + j + "/" + g).value
                    n_reg.append(n)
            datapv = np.asarray(n_reg).transpose()
            for col in list(range(np.shape(datapv)[1])):
              columndata = datapv[:,col]
              sampled = np.random.choice(len(columndata), 100,replace = False)
              colorindex = sampled[np.argwhere((columndata[sampled]<=0.8) | (columndata[sampled]>=1.2))]
              grayindex = sampled[np.argwhere((columndata[sampled]>0.8) & (columndata[sampled]<1.2))]
              axe.scatter(velocity[colorindex],columndata[colorindex], color = regcolors[Reg], alpha = 0.5, s = 0.5, edgecolor = None)
              axe.scatter(velocity[grayindex],columndata[grayindex], color = "gainsboro", alpha = 0.2,s = 0.5, edgecolor = "white")
    axe.tick_params(labelsize = 14)
    plt.yscale("log")
    if axeidx>=(len(nosulfate)-3):
        axe.set_xlabel("Velocity (-)", fontsize = 14)
    plt.xticks ((0,0.5,1,1.5,2),(0,0.5,1,1.5,2),fontsize = 14)
for a in [0,1,2,3]:
    ax[a,0].set_ylabel(position[a%2], fontsize = 14)
for a in [0,1]:
    ix = [a*2,a*2 + 1]
    for i in ix:
        ax[i,-1].annotate(state[a],
                          xy=(0, 0.5),
                          xytext=(180, 0),
                          xycoords="axes fraction",
                          textcoords="offset points",
                          size="large",
                          ha="left",
                          va="center",
                          rotation="vertical",
                          fontsize=14)
plt.annotate("Normalized biomass",
             xy=(-0.1, 2.5),
             xytext=(-450, 0),
             xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        rotation="vertical",
        fontsize=16)
plt.yticks ((0.02,0.1,1,10),(0.02,0.1,1,10),fontsize = 14)
plt.savefig(os.path.join(hdf5directory,"Paper2","Figurecodes",r+"_Temp_variation_sampled.png"),
            dpi=300,
            bbox_inches="tight",
            pad=0.1)

### Distribute biomass subpopulations between just state/activity
States = ["Active", "Inactive"]
Locations = ["Immobile", "Mobile"]
gvarnames = list(t for t in vardict.keys() if vardict[t]["State"] in (States))
actgvar = list(t for t in gvarnames if vardict[t]["State"]=="Active")
inactgvar = list(t for t in gvarnames if vardict[t]["State"]=="Inactive")
immgvar = list(t for t in gvarnames if vardict[t]["Location"]=="Immobile")
mobgvar = list(t for t in gvarnames if vardict[t]["Location"]=="Mobile")

h5subpops = h5py.File(os.path.join(hdf5directory, "Paper2", "Figurecodes","Temporal_analysis_biomass_subpops.h5"), mode = 'w')
h5ratiopops = h5py.File(os.path.join(hdf5directory, "Paper2", "Figurecodes","Temporal_analysis_biomass_ratiopops.h5"), mode = 'w')
for Reg in Regimes:
    basedata = np.load(os.path.join(sourcedatadir, Reg + "AR_0/NS-AH/NS-AH_df.npy"))
    basevelocity = np.mean(basedata[2, -1, :, :])
    for t in imposedtimeseries:
        directory = os.path.join(sourcedatadir,Reg + "AR_" + t)
        print (Reg, t)
        for j in Trial:
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    basefile = os.path.join(Reg + "AR_0","NS-A"+j,"NS-A"+j+"_df.npy")
                    basedata = np.load(os.path.join(sourcedatadir, basefile))
                    bact = ta.biomass_time (basedata,0,-1,0,-1, actgvar, "Saturated")[-1,:]
                    binact = ta.biomass_time (basedata,0,-1,0,-1, inactgvar, "Saturated")[-1,:]
                    bimm = ta.biomass_time (basedata,0,-1,0,-1, immgvar, "Saturated")[-1,:]
                    bmob = ta.biomass_time (basedata,0,-1,0,-1, mobgvar, "Saturated")[-1,:]

                    data = np.load(os.path.join(directory, "NS-A"+j,"NS-A"+j+"_df.npy"))
                    act = ta.biomass_time (data,0,-1,0,-1, actgvar, "Saturated")
                    inact = ta.biomass_time (data,0,-1,0,-1, inactgvar, "Saturated")
                    imm = ta.biomass_time (data,0,-1,0,-1, immgvar, "Saturated")
                    mob = ta.biomass_time (data,0,-1,0,-1, mobgvar, "Saturated")

                    for subpop,data_array, base in zip(["Active", "Inactive", "Immobile", "Mobile"],[act,inact,imm,mob], [bact,binact,bimm,bmob]):
                        dataset_name = t + "/" + Reg + "/" + j + "/" + subpop
                        print(dataset_name)
                        subconcs = np.sum(data_array, axis = 1)
                        basesum = np.sum(base)
                        h5subpops.create_dataset(dataset_name, data=subconcs)
                        h5subpops.create_dataset("Norm_"+dataset_name, data=subconcs/basesum)

                    actinact = np.sum(act, axis = 1)/np.sum(inact, axis = 1)
                    immmob = np.sum(imm, axis = 1)/np.sum(mob, axis = 1)
                    bactinact = np.sum(bact)/np.sum(binact)
                    bimmmob = np.sum(bimm)/np.sum(bmob)
                    
                    for subpop,data_array, base in zip(["State_Ratio", "Location_Ratio"],[actinact,immmob], [bactinact, bimmmob]):
                        dataset_name = t + "/" + Reg + "/" + j + "/" + subpop
                        print(dataset_name)
                        h5ratiopops.create_dataset(dataset_name, data=data_array)
                        h5ratiopops.create_dataset("Norm_"+dataset_name, data=data_array/base)
h5subpops.close()
h5ratiopops.close()

basedata = np.load(os.path.join(sourcedatadir, "SlowAR_0/NS-AH/NS-AH_df.npy"))
basevelocity = np.mean(basedata[2, -1, :, :])
ratio_plot = ["State_Ratio", "Location_Ratio","Active", "Inactive", "Immobile", "Mobile"]
hrsubpops = h5py.File(os.path.join(hdf5directory, "Paper2", "Figurecodes","Temporal_analysis_biomass_subpops.h5"), mode = 'r')
hrsubratio = h5py.File(os.path.join(hdf5directory, "Paper2", "Figurecodes","Temporal_analysis_biomass_ratiopops.h5"), mode = 'r')
#hr = hrsubpops
fig,ax = plt.subplots(nrows = 6, ncols = 3, figsize = (8,12), sharex = True, sharey = True)
ax.flat[0].set_ylabel(r"$\frac{Active}{Inactive}^\ast$", fontsize = 14)
ax.flat[3].set_ylabel(r"$\frac{Immobile}{Mobile}^\ast$", fontsize = 14)
for g in ratio_plot:
    grow = ratio_plot.index(g)
    if grow>1:
        hr = hrsubpops
    else:
        hr = hrsubratio
    for Reg in Regimes:
        if Reg == "Equal":
            r = "Medium"
        else:
            r = Reg 
        data_toplot = pd.DataFrame(columns = ["x","y"])
        for t in imposedtimeseries:
            tdirectory = os.path.join(sourcedatadir, "SlowAR_" + t)
            n_reg = []
            for j in ["37"]:
                data = np.load(os.path.join(tdirectory,"NS-A"+j,"NS-A"+j+"_df.npy"))
                velocity = np.mean(data[2, :, :, :], axis = (-1,-2))/basevelocity    
            for j in Trial:
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    n = hr.get("Norm_"+t + "/" + Reg + "/" + j + "/" + g).value
                    n_reg.append(n)
            datapv = np.asarray(n_reg).transpose()
            x_plot = np.tile(velocity[1:], np.shape(datapv)[1])
            y_plot = datapv.flatten('F')
            data_reg = pd.DataFrame({"x":x_plot, "y":y_plot})
            print(g,t,r)
            data_toplot = data_toplot.append(data_reg)
        axeidx = grow*len(Regimes) + Regimes.index(Reg)
        axe = ax.flat[axeidx]
        if axeidx<3:
            axe.set_title(r+"\nFlow", fontsize = 14)
        else:
            axe.set_xlabel(r"$Velocity^\ast$", fontsize =14)
        axe.hist2d(data_toplot["x"],data_toplot["y"], bins = (50,50), cmap = regmapcolors[Reg])
        axe.tick_params(labelsize = 14)
plt.yticks ((0,0.5,1,1.5,2),(0,0.5,1,1.5,2),fontsize = 14)
plt.xticks ((0,0.5,1,1.5,2),(0,0.5,1,1.5,2),fontsize = 14)
plt.savefig(os.path.join(hdf5directory, "Paper2", "Figurecodes","FigureXX_Temp_variation_biomass_ratios.png"),
            dpi=300,
            bbox_inches="tight", pad=0.1)
hr.close()

fig,ax = plt.subplots(nrows = 1, ncols = 2, figsize = (8,4), sharex = True, sharey = True)
ax.flat[0].set_ylabel(r"$\frac{Active}{Inactive}^\ast$", fontsize = 14)
#ax.flat[3].set_ylabel(r"$\frac{Immobile}{Mobile}^\ast$", fontsize = 14)
plt.xscale("log")
for g in ratio_plot[:2]:
    grow = ratio_plot.index(g)
    if grow>1:
        hr = hrsubpops
    else:
        hr = hrsubratio
    data_toplot = pd.DataFrame(columns = ["x","y"])
    n_reg = []
    for Reg in Regimes:
        if Reg == "Equal":
            r = "Medium"
        else:
            r = Reg 
        for t in imposedtimeseries:
            tdirectory = os.path.join(sourcedatadir, Reg+"AR_" + t)
            for j in ["37"]:
                data = np.load(os.path.join(tdirectory,"NS-A"+j,"NS-A"+j+"_df.npy"))
                velocity = np.abs(np.mean(data[2, :, :, :], axis = (-1,-2)))#/basevelocity
            for j in Trial:
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    n = hr.get("Norm_"+t + "/" + Reg + "/" + j + "/" + g).value
                    n_reg.append(n)
            datapv = np.asarray(n_reg).transpose()
            x_plot = np.tile(velocity[1:], np.shape(datapv)[1])
            y_plot = datapv.flatten('F')
            data_reg = pd.DataFrame({"x":x_plot, "y":y_plot})
            print(g,t,r)
            data_toplot = data_toplot.append(data_reg)
        axeidx = grow#*len(Regimes) + Regimes.index(Reg)
        axe = ax.flat[axeidx]
    axe.hist2d(data_toplot["x"],data_toplot["y"], bins = (50,50), cmap = regmapcolors[Reg])
    axe.tick_params(labelsize = 14)
plt.yticks ((0,0.5,1,1.5,2),(0,0.5,1,1.5,2),fontsize = 14)
#plt.xticks ((0,0.5,1,1.5,2),(0,0.5,1,1.5,2),fontsize = 14)
plt.savefig(os.path.join(hdf5directory, "Paper2", "Figurecodes","FigureXX_Temp_variation_biomass_ratios.png"),
            dpi=300,
            bbox_inches="tight", pad=0.1)
hr.close()

