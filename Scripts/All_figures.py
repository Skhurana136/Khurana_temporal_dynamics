# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 17:36:05 2020

@author: khurana
"""
#%%
# loading required libraries
import os
import pandas as pd
import numpy as np
import h5py
from statsmodels.tsa import stattools

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import DS.plots.saturated_transient as stp
import DS.data_reader.data_processing as proc

#%%
# Set up directories
DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
results_dir = os.path.join(DIR,"Results")
sourcedatadir = "D:/Data/Zenodo_temporal_dynamics"
os.chdir(results_dir)

legendkw = {'fontsize' : 14}
labelkw = {'labelsize' : 14}
secondlabelkw = {'labelsize' : 16}
suptitlekw = {'fontsize' : 18}
titlekw = {'fontsize' : 16}
mpl.rc('font',family='Arial')

#---- FIGURES ----

#%%
#----  Figure 1 Temporal variation: Chemical species ----

# Define commonly used scenarios
Regimes = ["Slow", "Equal", "Fast"]
reglist = ["Slow", "Medium", "Fast"]
gvarnames = ["DO", "DOC", "Nitrate", "Ammonium"]#, "TOC", "Nitrogen"]
imposedtimeseries = ["1","2","5"]
Trial = proc.masterscenarios("Saturated").keys()
blue_patch = mpatches.Patch(color="steelblue", label= 'log$_{10}$Da < -1', alpha = 0.8)
orange_patch = mpatches.Patch(color = "orange", label =  '-1 < log$_{10}$Da < 0', alpha = 0.8)
green_patch = mpatches.Patch(color="g", label='0 < log$_{10}$Da < 0.5', alpha = 0.8)
red_patch = mpatches.Patch(color="indianred", label='log$_{10}$Da > 0.5', alpha = 0.8)
patchlist = [blue_patch, orange_patch, green_patch, red_patch]

# Define plotting related parameters
colorlist = ["indianred", "g", "steelblue"]
PeDapalette = {0: "steelblue", 1: "orange", 2: "g", 3: "indianred"}

# Generate Dat from Da and Pe
data  = pd.read_csv("Da_29012021_95pcloss.csv", sep=",")
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
basedata = np.load(os.path.join(sourcedatadir, "SlowAR_0_NS-AH_df.npy"))
basevelocity = np.mean(basedata[2, -1, :, :])
hr = h5py.File("Temporal_analysis_full_data.h5", mode = 'r')
fig,ax = plt.subplots(2,2, figsize = (7,7), sharex = True, sharey = True)
for a,b in zip(ax.flat[:], ["A", "B", "C", "D"]):
    a.text(0.1, 10, b, c = "black", fontsize = 'xx-large')
for danum in [0,1,2,3]:
    n_danum = []
    vel_all3 = []
    axe = ax.flat[danum]
    axe.set_title(labels[danum], fontsize = 16)
    axe.axhline(y=1.2, color = "black", linestyle = "dashed")
    axe.axhline(y=0.80, color = "black", linestyle = "dashed")
    axe.axvline(x=1.2, color = "black", linestyle = "dashed")
    axe.axvline(x=0.80, color = "black", linestyle = "dashed")
    for t in imposedtimeseries:
        n_danumt = []
        for j in ["37"]:
            data = np.load(os.path.join(sourcedatadir,"SlowAR_"+t+"_NS-A"+j+"_df.npy"))
            velocity = np.mean(data[2, :, :, :], axis = (-1,-2))/basevelocity
        vel_all3.append(velocity)
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
                        n = hr.get(t + "/" + Reg + "/" + j + "/" + g)[:]
                        n_danum.append(n)
                        n_danumt.append(n)
                    else:
                        pass
        datapv = np.asarray(n_danumt).transpose()
        for col in list(range(np.shape(datapv)[1])):
            columndata = datapv[:,col]
            colorindex = np.argwhere((columndata<=0.8) | (columndata>=1.2))
            grayindex = np.argwhere((columndata>0.8) & (columndata<1.2))
            axe.scatter(velocity[colorindex],columndata[colorindex], color = PeDapalette[danum], alpha = 0.5, s = 0.5, edgecolor = None)
            axe.scatter(velocity[grayindex],columndata[grayindex], color = "gainsboro", alpha = 0.4,s = 0.3, facecolor = None, edgecolor = "grey")
    axe.tick_params(labelsize = 14)
    axe.set_yscale("log")
    if danum %2==0:
        axe.set_ylabel("Concentration (-)", fontsize = 14)
    if danum >1:
        axe.set_xlabel("Velocity (-)", fontsize = 14)
    plt.xticks ((0,0.5,1,1.5,2),(0,0.5,1,1.5,2),fontsize = 14)
    plt.yticks ((0.01,0.1,1,10),(0.01,0.1,1,10),fontsize = 14)

plt.legend(handles=patchlist, ncol = 2,
        bbox_to_anchor=(0.01, -0.48),
        loc="center",
        title="Reactive systems",
        fontsize=14, title_fontsize = 14)
picname = "Figure1_temp_variation_chem.jpg"
plt.savefig(picname, dpi=300, bbox_inches="tight", pad=0.1)

#%%
# Figure 2: Time series: chemical
#Characteristic to plot:
features = ["median"]
featurestyle = ["solid"]
imposedtimeseries = ["1","2","5"]
veline = mlines.Line2D([], [], linestyle = 'solid', color='grey', markersize=15, label='Velocity')
da0 = mlines.Line2D([], [], linestyle = 'solid', color='steelblue', markersize=15, label='log$_{10}$Da < -1')
da1 = mlines.Line2D([], [], linestyle = 'solid', color='orange', markersize=15, label='-1 < log$_{10}$Da < 0')
da2 = mlines.Line2D([], [], linestyle = 'solid', color='g', markersize=15, label='0 < log$_{10}$Da < 0.5')
da3 = mlines.Line2D([], [], linestyle = 'solid', color='indianred', markersize=15, label='log$_{10}$Da > 0.5')

#Load dataset for time series in terms of Dat
hr = h5py.File("Temporal_analysis_full_Dat.h5", mode = 'r')
#Load dataset for base velocity ratio values which will be plotted in the graph
basedata = np.load(os.path.join(sourcedatadir, "SlowAR_0_NS-AH_df.npy"))
basevelocity = np.mean(basedata[2, -1, :, :])
fig, ax = plt.subplots(4, 3, figsize = (16,8), sharey = True, sharex = True)
for col,num in zip ([0,1,2],imposedtimeseries):
    for a,b in zip(ax[:,col], ["A", "B", "C", "D"]):
        a.text(0.2, 2.0, b+num, c = "black", fontsize = 'large')

for t in imposedtimeseries:
    for j in ["37"]:
        data = np.load(os.path.join(sourcedatadir,"SlowAR_" + t+"_NS-A"+j+"_df.npy"))
        velocity = np.mean(data[2, :, :, :], axis = (-1,-2))/basevelocity
        for a in ax[:,imposedtimeseries.index(t)]:
            a.plot(np.abs(velocity[1:]), color= "gray", alpha = 0.7)
    for Dat,pedacol in zip([0,1,2,3],PeDapalette):
        axidx = Dat*len(imposedtimeseries) + imposedtimeseries.index(t)
        a = ax.flat[axidx]
        for datafeature,featureline in zip(features,featurestyle):
            n = hr.get(t + "/Dat" + str(Dat) + "/"+datafeature+"/")[:]
            nmin = hr.get(t + "/Dat" + str(Dat) + "/q25/")[:]
            nmax = hr.get(t + "/Dat" + str(Dat) + "/q75/")[:]
            print(np.shape(n))
            a.plot(n[1:], label = labels[Dat], color = PeDapalette[Dat], linestyle=featureline)
            a.fill_between(np.asarray(list(range(1095))),nmin, nmax, color = PeDapalette[Dat], alpha = 0.4)
            a.tick_params(labelsize = 14)
        if axidx <3:
            a.set_title("Time series: T"+ t, fontsize = 16)
plt.xticks((0,365,730,1095), (0,5,10,15))
legend_all = fig.legend(handles=[da0, da1, da2, da3, veline],
                        fontsize = 14, bbox_to_anchor=(0.5, 0.0),
                         ncol=5,loc="center")
plt.gca().add_artist(legend_all)
plt.annotate("Ratio with respect to steady state conditions",
        xy=(-2.6, 2.25),
        xytext=(0, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=18)
plt.annotate("Time (years)",
        xy=(-1.0,-0.3),
        xytext=(0, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        fontsize=18)
picname = "Figure2_Chem_time_series.jpg"
plt.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)
hr.close()

#%%
# Figure 3: Temporal variation: Microbial species
regmapcolors = {"Slow":plt.cm.Reds, "Equal":plt.cm.Greens, "Fast":plt.cm.Blues}
basedata = np.load(os.path.join(sourcedatadir, "SlowAR_0_NS-AH_df.npy"))
basevelocity = np.mean(basedata[2, -1, :, :])
ratio_plot = ["State_Ratio", "Location_Ratio"]
hrsubratio = h5py.File("Temporal_analysis_biomass_ratiopops.h5", mode = 'r')
#hr = hrsubpops
fig,ax = plt.subplots(nrows = 2, ncols = 3, figsize = (6,4), sharex = True, sharey = True)
for row,num in zip ([0,1],["1", "2"]):
    for a,b in zip(ax[row,:], ["A", "B", "C"]):
        a.text(0.2, 1.7, b+num, c = "black", fontsize = 'large')
ax.flat[0].set_ylabel(r"$\frac{Active}{Inactive}^\ast$", fontsize = 14)
ax.flat[3].set_ylabel(r"$\frac{Immobile}{Mobile}^\ast$", fontsize = 14)

for g in ratio_plot:
    grow = ratio_plot.index(g)
    hr = hrsubratio
    for Reg in Regimes:
        if Reg == "Equal":
            r = "Medium"
        else:
            r = Reg 
        data_toplot = pd.DataFrame(columns = ["x","y"])
        for t in imposedtimeseries:
            n_reg = []
            for j in ["37"]:
                data = np.load(os.path.join(sourcedatadir, "SlowAR_" + t+"_NS-A"+j+"_df.npy"))
                velocity = np.mean(data[2, :, :, :], axis = (-1,-2))/basevelocity    
            for j in Trial:
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    n = hr.get("Norm_"+t + "/" + Reg + "/" + j + "/" + g)[:]
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
plt.savefig("Figure3_Temp_variation_biomass_ratios.jpg",
            dpi=300, bbox_inches="tight", pad=0.1)
hr.close()

#%%
#Figure 4 Amplitude - Heterogeneous case
mymarklist = ["o", "^", "s", "d"]
colorlist = ["steelblue", "orange","g", "indianred"]

data  = pd.read_csv("mass_flux_sensitivity_generalized_19022021.csv")
gvarnames = ["DO", "Nitrogen", "TOC"]
finaldata = data[data['Chem'].isin (gvarnames)]

finaldata["logDa"] = np.log10(finaldata.Da63)
finaldata.loc[finaldata["logDa"] < -1, "PeDamark"] = "0"
finaldata.loc[(finaldata["logDa"] > -1) & (finaldata["logDa"] < 0), "PeDamark"] = "1"
finaldata.loc[(finaldata["logDa"] > 0) & (finaldata["logDa"] <0.5), "PeDamark"] = "2"
finaldata.loc[(finaldata["logDa"] > 0.5), "PeDamark"] = "3"
labels = {"3" : "log$_{10}$Da > 0.5",
          "2" : "0 < log$_{10}$Da < 0.5",
          "1" : "-1 < log$_{10}$Da < 0",
         "0" : "log$_{10}$Da < -1"}

plt.figure(figsize = [7,7])
for frac in [0, 1, 2, 3]:
    sub = finaldata[finaldata["PeDamark"]==str(frac)]
    plt.scatter(sub["Time"], sub["Senssquared"], s = 80,
                c = colorlist[frac], marker = mymarklist[frac],
                label = labels[str(frac)],
                alpha = 0.5)
plt.ylabel("Normalised Responsiveness", fontsize = 18)
plt.xlabel("Residence time of solutes (days)", fontsize = 18)
plt.yticks(fontsize = 18)
plt.xticks(fontsize = 18)
plt.xscale("log")
plt.yscale("log")
plt.yticks((40,10,1,0.1), (40,10,1,0.1),fontsize = 18)
plt.xticks((100,10,1,0.1), (100,10,1,0.1),fontsize = 18)
plt.legend(fontsize = 14)
picname = "Figure4_logDa_Normalized_Amplitude_chem.jpg"
plt.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)

#%%
#Figure 5: Normalised Amplitude: Biomass
filename = "biomass_sensitivity_generalized_19022021.csv"

data = pd.read_csv(filename)
print(data.columns)
print(data.shape)
print(data.Chem.unique())
data["Regime"] = data["Regime"].replace(["Equal"], "Medium")

Chemseries = ["Immobile active aerobic degraders", "Immobile active ammonia oxidizers", "Immobile active nitrate reducers"]
Regimestoplot = ["Slow", "Medium", "Fast"]
species = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
markerseries = ["d", "^", "o"]
colorlist = ["indianred", "g", "steelblue"]

fig, axes = plt.subplots(nrows=3, ncols=1, figsize=[8, 8], sharex=True, sharey = True)
for a,b in zip(axes.flat[:], ["A", "B", "C"]):
    a.text(0.12, 1.0, b, c = "black", fontsize = 'xx-large')
for k in Chemseries:
    dfc = data[(data["Chem"] == k)]
    colidx1 = Chemseries.index(k)
    a = axes.flat[colidx1]
    for i in Regimestoplot:
        dfctemp = dfc
        dfcr = dfctemp[dfctemp["Regime"] == i]
        a.scatter("Time","Senssquared",s = 80,
                                   c=colorlist[Regimestoplot.index(i)],
                                   alpha = 0.7,
                                   marker = markerseries[Regimestoplot.index(i)],
                 data=dfcr,label = i + " flow"
        )
        a.tick_params(axis="y", labelsize=15)
        a.tick_params(axis="x", labelsize=15)
        if colidx1==2:
            a.legend(title = "Flow regime", title_fontsize = 12, fontsize = 12,
            bbox_to_anchor=(0.9, -0.35), ncol = 3)# loc = "best")
            a.set_xlabel("Residence time of solutes (days)",fontsize =18)
        else:
            le2 = a.legend([])
            le2.remove()
plt.xscale("log")
plt.xticks((0.1,1,10,100), (0.1,1,10,100))

for ax, typsp in zip(axes, species):
    ax.set_title(typsp, fontsize=15)
axes.flat[1].set_ylabel("Normalised Responsiveness", fontsize = 18)
plt.savefig("Figure5_biomass_normalized_amplitude.jpg",
    dpi=300,
    bbox_inches="tight",
    pad=0.1)
#%%
#Figure S1: Autocorrelation discussion for imposed temporal dynamics

head_path = "headatinlet.csv"
head = pd.read_csv(head_path, sep = ",")
cov1 = np.round(np.cov(head["H1"]),2)
cov2 = np.round(np.cov(head["H2"]),2)
cov3 = np.round(np.cov(head["H3"]),2)

def corrfunc (data):
    normdata = data - np.mean(data)
    autocorr_f = (np.correlate(normdata, normdata, mode='full'))/(np.std(data)**2 * np.shape(data)[0])
    return autocorr_f

cor1 = corrfunc(head["H1"])[5474:]
cor2 = corrfunc(head["H2"])[5474:]
cor3 = corrfunc(head["H3"])[5474:]

cor1tim = np.where(cor1 < 0.75)[0][0]
cor2tim = np.where(cor2 < 0.75)[0][0]
cor3tim = np.where(cor3 < 0.75)[0][0]
mincor = cor2.min()

plt.figure(figsize = (8,8))
plt.title("Autocorrelation of imposed temporal dynamics\nat the inlet of the domain", fontsize = 18)
plt.plot(cor1, label = 'T1: covariance= ' + str(cov1))
plt.plot(cor2, label = 'T2: covariance= ' + str(cov2))
plt.plot(cor3, label = 'T5: covariance= ' + str(cov3))
plt.vlines(cor1tim, mincor-0.1, cor1[cor1tim], colors = 'grey', linestyles = 'dashed')
plt.vlines(cor2tim, mincor-0.1, cor2[cor2tim], colors = 'grey', linestyles = 'dashed')
plt.vlines(cor3tim, mincor-0.1, cor3[cor3tim], colors = 'grey', linestyles = 'dashed')
plt.xscale("log")
plt.xlabel ("Time (days)", fontsize = 18)
plt.ylabel ("Normalized autocorrelation (-)", fontsize = 18)
plt.yticks (fontsize = 18)
plt.xticks(fontsize = 18)
plt.legend(fontsize = 14, title = "Time series")
plt.ylim(bottom = mincor - 0.1 )
picname = "FigureS1_Autoccorelation_head.png"
plt.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)

# Figure S2: Variation in concentration profile between max, average, and min in homogeneous domain
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

import analyses.transient as ta
import matplotlib.patches as mpatches
dashedline = mlines.Line2D([], [], linestyle = '--', color='grey', markersize=15, label='Maximum')
dotline = mlines.Line2D([], [], linestyle = 'dotted', color='grey', markersize=15, label='Minimum')
solidline = mlines.Line2D([], [], linestyle = 'solid', color='grey', markersize=15, label='Average')
blue_patch = mpatches.Patch(color="blue", label= 'Ammonium', alpha = 0.5)
black_patch = mpatches.Patch(color="black", label= 'DOC', alpha = 0.5)
green_patch = mpatches.Patch(color="darkgreen", label='Nitrate', alpha = 0.5)
red_patch = mpatches.Patch(color="red", label='DO', alpha = 0.5)
patchlist = [blue_patch, green_patch, black_patch, red_patch, solidline, dashedline, dotline]
Regimes = ["Slow", "Equal", "Fast"]
trialist = proc.masterscenarios("Saturated")
Trial = ["H"]
species = proc.speciesdict("Saturated")
gvarnames = ["DO", "DOC", "Ammonium", "Nitrate"]
cvars = list(species[g]['TecIndex'] for g in gvarnames)
velindex = 2
colors = ["red", "black", "blue", "darkgreen"]
columntitles = ["Slow flow", "Medium flow", "Fast flow"]
raw_directory = "E:/Zenodo_temporal_dynamics"
pad = 230
figbig, axes = plt.subplots(1,3, figsize=(10, 3), sharey = True, sharex = True)
for t in Trial:
    for r in Regimes:
        if r == "Equal":
            rtitle = "Medium"
        else:
            rtitle = r
        i = Trial.index(t)*len(Regimes) + Regimes.index(r)
        host = axes.flat[i]
        file = os.path.join(raw_directory, r + "AR_5_NS-A"+str(t)+"_df.npy")
        data = np.load(file)
        conctime, TotalFlow, Headinlettime = ta.conc_time (data,0,50,0,30, 51, gvarnames,"Saturated")
        yindex = list(range(51))
        #fig, host = axe.subplots()
        host.plot(np.max(conctime[1:, :, 0], axis = 0),yindex,label=gvarnames[0],color=colors[0],linestyle="--")
        host.plot(np.min(conctime[1:, :, 0], axis = 0),yindex,label=gvarnames[0],color=colors[0],linestyle="dotted")
        host.plot(np.mean(conctime[1:, :, 0], axis = 0),yindex,label=gvarnames[0],color=colors[0],linestyle="-")
        host.plot(np.max(conctime[1:, :, 1], axis = 0),yindex,label=gvarnames[1],color=colors[1],linestyle="--")
        host.plot(np.min(conctime[1:, :, 1], axis = 0),yindex,label=gvarnames[1],color=colors[1],linestyle="dotted")
        host.plot(np.mean(conctime[1:, :, 1], axis = 0),yindex,label=gvarnames[1],color=colors[1],linestyle="-")
        par1 = host.twiny()
        par2 = host.twiny()

        # Offset the top spine of par2.  The ticks and label have already been
        # placed on the top by twiny above.
        par2.spines["top"].set_position(("axes", 1.2))
        # Having been created by twinx, par2 has its frame off, so the line of its
        # detached spine is invisible.  First, activate the frame but make the patch
        # and spines invisible.
        make_patch_spines_invisible(par2)
        # Second, show the right spine.

        par1.plot(np.max(conctime[1:, :, 2], axis = 0),yindex,label=gvarnames[2],color=colors[2],linestyle="--")
        par1.plot(np.min(conctime[1:, :, 2], axis = 0),yindex,label=gvarnames[2],color=colors[2],linestyle="dotted")
        par1.plot(np.mean(conctime[1:, :, 2], axis = 0),yindex,label=gvarnames[2],color=colors[2],linestyle="-")
        par2.plot(np.max(conctime[1:, :, 3], axis = 0),yindex,label=gvarnames[3],color=colors[3],linestyle="--")
        par2.plot(np.min(conctime[1:, :, 3], axis = 0),yindex,label=gvarnames[3],color=colors[3],linestyle="dotted")
        par2.plot(np.mean(conctime[1:, :, 3], axis = 0),yindex,label=gvarnames[3],color=colors[3],linestyle="-")
        

        host.set_ylim(0, 51)
        host.set_xlim(-1, 800)
        par1.set_xlim(0, 65)
        par2.set_xlim(0, 260)
        host.xaxis.label.set_color("black")
        tkw = dict(size=4, width=1.5, labelsize=14)
        host.tick_params(axis="x", colors="black", **tkw)
        host.tick_params(axis="y", **tkw)
        if i < 3:
            host.set_title (rtitle + " flow", **titlekw)
            par2.spines["top"].set_visible(True)
            par1.xaxis.label.set_color("blue")
            par2.xaxis.label.set_color("darkgreen")
            par1.tick_params(axis="x", colors="blue", **tkw)
            par2.tick_params(axis="x", colors="darkgreen", **tkw)
            par1.set_xlabel(str(gvarnames[2]) + ' ($\mu$M)', **legendkw)
            par2.set_xlabel(str(gvarnames[3]) + ' ($\mu$M)', **legendkw)
        elif i > 5:
            host.set_xlabel("DOC, DO ($\mu$M)", **legendkw)
            par1.set_xticks([])
            par2.set_xticks([])
        else:
            par1.set_xticks([])
            par2.set_xticks([])
figbig.gca().invert_yaxis()
figbig.subplots_adjust(top=1.0, hspace = 0.2, wspace = 0.2)
for t,a in zip(Trial[::-1],range(3)):
    #plt.annotate("Variance: " + str(trialist[t]["Het"])+ " &\nAnisotropy: " + str(trialist[t]["Anis"]),
    #             xy=(0.1, 0.17), xytext=(-50, 0.7 + pad*a),
    #            xycoords='figure fraction', textcoords='offset points',
    #            rotation = "vertical",
    #            size='large', ha='center', va='baseline',
    #            fontsize = 16)
    axes.flat[3*a].set_ylabel("Y (cm)", **legendkw)
plt.legend(handles = patchlist, ncol = 3, **legendkw,
           bbox_to_anchor = (0.5,-0.6),
           loc = 'lower right')
plt.grid(False)
picname = "FigureS2_dissolved_species_1D_homogeneous_5.png"
plt.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)

# Figure S3: Variation in concentration profile between max, average, and min at
#each y coordinate for select scenarios in medium flow regime

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

import analyses.transient as ta
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
dashedline = mlines.Line2D([], [], linestyle = '--', color='grey', markersize=15, label='Maximum')
dotline = mlines.Line2D([], [], linestyle = 'dotted', color='grey', markersize=15, label='Minimum')
solidline = mlines.Line2D([], [], linestyle = 'solid', color='grey', markersize=15, label='Average')
blue_patch = mpatches.Patch(color="blue", label= 'Ammonium', alpha = 0.5)
black_patch = mpatches.Patch(color="black", label= 'DOC', alpha = 0.5)
green_patch = mpatches.Patch(color="darkgreen", label='Nitrate', alpha = 0.5)
red_patch = mpatches.Patch(color="red", label='DO', alpha = 0.5)
patchlist = [blue_patch, green_patch, black_patch, red_patch, solidline, dashedline, dotline]
Regimes = ["Slow", "Equal", "Fast"]
trialist = proc.masterscenarios("Saturated")
Trial = ["50", "73", "63"]
species = proc.speciesdict("Saturated")
gvarnames = ["DO", "DOC", "Ammonium", "Nitrate"]
cvars = list(species[g]['TecIndex'] for g in gvarnames)
velindex = 2
colors = ["red", "black", "blue", "darkgreen"]
columntitles = ["Slow flow", "Medium flow", "Fast flow"]
raw_directory = "E:/Zenodo_temporal_dynamics"
pad = 230
figbig, axes = plt.subplots(3,3, figsize=(13, 10), sharey = True, sharex = True)
for t in Trial:
    for r in Regimes:
        if r == "Equal":
            rtitle = "Medium"
        else:
            rtitle = r
        i = Trial.index(t)*len(Regimes) + Regimes.index(r)
        host = axes.flat[i]
        file = os.path.join(raw_directory, r + "AR_5_NS-A"+str(t)+"_df.npy")
        data = np.load(file)
        conctime, TotalFlow, Headinlettime = ta.conc_time (data,0,50,0,30, 51, gvarnames,"Saturated")
        yindex = list(range(51))
        #fig, host = axe.subplots()
        host.plot(np.max(conctime[1:, :, 0], axis = 0),yindex,label=gvarnames[0],color=colors[0],linestyle="--")
        host.plot(np.min(conctime[1:, :, 0], axis = 0),yindex,label=gvarnames[0],color=colors[0],linestyle="dotted")
        host.plot(np.mean(conctime[1:, :, 0], axis = 0),yindex,label=gvarnames[0],color=colors[0],linestyle="-")
        host.plot(np.max(conctime[1:, :, 1], axis = 0),yindex,label=gvarnames[1],color=colors[1],linestyle="--")
        host.plot(np.min(conctime[1:, :, 1], axis = 0),yindex,label=gvarnames[1],color=colors[1],linestyle="dotted")
        host.plot(np.mean(conctime[1:, :, 1], axis = 0),yindex,label=gvarnames[1],color=colors[1],linestyle="-")
        par1 = host.twiny()
        par2 = host.twiny()

        # Offset the top spine of par2.  The ticks and label have already been
        # placed on the top by twiny above.
        par2.spines["top"].set_position(("axes", 1.2))
        # Having been created by twinx, par2 has its frame off, so the line of its
        # detached spine is invisible.  First, activate the frame but make the patch
        # and spines invisible.
        make_patch_spines_invisible(par2)
        # Second, show the right spine.

        par1.plot(np.max(conctime[1:, :, 2], axis = 0),yindex,label=gvarnames[2],color=colors[2],linestyle="--")
        par1.plot(np.min(conctime[1:, :, 2], axis = 0),yindex,label=gvarnames[2],color=colors[2],linestyle="dotted")
        par1.plot(np.mean(conctime[1:, :, 2], axis = 0),yindex,label=gvarnames[2],color=colors[2],linestyle="-")
        par2.plot(np.max(conctime[1:, :, 3], axis = 0),yindex,label=gvarnames[3],color=colors[3],linestyle="--")
        par2.plot(np.min(conctime[1:, :, 3], axis = 0),yindex,label=gvarnames[3],color=colors[3],linestyle="dotted")
        par2.plot(np.mean(conctime[1:, :, 3], axis = 0),yindex,label=gvarnames[3],color=colors[3],linestyle="-")
        

        host.set_ylim(0, 51)
        host.set_xlim(-1, 800)
        par1.set_xlim(0, 65)
        par2.set_xlim(0, 260)
        host.xaxis.label.set_color("black")
        tkw = dict(size=4, width=1.5, labelsize=14)
        host.tick_params(axis="x", colors="black", **tkw)
        host.tick_params(axis="y", **tkw)
        if i < 3:
            host.set_title (rtitle + " flow", **titlekw)
            par2.spines["top"].set_visible(True)
            par1.xaxis.label.set_color("blue")
            par2.xaxis.label.set_color("darkgreen")
            par1.tick_params(axis="x", colors="blue", **tkw)
            par2.tick_params(axis="x", colors="darkgreen", **tkw)
            par1.set_xlabel(str(gvarnames[2]) + ' ($\mu$M)', **legendkw)
            par2.set_xlabel(str(gvarnames[3]) + ' ($\mu$M)', **legendkw)
        elif i > 5:
            host.set_xlabel("DOC, DO ($\mu$M)", **legendkw)
            par1.set_xticks([])
            par2.set_xticks([])
        else:
            par1.set_xticks([])
            par2.set_xticks([])
figbig.gca().invert_yaxis()
figbig.subplots_adjust(top=1.0, hspace = 0.2, wspace = 0.2)
for t,a in zip(Trial[::-1],range(3)):
    plt.annotate("Variance: " + str(trialist[t]["Het"])+ " &\nAnisotropy: " + str(trialist[t]["Anis"]),
                 xy=(0.1, 0.17), xytext=(-50, 0.7 + pad*a),
                xycoords='figure fraction', textcoords='offset points',
                rotation = "vertical",
                size='large', ha='center', va='baseline',
                fontsize = 16)
    axes.flat[3*a].set_ylabel("Y (cm)", **legendkw)
plt.legend(handles = patchlist, ncol = 4, **legendkw,
           bbox_to_anchor = (0.0,-0.6),
           loc = 'lower right')
plt.grid(False)
picname = "FigureS3_dissolved_species_1D_5.png"
plt.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)

#Figure S4: Variation in 1D profile biomass between max, average, and min in homogeneous domain
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

import analyses.transient as ta
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
dashedline = mlines.Line2D([], [], linestyle = '--', color='grey', markersize=15, label='Maximum')
dotline = mlines.Line2D([], [], linestyle = 'dotted', color='grey', markersize=15, label='Minimum')
solidline = mlines.Line2D([], [], linestyle = 'solid', color='grey', markersize=15, label='Average')

blue_patch = mpatches.Patch(color="blue", label= 'Ammonia oxidizers', alpha = 0.5)
black_patch = mpatches.Patch(color="black", label= 'Aerobes', alpha = 0.5)
green_patch = mpatches.Patch(color="darkgreen", label='Nitrate reducers', alpha = 0.5)
patchlist = [blue_patch, dashedline, black_patch, solidline, green_patch, dotline]
Regimes = ["Slow", "Equal", "Fast"]
trialist = proc.masterscenarios("Saturated")
Trial = ["H"]
species = proc.speciesdict("Saturated")
gvarnames = list(g for g in species.keys() if (species[g]["State"] == "Active") and (species[g]["Location"] == "Immobile"))
gvarnames.remove('Immobile active sulphate reducers')
cvars = list(species[g]['TecIndex'] for g in gvarnames)
velindex = 2
colors = ["black", "darkgreen","blue"]
columntitles = ["Slow flow", "Medium flow", "Fast flow"]
pad = 230

raw_directory = "E:/Zenodo_temporal_dynamics"

figbig, axes = plt.subplots(1,3, figsize=(10, 3), sharey = True, sharex = True)
for t in Trial:
    for r in Regimes:
        if r == "Equal":
            rtitle = "Medium"
        else:
            rtitle = r
        i = Trial.index(t)*len(Regimes) + Regimes.index(r)
        host = axes.flat[i]
        file = os.path.join(raw_directory, r + "AR_5_NS-A"+str(t)+"_df.npy")
        data = np.load(file)
        masstime, conctime = ta.biomasstimefunc(data,0,50,0,30, 51, gvarnames,"Saturated")
        yindex = list(range(51))
        #fig, host = axe.subplots()
        host.plot(np.max(conctime[1:, :, 0], axis = 0),yindex,label=gvarnames[0],color=colors[0],linestyle="--")
        host.plot(np.min(conctime[1:, :, 0], axis = 0),yindex,label=gvarnames[0],color=colors[0],linestyle="dotted")
        host.plot(np.mean(conctime[1:, :, 0], axis = 0),yindex,label=gvarnames[0],color=colors[0],linestyle="-")
        par1 = host.twiny()
        par2 = host.twiny()

        # Offset the top spine of par2.  The ticks and label have already been
        # placed on the top by twiny above.
        par2.spines["top"].set_position(("axes", 1.2))
        # Having been created by twinx, par2 has its frame off, so the line of its
        # detached spine is invisible.  First, activate the frame but make the patch
        # and spines invisible.
        make_patch_spines_invisible(par2)
        # Second, show the right spine.
        par1.plot(np.max(conctime[1:, :, 1], axis = 0),yindex,label=gvarnames[1],color=colors[1],linestyle="--")
        par1.plot(np.min(conctime[1:, :, 1], axis = 0),yindex,label=gvarnames[1],color=colors[1],linestyle="dotted")
        par1.plot(np.mean(conctime[1:, :, 1], axis = 0),yindex,label=gvarnames[1],color=colors[1],linestyle="-")
        par2.plot(np.max(conctime[1:, :, 2], axis = 0),yindex,label=gvarnames[2],color=colors[2],linestyle="--")
        par2.plot(np.min(conctime[1:, :, 2], axis = 0),yindex,label=gvarnames[2],color=colors[2],linestyle="dotted")
        par2.plot(np.mean(conctime[1:, :, 2], axis = 0),yindex,label=gvarnames[2],color=colors[2],linestyle="-")
        host.set_ylim(0, 51)
        host.set_xlim(0, 500)
        par1.set_xlim(0, 200)
        par2.set_xlim(0, 30)
        host.xaxis.label.set_color("black")
        tkw = dict(size=4, width=1.5, labelsize=14)
        host.tick_params(axis="x", colors="black", **tkw)
        host.tick_params(axis="y", **tkw)
        if i < 3:
            host.set_title (rtitle + " flow", **titlekw)
            par2.spines["top"].set_visible(True)
            par1.xaxis.label.set_color("darkgreen")
            par2.xaxis.label.set_color("blue")
            par1.tick_params(axis="x", colors="darkgreen", **tkw)
            par2.tick_params(axis="x", colors="blue", **tkw)
            par1.set_xlabel(species[gvarnames[1]]["Graphname"] + " ($\mu$M)", **legendkw)
            par2.set_xlabel(species[gvarnames[2]]["Graphname"] + " ($\mu$M)", **legendkw)
            host.set_xlabel(species[gvarnames[0]]["Graphname"] + " ($\mu$M)", **legendkw)
        elif i > 5:
            host.set_xlabel(species[gvarnames[0]]["Graphname"], **legendkw)
            par1.set_xticks([])
            par2.set_xticks([])
        else:
            par1.set_xticks([])
            par2.set_xticks([])
figbig.gca().invert_yaxis()
figbig.subplots_adjust(top=1.0, hspace = 0.2, wspace = 0.2)
for t,a in zip(Trial[::-1],range(3)):
#    plt.annotate("Variance: " + str(trialist[t]["Het"])+ " &\nAnisotropy: " + str(trialist[t]["Anis"]),
#                 xy=(0.1, 0.17), xytext=(-50, 0.7 + pad*a),xycoords='figure fraction',textcoords='offset points',rotation = "vertical",size='large', ha='center', va='baseline',fontsize = 16)
    axes.flat[3*a].set_ylabel("Y (cm)", **titlekw)
plt.legend(handles = patchlist, ncol = 3, bbox_to_anchor = (-0.6,-0.6),loc = 'lower center', **legendkw)
plt.grid(False)
picname = "FigureS4_biomass_1D.png"
plt.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)

#Figure S5: Variation in 1D profile biomass between max, average, and min in heterogeneous domains
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

import analyses.transient as ta
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
dashedline = mlines.Line2D([], [], linestyle = '--', color='grey', markersize=15, label='Maximum')
dotline = mlines.Line2D([], [], linestyle = 'dotted', color='grey', markersize=15, label='Minimum')
solidline = mlines.Line2D([], [], linestyle = 'solid', color='grey', markersize=15, label='Average')

blue_patch = mpatches.Patch(color="blue", label= 'Ammonia oxidizers', alpha = 0.5)
black_patch = mpatches.Patch(color="black", label= 'Aerobes', alpha = 0.5)
green_patch = mpatches.Patch(color="darkgreen", label='Nitrate reducers', alpha = 0.5)
patchlist = [blue_patch, dashedline, black_patch, solidline, green_patch, dotline]
Regimes = ["Slow", "Equal", "Fast"]
trialist = proc.masterscenarios("Saturated")
Trial = ["50", "73", "63"]
species = proc.speciesdict("Saturated")
gvarnames = list(g for g in species.keys() if (species[g]["State"] == "Active") and (species[g]["Location"] == "Immobile"))
gvarnames.remove('Immobile active sulphate reducers')
cvars = list(species[g]['TecIndex'] for g in gvarnames)
velindex = 2
colors = ["black", "darkgreen","blue"]
columntitles = ["Slow flow", "Medium flow", "Fast flow"]
pad = 230

raw_directory = "E:/Zenodo_temporal_dynamics"

figbig, axes = plt.subplots(3,3, figsize=(13, 10), sharey = True, sharex = True)
for t in Trial:
    for r in Regimes:
        if r == "Equal":
            rtitle = "Medium"
        else:
            rtitle = r
        i = Trial.index(t)*len(Regimes) + Regimes.index(r)
        host = axes.flat[i]
        file = os.path.join(raw_directory, r + "AR_5_NS-A"+str(t)+"_df.npy")
        data = np.load(file)
        masstime, conctime = ta.biomasstimefunc(data,0,50,0,30, 51, gvarnames,"Saturated")
        yindex = list(range(51))
        #fig, host = axe.subplots()
        host.plot(np.max(conctime[1:, :, 0], axis = 0),yindex,label=gvarnames[0],color=colors[0],linestyle="--")
        host.plot(np.min(conctime[1:, :, 0], axis = 0),yindex,label=gvarnames[0],color=colors[0],linestyle="dotted")
        host.plot(np.mean(conctime[1:, :, 0], axis = 0),yindex,label=gvarnames[0],color=colors[0],linestyle="-")
        par1 = host.twiny()
        par2 = host.twiny()

        # Offset the top spine of par2.  The ticks and label have already been
        # placed on the top by twiny above.
        par2.spines["top"].set_position(("axes", 1.2))
        # Having been created by twinx, par2 has its frame off, so the line of its
        # detached spine is invisible.  First, activate the frame but make the patch
        # and spines invisible.
        make_patch_spines_invisible(par2)
        # Second, show the right spine.
        par1.plot(np.max(conctime[1:, :, 1], axis = 0),yindex,label=gvarnames[1],color=colors[1],linestyle="--")
        par1.plot(np.min(conctime[1:, :, 1], axis = 0),yindex,label=gvarnames[1],color=colors[1],linestyle="dotted")
        par1.plot(np.mean(conctime[1:, :, 1], axis = 0),yindex,label=gvarnames[1],color=colors[1],linestyle="-")
        par2.plot(np.max(conctime[1:, :, 2], axis = 0),yindex,label=gvarnames[2],color=colors[2],linestyle="--")
        par2.plot(np.min(conctime[1:, :, 2], axis = 0),yindex,label=gvarnames[2],color=colors[2],linestyle="dotted")
        par2.plot(np.mean(conctime[1:, :, 2], axis = 0),yindex,label=gvarnames[2],color=colors[2],linestyle="-")
        host.set_ylim(0, 51)
        host.set_xlim(0, 500)
        par1.set_xlim(0, 200)
        par2.set_xlim(0, 30)
        host.xaxis.label.set_color("black")
        tkw = dict(size=4, width=1.5, labelsize=14)
        host.tick_params(axis="x", colors="black", **tkw)
        host.tick_params(axis="y", **tkw)
        if i < 3:
            host.set_title (rtitle + " flow", **titlekw)
            par2.spines["top"].set_visible(True)
            par1.xaxis.label.set_color("darkgreen")
            par2.xaxis.label.set_color("blue")
            par1.tick_params(axis="x", colors="darkgreen", **tkw)
            par2.tick_params(axis="x", colors="blue", **tkw)
            par1.set_xlabel(species[gvarnames[1]]["Graphname"] + " ($\mu$M)", **legendkw)
            par2.set_xlabel(species[gvarnames[2]]["Graphname"] + " ($\mu$M)", **legendkw)
        elif i > 5:
            host.set_xlabel(species[gvarnames[0]]["Graphname"]+ " ($\mu$M)", **legendkw)
            par1.set_xticks([])
            par2.set_xticks([])
        else:
            par1.set_xticks([])
            par2.set_xticks([])
figbig.gca().invert_yaxis()
figbig.subplots_adjust(top=1.0, hspace = 0.2, wspace = 0.2)
for t,a in zip(Trial[::-1],range(3)):
    plt.annotate("Variance: " + str(trialist[t]["Het"])+ " &\nAnisotropy: " + str(trialist[t]["Anis"]),
                 xy=(0.1, 0.17), xytext=(-50, 0.7 + pad*a),xycoords='figure fraction',textcoords='offset points',rotation = "vertical",size='large', ha='center', va='baseline',fontsize = 16)
    axes.flat[3*a].set_ylabel("Y (cm)", **titlekw)
plt.legend(handles = patchlist, ncol = 3, bbox_to_anchor = (-0.6,-0.6),loc = 'lower center', **legendkw)
plt.grid(False)
picname = "FigureS5_biomass_1D.png"
plt.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)

#Figure S6 Amplitude - Homogeneous case - chemical species
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
grey_dot = mlines.Line2D([], [], linestyle = '', marker = "o", markerfacecolor = "grey", markeredgecolor = "grey", markersize=10, label='DO', alpha = 0.5)
grey_triangle = mlines.Line2D([], [], linestyle = '', marker = "^", markerfacecolor = "grey", markeredgecolor = "grey",markersize=10, label='Nitrogen', alpha = 0.5)
grey_square = mlines.Line2D([], [], linestyle = '', marker = "s", markerfacecolor = "grey", markeredgecolor = "grey",markersize=10, label='DOC', alpha = 0.5)
my_pal = {3:"indianred", 2: "g", 0:"steelblue", 1 :"orange"}
blue_patch = mpatches.Patch(color="steelblue", label= "Fast flow", alpha = 0.5)
green_patch = mpatches.Patch(color="g", label="Medium flow", alpha = 0.5)
red_patch = mpatches.Patch(color="indianred", label="Slow flow", alpha = 0.5)
patchlist = [blue_patch, green_patch, red_patch, grey_square, grey_dot, grey_triangle]
mymarklist = ["o", "^", "s", "d"]
reglist = ["Slow", "Medium", "Fast"]
colorlist = ["indianred", "g", "steelblue"]

data  = pd.read_csv("mass_flux_sensitivity_generalized_19022021.csv")
gvarnames = ["DO", "Nitrogen", "DOC"]
finaldata = data[data['Chem'].isin (gvarnames)]

finaldata["xaxis"] = finaldata["Time"]

base= finaldata[finaldata["Trial"]=='H']

plt.figure(figsize = [7,7])
for r in reglist:
    for g in gvarnames:
        sub = base[(base["Regime"]==r) & (base["Chem"]==g)]
        plt.scatter(sub["xaxis"], sub["Sensitivity%"], s = 80, c = colorlist[reglist.index(r)], marker = mymarklist[gvarnames.index(g)], alpha = 0.5, label = g + " in " + r + " flow regime")
plt.ylabel("Responsiveness (%)", fontsize = 18)
plt.xlabel("Residence time (days)", fontsize = 18)
plt.xscale("log")
plt.yscale("log")
plt.yticks((100,10,1,0.1), (100,10,1,0.1),fontsize = 18)
plt.xticks((100,10,1,0.1), (100,10,1,0.1),fontsize = 18)
plt.legend(handles = patchlist, fontsize = 14)
picname = "FigureS6_amplitude_basecase_homogeneous.png"
plt.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)

#Figure S7 Amplitude - Homogeneous case - biomass

filename = "Normalized_RMSamplitude_biomass.csv"
import matplotlib.patches as mpatches
import matplotlib as mpl
Redscmap = mpl.cm.Reds(np.linspace(0, 1, 30))
Greenscmap = mpl.cm.Greens(np.linspace(0, 1, 30))
Bluescmap = mpl.cm.Blues(np.linspace(0, 1, 30))
Redscmap = mpl.colors.ListedColormap(Redscmap[10:, :-1])
Greenscmap = mpl.colors.ListedColormap(Greenscmap[10:, :-1])
Bluescmap = mpl.colors.ListedColormap(Bluescmap[10:, :-1])
#colseries = [Redscmap, Greenscmap, Bluescmap]
colseries = ["indianred", "g", "steelblue"]

data = pd.read_csv(filename, sep="\t")
print(data.columns)
print(data.shape)
print(data.Chem.unique())
data["Regime"] = data["Regime"].replace(["Equal"], "Medium")

master = data

data = master[master["Trial"]=='H']
Chemseries = ["Immobile active aerobic degraders", "Immobile active ammonia oxidizers", "Immobile active nitrate reducers"]
Regimestoplot = ["Slow", "Medium", "Fast"]
species = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
markerseries = ["d", "^", "o"]
data["xaxis"] = data["Time"]

fig, axes = plt.subplots(nrows=3, ncols=1, figsize=[8, 8], sharex=True, sharey = True)
for k in Chemseries:
    dfc = data[(data["Chem"] == k)]
    colidx1 = Chemseries.index(k)
    for i in Regimestoplot:
        dfctemp = dfc
        dfcr = dfctemp[dfctemp["Regime"] == i]
        axes.flat[colidx1].scatter("xaxis","Sensitivity%",s = 80, c=colseries[Regimestoplot.index(i)], alpha = 0.7,
                 data=dfcr,label = i + " flow"
        )
        axes.flat[colidx1].tick_params(axis="y", labelsize=18)
        axes.flat[colidx1].tick_params(axis="x", labelsize=18)
plt.xscale("log")
plt.legend(fontsize = 14)
for ax, typsp in zip(axes, species):
    ax.set_title(typsp, fontsize=18)
    ax.set_xticks((1,10,100),(1,10,100))
axes.flat[1].set_ylabel("Responsiveness (%)", fontsize = 18)
plt.annotate(
        "Residence time (days)",
        xy=(0.5, 0),
        xytext=(0, -50),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=18,
    )
plt.savefig("FigureS7_biomass_amplitude_homogeneous.png",
    dpi=300,
    bbox_inches="tight",
    pad=0.1)

#Figure S8 Peak cross-correlation - chemical species

directory = r"Y:\Home\khurana\4. Publications\Restructuring\Paper2\Figurecodes\/"
filename = "crosschor_memory_chem.csv"

data = pd.read_csv(filename, sep="\t")
print(data.columns)
print(data.shape)
print(data.Chem.unique())
data["Regime"] = data["Regime"].replace(["Equal"], "Medium")
chemstoplot = ["DOC", "DO", "TOC", "Nitrogen"]
dummy = stp.correlationdistributionchem(data[(data["Trial"] != "52") & (data["Trial"] != "43")], chemstoplot)
dummy.savefig("FigureS8_cross_correlation_distribution_chem.png",
    dpi=300,
    bbox_inches="tight",
    pad=0.1,
)

#Figure S9 Peak cross-correlation -biomass

filename = "crosschor_memory_biomass.csv"

data = pd.read_csv(filename, sep="\t")
print(data.columns)
print(data.shape)
print(data.Chem.unique())
data["Regime"] = data["Regime"].replace(["Equal"], "Medium")

biomasstoplot = ["Immobile active aerobic degraders", "Immobile active ammonia oxidizers", "Immobile active nitrate reducers",
                 "Mobile active aerobic degraders", "Mobile active ammonia oxidizers", "Mobile active nitrate reducers"]

dummy = stp.correlationdistributionbiomass(data[(data["Trial"] != "52") & (data["Trial"] != "43")], biomasstoplot)
dummy.savefig("FigureS9_cross_correlation_distribution_biomass.png",
    dpi=300,
    bbox_inches="tight",
    pad=0.1,
)

#Figure S10 System memory: Chemical species

low_var = mpatches.Patch(color="silver", label="Low")
mid_var = mpatches.Patch(color="darkgray", label="Medium")
high_var = mpatches.Patch(color="black", label="High")
my_pal = {3:"indianred", 2: "g", 0:"steelblue", 1 :"orange"}
blue_patch = mpatches.Patch(color="steelblue", label= "Fast flow", alpha = 0.5)
green_patch = mpatches.Patch(color="g", label="Medium flow", alpha = 0.5)
red_patch = mpatches.Patch(color="indianred", label="Slow flow", alpha = 0.5)
patchlist = [blue_patch, green_patch, red_patch, low_var, mid_var, high_var]
mymarklist = ["o", "^", "s", "d"]
reglist = ["Slow", "Medium", "Fast"]
colorlist = ["indianred", "g", "steelblue"]
Redscmap = mpl.cm.Reds(np.linspace(0, 1, 30))
Greenscmap = mpl.cm.Greens(np.linspace(0, 1, 30))
Bluescmap = mpl.cm.Blues(np.linspace(0, 1, 30))
Redscmap = mpl.colors.ListedColormap(Redscmap[10:, :-1])
Greenscmap = mpl.colors.ListedColormap(Greenscmap[10:, :-1])
Bluescmap = mpl.colors.ListedColormap(Bluescmap[10:, :-1])
colseries = [Redscmap, Greenscmap, Bluescmap]

filename = "crosschor_memory_chem.csv"

data = pd.read_csv(filename, sep="\t")
print(data.columns)
print(data.shape)
print(data.Chem.unique())

head_path = "headatinlet.csv"
head = pd.read_csv(head_path, sep= ",")
cor1tim = np.where(stattools.acf(head["H1"], fft = True, nlags = 5476) < 0.75)[0][0]
cor2tim = np.where(stattools.acf(head["H2"], fft = True, nlags = 5476) < 0.75)[0][0]
cor3tim = np.where(stattools.acf(head["H3"], fft = True, nlags = 5476) < 0.75)[0][0]

cov1 = np.round(stattools.acovf(head["H1"], fft = True),2)[0]
cov2 = np.round(stattools.acovf(head["H2"], fft = True),2)[0]
cov3 = np.round(stattools.acovf(head["H3"], fft = True),2)[0]

data["timtrace"]=data["Time_series"]
data["timtrace"]=data["timtrace"].replace([1], cor1tim)
data["timtrace"]=data["timtrace"].replace([2], cor2tim)
data["timtrace"]=data["timtrace"].replace([5], cor3tim)
data['fraction%'] = data['fraction']*100

data["normmem"] = data["Memory"]/data["timtrace"]
    
data["Regime"] = data["Regime"].replace(["Equal"], "Medium")
gvarnames = ["DOC", "DO", "TOC", "Nitrogen"]

indices = data.Regime.unique().tolist()
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=[11, 8], sharex=True, sharey = True)
for k in gvarnames:
    dfc = data[(data["Chem"] == k)]  # &(data['Time_series']==3)]
    colidx1 = gvarnames.index(k)
    for i in indices:
        dfctemp = dfc
        dfcr = dfctemp[dfctemp["Regime"] == i]
        axes.flat[colidx1].scatter(
            "fraction%",
            "normmem",
            c="Time_series",
            cmap=colseries[indices.index(i)],
            data=dfcr,
            label="Minimum",
            marker=mymarklist[indices.index(i)],
        )
        axes.flat[colidx1].tick_params(axis="y", labelsize=15)
        axes.flat[colidx1].tick_params(axis="x", labelsize=15)
        axes.flat[colidx1].set_title(k, fontsize=15)
plt.annotate("Ratio of backward traceability and\nmemory of forcing",
        xy=(-1.1, 1.1),
        xytext=(-100, 5),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=18)
plt.annotate("Residence time of solutes (%)",#"Temporal dynamics factor",
        xy=(-0.1, 0),
        xytext=(0, -40),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=18,
)
legend_flow = plt.legend(handles=patchlist[:3], ncol = 2,
        bbox_to_anchor=(-0.6, -0.5),
        loc="center",
        title="Flow regime",
        fontsize=14, title_fontsize = 14)
plt.legend(handles=patchlist[-3:], ncol = 2,
        bbox_to_anchor=(0.3, -0.5),
        loc="center",
        title="Memory of forcing",
        fontsize=14, title_fontsize = 14)
plt.gca().add_artist(legend_flow)
picname = "FigureS10_Memory_chem.png"
plt.savefig(picname, dpi=300, bbox_inches="tight", pad=0.1)

#Figure S11 System memory: Microbial species

filename = "crosschor_memory_biomass.csv"

data = pd.read_csv(filename, sep="\t")
print(data.columns)
print(data.shape)
print(data.Chem.unique())

data["timtrace"]=data["Time_series"]
data["timtrace"]=data["timtrace"].replace([1], cor1tim)
data["timtrace"]=data["timtrace"].replace([2], cor2tim)
data["timtrace"]=data["timtrace"].replace([5], cor3tim)

data["normmem"] = data["Memory"]/data['timtrace']#data["Time"]
data['fraction%'] = data['fraction']*100
data["Regime"] = data["Regime"].replace(["Equal"], "Medium")

biomasstoplot = ["Immobile active aerobic degraders", "Immobile active ammonia oxidizers", "Immobile active nitrate reducers",
                 "Mobile active aerobic degraders", "Mobile active ammonia oxidizers", "Mobile active nitrate reducers"]
Regimestoplot = data.Regime.unique().tolist()

species = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
position = ["Immobile", "Mobile"]
markerseries = ["d", "^", "o"]
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=[10, 6], sharex=True, sharey = True)
for k in biomasstoplot:
    dfc = data[(data["Chem"] == k)]
    colidx1 = biomasstoplot.index(k)
    for i in Regimestoplot:
        dfctemp = dfc
        dfcr = dfctemp[dfctemp["Regime"] == i]
        axes.flat[colidx1].scatter(
            "fraction%",
            "normmem",#Memory",#"normmem",
            c="Time_series",
            cmap=colseries[Regimestoplot.index(i)],
            data=dfcr,
            marker=markerseries[Regimestoplot.index(i)])
        axes.flat[colidx1].tick_params(axis="y", labelsize=15)
        axes.flat[colidx1].tick_params(axis="x", labelsize=15)
#plt.yscale("log")
plt.xticks((100,80,60,40,20),(100,80,60,40,20),fontsize = 15)
#plt.yticks((1,10,100,1000),(1,10,100,1000),fontsize = 15)
for ax, typsp in zip(axes[0, :], species):
    ax.set_title(typsp, fontsize=15)
axes[0, -1].annotate(position[0],
        xy=(0, 0.5),
        xytext=(180, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15)
axes[1, -1].annotate(position[1],
        xy=(0, 0.5),
        xytext=(180, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15)
plt.annotate("Ratio of backward traceability and\nmemory of the forcing",
        xy=(-0.1, 2.0),
        xytext=(-450, -150),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        rotation="vertical",
        fontsize=18)
plt.annotate("Residence time of solutes (%)",
        xy=(-0.4, -0.3),
        xytext=(-50, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=18)
legend_flow = plt.legend(handles=patchlist[:3], ncol = 2,
        bbox_to_anchor=(-1.5, -0.6),
        loc="center",
        title="Flow regime",
        fontsize=14, title_fontsize = 14)
plt.legend(handles=patchlist[-3:], ncol = 2,
        bbox_to_anchor=(0.3, -0.6),
        loc="center",
        title="Memory of forcing",
        fontsize=14, title_fontsize = 14)
plt.gca().add_artist(legend_flow)
plt.savefig("FigureS11_Memory_biomass.png",
    dpi=300,
    bbox_inches="tight",
    pad=0.1)

#Figure S12 Normalized amplitude - Heterogeneous case
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
grey_dot = mlines.Line2D([], [], linestyle = '', marker = "o", markerfacecolor = "grey", markeredgecolor = "grey", markersize=10, label='DO', alpha = 0.5)
grey_triangle = mlines.Line2D([], [], linestyle = '', marker = "^", markerfacecolor = "grey", markeredgecolor = "grey",markersize=10, label='Nitrogen', alpha = 0.5)
grey_square = mlines.Line2D([], [], linestyle = '', marker = "s", markerfacecolor = "grey", markeredgecolor = "grey",markersize=10, label='TOC', alpha = 0.5)
my_pal = {3:"indianred", 2: "g", 0:"steelblue", 1 :"orange"}
blue_patch = mpatches.Patch(color="steelblue", label= "Fast flow", alpha = 0.5)
green_patch = mpatches.Patch(color="g", label="Medium flow", alpha = 0.5)
red_patch = mpatches.Patch(color="indianred", label="Slow flow", alpha = 0.5)
patchlist = [blue_patch, green_patch, red_patch, grey_square, grey_dot, grey_triangle]
mymarklist = ["o", "^", "s", "d"]
reglist = ["Slow", "Medium", "Fast"]
colorlist = ["indianred", "g", "steelblue"]

data  = pd.read_csv("mass_flux_sensitivity_generalized_19022021.csv")
gvarnames = ["DO", "Nitrogen", "TOC"]
finaldata = data[data['Chem'].isin (gvarnames)]

plt.figure(figsize = [10,7])
for r in reglist:
    for g in gvarnames:
        sub = finaldata[(finaldata["Regime"]==r) & (finaldata["Chem"]==g)]
        plt.scatter(sub["Time"], sub["Senssquared"], s = 80, c = colorlist[reglist.index(r)], marker = mymarklist[gvarnames.index(g)], alpha = 0.5, label = g + " in " + r + " flow regime")
plt.ylabel("Normalised responsiveness", fontsize = 18)
plt.xlabel("Residence time of solutes (days)", fontsize = 18)
plt.xscale("log")
plt.yscale("log")
plt.yticks(fontsize = 18)
plt.xticks((0.1,1,10,100), (0.1,1,10,100),fontsize = 18)
plt.legend(handles = patchlist, fontsize = 14)
picname = "FigureS12_normalised_sensitivity_chem.png"
plt.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)

#Figure S13 Normalised amplitude: Microbial species
low_var = mpatches.Patch(color="silver", label="Low")
mid_var = mpatches.Patch(color="darkgray", label="Medium")
high_var = mpatches.Patch(color="black", label="High")
blue_patch = mpatches.Patch(color="steelblue", label= "Fast flow", alpha = 0.5)
green_patch = mpatches.Patch(color="g", label="Medium flow", alpha = 0.5)
red_patch = mpatches.Patch(color="indianred", label="Slow flow", alpha = 0.5)
patchlist = [blue_patch, green_patch, red_patch, low_var, mid_var, high_var]
filename = "Normalized_RMSamplitude_biomass_190022021.csv"

data = pd.read_csv(filename)
print(data.columns)
print(data.shape)
print(data.Chem.unique())
data["Regime"] = data["Regime"].replace(["Equal"], "Medium")

biomasstoplot = ["Immobile active aerobic degraders", "Immobile active ammonia oxidizers", "Immobile active nitrate reducers",
                 "Mobile active aerobic degraders", "Mobile active ammonia oxidizers", "Mobile active nitrate reducers"]

amplitude = stp.RMSamp_biomass(
    data[data["Trial"] != "52"],
    biomasstoplot,
    "Sensitivity%",
    "Normalized responsiveness (%)",
    ["Slow", "Medium", "Fast"],
)
legend_flow = plt.legend(handles=patchlist[:3], ncol = 2,
        bbox_to_anchor=(-1.5, -0.6),
        loc="center",
        title="Flow regime",
        fontsize=14, title_fontsize=14)
plt.legend(handles=patchlist[-3:], ncol = 2,
        bbox_to_anchor=(0.3, -0.6),
        loc="center",
        title="Memory of forcing",
        fontsize=14, title_fontsize=14)
plt.gca().add_artist(legend_flow)
amplitude.savefig("FigureS13_Normalized_Maximum%_RMS_amp_biomass.png",
    dpi=300,
    bbox_inches="tight",
    pad=0.1,
)

##Figure S14: Microbial subpopulation variation with changing velocity

basedata = np.load(os.path.join(sourcedatadir, "SlowAR_0_NS-AH_df.npy"))
basevelocity = np.mean(basedata[2, -1, :, :])
vardict = proc.speciesdict("Saturated")
position = ["Immobile", "Mobile"]
state = ["Active", "Inactive"]
gvarnames = list(t for t in vardict.keys() if vardict[t]["State"] in (state))
nosulfate = list(g for g in gvarnames if "sulphate" not in g)
species = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]

blue_patch = mpatches.Patch(color="steelblue", label= 'Fast flow', alpha = 0.5)
green_patch = mpatches.Patch(color="g", label='Medium flow', alpha = 0.5)
red_patch = mpatches.Patch(color="indianred", label='Slow flow', alpha = 0.5)
patchlist = [blue_patch, green_patch, red_patch]
regcolors = {"Slow":"indianred", "Equal":"g", "Fast":"steelblue"}

hr = h5py.File("Temporal_analysis_biomass.h5", mode = 'r')
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
        for j in ["37"]:
            data = np.load(os.path.join(sourcedatadir, "SlowAR_" + t+"_NS-A"+j+"_df.npy"))
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
plt.legend(handles=patchlist, ncol = 3,
        bbox_to_anchor=(-0.6, -0.55),
        loc="center",
        title="Flow regime",
        fontsize=14, title_fontsize = 14)
plt.savefig("FigureS14_Temp_variation_biomass_sampled.jpg",
            dpi=300,
            bbox_inches="tight",
            pad=0.1)