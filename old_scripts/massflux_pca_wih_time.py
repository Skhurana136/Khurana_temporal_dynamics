# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 17:21:42 2020

@author: khurana
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns

from sklearn.decomposition import PCA

from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()

raw_dir = "E:/Zenodo_temporal_dynamics"
out_dir = "Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes"

Regimes = ["Slow", "Equal", "Fast"]

import data_reader.data_processing as proc
import analyses.transient as ta

scdict = proc.masterscenarios("Saturated") #master dictionary of all spatially heterogeneous scenarios that were run

# Scenarios to investigate:
Trial = list(scdict.keys())

vardict = proc.speciesdict("Saturated")
States = ["Active", "Inactive"]
gvarnames = ["DOC", "DO", "TOC", "Ammonium","Nitrogen"]

timrange = list(range(1095))
#Distribution of biomass ratio
rtji = 0
for Reg in Regimes:
    reglist = [Reg]*1095
    basedata = np.load(os.path.join(raw_dir, Reg + "AR_0_NS-AH_df.npy"))
    basevelocity = np.mean(basedata[2, -1, :, :])
    for t in ["1", "2", "5"]:
        print (Reg, t)
        serieslist = [t]*1095
        for j in Trial:
            jlist = [j]*1095
            vlist = [scdict[j]["Het"]]*1095
            alist = [scdict[j]["Anis"]]*1095
            if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                pass
            else:
                data_path = os.path.join(raw_dir, Reg+"AR_"+t+"_NS-A"+j+"_df.npy")
                data = np.load(data_path)
                conc,v,h = ta.conc_time(data, 0, -1, 0, -1, 51, gvarnames, "Saturated")
                conc_df = pd.DataFrame(conc[1:,-1,:], columns = gvarnames)
                velocity = np.mean(data[2, 1:, :, :], axis = (-1,-2))/basevelocity
                g_df = pd.DataFrame(np.column_stack([reglist, serieslist, jlist, vlist, alist, velocity]), 
                               columns=['Regime','Time_series', 'Trial', 'Variance', 'Anisotropy', "Normalised_vel"])
                g_df = pd.concat([g_df,conc_df], axis = 1)
            if rtji == 0:
                results_df = g_df
            else:
                results_df = results_df.append(g_df,ignore_index=True)
            rtji += 1

head_path = os.path.join(out_dir, "headatinlet.csv")
head = pd.read_csv(head_path, sep= ",")
cor1tim = 39
cor2tim = 224
cor3tim = 250
results_df["timtrace"]=results_df["Time_series"]
results_df["timtrace"]=results_df["timtrace"].replace([1], cor1tim)
results_df["timtrace"]=results_df["timtrace"].replace([2], cor2tim)
results_df["timtrace"]=results_df["timtrace"].replace([5], cor3tim)
results_df["uid_v"] = results_df["Regime"]+"_"+results_df["Time_series"]+"_"+results_df["Variance"]
results_df["uid"] = results_df["Regime"]+"_"+results_df["Time_series"]+"_"+results_df["Trial"]

filename = os.path.join(raw_dir, "pca_chem.csv")
results_df.to_csv(filename, index=False)

results_df = pd.read_csv(filename)
X = results_df[gvarnames]
y = results_df["uid"]

X1 = scaler.fit_transform(X)
pca = PCA(n_components=2)
pca.fit(X1)
X2 = pca.transform(X1)

pca.explained_variance_ratio_
pca.singular_values_
pca.components_

x1_df = pd.DataFrame(X2, columns = ["PCA1", "PCA2"])

res_df = pd.concat([results_df, x1_df], axis = 1)

filename = os.path.join(raw_dir, "pca_chem_post.csv")
res_df.to_csv(filename, index = False)

g = sns.FacetGrid(res_df, col = "Variance", row = "Regime")
g.map_dataframe(sns.scatterplot, x = "PCA1",
                y = "PCA2", hue = "Normalised_vel",
                #size = "Variance",
                cmap = "flare",
                alpha = 0.6,
                edgecolor = "None")
g.add_legend()

sns.scatterplot(x = "PCA1",y= "PCA2", data = res_df, style = "Regime", 
                size = "Normalised_vel",
                hue = "Normalised_vel",cmap = "flare_r",
                edgecolor = "None")