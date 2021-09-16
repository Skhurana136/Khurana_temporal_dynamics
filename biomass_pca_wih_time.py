# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 17:21:42 2020

@author: khurana
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


from sklearn.decomposition import PCA

from sklearn.preprocessing import StandardScaler

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
gvarnames = list(t for t in vardict.keys() if vardict[t]["State"] in (States))

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
                mass = ta.biomass_time(data, 0, -1, 0, -1, gvarnames, "Saturated")
                summass = np.sum(mass, axis = 1)
                velocity = np.mean(data[2, 1:, :, :], axis = (-1,-2))/basevelocity
                g_df = pd.DataFrame(np.column_stack([reglist, serieslist, jlist, vlist, alist, velocity]), 
                               columns=['Regime','Time_series', 'Trial', 'Variance', 'Anisotropy', "Normalised_vel"])
                for g in gvarnames:
                    ratio = np.divide(mass[:,gvarnames.index(g)],summass)
                    ratio_df = pd.DataFrame(ratio, columns = [g])
                    g_df = pd.concat([g_df,ratio_df], axis = 1)
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

filename = os.path.join(raw_dir, "pca.csv")
results_df.to_csv(filename, index=False)


filename = os.path.join(raw_dir, "pca.csv")
results_df = pd.read_csv(filename)

gvarnames_short = list(t for t in gvarnames if 'sulphate' not in t)
X = results_df[gvarnames_short+['Variance']]
y = results_df["uid"]

pca = PCA(n_components=2)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
X1 = pca.fit_transform(X_scaled)

pca.explained_variance_ratio_
pca.singular_values_
pca.components_

weights = pca.components_

pca_weightst = pd.DataFrame(weights, columns = gvarnames_short+['Variance'])
pca_weights = pca_weightst.T

plt.figure(figsize = (2,4))
b = sns.heatmap(pca_weights, cmap = "vlag",
                xticklabels = ["Component\n1", "Component\n2"])
b.tick_params(labelsize = 12)
plt.xticks(rotation = 60)

x1_df = pd.DataFrame(X1, columns = ["PCA1", "PCA2"])

res_df = pd.concat([results_df, x1_df], axis = 1)

filename = os.path.join(raw_dir, "pca_post_no_sulphate.csv")
res_df.to_csv(filename, index = False)

res_df=pd.read_csv(filename)
res_df_het = res_df[res_df['Variance'] != 0]

r = "Equal"
res_df_r = res_df_het[res_df_het['Regime'] == r]
g = sns.FacetGrid(res_df_r, col = "Variance", row = "Anisotropy")
g.map_dataframe(sns.scatterplot, x = "PCA1",
                y = "PCA2", hue = "Normalised_vel",
                style = "Regime",
                cmap = "flare",
                alpha = 0.6,
                edgecolor = "None")
g.add_legend()
picname = os.path.join(out_dir, r+"_biomass_pca.png")
g.savefig(picname, dpi = 300, pad = 0.01, bbox_to_inches = 'tight')

