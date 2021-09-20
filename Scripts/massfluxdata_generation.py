# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 12:05:06 2020

@author: khurana
"""
import os
import numpy as np
import data_reader.data_processing as proc
import pandas as pd
import analyses.steady_state as sssa
import analyses.transient as sta

#set up basic constants 
Regimes = ["Slow", "Equal", "Fast"]

fpre = "NS-A"
scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run
fsuf = r"/"
filename = "model_domain_quad.tec" #same filename in each subfolder
horiznodes = 31
raw_data = "E:/Saturated_flow/EGUGoldschmidtdataset6"
output_dir = "Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes"
# Scenarios to investigate:
Trial = list(scdict.keys())

vardict = proc.speciesdict("Saturated")
gvarnames = list(vardict.keys()) + ["Nitrogen", "TOC"]

row = []
for Reg in Regimes:
    for t in ["0", "2", "5", "1"]:
        directory = os.path.join(raw_data, Reg + "AR_" + t)
        print (Reg, t)
        for j in Trial:
            if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                pass
            else:
                file = os.path.join(directory, "NS-A"+j, "NS-A"+j+"_df.npy")
                data = np.load(file)
                if t == "0":
                    massfluxin, massfluxout = sssa.massflux(data, 0, -1, 0, -1, gvarnames, "Saturated")
                else:
                    massfluxin, massfluxout = sta.calcmft_temp(data, 0, -1, 0, -1, gvarnames, "Saturated")
                delmassflux = massfluxin - massfluxout
                reldelmassflux = 100*delmassflux/massfluxin
                normmassflux = massfluxout/massfluxin
                for g in gvarnames:
                    row.append([j,scdict[j]['Het'], scdict[j]['Anis'], Reg, t, g, delmassflux[gvarnames.index(g)], reldelmassflux[gvarnames.index(g)], normmassflux[gvarnames.index(g)]])

massfluxdata = pd.DataFrame.from_records (row, columns = ["Trial", "Variance", "Anisotropy", "Domain", "Regime", "Time_series", "Chem", "delmassflux", "reldelmassflux", "normmassflux"])

#Load tracer data
path_tr_data = os.path.join(output_dir, "tracer_combined_05032020.csv")
tr_data = pd.read_csv(path_tr_data, sep = "\t")
tr_data.columns

#Merge the datasets and save
cdata = pd.merge(massfluxdata, tr_data[["Trial", "Regime", "Time", "fraction"]], on = ["Regime", "Trial"])

cdata.to_csv(os.path.join(output_dir, "massflux_complete_28022021.csv"))