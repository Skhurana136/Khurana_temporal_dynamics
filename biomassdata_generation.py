# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 12:05:06 2020

@author: khurana
"""

import numpy as np
import data_reader.data_processing as proc
import pandas as pd
import analyses.saturated_steady_state as sssa
import analyses.saturated_transient as sta

#set up basic constants 
Regimes = ["Slow", "Equal", "Fast"]
domains = ["Original", "Half", "Double", "Big"]
domainodes = {"Original": {'ynodes' : 51},
              "Big" : {'ynodes' : 126},
              "Double" : {'ynodes' : 101},
              "Half" : {'ynodes' : 26}}
fpre = "NS-A"
scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run
fsuf = r"/"
filename = "model_domain_quad.tec" #same filename in each subfolder
horiznodes = 31

# Scenarios to investigate:
Trial = list(scdict.keys())
#Trial = ["H", "37", "38", "39", "40", "41", "42", "43", "44", "45"]
reginvest = Regimes
domaininvest = list(domainodes.keys())[:1]

vardict = proc.speciesdict("Saturated")
states = ["Active", "Inactive"]
gvarnames = list(t for t in vardict.keys() if vardict[t]["State"] in (states))

row = []
for Reg in reginvest:
    for domain in domaininvest:
        if domain != "Original":
            domadd = domain + "_"
        else:
            domadd = ""
        for t in ["1"]:
            directory = "D:/Saturated_flow/EGUGoldschmidtdataset6/" + domadd + Reg + "AR_" + t + "/"
            #directory = "X:/Saturated_flow/changedkindox_transient/" + domadd + Reg + "AR_" + t + "/"#change directory as per flow regime
            print (Reg, domain, t)
            for j in Trial[8:]:
                data = np.load(directory + "NS-A"+j+"/NS-A"+j+"_df.npy")
                if t == "0":
                    meanmass = sssa.calcsum(data, 0, -1, 0, -1, gvarnames, "Saturated")
                else:
                    if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                        pass
                    else:
                        mass = sta.calcsum_temp(data, 0, -1, 0, -1, gvarnames, "Saturated")
#                        mass = sta.biomasstimefunc(data, 0, -1, 0, -1, domainodes[domain]['ynodes'], gvarnames, "Saturated")
                        meanmass = np.mean(mass, axis = 0)
                summass = sum(meanmass)
                masscontribution = meanmass/summass
                for g in gvarnames + ["Total"]:
                    if g == "Total":
                        row.append([j,scdict[j]['Het'], scdict[j]['Anis'], domain, Reg, t, g, summass, 1])
                    else:
                        row.append([j,scdict[j]['Het'], scdict[j]['Anis'], domain, Reg, t, g, meanmass[gvarnames.index(g)], masscontribution[gvarnames.index(g)]])

massdata = pd.DataFrame.from_records (row, columns = ["Trial", "Variance", "Anisotropy", "Domain", "Regime", "Time_series", "Chem", "Mass", "Contribution"])

#Load tracer data
path_tr_data = "Z:/tracer_combined_05032020.csv"
tr_data = pd.read_csv(path_tr_data, sep = "\t")
tr_data.columns

#Merge the datasets and save
cdata = pd.merge(massdata, tr_data[["Trial", "Regime", "Time", "fraction"]], on = ["Regime", "Trial"])

cdata.to_csv("//msg-filer2/scratch_60_days/khurana/biomass_Original_complete.csv", sep = "\t")