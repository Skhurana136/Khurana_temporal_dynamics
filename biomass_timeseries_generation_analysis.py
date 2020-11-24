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
import numpy as np
import pandas as pd
import data_reader.data_processing as proc
import analyses.saturated_transient as sta

#set up basic constants 
Regimes = ["Slow", "Equal", "Fast"]
domains = ["Original", "Half", "Double", "Big"]
domainodes = {"Original": {'ynodes' : 51},
              "Big" : {'ynodes' : 126},
              "Double" : {'ynodes' : 101},
              "Half" : {'ynodes' : 26}}
scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run

# Scenarios to investigate:
Trial = list(scdict.keys())
#Trial = ["H", "37", "38", "39", "40", "41", "42", "43", "44", "45"]
reginvest = Regimes
domaininvest = list(domainodes.keys())[:1]

vardict = proc.speciesdict("Saturated")
States = ["Active", "Inactive"]
gvarnames = list(t for t in vardict.keys() if vardict[t]["State"] in (States))

#Sensitivity
row = []
for Reg in reginvest:
    for domain in domaininvest:
        if domain != "Original":
            domadd = domain + "_"
        else:
            domadd = ""
        benchmark = np.load("D:/Saturated_flow/EGUGoldschmidtdataset6/" + domadd + Reg + "AR_0/NS-AH/NS-AH_df.npy")
        for t in ["1", "2", "5"]:
#            directory = "//tsclient/D/Saturated_flow/EGUGoldschmidtdataset6/" + domadd + Reg + "AR_" + t + "/"
            directory = "D:/Saturated_flow/EGUGoldschmidtdataset6/" + domadd + Reg + "AR_" + t + "/"
            #directory = "X:/Saturated_flow/changedkindox_transient/" + domadd + Reg + "AR_" + t + "/"#change directory as per flow regime
            print (Reg, domain, t)
            for j in Trial:
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    data = np.load(directory + "NS-A"+j+"/NS-A"+j+"_df.npy")
                    amplitude, baseamp = sta.mass_norm_amplitude(data, benchmark, 0, -1, 0, -1, domainodes[domain]['ynodes'], gvarnames, "Saturated")
                    for g in gvarnames:
                        row.append([j,scdict[j]['Het'], scdict[j]['Anis'], domain, Reg, t, g, amplitude[gvarnames.index(g)], baseamp[gvarnames.index(g)]])

sensitivitydata = pd.DataFrame.from_records (row, columns = ["Trial", "Variance", "Anisotropy", "Domain", "Regime", "Time_series", "Chem", "Sensitivity", "Sensitivitybase"])

tracerdata = pd.read_csv("Z:/tracer_combined_05032020.csv", sep = "\t")

ampbth = pd.merge(sensitivitydata, tracerdata[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"])

ampbth["Sensitivity%"] = ampbth["Sensitivity"] * 100
ampbth["Sensitivitybase%"] = ampbth["Sensitivitybase"] * 100

ampbth.to_csv("Z:/Normalized_RMSamplitude_biomass.csv", sep="\t")

#Cross-correlation

criteria = 0.7
datafreq = 5
row = []
for Reg in reginvest:
    for domain in domaininvest:
        if domain != "Original":
            domass = domain + "_"
        else:
            domadd = ""
        for t in ["1", "2", "5"]:
#            directory = "//tsclient/D/Saturated_flow/EGUGoldschmidtdataset6/" + domadd + Reg + "AR_" + t + "/"
            directory = "D:/Saturated_flow/EGUGoldschmidtdataset6/" + domadd + Reg + "AR_" + t + "/"
            #directory = "X:/Saturated_flow/changedkindox_transient/" + domadd + Reg + "AR_" + t + "/"#change directory as per flow regime
            print (Reg, domain, t)
            for j in Trial:
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    data = np.load(directory + "NS-A"+j+"/NS-A"+j+"_df.npy")
                    acfchem, Headinlettime = sta.mass_correlation(data, 0, -1, 0, -1, domainodes[domain]['ynodes'], gvarnames, "Saturated")
                    for g in gvarnames:
                        k = gvarnames.index(g)
                        maxchem = np.argmax(np.abs(acfchem[np.shape(Headinlettime)[0] - 1 :, k]))
                        ychem = acfchem[np.shape(Headinlettime)[0] - 1 + maxchem :, k]
                        memorychem = np.where((ychem <= criteria))[0][0]
                        val = acfchem[np.shape(Headinlettime)[0] - 1 + maxchem, k]
                        row.append([j,scdict[j]['Het'], scdict[j]['Anis'], domain, Reg, t, g, maxchem*datafreq, memorychem*datafreq, val])

Ampchem = pd.DataFrame.from_records(row,
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Domain",
        "Regime",
        "Time_series",
        "Chem",
        "Delay",
        "Memory",
        "Crosscorrelation"])

tracerdata = pd.read_csv("Z:/tracer_combined_05032020.csv", sep = "\t")

dfall2 = pd.merge(Ampchem, tracerdata[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"])
dfall2.to_csv("Z:/crosschor_memory_biomass.csv", sep="\t")