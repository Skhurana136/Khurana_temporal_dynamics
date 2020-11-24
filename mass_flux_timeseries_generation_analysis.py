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
gvarnames = list(t for t in vardict.keys() if vardict[t]["Location"]=="Mobile") + ["Nitrogen", "TOC"]

#Sensitivity
row = []
for Reg in reginvest:
    for domain in domaininvest:
        if domain != "Original":
            domadd = domain + "_"
        else:
            domadd = ""
        basedata = np.load("D:/Saturated_flow/EGUGoldschmidtdataset6/" + domadd + Reg + "AR_0/NS-AH/NS-AH_df.npy")
        for t in ["1", "2", "5"]:
            directory = "D:/Saturated_flow/EGUGoldschmidtdataset6/" + domadd + Reg + "AR_" + t + "/"
            #directory = "X:/Saturated_flow/changedkindox_transient/" + domadd + Reg + "AR_" + t + "/"#change directory as per flow regime
            print (Reg, domain, t)
            for j in Trial:
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    data = np.load(directory + "NS-A"+j+"/NS-A"+j+"_df.npy")
                    amplitude, ampmax, amplitudebase, basemax = sta.conc_norm_amplitude(data, basedata, 0, -1, 0, -1, domainodes[domain]['ynodes'], gvarnames, "Saturated")
                    for g in gvarnames:
                        row.append([j,scdict[j]['Het'], scdict[j]['Anis'], domain, Reg, t, g, amplitude[gvarnames.index(g)], ampmax[gvarnames.index(g)], amplitudebase[gvarnames.index(g)], basemax[gvarnames.index(g)]])

sensitivitydata = pd.DataFrame.from_records (row, columns = ["Trial", "Variance", "Anisotropy", "Domain", "Regime", "Time_series", "Chem", "Sensitivity", "Timloc_max", "Sensitivitybase", "Timloc_maxbase"])

tracerdata = pd.read_csv("Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/tracer_combined_05032020.csv", sep = "\t")

ampbth = pd.merge(sensitivitydata, tracerdata[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"])

ampbth["Sensitivity%"] = ampbth["Sensitivity"] * 100
ampbth["Sensitivitybase%"] = ampbth["Sensitivitybase"] * 100

ampbth.to_csv("Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/Normalized_RMSamplitude_chem_withloc.csv", sep="\t")

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
                    amplitude, ampmax, baseamp, basemax = sta.mass_norm_amplitude(data, benchmark, 0, -1, 0, -1, domainodes[domain]['ynodes'], gvarnames, "Saturated")
                    for g in gvarnames:
                        row.append([j,scdict[j]['Het'], scdict[j]['Anis'], domain, Reg, t, g, amplitude[gvarnames.index(g)], ampmax[gvarnames.index(g)], baseamp[gvarnames.index(g)], baseamp[gvarnames.index(g)]])

sensitivitydata = pd.DataFrame.from_records (row, columns = ["Trial", "Variance", "Anisotropy", "Domain", "Regime", "Time_series", "Chem", "Sensitivity", "Sensitivitybase"])

tracerdata = pd.read_csv("Z:/tracer_combined_05032020.csv", sep = "\t")

ampbth = pd.merge(sensitivitydata, tracerdata[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"])

ampbth["Sensitivity%"] = ampbth["Sensitivity"] * 100
ampbth["Sensitivitybase%"] = ampbth["Sensitivitybase"] * 100

ampbth.to_csv("Z:/Normalized_RMSamplitude_biomass.csv", sep="\t")

#Sensitivity comparison
head_path = "Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/headatinlet.csv"
head = pd.read_csv(head_path, sep = ",")
cov1 = np.round(np.cov(head["H1"]),2)
cov2 = np.round(np.cov(head["H2"]),2)
cov3 = np.round(np.cov(head["H3"]),2)

from statsmodels.tsa import stattools
cov1 = np.round(stattools.acovf(head["H1"], fft = True),2)
cov2 = np.round(stattools.acovf(head["H2"], fft = True),2)
cov3 = np.round(stattools.acovf(head["H3"], fft = True),2)

cov1tim = np.where(stattools.acf(head["H1"], fft = True, nlags = 5476) < 0.75)[0][0]
cov2tim = np.where(stattools.acf(head["H2"], fft = True, nlags = 5476) < 0.75)[0][0]
cov3tim = np.where(stattools.acf(head["H3"], fft = True, nlags = 5476) < 0.75)[0][0]

path_da_data = "Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/Conc_da_ss.csv"
da = pd.read_csv(path_da_data, sep = "\t")

directory = r"Y:\Home\khurana\4. Publications\Restructuring\Paper2\Figurecodes\/"
filename = "Normalized_RMSamplitude_chem.csv"
sens = pd.read_csv(directory + filename, sep="\t")
sens["Regime"] = sens["Regime"].replace(["Equal"], "Medium")

print(da.columns)
print(sens.columns)

data = pd.merge(sens, da[["PeDa", "conc_Da","Pe", "Regime", "Trial", "Chem"]], on = ["Regime", "Trial", "Chem"])
data["cov"]=data["Time_series"]
data["cov"]=data["cov"].replace([1], cov1)
data["cov"]=data["cov"].replace([2], cov2)
data["cov"]=data["cov"].replace([5], cov3)

gvarnames = data.Chem.unique().tolist()
reglist = data.Regime.unique().tolist()

for t in [1,2,5]:
    for r in reglist:
        for g in gvarnames:
            base = data[(data["Time_series"]==t) & (data["Regime"]==r) & (data["Chem"]==g) & (data["Trial"]=='H')]["Sensitivitybase%"].values[0]
            data.loc[(data.Regime == r) & (data.Chem == g) & (data.Time_series == t), 'sensbase'] = base

data["Senssquared"] = data["Sensitivitybase%"]/data["sensbase"]

data.to_csv("Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/mass_flux_sensitivity_generalized.csv", sep="\t")

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
            directory = "D:/Saturated_flow/EGUGoldschmidtdataset6/" + domadd + Reg + "AR_" + t + "/"
            #directory = "X:/Saturated_flow/changedkindox_transient/" + domadd + Reg + "AR_" + t + "/"#change directory as per flow regime
            print (Reg, domain, t)
            for j in Trial:
                if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                    pass
                else:
                    data = np.load(directory + "NS-A"+j+"/NS-A"+j+"_df.npy")
                    acfchem, Headinlettime = sta.correlation(data, 0, -1, 0, -1, domainodes[domain]['ynodes'], gvarnames, "Saturated")
                    for g in gvarnames:
                        k = gvarnames.index(g)
                        maxchem = np.argmax(np.abs(acfchem[np.shape(Headinlettime)[0]-1:, k]))
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

tracerdata = pd.read_csv("Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/tracer_combined_05032020.csv", sep = "\t")

dfall2 = pd.merge(Ampchem, tracerdata[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"])
dfall2.to_csv("Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/crosschor_memory_chem.csv", sep="\t")