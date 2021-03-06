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
import pandas as pd
import data_reader.data_processing as proc
import analyses.transient as ta

#set up basic constants 
Regimes = ["Slow", "Equal", "Fast"]
scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run

# Scenarios to investigate:
Trial = list(scdict.keys())
reginvest = Regimes

vardict = proc.speciesdict("Saturated")
gvarnames = list(t for t in vardict.keys() if vardict[t]["Location"]=="Mobile") + ["Nitrogen", "TOC"]

parent_directory = "E:/Zenodo_temporal_dynamics"
output_directory = "Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes"
#Sensitivity
row = []
for Reg in reginvest:
    basefile = os.path.join(parent_directory, Reg+"AR_0_NS-AH_df.npy")
    basedata = np.load(basefile)
    for t in ["1", "2", "5"]:
       print (Reg, t)
       for j in Trial:
           if ((j == '52' and t == "5") or (j == '43' and t == "1")):
               pass
           else:
               file = os.path.join(parent_directory, Reg+"AR_"+str(t)+"_NS-A"+str(j)+"_df.npy")
               data = np.load(file)
               amplitude, ampmax, amplitudebase, basemax = ta.conc_norm_amplitude(data, basedata, 0, -1, 0, -1, 51, gvarnames, "Saturated")
               for g in gvarnames:
                   row.append([j,scdict[j]['Het'], scdict[j]['Anis'], Reg, t, g, amplitude[gvarnames.index(g)], ampmax[gvarnames.index(g)], amplitudebase[gvarnames.index(g)], basemax[gvarnames.index(g)]])

sensitivitydata = pd.DataFrame.from_records (row, columns = ["Trial", "Variance", "Anisotropy", "Regime", "Time_series", "Chem", "Sensitivity", "Timloc_max", "Sensitivitybase", "Timloc_maxbase"])

tracerdata = pd.read_csv("Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/tracer_combined_05032020.csv", sep = "\t")

ampbth = pd.merge(sensitivitydata, tracerdata[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"])

ampbth["Sensitivity%"] = ampbth["Sensitivity"] * 100
ampbth["Sensitivitybase%"] = ampbth["Sensitivitybase"] * 100

ampbth.to_csv(os.path.join(output_directory, "normalized_sensivity_19022021.csv"))

import numpy as np
import pandas as pd
import data_reader.data_processing as proc
import analyses.transient as ta

#set up basic constants 
Regimes = ["Slow", "Equal", "Fast"]
scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run

# Scenarios to investigate:
Trial = list(scdict.keys())
#Trial = ["H", "37", "38", "39", "40", "41", "42", "43", "44", "45"]
reginvest = Regimes

vardict = proc.speciesdict("Saturated")
States = ["Active", "Inactive"]
gvarnames = list(t for t in vardict.keys() if vardict[t]["State"] in (States))

#Sensitivity
row = []
for Reg in reginvest:
    benchmark = np.load(os.path.join(parent_directory, Reg + "AR_0_NS-AH_df.npy"))
    for t in ["1", "2", "5"]:
        print (Reg, t)
        for j in Trial:
            if ((j == '52' and t == "5") or (j == '43' and t == "1")):
                pass
            else:
                data = np.load(os.path.join(parent_directory, Reg + "AR_"+str(t)+"_NS-A"+str(j)+"_df.npy"))
                amplitude, ampmax, baseamp, basemax = ta.mass_norm_amplitude(data, benchmark, 0, -1, 0, -1, 51, gvarnames, "Saturated")
                for g in gvarnames:
                    row.append([j,scdict[j]['Het'], scdict[j]['Anis'], Reg, t, g, amplitude[gvarnames.index(g)], ampmax[gvarnames.index(g)], baseamp[gvarnames.index(g)], baseamp[gvarnames.index(g)]])

sensitivitydata = pd.DataFrame.from_records (row, columns = ["Trial", "Variance", "Anisotropy", "Regime", "Time_series", "Chem", "Sensitivity", "Timloc_max", "Sensitivitybase", "Timloc_maxbase"])

tracerdata = pd.read_csv("Z:/tracer_combined_05032020.csv", sep = "\t")

ampbth = pd.merge(sensitivitydata, tracerdata[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"])

ampbth["Sensitivity%"] = ampbth["Sensitivity"] * 100
ampbth["Sensitivitybase%"] = ampbth["Sensitivitybase"] * 100

ampbth.to_csv(os.path.join(output_directory, "Normalized_RMSamplitude_biomass_190022021.csv"))

#Sensitivity comparison
head_path = os.path.join(output_directory, "headatinlet.csv")
head = pd.read_csv(head_path, sep = ",")
cov1 = np.round(np.cov(head["H1"]),2)
cov2 = np.round(np.cov(head["H2"]),2)
cov3 = np.round(np.cov(head["H3"]),2)

from statsmodels.tsa import stattools
cov1 = np.round(stattools.acovf(head["H1"], fft = True),2)[0]
cov2 = np.round(stattools.acovf(head["H2"], fft = True),2)[0]
cov3 = np.round(stattools.acovf(head["H3"], fft = True),2)[0]

cov1tim = np.where(stattools.acf(head["H1"], fft = True, nlags = 5476) < 0.75)[0][0]
cov2tim = np.where(stattools.acf(head["H2"], fft = True, nlags = 5476) < 0.75)[0][0]
cov3tim = np.where(stattools.acf(head["H3"], fft = True, nlags = 5476) < 0.75)[0][0]

path_da_data = "Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/Da_29012021_95pcloss.csv"
da = pd.read_csv(path_da_data, sep = ",")

filename = "normalized_sensivity_19022021.csv"
sens = pd.read_csv(os.path.join(output_directory, filename))
sens["Regime"] = sens["Regime"].replace(["Equal"], "Medium")

print(da.columns)
print(sens.columns)

data = pd.merge(sens, da[["Da63", "Pe", "Regime", "Trial", "Chem"]], on = ["Regime", "Trial", "Chem"])
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

data.to_csv(os.path.join(output_directory, "mass_flux_sensitivity_generalized_19022021.csv"))

head_path = os.path.join(output_directory, "headatinlet.csv")
head = pd.read_csv(head_path, sep = ",")
cov1 = np.round(np.cov(head["H1"]),2)
cov2 = np.round(np.cov(head["H2"]),2)
cov3 = np.round(np.cov(head["H3"]),2)

from statsmodels.tsa import stattools
cov1 = np.round(stattools.acovf(head["H1"], fft = True),2)[0]
cov2 = np.round(stattools.acovf(head["H2"], fft = True),2)[0]
cov3 = np.round(stattools.acovf(head["H3"], fft = True),2)[0]

cov1tim = np.where(stattools.acf(head["H1"], fft = True, nlags = 5476) < 0.75)[0][0]
cov2tim = np.where(stattools.acf(head["H2"], fft = True, nlags = 5476) < 0.75)[0][0]
cov3tim = np.where(stattools.acf(head["H3"], fft = True, nlags = 5476) < 0.75)[0][0]

filename = "Normalized_RMSamplitude_biomass_190022021.csv"
data = pd.read_csv(os.path.join(output_directory, filename))
data["Regime"] = data["Regime"].replace(["Equal"], "Medium")

print(sens.columns)

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

data.to_csv(os.path.join(output_directory, "biomass_sensitivity_generalized_19022021.csv"))

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
            directory = "E:/Saturated_flow/EGUGoldschmidtdataset6/" + domadd + Reg + "AR_" + t + "/"
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