# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 17:21:42 2020

@author: khurana
"""

#temporal_variation_factor_analysis
#!pip install factor_analyzer
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
import seaborn as sns
import statsmodels.formula.api as smf

#Figure S1: Autocorrelation discussion for imposed temporal dynamics

head_path = "Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes/headatinlet.csv"
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

cleanup_tim = {"cov" : {1 : cov1, 2 : cov2, 5 : cov3},
               "timtrace" : {1 : cor1tim, 2: cor2tim, 5 : cor3tim},
               "chem_factors" : {'DOC':0,'DO':1, 'Ammonium':2, 'Nitrate':3, 'Sulphate':4, 'Particulate organic matter':5,
                                 'Mobile active aerobic degraders':6, 'Mobile active nitrate reducers':7,
                                 'Mobile active sulphate reducers':8, 'Mobile active ammonia oxidizers':9,
                                 'Mobile inactive aerobic degraders':10, 'Mobile inactive nitrate reducers':11,
                                 'Mobile inactive sulphate reducers':12, 'Mobile inactive ammonia oxidizers':13,
                                 'Nitrogen':14, 'TOC':15},
                   "Pe" : {'Slow':2,'Medium':11, 'Fast':22}}

directory = r"Y:\Home\khurana\4. Publications\Restructuring\Paper2\Figurecodes\/"
filename = "Normalized_RMSamplitude_chem.csv"
data = pd.read_csv(directory + filename, sep="\t")
print(data.columns)
print(data.shape)
print(data.Chem.unique())
data["Regime"] = data["Regime"].replace(["Equal"], "Medium")
chemstoplot = ["DOC", "DO", "TOC", "Nitrogen"]
data["cov"] = data ["Time_series"]

data  = pd.read_csv("Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/mass_flux_sensitivity_generalized.csv", sep="\t")

data["timtrace"] = data ["Time_series"]
data["chem_factors"] = data ["Chem"]
data["Pe"] = data ["Regime"]
data.replace(cleanup_tim, inplace = True)

dadata = pd.read_csv("Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/Conc_da_ss.csv", sep = "\t")

mdata = pd.merge(data, dadata[["Regime", "Trial", "Chem", "%reldelmassflux", "%fraction_rel_delmf"]], on = ["Regime", "Trial", "Chem"])
mdata.loc[mdata["PeDa"] > 40, "PeDamark"] = 3
mdata.loc[(mdata["PeDa"] > 20) & (mdata["PeDa"] < 40), "PeDamark"] = 2
mdata.loc[(mdata["PeDa"] > 1) & (mdata["PeDa"] < 20), "PeDamark"] = 1
mdata.loc[mdata["PeDa"] < 1, "PeDamark"] = 0
labels = {0 : "Da/Pe < 1",
          1 : "1 < Da/Pe < 15",
          2 : "15 < Da/Pe < 40",
          3 : "Da/Pe > 40"}

corrMatrix = mdata.drop(['Unnamed: 0', 'Unnamed: 0.1', 'Sensitivity%', 'Sensitivitybase%', 'PeDamark', 'sensbase', 'timtrace', 'chem_factors', 'Time_series'], axis = 1).corr()
sns.heatmap(corrMatrix, annot = True)
mask = np.zeros_like(corrMatrix)
mask[np.triu_indices_from(mask)] = True
with sns.axes_style("white"):
    p2 = sns.heatmap(corrMatrix, mask=mask, square=True, cmap = "vlag")

for k in ['pearson', 'kendall', 'spearman']:
    corrMatrix = mdata.drop(['Unnamed: 0', 'Unnamed: 0.1', 'Sensitivity%', 'Sensitivitybase%', 'PeDamark', 'sensbase', 'timtrace', 'chem_factors', 'Time_series'], axis = 1).corr(method = k)
    plt.figure()
    sns.heatmap(corrMatrix, annot = True)