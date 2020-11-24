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
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
import statsmodels.formula.api as smf
import statsmodels.api as sm

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

allfeatures = dadata.columns.tolist()
features = ['Trial', 'Variance', 'Anisotropy', 'Pe', 'Breakthroughtime', 'fraction', 'avgconc_in', 'avgconc_out']
#yfeature = ['%reldelmassflux', 'cov', 'timtrace','Time', 'PeDa', '%fraction_rel_delmf']
yfeature = ['%reldelmassflux']

mydata = dadata[features+yfeature]
mydata.rename(columns = {'%reldelmassflux':'perreldelmassflux'}, inplace = True)
mydata.Trial = mydata.Trial.replace(['H'], 0).astype(int)
scaler.fit(mydata)
scaled = scaler.transform(mydata)
mydata_scaled = pd.DataFrame(data = scaled, columns = features + ['perreldelmassflux'])

#OLS
plt.figure()
symbols = ['D', '^', 's', 'o']
colors = ['r', 'g', 'blue', 'orange']
factor_groups = mydata.groupby(['Variance', 'Anisotropy'])
for values, group in factor_groups:
    i,j = values
    plt.scatter(group['Trial'], group['perreldelmassflux'], marker=symbols[j], color=colors[i-1],
               s=144)
plt.xlabel('Experience');
plt.ylabel('Salary');

from urllib.request import urlopen
url = 'http://stats191.stanford.edu/data/salary.table'
fh = urlopen(url)
salary_table = pd.read_table(fh)
salary_table.to_csv('salary.table')

E = salary_table.E
M = salary_table.M
X = salary_table.X
S = salary_table.S

symbols = ['D', '^']
colors = ['r', 'g', 'blue']
factor_groups = salary_table.groupby(['E','M'])
for values, group in factor_groups:
    i,j = values
    plt.scatter(group['X'], group['S'], marker=symbols[j], color=colors[i-1],
               s=144)
plt.xlabel('Experience');
plt.ylabel('Salary');

anovasummaries = pd.DataFrame(columns = summarycols)
for g in groups:
    f1 = 'perreldelmassflux ~ '
    for c in reversed(indepvar):
        f1 += ' + ' + c
        print(f1)
        mymd = smf.mixedlm(formula = f1, data = mydata_scaled, groups = mydata_scaled[g])
        summary = sm.stats.anova_lm(mymd, typ = 2)
        mdf = mymd.fit()
        summary = pd.concat([mdf.summary().tables[0],mdf.summary().tables[1]], axis = 0)
        summary["Group"] = g
        allsummaries = pd.concat([allsummaries, summary[summarycols]],axis = 0) 

anovasummaries.to_csv("X:/Grouped_summary_ANOVA.csv", sep = ',')

#LMEM on steady state results
summarycols = ["Group", 0, 1, 2, 3, 'Coef.', 'Std.Err.', 'z', 'P>|z|', '[0.025', '0.975]']
indepvar = ['Breakthroughtime', 'fraction', 'avgconc_in', 'avgconc_out']
#yfeature = ['%reldelmassflux', 'cov', 'timtrace','Time', 'PeDa', '%fraction_rel_delmf']
yfeature = ['%reldelmassflux']
groups = ["Pe", "Trial", "Variance", "Anisotropy"]
allsummaries = pd.DataFrame(columns = summarycols)
for g in groups:
    f1 = 'perreldelmassflux ~ '
    for c in reversed(indepvar):
        f1 += ' + ' + c
        print(f1)
        mymd = smf.mixedlm(formula = f1, data = mydata_scaled, groups = mydata_scaled[g])
        mdf = mymd.fit()
        summary = pd.concat([mdf.summary().tables[0],mdf.summary().tables[1]], axis = 0)
        summary["Group"] = g
        allsummaries = pd.concat([allsummaries, summary[summarycols]],axis = 0) 

allsummaries.to_csv("X:/Grouped_summary_.csv", sep = ',')

#ANOVA



f1 = 'perreldelmassflux ~ Variance + Anisotropy + Pe + Breakthroughtime + fraction + avgconc_in + avgconc_out'

import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    likev = mdf.profile_re(0, dist_low=0.00000000001, dist_high=0.1, vtype = 're')
    
plt.figure(figsize=(10,8))
plt.plot(likev[:,0], 2*likev[:,1])
plt.xlabel("Variance of random slope", size=17)
plt.ylabel("-2 times profile log likelihood", size=17)

re = mdf.cov_re.iloc[1, 1]
likev = mdf.profile_re(1, dist_low=.5*re, dist_high=0.8*re)

plt.figure(figsize=(10, 8))
plt.plot(likev[:,0], 2*likev[:,1])
plt.xlabel("Variance of random slope", size=17)
plt.ylabel("-2 times profile log likelihood", size=17)