# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 12:05:06 2020

@author: khurana
"""

import pandas as pd

#Load data
path_data = "//msg-filer2/scratch_60_days/khurana/massflux_Original_complete.csv"
data = pd.read_csv(path_data, sep = "\t")
data.columns
data.dtypes

regimes = data.Regime.unique().tolist()
Time_series = data.Time_series.unique().tolist()
chem_series = data.Chem.unique().tolist()
trial_series = data.Trial.unique().tolist()
spatial_base = 'H'
temp_base = 0
data['normmassflux'] = 1/data['normmassflux']

for r in regimes:
    for t in trial_series:
        for c in chem_series:
            rel_base = data.loc[(data.Regime == r) & (data.Chem == c) & (data.Trial == t) & (data.Time_series == 0)]['normmassflux'].values[0]
            print(rel_base)
            data.loc[(data.Regime == r) & (data.Chem == c) & (data.Trial == t), 'relbase'] = rel_base
        
data['normmassflux_fraction'] = data['normmassflux']/data['relbase']

data.to_csv("//msg-filer2/scratch_60_days/khurana/massflux_comparison_Original_complete.csv", sep ="\t")