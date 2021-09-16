# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 12:05:06 2020

@author: khurana
"""
import os
import pandas as pd

#Load data
path_dir = "Y:/Home/khurana/4. Publications/Restructuring/Paper2/Figurecodes"
data_filename = "massflux_Original_complete.csv"
data = pd.read_csv(os.path.join(path_dir, data_filename), sep = "\t")
data.columns
data.dtypes

regimes = data.Regime.unique().tolist()
Time_series = data.Time_series.unique().tolist()
chem_series = data.Chem.unique().tolist()
trial_series = data.Trial.unique().tolist()

for r in regimes:
    for t in trial_series:
        for c in chem_series:
            reldelmf_h_base = data.loc[(data.Regime == r) & (data.Chem == c) & (data.Trial == "H") & (data.Time_series == 0)]['reldelmassflux'].values[0]
            reldelmf_t_base = data.loc[(data.Regime == r) & (data.Chem == c) & (data.Trial == t) & (data.Time_series == 0)]['reldelmassflux'].values[0]
            data.loc[(data.Regime == r) & (data.Chem == c) & (data.Trial == t), 'reldelmf_h_base'] = reldelmf_h_base
            data.loc[(data.Regime == r) & (data.Chem == c) & (data.Trial == t), 'reldelmf_t_base'] = reldelmf_t_base
        
data['reldelmf_h_fraction'] = data['reldelmassflux']/data['reldelmf_h_base']
data['reldelmf_t_fraction'] = data['reldelmassflux']/data['reldelmf_t_base']

data.to_csv(os.path.join(path_dir, "massflux_comparison_Original_complete_16062021.csv"), index = False)