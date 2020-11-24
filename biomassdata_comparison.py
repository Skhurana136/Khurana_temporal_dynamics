# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 12:05:06 2020

@author: khurana
"""

import pandas as pd

#Load data
path_data = "//msg-filer2/scratch_60_days/khurana/biomass_Original_complete.csv"
data = pd.read_csv(path_data, sep = "\t")
data.columns
data.dtypes

regimes = data.Regime.unique().tolist()
Time_series = data.Time_series.unique().tolist()
chem_series = data.Chem.unique().tolist()
trial_series = data.Trial.unique().tolist()

for r in regimes:
    for t in trial_series:
        for c in chem_series:
            mass_base = data.loc[(data.Regime == r) & (data.Chem == c) & (data.Trial == t) & (data.Time_series == 0)]['Mass'].values[0]
            cont_base = data.loc[(data.Regime == r) & (data.Chem == c) & (data.Trial == t) & (data.Time_series == 0)]['Mass'].values[0]
            data.loc[(data.Regime == r) & (data.Chem == c) & (data.Trial == t), 'Mass_base'] = mass_base
            data.loc[(data.Regime == r) & (data.Chem == c) & (data.Trial == t), 'Cont_base'] = cont_base
        
data['meanmass_fraction'] = data['Mass']/data['Mass_base']
data['contribution_fraction'] = data['Contribution']/data['Cont_base']

data.to_csv("//msg-filer2/scratch_60_days/khurana/biomass_comparison_Original_complete.csv", sep ="\t")