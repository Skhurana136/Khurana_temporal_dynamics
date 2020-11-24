# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 21:20:34 2020

@author: khurana
"""

#Creating master file of all the outputs
import pandas as pd

# Saturated flow regime
directory = r"Y:\Home\khurana\4. Publications\Restructuring\Paper2\Figurecodes\/"
filename = "mass_flux_temporal_impact.csv"
data = pd.read_csv(directory + filename, sep = "\t")

print(data.columns)
print(data.shape)

trial = data.Trial.unique().tolist()
het = data.Variance.unique().tolist()
anis = data.Anisotropy.unique().tolist()
time_series = data.Time_series.unique().tolist()
regime_series = data.Regime.unique().tolist()

chemstoplot = ["DOC", "DO", "Ammonium", "Nitrate", "Nitrogen", "TOC"]
biomasstoplot = [
    "Active fixed Aerobes",
    "Active fixed Ammonia oxidizers",
    "Active fixed Nitrate reducers",
    "Active mobile Aerobes",
    "Active mobile Ammonia oxidizers",
    "Active mobile Nitrate reducers",
    "Inactive fixed Aerobes",
    "Inactive fixed Ammonia oxidizers",
    "Inactive fixed Nitrate reducers",
    "Inactive mobile Aerobes",
    "Inactive mobile Ammonia oxidizers",
    "Inactive mobile Nitrate reducers",
]

data["Relative_Removal"] = 100 * (data["Massflux_in"] - data["Massflux_out"])/data["Massflux_in"]

subset = data[data["Chem"].isin (chemstoplot)]

basecase = subset[subset["Trial"] == "H"]

print(basecase.set_index(['Regime', 'Time_series','Chem']))

basecase[basecase['Chem']=="DO"]