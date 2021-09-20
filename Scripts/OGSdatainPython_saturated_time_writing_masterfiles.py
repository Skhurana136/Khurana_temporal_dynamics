# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 21:20:34 2020

@author: khurana
"""

#Creating master file of all the outputs
import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt
import data_reader.data_processing as proc
import analyses.saturated_transient as sta
import plots.saturated_transient as stp
import data_reader.reader as rdr
import os

# Saturated flow regime
Reg = "Fast"
directory = r"Z:/Saturated_flow/diffusion_transient/"
fpre = "/NS-A"
fsuf = r"/"
gw = 1

scdict = proc.masterscenarios() #master dictionary of all spatially heterogeneous scenarios that were run

# Default:
Trial = list(t for t,values in scdict.items() if t not in ['43', '52'])
Het = list(values['Het'] for t,values in scdict.items() if t not in ['43', '52'])
Anis = list(values['Anis'] for t,values in scdict.items() if t not in ['43', '52'])

Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
nScenarios = 3

# Constants
yout = 50
yin = 0
xleft = 0
xright = 30
# Assign index to Variable
doc1 = 10 - gw
dox1 = 11 - gw
Amm1 = 12 - gw
nitra1 = 17 - gw
sulpha1 = 22 - gw
tr1 = 29 - gw
Bfo1 = 8 - gw
Bfn1 = 15 - gw
Bfs1 = 20 - gw
Bfa1 = 25 - gw
Bmo1 = 9 - gw
Bmn1 = 16 - gw
Bms1 = 21 - gw
Bma1 = 26 - gw
Bifo1 = 13 - gw
Bifn1 = 18 - gw
Bifs1 = 23 - gw
Bifa1 = 27 - gw
Bimo1 = 14 - gw
Bimn1 = 19 - gw
Bims1 = 24 - gw
Bima1 = 28 - gw
head1 = 1 - gw
vely = 5
velx = 4
vars = [doc1, dox1, Amm1, nitra1, sulpha1, Bmo1, Bma1, Bmn1, Bimo1, Bima1, Bimn1, head1]
gvarnames = [
    "DOC",
    "DO",
    "Ammonium",
    "Nitrate",
    "Sulphate",
    "Active mobile Aerobes",
    "Active mobile Ammonia oxidizers",
    "Active mobile Nitrate reducers",
    "Inactive mobile Aerobes",
    "Inactive mobile Ammonia oxidizers",
    "Inactive mobile Nitrate reducers",
    "Head",
    "Nitrogen",
    "TOC",
]
AFbiomassvars = [
    Bfo1,
    Bfa1,
    Bfn1,
    Bmo1,
    Bma1,
    Bmn1,
    Bifo1,
    Bifa1,
    Bifn1,
    Bimo1,
    Bima1,
    Bimn1,
]
AFbiomassgvarnames = [
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

mastermf = pd.DataFrame(
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Chem",
        "Massflux_in",
        "Massflux_out",
        "Removal",
        "Normalized_removal",
        "UniformHet_normalized_removal",
        "Time_series",
        "Varyinghom_normalized_removal",
        "Regime",
    ]
)
masterbiomass = pd.DataFrame(
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Chem",
        "Total_biomass_umoles",
        "Hom_normalized_change",
        "Unihet_normalized_change",
        "Varyinghom_normalized_change",
        "Time_series",
        "species_fraction",
        "Hom_normalized_fractionchange",
        "Unihet_normalized_fractionchange",
        "Varyinghom_normalized_fractionchange",
        "Regime",
    ]
)

Regimes = ["Slow", "Equal", "Fast"]
for Reg in Regimes:
    Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
    mft, calcsum = sta.calcmft_temp(
        Tforfpre,
        Trial,
        vars,
        gvarnames,
        AFbiomassvars,
        AFbiomassgvarnames,
        directory,
        fpre,
        fsuf,
        Het,
        Anis,
        gw,
    )
    
    dfmft = pd.DataFrame(
        mft,
        columns=[
            "Trial",
            "Variance",
            "Anisotropy",
            "Chem",
            "Massflux_in",
            "Massflux_out",
            "Removal",
            "Normalized_removal",
            "UniformHet_normalized_removal",
            "Time_series",
            "Varyinghom_normalized_removal",
            "Regime",
        ],
    )
    
    dfmft["Regime"] = Regimes.index(Reg)
    
    dfbiomasssum = pd.DataFrame(
        calcsum,
        columns=[
            "Trial",
            "Variance",
            "Anisotropy",
            "Chem",
            "Total_biomass_umoles",
            "Hom_normalized_change",
            "Unihet_normalized_change",
            "Varyinghom_normalized_change",
            "Time_series",
            "species_fraction",
            "Hom_normalized_fractionchange",
            "Unihet_normalized_fractionchange",
            "Varyinghom_normalized_fractionchange",
            "Regime",
        ],
    )
    
    dfbiomasssum["Regime"] = Regimes.index(Reg)
    
    dfmft2 = proc.processdataframe(dfmft, gvarnames)
    dfbiomass2 = proc.processdataframe(dfbiomasssum, AFbiomassgvarnames)
    
    fnamemf = (
        "Z:/Saturated_flow/diffusion_transient/mass_flux_temporal_impact_"
        + Reg
        + ".csv"
    )
    fnamebiomasssum = (
        "Z:/Saturated_flow/diffusion_transient/biomass_temporal_impact_" + Reg + ".csv"
    )
    
    dfmft2.to_csv(fnamemf, sep="\t")
    
    dfbiomass2.to_csv(fnamebiomasssum, sep="\t")
    
    mastermf = pd.concat([mastermf, dfmft])
    masterbiomass = pd.concat([masterbiomass, dfbiomasssum])

    print("Files written, processing as dataframes ...")

mastermf2 = proc.processdataframe(mastermf, gvarnames)
masterbiomass2 = proc.processdataframe(masterbiomass, AFbiomassgvarnames)
mastermf.to_csv(
    "Z:/Saturated_flow/diffusion_transient/mass_flux_temporal_impact.csv", sep="\t"
)
masterbiomass.to_csv(
    "Z:/Saturated_flow/diffusion_transient/biomass_temporal_impact.csv", sep="\t"
)