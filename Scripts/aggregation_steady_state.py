# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 16:21:55 2020

@author: khurana
"""
#Aggregation/Summarization of results

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plots.saturated_steady_state as sssp
import plots.general_plots as gp
my_pal = {2:"indianred", 11:"g", 22:"steelblue", "DO":"indianred", "Nitrate":"g", "Ammonium":"steelblue",
          "Slow":"indianred", "Medium":"g", "Fast":"steelblue"}

#path_da_data = "Z:/Saturated_flow/diffusion_transient/Complex_Da_mean_ss.csv" #old_file_lost
path_da_data = "Y:/Home/khurana/4. Publications/Restructuring/Paper1/Figurecodes/Conc_da_ss.csv"
da = pd.read_csv(path_da_data, sep = "\t")
da.columns
da.shape

gvarnames = ["DO", "Nitrogen", "TOC"]
finaldata = da[da['Chem'].isin (gvarnames)]

#Show distribution of Da numbers differentiated with Pe numbers investigated
#Kerneldensityplot

def plotdistribution(data, variable):
    sns.distplot(finaldata[data["Regime"] == "Slow"][variable], color = "indianred", label = "Slow flow", kde = False, bins = 20)
    sns.distplot(finaldata[data["Regime"] == "Medium"][variable], color = "g", label = "Medium flow", kde = False, bins = 20)
    sns.distplot(finaldata[data["Regime"] == "Fast"][variable], color = "steelblue", label = "Fast flow", kde = False, bins = 20)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.xlabel ("Damk$\ddot{o}$hler number", fontsize = 15)
    plt.ylabel ("Number of scenarios", fontsize = 15)
    plt.title ("Distribution of Damk$\ddot{o}$hler numbers", fontsize = 15)
    #plt.xscale ("log")
    plt.legend(fontsize = 12)

plotdistribution(da, "PeDa")

#Show variation of Da with spatial heterogeneity

sssp.scatter_chem_regime (finaldata, "fraction", "Spatial heterogeneity", "del2massflux", "Impact on O, N removal", "None")
plt.xlim(left = 0.01)
plt.title("Impact on chemical removal (normalized by that in base case) with spatial heterogeneity")

sssp.scatter_chem_regime(finaldata, "fraction", "Spatial heterogeneity", "fraction_da", "Da", "Logy")
#plt.ylim(bottom = 1)
#plt.ylim(bottom = 10)
plt.title("Concentration derived Da# with spatial heterogeneity")

sssp.scatter_chem_regime(finaldata, "fraction", "Spatial heterogeneity", "fraction_PeDa", "Da", "Logy")
#plt.ylim(bottom = 1)
#plt.ylim(bottom = 10)
plt.title("Concentration derived Da# with spatial heterogeneity")
          
sssp.scatter_chem_regime (finaldata, "PeDa", "PeDa", "%del2massflux", "Impact on O, N removal", "True")
plt.xlim(left = 0.001)
plt.title("Impact on chemical removal (normalized by that in base case) with concentration derived Da")

sssp.scatter_chem_regime (finaldata, "conc_Da", "Da", "del2massflux", "Impact on O, N removal", "Logx")
plt.xlim(left = 1)
plt.title("Impact on chemical removal (normalized by that in base case) with spatial heterogeneity")

#!!!PROMISING!!!!!
sns.boxplot(data = finaldata, x = 'Chem', y = 'PeDa', hue = "Regime", palette = my_pal)
plt.xticks(rotation = 45)
plt.title("Distribution of concentration derived Da")

sssp.scatter_chem_regime(finaldata, "fraction", "Spatial heterogeneity", "fraction_PeDa", "Relative_PeDa", "None")
#plt.ylim(bottom = 1)
#plt.ylim(bottom = 10)
plt.title("Concentration derived Da# with spatial heterogeneity")

import matplotlib.pyplot as plt          
import numpy as  np
plt.scatter(finaldata["fraction"], finaldata["fraction_PeDa"], s = 10*finaldata["PeDabase"], alpha = 0.5)
#plt.ylim(bottom = 1)
#plt.ylim(bottom = 10)
plt.title("Concentration derived Da# with spatial heterogeneity")          

sssp.scatter_chem_regime (finaldata, "%fraction_da", "Concentration derived Da", "%fraction_rel_delmf", "O, N removal", "None")
#plt.xlim(left = 0.001)
plt.title("Impact on chemical removal (normalized by that in base case) with Da")

#Generic plots without distinguishing between chemicals and flow regimes
gp.scatter(finaldata, "fraction", "Spatial heterogeneity", "del2massflux", "Impact on O, N removal", "None")
plt.xlim(left = 0.01)
plt.title("Impact on chemical removal (normalized by that in base case) with spatial heterogeneity")

gp.scatter (finaldata, "fraction", "Spatial heterogeneity", "conc_Da", "Da", "Logy")
plt.ylim(bottom = 0.05)
plt.title("Concentration derived Da with spatial heterogeneity")

gp.scatter (finaldata, "conc_Da", "Da", "del2massflux", "Impact on O, N removal", "Logx")
plt.xlim(left = 1)
plt.title("Impact on chemical removal (normalized by that in base case) with Da")

gp.scatter (finaldata, "PeDa", "PeDa", "%del2massflux", "Impact on O, N removal", "True")
plt.xlim(left = 0.1)
plt.title("Impact on chemical removal (normalized by that in base case) with impact on normalized Da")

gp.scatter (da, "PeDa", "Concentration derived Da", "%fraction_rel_delmf", "O, N removal", "True")
plt.xlim(left = 0.001)
#plt.ylim(top = 1000)
plt.title("Impact on chemical removal (normalized by that in base case) with Da")
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure()
sns.lmplot(data = da, x = "PeDa", y = "%fraction_rel_delmf", hue = "Chem", fit_reg = False)
plt.yscale ("log")
plt.xscale("log")