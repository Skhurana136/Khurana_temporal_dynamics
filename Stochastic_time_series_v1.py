# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 16:02:44 2019

@author: khurana
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
import seaborn as sns
import scipy
from gstools import SRF, Exponential

directory = r"C:\Users\khurana\Documents\AquaDiva\DatafromSAgroup\TLUG_GW\/"
file = "20141113_Unstrut_GWM_Menge_Auswahl_neu_sortiert.csv"
csv.register_dialect(
    "myDialect", delimiter=";", quoting=csv.QUOTE_ALL, skipinitialspace=True
)
site = [
    "Großengottern",
    "Grüningen",
    "GWM Frohndorf (1/2004)",
    "Hy Stotternheim 1/1977",
    "Nägelstedt (E13/1975)",
]
site = ["GWM Frohndorf (1/2004)"]
for s in site:
    values = []
    with open(directory + file, "r") as f:
        reader = csv.reader(f, dialect="myDialect")
        for row in reader:
            if row[2] == s:
                values.append(row[3:6])
    C = []
    A = []
    for line in range(len(values) - 1):
        C.append(float(values[line][2].replace(",", ".")))
        A.append(values[line][0])
    #    C = np.log(C)
    #    print("Mean :",np.mean(C), "Variance :", np.var(C), "Median :", scipy.median(C))
    #    Cfiltered = [i for i in C if i>152.5]
    #    C3 = np.log(C3)
    #    print("Mean :",np.mean(C3), "Variance :", np.var(C3), "Median :", scipy.median(C3))
    plt.figure()
    sns.distplot(C, fit=scipy.stats.norm, kde=True)
    #    plt.plot(C)
    plt.title(s)
    plt.xlabel("GW Elevation (m amsl)")
    plt.ylabel("Density")
#    plt.figure()
#    plt.hist(C, bins = 20)
#    plt.title(s)

Cfiltered = [i for i in C if i < 185]
plt.figure()
sns.distplot(Cfiltered, fit=scipy.stats.norm, kde=True)
plt.title(s)
plt.xlabel("GW Elevation (m amsl)")
plt.ylabel("Density")

Csq = np.sqrt(Cfiltered)
plt.figure()
sns.distplot(Csq, fit=scipy.stats.norm, kde=True)
plt.title(s)
plt.xlabel("GW Elevation (m amsl)")
plt.ylabel("Density")

mu = np.mean(C)
median = np.median(C)
V = np.var(C)
S = scipy.stats.skew(C)
K = scipy.stats.kurtosis(C)

muratio = mu / mu
medianratio = median / mu
Vratio = V / mu
Sratio = S / mu
Kratio = K / mu

Newmu = muratio * 0.000123
Newmedian = medianratio * Newmu
NewV = Vratio * Newmu
NewS = Sratio * Newmu
NewK = Kratio * Newmu

mu = 1
Vratio = [1, 2, 5]
V = Vratio * mu
#T = range(365 * 15)
T = range(365 * 1000) #trying a longer time series to check FFT results
counter = 0
figs = [30, 20]
for uservar in V:
    model1 = Exponential(dim=1, var=uservar, len_scale=365 * 2)
    model2 = Exponential(dim=1, var=uservar, len_scale=365)
    #    model3 = Exponential(dim = 1, var = uservar, len_scale = 360/2)
    srf1 = SRF(model1, mean=mu, force_moments=True, seed=counter)
    srf2 = SRF(model2, mean=mu, force_moments=True, seed=counter)
    #    srf3 = SRF(model3, mean = mu, force_moments = True, seed = counter)
    field1 = srf1.unstructured([T])
    field2 = srf2.unstructured([T])
    #    field3 = srf3.unstructured([T])
    field = field1 + field2  # + field3
    # plt.plot(field)
    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=figs)
    #    ax1 = fig.add_subplot(231)
    #    ax2 = fig.add_subplot(234)
    #    ax3 = fig.add_subplot(232)
    #    ax4 = fig.add_subplot(235)
    #    ax5 = fig.add_subplot(233)
    #    ax6 = fig.add_subplot(236)
    #    ax7 = fig.add_subplot(244)
    #    ax8 = fig.add_subplot(248)
    sns.distplot(field1, fit=scipy.stats.norm, kde=True, ax=ax.flat[0])
    ax.flat[0].set_title("Density plot_Superannualcycle")
    ax.flat[0].set_ylabel("Density")
    ax.flat[0].set_xlabel("Elevation_amasl")
    ax.flat[3].plot(field1)
    ax.flat[3].set_title("Time Series")
    ax.flat[3].set_ylabel("Elevation amasl")
    ax.flat[3].set_xlabel("Time")
    sns.distplot(field2, fit=scipy.stats.norm, kde=True, ax=ax.flat[1])
    ax.flat[1].set_title("Density plot_Annualcycle")
    ax.flat[1].set_ylabel("Density")
    ax.flat[1].set_xlabel("Elevation_amasl")
    ax.flat[4].plot(field2)
    ax.flat[4].set_title("Time Series")
    ax.flat[4].set_ylabel("Elevation amasl")
    ax.flat[4].set_xlabel("Time")
    #    sns.distplot(field3, fit = scipy.stats.norm, kde = True, ax = ax5)
    #    ax5.set_title("Density plot_Biannualcycle")
    #    ax5.set_ylabel("Density")
    #    ax5.set_xlabel("Elevation_amasl")
    #    ax6.plot(field3)
    #    ax6.set_title("Time Series")
    #    ax6.set_ylabel("Elevation amasl")
    #    ax6.set_xlabel("Time")
    sns.distplot(field, fit=scipy.stats.norm, kde=True, ax=ax.flat[2])
    ax.flat[2].set_title("Density plot_Superimposedcycle")
    ax.flat[2].set_ylabel("Density")
    ax.flat[2].set_xlabel("Elevation_amasl")
    ax.flat[5].plot(field)
    ax.flat[5].set_title("Time Series")
    ax.flat[5].set_ylabel("Elevation amasl")
    ax.flat[5].set_xlabel("Time")
    fig.suptitle("Variance :" + str(uservar))
    plt.savefig("timeseries" + str(uservar) + ".jpg")
    #    fname = "try1"+str(uservar)+".csv"
    #    csvfile= open(fname, "w")
    #    writer = csv.writer(csvfile, delimiter=' ', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
    #    writer.writerow(["Seed", "Variance", "Mean_1","Variance_1","Ratio_1","Mean_2","Variance_2","Ratio_2","Mean_3","Variance_3","Ratio_3","Mean_4","Variance_4","Ratio_4"])
    #    writer.writerow([counter, uservar, np.mean(field1), np.var(field1), np.var(field1)/np.mean(field1), np.mean(field2), np.var(field2), np.var(field2)/np.mean(field2),np.mean(field3), np.var(field3), np.var(field3)/np.mean(field3), np.mean(field), np.var(field), np.var(field)/np.mean(field)])
    #    csvfile.close()
    field = field + abs(np.min(field))
    print(np.mean(field1), np.var(field1), np.var(field1) / np.mean(field1))
    print(np.mean(field2), np.var(field2), np.var(field2) / np.mean(field2))
    print(np.mean(field), np.var(field), np.var(field) / np.mean(field))
    fname = "Exponential_" + str(uservar) + "_1000years.csv"
    csvfile = open(fname, "w")
    writer = csv.writer(
        csvfile,
        delimiter=",",
        quotechar="'",
        quoting=csv.QUOTE_MINIMAL,
        lineterminator="\n",
    )
    writer.writerow(["time", "Ratio", "Originalvalue"])
    for x in range(len(T)):
        writer.writerow([T[x], field[x] / np.mean(field), field[x]])
    csvfile.close()
    counter = counter + 1
