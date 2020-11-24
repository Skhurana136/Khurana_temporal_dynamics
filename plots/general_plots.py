# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 11:17:20 2020

@author: khurana
"""
import seaborn as sns
import matplotlib.pyplot as plt
from data_reader.data_processing import tracerstudies

def plot_tracer(filename):
    combined_tracer = tracerstudies(filename)
    combined_tracer["%fraction"] = combined_tracer["fraction"] * 100
    combined_tracer["%fraction_withslow"] = combined_tracer["fraction_withslow"] * 100
    sns.set(rc={"figure.figsize": (7, 4)})
    sns.set_style("whitegrid")
    sns.boxplot(
        x="Xlabels",
        y="%fraction",
        hue="Regime",
        data=combined_tracer,
        hue_order=["Slow", "Medium", "Fast"],
        palette=["coral", "mediumseagreen", "steelblue"],
    )
    plt.xlabel("Variance:Anisotropy")
    plt.ylabel("% of homogeneous scenario")
    plt.title("Time taken for tracer breakthrough")

    return plt


def scatter (data, Xvariable, Xlabel, Yvariable, Ylabel, scaling):
    plt.scatter(Xvariable, Yvariable, s = 100, linewidths = 1, alpha = .7, edgecolor = 'k', data = data)
    if scaling == "Logx":
        plt.xscale("log")
    elif scaling == "Logy":
        plt.yscale("log")
    elif scaling == "True":
        plt.xscale("log")
        plt.yscale("log")
        plt.axhline(10, linestyle = '--', color = 'gray')
        plt.axhline(100, linestyle = '--', color = 'gray')
    plt.ylabel (Ylabel)
    plt.xlabel (Xlabel)
    
    return plt