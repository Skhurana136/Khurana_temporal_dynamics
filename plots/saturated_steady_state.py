# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:57:12 2020

@author: khurana
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import pandas as pd
from matplotlib.colors import LogNorm
from analyses.saturated_steady_state import calcoxiccells
from analyses.saturated_transient import (
    calcconcmasstime,
    biomasstimefunc,
    calcconcmasstimeX,
)
from data_reader.data_processing import tracerstudies
import matplotlib.patches as mpatches
import matplotlib as mpl

FlowRegimes = ["Slow", "Medium", "Fast"]
markerseries = ["d", "^", "o"]
defaultcolors = ["indianred", "g", "steelblue"]
Redscmap = mpl.cm.Reds(np.linspace(0, 1, 30))
Greenscmap = mpl.cm.Greens(np.linspace(0, 1, 30))
Bluescmap = mpl.cm.Blues(np.linspace(0, 1, 30))
Redscmap = mpl.colors.ListedColormap(Redscmap[10:, :-1])
Greenscmap = mpl.colors.ListedColormap(Greenscmap[10:, :-1])
Bluescmap = mpl.colors.ListedColormap(Bluescmap[10:, :-1])
colseries = [Redscmap, Greenscmap, Bluescmap]
red_patch = mpatches.Patch(color="indianred", label="Slow flow")
green_patch = mpatches.Patch(color="g", label="Medium flow")
blue_patch = mpatches.Patch(color="steelblue", label="Fast flow")
patchlist = [red_patch, green_patch, blue_patch]

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


def norm_mf_1col(data, variablenames):
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = ["DOC", "DO", "Nitrogen", "TOC"]
    colseries = ["indianred", "g", "steelblue"]
    nrows = len(Chems)
    fig, axes = plt.subplots(nrows=nrows, figsize=[4, 10])
    # plt.suptitle("Change in removal of carbon and nitrogen with respect to homogeneous scenario at steady state", fontsize = 20)
    for k in Chems:
        dfc = data[data["Chem"] == k]
        colidx1 = Chems.index(k)
        for i in Regimes:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp["Regime"] == i]
            print(i)
            axes.flat[colidx1].scatter(
                "fraction",
                "del2massflux",
                color=colseries[Regimes.index(i)],
                data=dfcr,
                label=i + " flow",
            )
            axes.flat[colidx1].set_ylabel(k, fontsize=15)
            axes.flat[colidx1].tick_params(axis="y", labelsize=15)
            if Chems.index(k) != len(Chems) - 1:
                axes.flat[colidx1].set_xlabel("")
                axes.flat[colidx1].set_xticklabels([])
            else:
                axes.flat[colidx1].tick_params(axis="x", labelsize=15)
                axes.flat[colidx1].set_xlabel(
                    "Fraction of breakthrough time in base case", fontsize=15
                )
    plt.legend(loc="best", fontsize=12)
    plt.annotate(
        "Normalized reduction in mass flux in the domain",
        xy=(0, 2.2),
        xytext=(-80, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )

    return plt


def norm_mf_2x2(data, variablenames):
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = ["DOC", "DO", "Nitrogen", "TOC"]
    colseries = ["indianred", "g", "steelblue"]
    nrows = 2
    ncols = 2
    fig, axes = plt.subplots(nrows=nrows, ncols=2, figsize=[11, 8], sharex=True)
    # plt.suptitle("Change in removal of carbon and nitrogen with respect to homogeneous scenario at steady state", fontsize = 20)
    for k in Chems:
        dfc = data[data["Chem"] == k]
        colidx1 = Chems.index(k)
        for i in Regimes:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp["Regime"] == i]
            print(i)
            axes.flat[colidx1].scatter(
                "fraction",
                "del2massflux",
                color=colseries[Regimes.index(i)],
                data=dfcr,
                label=i + " flow",
            )
            axes.flat[colidx1].set_ylabel(k, fontsize=15)
            axes.flat[colidx1].tick_params(axis="y", labelsize=15)
            axes.flat[colidx1].tick_params(axis="x", labelsize=15)
    #            if(Chems.index(k) != len(Chems)-1):
    #                axes.flat[colidx1].set_xlabel('')
    #                axes.flat[colidx1].set_xticklabels([])
    #            else:
    #                axes.flat[colidx1].tick_params(axis = 'x', labelsize = 15)
    #                axes.flat[colidx1].set_xlabel('Fraction of breakthrough time in base case', fontsize = 15)
    plt.legend(loc="best", fontsize=12)
    plt.annotate(
        "Normalized reduction in mass flux in the domain",
        xy=(-1.2, 1),
        xytext=(-80, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )
    plt.annotate(
        "Fraction of breakthrough time in base case",
        xy=(-0.7, 0),
        xytext=(0, -50),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        fontsize=15,
    )

    return plt


def steadystate_biomass(data, biomassvariablenames):
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = biomassvariablenames
    species = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    position = ["Immobilized", "Mobile"]
    activity = ["Active", "Inactive"]
    colseries = ["indianred", "g", "steelblue"]
    fig, axes = plt.subplots(
        nrows=4, ncols=3, figsize=[len(Chems) * 1.5, len(Chems) * 1.5]
    )
    for k in Chems:
        dfc = data[data["Chem"] == k]
        colidx1 = Chems.index(k)
        for i in Regimes:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp["Regime"] == i]
            print(i)
            axes.flat[colidx1].scatter(
                "fraction",
                "Change_umoles",
                color=colseries[Regimes.index(i)],
                data=dfcr,
                label=i + " flow",
            )
            axes.flat[colidx1].tick_params(axis="y", labelsize=15)
            if Chems.index(k) < len(Chems) - 3:
                axes.flat[colidx1].set_xlabel("")
                axes.flat[colidx1].set_xticklabels([])
            else:
                axes.flat[colidx1].tick_params(axis="x", labelsize=15)
    plt.legend(loc="best", fontsize=12)
    for ax, typsp in zip(axes[0, :], species):
        ax.set_title(typsp, fontsize=15)
    axes[0, -1].annotate(
        position[0],
        xy=(0, 0.5),
        xytext=(300, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )
    axes[1, -1].annotate(
        position[0],
        xy=(0, 0.5),
        xytext=(300, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )
    axes[2, -1].annotate(
        position[1],
        xy=(0, 0.5),
        xytext=(300, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )
    axes[3, -1].annotate(
        position[1],
        xy=(0, 0.5),
        xytext=(300, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )
    axes[0, 0].set_ylabel(activity[0], fontsize=15)
    axes[1, 0].set_ylabel(activity[1], fontsize=15)
    axes[2, 0].set_ylabel(activity[0], fontsize=15)
    axes[3, 0].set_ylabel(activity[1], fontsize=15)
    plt.annotate(
        "Total biomass in the domain (uM C)",
        xy=(0, 2.2),
        xytext=(-800, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )
    plt.annotate(
        "Fraction of breakthrough time in base case",
        xy=(-0.4, -0.3),
        xytext=(-100, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=15,
    )

    return plt


def steadystate_active_biomass(data, biomassvariablenames):
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = biomassvariablenames
    species = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    position = ["Immobile", "Mobile"]
    colseries = ["indianred", "g", "steelblue"]
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=[len(Chems) * 2, len(Chems)])
    for k in Chems:
        dfc = data[data["Chem"] == k]
        colidx1 = Chems.index(k)
        for i in Regimes:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp["Regime"] == i]
            print(i)
            axes.flat[colidx1].scatter(
                "fraction",
                "Change_umoles",
                color=colseries[Regimes.index(i)],
                data=dfcr,
                label=i + " flow",
            )
            axes.flat[colidx1].tick_params(axis="y", labelsize=15)
            if Chems.index(k) < len(Chems) - 3:
                axes.flat[colidx1].set_xlabel("")
                axes.flat[colidx1].set_xticklabels([])
            else:
                axes.flat[colidx1].tick_params(axis="x", labelsize=15)
    plt.legend(loc="best", fontsize=12)
    for ax, typsp in zip(axes[0, :], species):
        ax.set_title(typsp, fontsize=15)
    axes[0, -1].annotate(
        position[0],
        xy=(0, 0.5),
        xytext=(200, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )
    axes[1, -1].annotate(
        position[1],
        xy=(0, 0.5),
        xytext=(200, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )
    plt.annotate(
        "Normalized biomass in the domain",
        xy=(-2.5, 1),
        xytext=(-50, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )
    plt.annotate(
        "Fraction of breakthrough time in base case",
        xy=(-0.4, -0.3),
        xytext=(-50, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=15,
    )

    return plt


def heatmapvelocity(
    directory, fpre, fsuf, trialseries, varianceseries, anisotropyseries
):
    # Heatmaps
    from matplotlib.ticker import FuncFormatter

    fmtvel = lambda x, pos: "{:1.1e}".format(x)
    subtitlesize = 18
    rows = 3
    cols = 4
    size = [14, 12]
    fig, ax = plt.subplots(nrows=rows, ncols=cols, figsize=size)
    #    plt.suptitle("Velocity distribution in heterogeneous scenarions (Variance:Anisotropy)", fontsize = subtitlesize)
    for k in trialseries:
        data = np.load(directory + fpre + str(k) + fsuf + fpre + str(k) + "_df.npy")
        #    title = "Variance "+str(Het[Trial.index(k)])+" : Anisotropy "+str(Anis[Trial.index(k)])
        vel = data[2, -1, :, :] * -1
        sns.heatmap(
            vel,
            square=False,
            cmap="YlGnBu",
            xticklabels=False,
            yticklabels=False,
            ax=ax.flat[trialseries.index(k)],
            cbar_kws={"format": FuncFormatter(fmtvel)},
        )
        #        ax.flat[trialseries.index(k)].set_title("Variance: " + str(varianceseries[trialseries.index(k)]) + "& Anisotropy: "+str(anisotropyseries[trialseries.index(k)]), fontsize = subtitlesize)
        ax.flat[trialseries.index(k)].set_title(
            str(varianceseries[trialseries.index(k)])
            + ":"
            + str(anisotropyseries[trialseries.index(k)]),
            fontsize=subtitlesize,
        )
    return plt


def heatmapconcdist(df, vars, Trial, gvarnames, d, fpre, title):
    # Heatmaps
    from matplotlib.ticker import FuncFormatter

    fmt = lambda x, pos: "{:4.0f}".format(x)
    fmtvel = lambda x, pos: "{:1.1e}".format(x)
    titlesize = 20
    subtitlesize = 18
    vel = df[2, np.shape(df)[1] - 1, :, :] * -1
    if "DOC" in gvarnames:
        fig, ax = plt.subplots(
            nrows=2, ncols=3, figsize=(8, 8), sharex=True, sharey=True
        )
        #        plt.tight_layout()
        plt.suptitle(title, ha="center", va="center", fontsize=titlesize)
        for i in range(len(vars)):
            name = sns.heatmap(
                df[vars[i] - 3, np.shape(df)[1] - 1, :, :],
                square=False,
                cmap="YlGnBu",
                xticklabels=10,
                yticklabels=10,
                ax=ax.flat[i],
                cbar_kws={"format": FuncFormatter(fmt)},
            )
            ax.flat[i].set_title(gvarnames[i], fontsize=subtitlesize)
        name = sns.heatmap(
            vel,
            square=False,
            cmap="YlGnBu",
            xticklabels=False,
            yticklabels=False,
            ax=ax.flat[5],
            cbar_kws={"format": FuncFormatter(fmtvel)},
        )
        ax.flat[5].set_title("Velocity", fontsize=subtitlesize)
        picname = (
            d
            + fpre
            + str(Trial)
            + "_"
            + "heatmap"
            + datetime.now().strftime("%d%m%Y%H%M")
            + "_chem.png"
        )
    elif "Aerobes" in gvarnames:
        fig, ax = plt.subplots(
            nrows=2, ncols=2, figsize=(7, 9), sharex=True, sharey=True
        )
        #        plt.tight_layout()
        plt.suptitle(title, ha="center", va="center", fontsize=titlesize)
        for i in range(len(vars)):
            name = sns.heatmap(
                df[vars[i] - 3, np.shape(df)[1] - 1, :, :],
                square=False,
                cmap="YlGnBu",
                xticklabels=10,
                yticklabels=10,
                ax=ax.flat[i],
                cbar_kws={"format": FuncFormatter(fmt)},
            )
            ax.flat[i].set_title(gvarnames[i], fontsize=subtitlesize)
        name = sns.heatmap(
            vel,
            square=False,
            cmap="YlGnBu",
            xticklabels=False,
            yticklabels=False,
            ax=ax.flat[3],
            cbar_kws={"format": FuncFormatter(fmtvel)},
        )
        ax.flat[3].set_title("Velocity", fontsize=subtitlesize)
    return plt


def heatmapconcdistLOG(df, vars, Trial, gvarnames, d, fpre, title):
    # Heatmaps
    from matplotlib.ticker import FuncFormatter

    fmt = lambda x, pos: "{:4.0f}".format(x)
    fmtvel = lambda x, pos: "{:1.1e}".format(x)
    vel = df[2, np.shape(df)[1] - 1, :, :] * -1
    if "DOC" in gvarnames:
        fig, ax = plt.subplots(
            nrows=2, ncols=3, figsize=(10, 10), sharex=True, sharey=True
        )
        plt.suptitle(title)
        for i in range(len(vars)):
            log_norm = LogNorm(
                vmin=df[vars[i] - 3, np.shape(df)[1] - 1, :, :].min().min(),
                vmax=df[vars[i] - 3, np.shape(df)[1] - 1, :, :].max().max(),
            )
            cbar_ticks = [
                math.pow(10, x)
                for x in range(
                    int(
                        math.floor(
                            math.log10(
                                df[vars[i] - 3, np.shape(df)[1] - 1, :, :].min().min()
                            )
                        )
                    ),
                    1
                    + int(
                        math.ceil(
                            math.log10(
                                df[vars[i] - 3, np.shape(df)[1] - 1, :, :].max().max()
                            )
                        )
                    ),
                )
            ]
            name = sns.heatmap(
                df[vars[i] - 3, np.shape(df)[1] - 1, :, :],
                norm=log_norm,
                square=False,
                cmap="YlGnBu",
                xticklabels=10,
                yticklabels=10,
                ax=ax.flat[i],
                cbar_kws={"format": FuncFormatter(fmt), "ticks": cbar_ticks},
                vmin=df[vars[i] - 3, np.shape(df)[1] - 1, :, :].min().min(),
                vmax=df[vars[i] - 3, np.shape(df)[1] - 1, :, :].min().max().max(),
            )
            ax.flat[i].set_title(gvarnames[i])
        name = sns.heatmap(
            vel,
            square=False,
            norm=LogNorm(vmin=abs(vel).min().min(), vmax=abs(vel).max().max()),
            cmap="YlGnBu",
            xticklabels=False,
            yticklabels=False,
            ax=ax.flat[5],
            cbar_kws={
                "format": FuncFormatter(fmtvel),
                "ticks": [
                    math.pow(10, x)
                    for x in range(
                        int(math.floor(math.log10(np.min(abs(vel))))),
                        1 + int(math.ceil(math.log10(np.max(abs(vel))))),
                    )
                ],
            },
            vmin=abs(vel).min().min(),
            vmax=abs(vel).max().max(),
        )
        ax.flat[5].set_title("Velocity")
        picname = (
            d
            + fpre
            + str(Trial)
            + "_"
            + "heatmap"
            + datetime.now().strftime("%d%m%Y%H%M")
            + "_chem.png"
        )
    elif "Aerobes" in gvarnames:
        fig, ax = plt.subplots(
            nrows=2, ncols=2, figsize=(7, 10), sharex=True, sharey=True
        )
        plt.suptitle(title)
        for i in range(len(vars)):
            log_norm = LogNorm(
                vmin=df[vars[i] - 3, np.shape(df)[1] - 1, :, :].min().min(),
                vmax=df[vars[i] - 3, np.shape(df)[1] - 1, :, :].max().max(),
            )
            cbar_ticks = [
                math.pow(10, x)
                for x in range(
                    int(
                        math.floor(
                            math.log10(
                                df[vars[i] - 3, np.shape(df)[1] - 1, :, :].min().min()
                            )
                        )
                    ),
                    1
                    + int(
                        math.ceil(
                            math.log10(
                                df[vars[i] - 3, np.shape(df)[1] - 1, :, :].max().max()
                            )
                        )
                    ),
                )
            ]
            name = sns.heatmap(
                df[vars[i] - 3, np.shape(df)[1] - 1, :, :],
                norm=log_norm,
                square=False,
                cmap="YlGnBu",
                xticklabels=10,
                yticklabels=10,
                ax=ax.flat[i],
                cbar_kws={"format": FuncFormatter(fmt), "ticks": cbar_ticks},
                vmin=df[vars[i] - 3, np.shape(df)[1] - 1, :, :].min().min(),
                vmax=df[vars[i] - 3, np.shape(df)[1] - 1, :, :].min().max().max(),
            )
            ax.flat[i].set_title(gvarnames[i])
        name = sns.heatmap(
            vel,
            square=False,
            norm=LogNorm(vmin=abs(vel).min().min(), vmax=abs(vel).max().max()),
            cmap="YlGnBu",
            xticklabels=False,
            yticklabels=False,
            ax=ax.flat[3],
            cbar_kws={
                "format": FuncFormatter(fmtvel),
                "ticks": [
                    math.pow(10, x)
                    for x in range(
                        int(math.floor(math.log10(np.min(abs(vel))))),
                        1 + int(math.ceil(math.log10(np.max(abs(vel))))),
                    )
                ],
            },
            vmin=abs(vel).min().min(),
            vmax=abs(vel).max().max(),
        )
        ax.flat[3].set_title("Velocity")
        picname = (
            d
            + fpre
            + str(Trial)
            + "_"
            + "heatmap"
            + datetime.now().strftime("%d%m%Y%H%M")
            + "_biomass.png"
        )
    pictosave = name.get_figure()
    pictosave.savefig(picname)
    return plt


def heatmapratedist(df, vars, ratenames, nrows, ncols, figsize):
    # Heatmaps
    fig, ax = plt.subplots(
        nrows=nrows, ncols=ncols, figsize=figsize, sharex=True, sharey=True
    )
    plt.suptitle("Impact on respiration rate distribution at steady state (uM)")
    for i in range(len(vars)):
        #    plt.title(gvarnames[i], ax = ax.flat[i])
        name = sns.heatmap(
            df[i, :, :],
            square=False,
            cmap="YlGnBu",
            xticklabels=False,
            yticklabels=False,
            ax=ax.flat[i],
        )
        ax.flat[i].set_title(ratenames[i])
    return name


def heatmapratedistLOG(df, vars, ratenames, nrows, ncols, figsize):
    # Heatmaps
    fig, ax = plt.subplots(
        nrows=nrows, ncols=ncols, figsize=figsize, sharex=True, sharey=True
    )
    plt.suptitle("Impact on respiration rate distribution at steady state (uM)")
    for i in range(len(vars)):
        log_normrate = LogNorm(
            vmin=abs(df[i, :, :].min().min()), vmax=abs(df[i, :, :].max().max())
        )
        cbar_ticks = [
            math.pow(10, x)
            for x in range(
                int(math.floor(math.log10(abs(df[i, :, :].min().min())))),
                1 + int(math.ceil(math.log10(abs(df[i, :, :].max().max())))),
            )
        ]
        name = sns.heatmap(
            df[i, :, :],
            square=False,
            cmap="YlGnBu",
            norm=log_normrate,
            cbar_kws={"ticks": cbar_ticks},
            vmin=abs(df[i, :, :].min().min()),
            vmax=abs(df[i, :, :].max().max()),
            xticklabels=False,
            yticklabels=False,
            ax=ax.flat[i],
        )
        ax.flat[i].set_title(ratenames[i])
    return name


def concprofileatss(
    Trial,
    Het,
    Anis,
    gw,
    d,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    vars,
    gvarnames,
    intsce,
    colors,
):
    legendsize = 16
    axissize = 16
    ticksize = 14
    titlesize = 20
    for i in intsce:
        if "DOC" in gvarnames:
            df, massendtime, masstime, conctime, Velocity, head = calcconcmasstime(
                Trial[Trial.index(i)],
                Het[Trial.index(i)],
                Anis[Trial.index(i)],
                gw,
                d,
                fpre,
                fsuf,
                yin,
                yout,
                xleft,
                xright,
                vars,
            )
            #    print (np.mean(df[6, np.shape(df)[1]-1, yin, :]), np.mean(df[6, np.shape(df)[1]-1, yout-1, :]), np.mean(df[6, np.shape(df)[1]-1, yout, :]))
            yindex = list(range(51))
            fig, host = plt.subplots()
            fig.subplots_adjust(top=0.8)

            par1 = host.twiny()
            par2 = host.twiny()

            # Offset the top spine of par2.  The ticks and label have already been
            # placed on the top by twiny above.
            par2.spines["top"].set_position(("axes", 1.2))
            # Having been created by twinx, par2 has its frame off, so the line of its
            # detached spine is invisible.  First, activate the frame but make the patch
            # and spines invisible.
            make_patch_spines_invisible(par2)
            # Second, show the right spine.
            par2.spines["top"].set_visible(True)

            p1, = host.plot(
                conctime[np.shape(conctime)[0] - 1, 0:51, 0],
                yindex,
                label=gvarnames[0],
                color=colors[0],
            )
            p2, = host.plot(
                conctime[np.shape(conctime)[0] - 1, 0:51, 1],
                yindex,
                label=gvarnames[1],
                color=colors[1],
            )
            p3, = par1.plot(
                conctime[np.shape(conctime)[0] - 1, 0:51, 2],
                yindex,
                label=gvarnames[2],
                color=colors[2],
            )
            p4, = par2.plot(
                conctime[np.shape(conctime)[0] - 1, 0:51, 3],
                yindex,
                label=gvarnames[3],
                color=colors[3],
            )

            host.set_ylim(0, 51)
            host.set_xlim(0, 800)
            par1.set_xlim(30, 60)
            par2.set_xlim(50, 260)

            plt.gca().invert_yaxis()

            host.set_ylabel("Y (cm)", fontsize=axissize)
            host.set_xlabel("DOC, DO (uM)", fontsize=axissize)
            par1.set_xlabel(str(gvarnames[2]) + " (uM)", fontsize=axissize)
            par2.set_xlabel(str(gvarnames[3]) + " (uM)", fontsize=axissize)

            host.xaxis.label.set_color(p1.get_color())
            par1.xaxis.label.set_color(p3.get_color())
            par2.xaxis.label.set_color(p4.get_color())

            tkw = dict(size=4, width=1.5, labelsize=ticksize)
            host.tick_params(axis="x", colors=p1.get_color(), **tkw)
            par1.tick_params(axis="x", colors=p3.get_color(), **tkw)
            par2.tick_params(axis="x", colors=p4.get_color(), **tkw)
            host.tick_params(axis="y", **tkw)

            lines = [p1, p2, p3, p4]

            host.legend(
                lines,
                [l.get_label() for l in lines],
                fontsize=legendsize,
                loc="lower right",
            )
            plt.title(
                "Variance "
                + str(Het[Trial.index(i)])
                + " : Anisotropy "
                + str(Anis[Trial.index(i)]),
                fontsize=titlesize,
                pad=100,
            )
            picname = (
                d
                + fpre
                + str(Trial[Trial.index(i)])
                + "_"
                + datetime.now().strftime("%d%m%Y%H%M")
                + "_chem_concatss.png"
            )

        elif "Aerobes" in gvarnames:
            df, massendtime, masstime, conctime = biomasstimefunc(
                Trial[Trial.index(i)],
                Het[Trial.index(i)],
                Anis[Trial.index(i)],
                gw,
                d,
                fpre,
                fsuf,
                yin,
                yout,
                xleft,
                xright,
                vars,
            )
            yindex = list(range(51))
            fig, host = plt.subplots()
            fig.subplots_adjust(top=0.8)

            par1 = host.twiny()
            par2 = host.twiny()

            # Offset the right spine of par2.  The ticks and label have already been
            # placed on the right by twinx above.
            par2.spines["top"].set_position(("axes", 1.2))
            # Having been created by twinx, par2 has its frame off, so the line of its
            # detached spine is invisible.  First, activate the frame but make the patch
            # and spines invisible.
            make_patch_spines_invisible(par2)
            # Second, show the right spine.
            par2.spines["top"].set_visible(True)

            p1, = host.plot(
                conctime[np.shape(conctime)[0] - 1, 0:51, 0],
                yindex,
                label=gvarnames[0],
                color=colors[0],
            )
            p2, = par1.plot(
                conctime[np.shape(conctime)[0] - 1, 0:51, 1],
                yindex,
                label=gvarnames[1],
                color=colors[2],
            )
            p3, = par2.plot(
                conctime[np.shape(conctime)[0] - 1, 0:51, 2],
                yindex,
                label=gvarnames[2],
                color=colors[3],
            )

            host.set_ylim(0, 51)
            host.set_xlim(0, 500)
            par1.set_xlim(0, 30)
            par2.set_xlim(0, 60)

            plt.gca().invert_yaxis()

            host.set_ylabel("Y (cm)", fontsize=axissize)
            host.set_xlabel(str(gvarnames[0]) + " (uM)", fontsize=axissize)
            par1.set_xlabel(str(gvarnames[1]) + " (uM)", fontsize=axissize)
            par2.set_xlabel(str(gvarnames[2]) + " (uM)", fontsize=axissize)

            host.xaxis.label.set_color(p1.get_color())
            par1.xaxis.label.set_color(p2.get_color())
            par2.xaxis.label.set_color(p3.get_color())

            tkw = dict(size=4, width=1.5, labelsize=ticksize)
            host.tick_params(axis="x", colors=p1.get_color(), **tkw)
            par1.tick_params(axis="x", colors=p2.get_color(), **tkw)
            par2.tick_params(axis="x", colors=p3.get_color(), **tkw)
            host.tick_params(axis="y", **tkw)

            lines = [p1, p2, p3]

            host.legend(
                lines,
                [l.get_label() for l in lines],
                fontsize=legendsize,
                loc="lower right",
            )
            plt.title(
                "Variance "
                + str(Het[Trial.index(i)])
                + " : Anisotropy "
                + str(Anis[Trial.index(i)]),
                fontsize=titlesize,
                pad=100,
            )
            picname = (
                d
                + fpre
                + str(Trial[Trial.index(i)])
                + "_"
                + datetime.now().strftime("%d%m%Y%H%M")
                + "_biomass_concatss.png"
            )
        plt.savefig(picname, bbox_inches="tight")
    return None


def concprofilewithH(
    Trial,
    Het,
    Anis,
    gw,
    d,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    cvars,
    gvarnames,
    intsce,
    colors,
    i,
    averagingtype,
):
    legendsize = 16
    axissize = 16
    ticksize = 14
    titlesize = 20
    if ("DOC" in gvarnames) & (averagingtype == "flux"):
        df, massendtime, masstime, conctime, Velocity, head = calcconcmasstime(
            Trial[Trial.index(i)],
            Het[Trial.index(i)],
            Anis[Trial.index(i)],
            gw,
            d,
            fpre,
            fsuf,
            yin,
            yout,
            xleft,
            xright,
            cvars,
            gvarnames,
        )
        dfh, massendtimeh, masstimeh, conctimeh, Velocityh, headh = calcconcmasstime(
            Trial[Trial.index("H")],
            Het[Trial.index("H")],
            Anis[Trial.index("H")],
            gw,
            d,
            fpre,
            fsuf,
            yin,
            yout,
            xleft,
            xright,
            cvars,
            gvarnames,
        )
        yindex = list(range(51))
        fig, host = plt.subplots()
        fig.subplots_adjust(top=0.8)

        par1 = host.twiny()
        par2 = host.twiny()

        # Offset the top spine of par2.  The ticks and label have already been
        # placed on the top by twiny above.
        par2.spines["top"].set_position(("axes", 1.2))
        # Having been created by twinx, par2 has its frame off, so the line of its
        # detached spine is invisible.  First, activate the frame but make the patch
        # and spines invisible.
        make_patch_spines_invisible(par2)
        # Second, show the right spine.
        par2.spines["top"].set_visible(True)

        p1, = host.plot(
            conctimeh[-1, 0:51, 0],
            yindex,
            label=gvarnames[0],
            color=colors[0],
            linestyle="-.",
        )
        p1, = host.plot(
            conctime[-1, 0:51, 0],
            yindex,
            label=gvarnames[0],
            color=colors[0],
            linestyle="-",
        )
        p2, = host.plot(
            conctimeh[-1, 0:51, 1],
            yindex,
            label=gvarnames[1],
            color=colors[1],
            linestyle="-.",
        )
        p2, = host.plot(
            conctime[-1, 0:51, 1],
            yindex,
            label=gvarnames[1],
            color=colors[1],
            linestyle="-",
        )
        p3, = par1.plot(
            conctimeh[-1, 0:51, 2],
            yindex,
            label=gvarnames[2],
            color=colors[2],
            linestyle="-.",
        )
        p3, = par1.plot(
            conctime[-1, 0:51, 2],
            yindex,
            label=gvarnames[2],
            color=colors[2],
            linestyle="-",
        )
        p4, = par2.plot(
            conctimeh[-1, 0:51, 3],
            yindex,
            label=gvarnames[3],
            color=colors[3],
            linestyle="-.",
        )
        p4, = par2.plot(
            conctime[-1, 0:51, 3],
            yindex,
            label=gvarnames[3],
            color=colors[3],
            linestyle="-",
        )

        host.set_ylim(0, 51)
        host.set_xlim(0, 800)
        par1.set_xlim(30, 60)
        par2.set_xlim(50, 260)

        host.set_ylabel("Y (cm)", fontsize=axissize)
        host.set_xlabel("DOC, DO (uM)", fontsize=axissize)
        par1.set_xlabel(str(gvarnames[2]) + " (uM)", fontsize=axissize)
        par2.set_xlabel(str(gvarnames[3]) + " (uM)", fontsize=axissize)

        host.xaxis.label.set_color(p1.get_color())
        par1.xaxis.label.set_color(p3.get_color())
        par2.xaxis.label.set_color(p4.get_color())

        tkw = dict(size=4, width=1.5, labelsize=ticksize)
        host.tick_params(axis="x", colors=p1.get_color(), **tkw)
        par1.tick_params(axis="x", colors=p3.get_color(), **tkw)
        par2.tick_params(axis="x", colors=p4.get_color(), **tkw)
        host.tick_params(axis="y", **tkw)

        plt.gca().invert_yaxis()

        lines = [p1, p2, p3, p4]

        host.legend(
            lines,
            [l.get_label() for l in lines],
            bbox_to_anchor=(1.9, 0.5),
            fontsize=legendsize,
            loc="center",
        )
        plt.title(
            "Variance "
            + str(Het[Trial.index(i)])
            + " : Anisotropy "
            + str(Anis[Trial.index(i)]),
            fontsize=titlesize,
            pad=100,
        )
        plt.tight_layout()

    elif ("Aerobes" in gvarnames) & (averagingtype == "spatial"):
        df, massendtime, ma, bioconctime = biomasstimefunc(
            Trial[Trial.index(i)],
            Het[Trial.index(i)],
            Anis[Trial.index(i)],
            gw,
            d,
            fpre,
            fsuf,
            yin,
            yout,
            xleft,
            xright,
            cvars,
        )
        dfh, massendtimeh, mh, bioconctimeh = biomasstimefunc(
            "H",
            Het[Trial.index("H")],
            Anis[Trial.index("H")],
            gw,
            d,
            fpre,
            fsuf,
            yin,
            yout,
            xleft,
            xright,
            cvars,
        )
        yindex = list(range(51))
        fig, host = plt.subplots()
        fig.subplots_adjust(top=0.8)

        par1 = host.twiny()
        par2 = host.twiny()

        # Offset the right spine of par2.  The ticks and label have already been
        # placed on the right by twinx above.
        par2.spines["top"].set_position(("axes", 1.2))
        # Having been created by twinx, par2 has its frame off, so the line of its
        # detached spine is invisible.  First, activate the frame but make the patch
        # and spines invisible.
        make_patch_spines_invisible(par2)
        # Second, show the right spine.
        par2.spines["top"].set_visible(True)

        p1, = host.plot(
            bioconctimeh[-1, :, 0],
            yindex,
            label=gvarnames[0],
            color=colors[0],
            linestyle="-.",
        )
        p1, = host.plot(
            bioconctime[-1, :, 0],
            yindex,
            label=gvarnames[0],
            color=colors[0],
            linestyle="-",
        )
        p2, = par1.plot(
            bioconctimeh[-1, :, 1],
            yindex,
            label=gvarnames[1],
            color=colors[2],
            linestyle="-.",
        )
        p2, = par1.plot(
            bioconctime[-1, :, 1],
            yindex,
            label=gvarnames[1],
            color=colors[2],
            linestyle="-",
        )
        p3, = par2.plot(
            bioconctimeh[-1, :, 2],
            yindex,
            label=gvarnames[2],
            color=colors[3],
            linestyle="-.",
        )
        p3, = par2.plot(
            bioconctime[-1, :, 2],
            yindex,
            label=gvarnames[2],
            color=colors[3],
            linestyle="-",
        )

        host.set_ylim(0, 51)
        host.set_xlim(0, 500)
        par1.set_xlim(0, 30)
        par2.set_xlim(0, 200)

        lblw = dict(fontsize=axissize)
        host.set_ylabel("Y (cm)", **lblw)
        host.set_xlabel(str(gvarnames[0]) + " (uM)", **lblw)
        par1.set_xlabel(str(gvarnames[1]) + " (uM)", **lblw)
        par2.set_xlabel(str(gvarnames[2]) + " (uM)", **lblw)

        host.xaxis.label.set_color(p1.get_color())
        par1.xaxis.label.set_color(p2.get_color())
        par2.xaxis.label.set_color(p3.get_color())

        tkw = dict(size=4, width=1.5, labelsize=ticksize)
        host.tick_params(axis="x", colors=p1.get_color(), **tkw)
        par1.tick_params(axis="x", colors=p2.get_color(), **tkw)
        par2.tick_params(axis="x", colors=p3.get_color(), **tkw)
        host.tick_params(axis="y", **tkw)

        plt.gca().invert_yaxis()

        lines = [p1, p2, p3]

        host.legend(
            lines,
            [l.get_label() for l in lines],
            bbox_to_anchor=(1.9, 0.5),
            fontsize=legendsize,
            loc="center",
        )
        plt.title(
            "Variance "
            + str(Het[Trial.index(i)])
            + " : Anisotropy "
            + str(Anis[Trial.index(i)]),
            fontsize=titlesize,
            pad=100,
        )
        plt.tight_layout()

    elif ("Aerobes" in gvarnames) & (averagingtype == "flux"):
        df, massendtime, masstime, bioconctime, Velocity, head = calcconcmasstime(
            Trial[Trial.index(i)],
            Het[Trial.index(i)],
            Anis[Trial.index(i)],
            gw,
            d,
            fpre,
            fsuf,
            yin,
            yout,
            xleft,
            xright,
            cvars,
            gvarnames,
        )
        dfh, massendtimeh, masstimeh, bioconctimeh, Velocityh, headh = calcconcmasstime(
            "H",
            Het[Trial.index("H")],
            Anis[Trial.index("H")],
            gw,
            d,
            fpre,
            fsuf,
            yin,
            yout,
            xleft,
            xright,
            cvars,
            gvarnames,
        )
        yindex = list(range(51))
        fig, host = plt.subplots()
        fig.subplots_adjust(top=0.8)

        par1 = host.twiny()
        par2 = host.twiny()

        # Offset the right spine of par2.  The ticks and label have already been
        # placed on the right by twinx above.
        par2.spines["top"].set_position(("axes", 1.2))
        # Having been created by twinx, par2 has its frame off, so the line of its
        # detached spine is invisible.  First, activate the frame but make the patch
        # and spines invisible.
        make_patch_spines_invisible(par2)
        # Second, show the right spine.
        par2.spines["top"].set_visible(True)

        p1, = host.plot(
            bioconctimeh[-1, :, 0],
            yindex,
            label=gvarnames[0],
            color=colors[0],
            linestyle="-.",
        )
        p1, = host.plot(
            bioconctime[-1, :, 0],
            yindex,
            label=gvarnames[0],
            color=colors[0],
            linestyle="-",
        )
        p2, = par1.plot(
            bioconctimeh[-1, :, 1],
            yindex,
            label=gvarnames[1],
            color=colors[2],
            linestyle="-.",
        )
        p2, = par1.plot(
            bioconctime[-1, :, 1],
            yindex,
            label=gvarnames[1],
            color=colors[2],
            linestyle="-",
        )
        p3, = par2.plot(
            bioconctimeh[-1, :, 2],
            yindex,
            label=gvarnames[2],
            color=colors[3],
            linestyle="-.",
        )
        p3, = par2.plot(
            bioconctime[-1, :, 2],
            yindex,
            label=gvarnames[2],
            color=colors[3],
            linestyle="-",
        )

        host.set_ylim(0, 51)
        host.set_xlim(0, 100)
        par1.set_xlim(0, 5)
        par2.set_xlim(0, 50)

        lblw = dict(fontsize=axissize)
        host.set_ylabel("Y (cm)", **lblw)
        host.set_xlabel(str(gvarnames[0]) + " (uM)", **lblw)
        par1.set_xlabel(str(gvarnames[1]) + " (uM)", **lblw)
        par2.set_xlabel(str(gvarnames[2]) + " (uM)", **lblw)

        host.xaxis.label.set_color(p1.get_color())
        par1.xaxis.label.set_color(p2.get_color())
        par2.xaxis.label.set_color(p3.get_color())

        tkw = dict(size=4, width=1.5, labelsize=ticksize)
        host.tick_params(axis="x", colors=p1.get_color(), **tkw)
        par1.tick_params(axis="x", colors=p2.get_color(), **tkw)
        par2.tick_params(axis="x", colors=p3.get_color(), **tkw)
        host.tick_params(axis="y", **tkw)

        plt.gca().invert_yaxis()

        lines = [p1, p2, p3]

        host.legend(
            lines,
            [l.get_label() for l in lines],
            bbox_to_anchor=(1.9, 0.5),
            fontsize=legendsize,
            loc="center",
        )
        plt.title(
            "Variance "
            + str(Het[Trial.index(i)])
            + " : Anisotropy "
            + str(Anis[Trial.index(i)]),
            fontsize=titlesize,
            pad=100,
        )
        plt.tight_layout()

    return plt


def concprofileatssX(
    Trial,
    Het,
    Anis,
    gw,
    d,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    vars,
    gvarnames,
    intsce,
    colors,
):
    for i in intsce:
        if "DOC" in gvarnames:
            df, massendtime, masstime, conctime, Velocity, head = calcconcmasstimeX(
                Trial[Trial.index(i)],
                Het[Trial.index(i)],
                Anis[Trial.index(i)],
                gw,
                d,
                fpre,
                fsuf,
                yin,
                yout,
                xleft,
                xright,
                vars,
            )
            #    print (np.mean(df[6, np.shape(df)[1]-1, yin, :]), np.mean(df[6, np.shape(df)[1]-1, yout-1, :]), np.mean(df[6, np.shape(df)[1]-1, yout, :]))
            fig, host = plt.subplots()
            fig.subplots_adjust(right=0.75)

            par1 = host.twinx()
            par2 = host.twinx()

            # Offset the right spine of par2.  The ticks and label have already been
            # placed on the right by twinx above.
            par2.spines["right"].set_position(("axes", 1.2))
            # Having been created by twinx, par2 has its frame off, so the line of its
            # detached spine is invisible.  First, activate the frame but make the patch
            # and spines invisible.
            make_patch_spines_invisible(par2)
            # Second, show the right spine.
            par2.spines["right"].set_visible(True)

            p1, = host.plot(
                conctime[np.shape(conctime)[0] - 1, 0:31, 0],
                label=gvarnames[0],
                color=colors[0],
            )
            p2, = host.plot(
                conctime[np.shape(conctime)[0] - 1, 0:31, 1],
                label=gvarnames[1],
                color=colors[1],
            )
            p3, = par1.plot(
                conctime[np.shape(conctime)[0] - 1, 0:31, 2],
                label=gvarnames[2],
                color=colors[2],
            )
            p4, = par2.plot(
                conctime[np.shape(conctime)[0] - 1, 0:31, 3],
                label=gvarnames[3],
                color=colors[3],
            )

            host.set_xlim(0, 31)
            host.set_ylim(0, 800)
            par1.set_ylim(30, 60)
            par2.set_ylim(50, 260)

            host.set_xlabel("X (cm)")
            host.set_ylabel("DOC, DO (uM)")
            par1.set_ylabel(str(gvarnames[2]) + " (uM)")
            par2.set_ylabel(str(gvarnames[3]) + " (uM)")

            host.yaxis.label.set_color(p1.get_color())
            par1.yaxis.label.set_color(p3.get_color())
            par2.yaxis.label.set_color(p4.get_color())

            tkw = dict(size=4, width=1.5)
            host.tick_params(axis="y", colors=p1.get_color(), **tkw)
            par1.tick_params(axis="y", colors=p3.get_color(), **tkw)
            par2.tick_params(axis="y", colors=p4.get_color(), **tkw)
            host.tick_params(axis="x", **tkw)

            lines = [p1, p2, p3, p4]

            host.legend(lines, [l.get_label() for l in lines])
            plt.title(Trial[Trial.index(i)])
            picname = (
                d
                + fpre
                + str(Trial[Trial.index(i)])
                + "_"
                + datetime.now().strftime("%d%m%Y%H%M")
                + "_chem_concatss_X.png"
            )

        elif "Aerobes" in gvarnames:
            df, massendtime, masstime, conctime = biomasstimefunc(
                Trial[Trial.index(i)],
                Het[Trial.index(i)],
                Anis[Trial.index(i)],
                gw,
                d,
                fpre,
                fsuf,
                yin,
                yout,
                xleft,
                xright,
                vars,
            )
            fig, host = plt.subplots()
            fig.subplots_adjust(right=0.75)

            par1 = host.twinx()
            par2 = host.twinx()

            # Offset the right spine of par2.  The ticks and label have already been
            # placed on the right by twinx above.
            par2.spines["right"].set_position(("axes", 1.2))
            # Having been created by twinx, par2 has its frame off, so the line of its
            # detached spine is invisible.  First, activate the frame but make the patch
            # and spines invisible.
            make_patch_spines_invisible(par2)
            # Second, show the right spine.
            par2.spines["right"].set_visible(True)

            p1, = host.plot(
                conctime[np.shape(conctime)[0] - 1, 0:31, 0],
                label=gvarnames[0],
                color=colors[0],
            )
            p2, = par1.plot(
                conctime[np.shape(conctime)[0] - 1, 0:31, 1],
                label=gvarnames[1],
                color=colors[2],
            )
            p3, = par2.plot(
                conctime[np.shape(conctime)[0] - 1, 0:31, 2],
                label=gvarnames[2],
                color=colors[3],
            )

            host.set_xlim(0, 31)
            host.set_ylim(0, 500)
            par1.set_ylim(0, 30)
            par2.set_ylim(0, 60)

            host.set_xlabel("X (cm)")
            host.set_ylabel(str(gvarnames[0]) + " (uM)")
            par1.set_ylabel(str(gvarnames[1]) + " (uM)")
            par2.set_ylabel(str(gvarnames[2]) + " (uM)")

            host.yaxis.label.set_color(p1.get_color())
            par1.yaxis.label.set_color(p2.get_color())
            par2.yaxis.label.set_color(p3.get_color())

            tkw = dict(size=4, width=1.5)
            host.tick_params(axis="y", colors=p1.get_color(), **tkw)
            par1.tick_params(axis="y", colors=p2.get_color(), **tkw)
            par2.tick_params(axis="y", colors=p3.get_color(), **tkw)
            host.tick_params(axis="x", **tkw)

            lines = [p1, p2, p3]

            host.legend(lines, [l.get_label() for l in lines])
            plt.title(Trial[Trial.index(i)])
            picname = (
                d
                + fpre
                + str(Trial[Trial.index(i)])
                + "_"
                + datetime.now().strftime("%d%m%Y%H%M")
                + "_biomass_concatss_X.png"
            )
        plt.savefig(picname)
    return None


def plotconcallatss(
    Trial,
    Het,
    Anis,
    gw,
    d,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    vars,
    gvarnames,
    intsce,
    colors,
    nrows,
    ncols,
    figsize,
):
    col = 0
    figbig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    for i in intsce:
        df, massendtime, masstime, conctime, Velocity, head = calcconcmasstime(
            Trial[Trial.index(i)],
            Het[Trial.index(i)],
            Anis[Trial.index(i)],
            gw,
            d,
            fpre,
            fsuf,
            yin,
            yout,
            xleft,
            xright,
            vars,
        )
        for k in range(len(vars) - 2):
            ax[int(math.floor(col / ncols))][col % ncols].plot(
                conctime[np.shape(conctime)[0] - 1, 0:51, k],
                label=gvarnames[k],
                color=colors[k],
            )
            ax[int(math.floor(col / ncols))][col % ncols].set_title(
                plt.title(str(Het[Trial.index(i)]) + " : " + str(Anis[Trial.index(i)]))
            )
        col = col + 1
    plt.legend()
    picname = (
        d
        + fpre
        + str(intsce[0])
        + "_"
        + str(intsce[-1])
        + "_"
        + datetime.now().strftime("%d%m%Y%H%M")
        + "_concatss.png"
    )
    plt.savefig(picname)

    return None


def plotconcwithhet(
    Regimes,
    intsce,
    chem,
    colorfamily,
    Trial,
    Het,
    Anis,
    gw,
    d,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    vars,
    gvarnames,
):
    legendsize = 14
    axissize = 14
    ticksize = 12
    titlesize = 15
    yindex = list(range(51))
    fig, axes = plt.subplots(nrows=1, ncols=len(Regimes), figsize=[16, 5])
    for r in Regimes:
        colors = sns.color_palette(colorfamily[Regimes.index(r)], len(intsce))
        d = r"X:/Saturated_flow/diffusion/" + r + "AR_changedkindox/"
        lines = []
        for i in intsce:
            df, massendtime, masstime, conctime, Velocity, head = calcconcmasstime(
                Trial[Trial.index(i)],
                Het[Trial.index(i)],
                Anis[Trial.index(i)],
                gw,
                d,
                fpre,
                fsuf,
                yin,
                yout,
                xleft,
                xright,
                vars,
            )
            p1, = axes.flat[Regimes.index(r)].plot(
                conctime[np.shape(conctime)[0] - 1, 0:51, gvarnames.index(chem)],
                yindex,
                label=str(Het[Trial.index(i)]) + ":" + str(Anis[Trial.index(i)]),
                color=colors[intsce.index(i)],
            )
            lines.append(p1)
            axes.flat[Regimes.index(r)].legend(
                lines, [l.get_label() for l in lines], fontsize=legendsize
            )

        if r == "Equal":
            axes.flat[Regimes.index(r)].set_title("Medium flow", fontsize=titlesize)
        else:
            axes.flat[Regimes.index(r)].set_title(r + " flow", fontsize=titlesize)

        axes.flat[Regimes.index(r)].set_xlim(0, 250)
        axes.flat[Regimes.index(r)].set_ylabel("Y (cm)", fontsize=axissize)
        axes.flat[Regimes.index(r)].set_xlabel("DO (uM)", fontsize=axissize)
        axes.flat[Regimes.index(r)].tick_params(axis="y", labelsize=ticksize)
        axes.flat[Regimes.index(r)].tick_params(axis="x", labelsize=ticksize)
        axes.flat[Regimes.index(r)].invert_yaxis()
    return fig


def plotconcallatsstime(
    Trialhead,
    Trial,
    Het,
    Anis,
    gw,
    d,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    vars,
    gvarnames,
    intsce,
    colors,
    nrows,
    ncols,
    figsize,
):
    col = 0
    figbig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    for i in intsce:
        Trial2 = str(Trialhead) + str(i)
        df, massendtime, masstime, conctime, Velocity, head = calcconcmasstime(
            Trial2[Trial.index(i)],
            Het[Trial.index(i)],
            Anis[Trial.index(i)],
            gw,
            d,
            fpre,
            fsuf,
            yin,
            yout,
            xleft,
            xright,
            vars,
        )
        for k in range(len(vars) - 2):
            ax[int(math.floor(col / ncols))][col % ncols].plot(
                conctime[np.shape(conctime)[0] - 1, 0:51, k],
                label=gvarnames[k],
                color=colors[k],
            )
            ax[int(math.floor(col / ncols))][col % ncols].set_title(
                Trial[Trial.index(i)]
            )
        col = col + 1
    plt.legend()
    picname = (
        d
        + fpre
        + str(intsce[0])
        + "_"
        + str(intsce[-1])
        + "_"
        + datetime.now().strftime("%d%m%Y%H%M")
        + "_concatss.png"
    )
    plt.savefig(picname)

    return None


def boxrestime_flux(dataset1, dataset2, dataset3, chemseries):

    dfall = pd.concat([dataset1, dataset2, dataset3], axis=0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall["Variance"][i]) + ":" + str(dfall["Anisotropy"][i]))

    dfall["Xlabels"] = l
    dfall = dfall.sort_values(by=["Variance", "Anisotropy"])

    bth = tracerstudies()

    dfall2 = pd.merge(
        dfall,
        bth[["Trial", "Regime", "Firsthit", "%ofhomogeneous"]],
        on=["Trial", "Regime"],
    ).rename(columns={"Firsthit": "Residencetime"})
    bins = [-60, -40, -20, 0, 20]
    dfall2["binned"] = pd.cut(dfall2["%ofhomogeneous"].astype(float), bins)

    Regimes = ["Slow", "Medium", "Fast"]
    Chems = chemseries
    colseries = ["Reds", "Greens", "Blues"]
    ncols = len(Regimes)
    nrows = len(Chems)

    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=[15, 10])
    plt.suptitle(
        "Change in removal of carbon and nitrogen with respect to homogeneous scenario at steady state"
    )
    for i in Regimes:
        #    dfall0 = dfall2[dfall2['Time1']==0]
        df = dfall2[dfall2["Regime"] == i]
        col = 0
        for k in Chems:
            dfc = df[df["Chem"] == k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.boxplot(
                x="binned",
                y="del2massflux%",
                palette=colseries[Regimes.index(i)],
                data=dfc,
                ax=axes[colidx1][colidx2],
            )
            axes[colidx1][colidx2].set_xlabel("")
            axes[colidx1][colidx2].set_ylabel("")
            col = col + 1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:, 0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:, 2], Chems):
        ax.annotate(
            row,
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad + 250, 0),
            xycoords="axes fraction",
            textcoords="offset points",
            size="large",
            ha="left",
            va="center",
        )
    for ax in axes[2]:
        ax.set_xlabel("Relative difference in residence time (%)")
    fig.subplots_adjust(left=0.15, top=0.9)

    return fig


def boxrestime_biomass(dataset1, dataset2, dataset3):

    dfall = pd.concat([dataset1, dataset2, dataset3], axis=0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall["Variance_b"][i]) + ":" + str(dfall["Anisotropy_b"][i]))

    dfall["Xlabels"] = l
    dfall = dfall.sort_values(by=["Variance_b", "Anisotropy_b"])

    bth = tracerstudies()

    dfall2 = pd.merge(
        dfall,
        bth[["Trial", "Regime", "Firsthit", "%ofhomogeneous"]],
        on=["Trial", "Regime"],
    ).rename(columns={"Firsthit": "Residencetime"})
    bins = [-60, -40, -20, 0, 20]
    dfall2["binned"] = pd.cut(dfall2["%ofhomogeneous"].astype(float), bins)

    Regimes = ["Slow", "Medium", "Fast"]
    Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    colseries = ["Reds", "Greens", "Blues"]
    ncols = len(Regimes)
    nrows = len(Chems)

    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=[15, 10])
    plt.suptitle(
        "Change in biomass with respect to homogeneous scenario at steady state"
    )
    for i in Regimes:
        #    dfall0 = dfall2[dfall2['Time1']==0]
        df = dfall2[dfall2["Regime"] == i]
        col = 0
        for k in Chems:
            dfc = df[df["Chemb"] == k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.boxplot(
                x="binned",
                y="delbiomass%",
                palette=colseries[Regimes.index(i)],
                data=dfc,
                ax=axes[colidx1][colidx2],
            )
            axes[colidx1][colidx2].set_xlabel("")
            axes[colidx1][colidx2].set_ylabel("")
            col = col + 1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:, 0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:, 2], Chems):
        ax.annotate(
            row,
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad + 250, 0),
            xycoords="axes fraction",
            textcoords="offset points",
            size="large",
            ha="left",
            va="center",
        )
    for ax in axes[2]:
        ax.set_xlabel("Relative difference in residence time (%)")
    fig.subplots_adjust(left=0.15, top=0.9)

    return fig


def scatterrestime_flux(dataset1, dataset2, dataset3, chemseries):
    legendsize = 14
    axissize = 14
    ticksize = 12
    titlesize = 15
    dfall = pd.concat([dataset2, dataset3], axis=0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall["Variance"][i]) + ":" + str(dfall["Anisotropy"][i]))

    dfall["Xlabels"] = l
    dfall = dfall.sort_values(by=["Variance", "Anisotropy"])

    bth = tracerstudies()

    dfall2 = pd.merge(
        dfall,
        bth[["Trial", "Regime", "Firsthit", "%ofhomogeneous"]],
        on=["Trial", "Regime"],
    ).rename(columns={"Firsthit": "Residencetime"})
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = chemseries
    colseries = ["indianred", "g", "steelblue"]
    ncols = len(Regimes)
    nrows = len(Chems)
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=[15, 10])
    plt.suptitle(
        "Change in removal of carbon and nitrogen with respect to homogeneous scenario at steady state",
        fontsize=20,
    )
    for i in Regimes:
        #    dfall0 = dfall2[dfall2['Time1']==0]
        df = dfall2[dfall2["Regime"] == i]
        col = 0
        for k in Chems:
            dfc = df[df["Chem"] == k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            #            pred = dfc["%ofhomogeneous"].to_numpy()
            #            resp = dfc["del2massflux%"].to_numpy()
            #            intercept , slope= np.polyfit(pred, resp,1)
            #            print (i, k, intercept, slope)
            dum = sns.lmplot(
                x="%ofhomogeneous",
                y="del2massflux%",
                data=dfall2,
                palette=colseries[Regimes.index(i)],
                col="Regime",
                row="Chem",
            )
            #            axes[colidx1][colidx2].scatter(pred, resp, c = colseries[Regimes.index(i)])
            #            axes[colidx1][colidx2].plot(pred, slope + intercept*pred, '-', color = colseries[Regimes.index(i)])
            axes[colidx1][colidx2].set_xlabel("")
            axes[colidx1][colidx2].set_ylabel("")
            axes[colidx1][colidx2].set_xticks(
                np.arange(
                    round(min(np.unique(dfall2["%ofhomogeneous"])), -1),
                    round(max(np.unique(dfall2["%ofhomogeneous"])), -1),
                    10.0,
                )
            )
            axes[colidx1][colidx2].set_xticklabels([])
            col = col + 1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col + " flow", ha="center", fontsize=titlesize)
    for ax, row in zip(axes[:, 2], Chems):
        ax.annotate(
            row,
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad + 255, 0),
            xycoords="axes fraction",
            textcoords="offset points",
            size="large",
            ha="left",
            va="center",
            rotation="vertical",
            fontsize=15,
        )
    for ax in axes[-1]:
        ax.set_xticklabels(
            np.arange(
                round(min(np.unique(dfall2["%ofhomogeneous"])), -1),
                round(max(np.unique(dfall2["%ofhomogeneous"])), -1),
                10,
            ),
            size=ticksize,
        )
    plt.figtext(
        0.5,
        0.08,
        "Relative difference in breakthrough time (%)",
        ha="center",
        va="center",
        fontsize=20,
    )
    plt.figtext(
        0.08,
        0.5,
        "Relative difference (%)",
        ha="center",
        va="center",
        rotation="vertical",
        fontsize=20,
    )

    return dfall2, fig


def scatterrestime_biomass(dataset1, dataset2, dataset3):

    dfall = pd.concat([dataset1, dataset2, dataset3], axis=0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall["Variance"][i]) + ":" + str(dfall["Anisotropy"][i]))

    dfall["Xlabels"] = l
    dfall = dfall.sort_values(by=["Variance", "Anisotropy"])

    bth = tracerstudies()

    dfall2 = pd.merge(
        dfall,
        bth[["Trial", "Regime", "Firsthit", "%ofhomogeneous"]],
        on=["Trial", "Regime"],
    ).rename(columns={"Firsthit": "Residencetime"})

    Regimes = ["Slow", "Medium", "Fast"]
    Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    colseries = ["indianred", "g", "steelblue"]
    ncols = len(Regimes)
    nrows = len(Chems)

    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=[15, 10], sharex="col")
    plt.suptitle(
        "Change in biomass with respect to homogeneous scenario at steady state",
        fontsize=20,
    )
    for i in Regimes:
        #    dfall0 = dfall2[dfall2['Time1']==0]
        df = dfall2[dfall2["Regime"] == i]
        col = 0
        for k in Chems:
            dfc = df[df["Chemb"] == k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            pred = dfc["%ofhomogeneous"].to_numpy()
            resp = dfc["delbiomass%"].to_numpy()
            intercept, slope = np.polyfit(pred, resp, 1)
            print(i, k, intercept, slope)
            axes[colidx1][colidx2].scatter(
                pred,
                resp,
                c=colseries[Regimes.index(i)],
                cmap=colseries[Regimes.index(i)],
            )
            axes[colidx1][colidx2].plot(
                pred, slope + intercept * pred, "-", color=colseries[Regimes.index(i)]
            )
            axes[colidx1][colidx2].set_xlabel("")
            axes[colidx1][colidx2].set_ylabel("")
            axes[colidx1][colidx2].set_xticks(
                np.arange(
                    round(min(np.unique(dfall2["%ofhomogeneous"])), -1),
                    round(max(np.unique(dfall2["%ofhomogeneous"])), -1),
                    10.0,
                )
            )
            axes[colidx1][colidx2].set_xticklabels([])
            #            axes[colidx1][colidx2].set_xticklabels([20, 0, -20, -40, -60, -80])
            col = col + 1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    #    for ax in axes[:,0]:
    #        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:, 2], Chems):
        ax.annotate(
            row,
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad + 250, 0),
            xycoords="axes fraction",
            textcoords="offset points",
            ha="left",
            va="center",
            rotation="vertical",
            fontsize=15,
        )
    for ax in axes[-1]:
        #        ax.set_xlabel("Relative difference in breakthrough time (%)")
        ax.set_xticklabels(
            np.arange(
                round(min(np.unique(dfall2["%ofhomogeneous"])), -1),
                round(max(np.unique(dfall2["%ofhomogeneous"])), -1),
                10,
            )
        )
    plt.figtext(
        0.5,
        0.08,
        "Relative difference in breakthrough time (%)",
        ha="center",
        va="center",
        fontsize=20,
    )
    plt.figtext(
        0.08,
        0.5,
        "Relative difference (%)",
        ha="center",
        va="center",
        rotation="vertical",
        fontsize=20,
    )
    fig.subplots_adjust(left=0.15, top=0.9)

    return dfall2, fig


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


def plotoxiccellssdo(
    limit,
    Trial,
    Het,
    Anis,
    gw,
    d,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    vars,
    gvarnames,
):
    axissize = 15
    ticksize = 15
    titlesize = 15
    oxiccells = calcoxiccells(
        limit,
        Trial,
        Het,
        Anis,
        gw,
        d,
        fpre,
        fsuf,
        yin,
        yout,
        xleft,
        xright,
        vars,
        gvarnames,
    )
    nrows = 2
    ncols = 4
    figsize = [16, 10]
    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, figsize=figsize, sharex=True, sharey=True
    )
    count = 0
    for j in Trial:
        df, massendtime, masstime, conctime, Velocity, head = calcconcmasstime(
            j,
            Het[Trial.index(j)],
            Anis[Trial.index(j)],
            gw,
            d,
            fpre,
            fsuf,
            yin,
            yout,
            xleft,
            xright,
            vars,
            gvarnames,
        )
        axes.flat[count].plot(conctime[-1, :, gvarnames.index("DO")], "r-")
        axes.flat[count].set_ylim(0, 260)
        axes.flat[count].tick_params(axis="y", colors="r", labelsize=ticksize)
        axes.flat[count].set_title(
            "Variance: "
            + str(Het[Trial.index(j)])
            + "& Anisotropy: "
            + str(Anis[Trial.index(j)]),
            fontsize=titlesize,
        )
        ax2 = axes.flat[count].twinx()
        ax2.plot((oxiccells[Trial.index(j), :, 0] / 31) * 100, "b-")
        ax2.set_ylim(0, 105)
        ax2.tick_params(axis="y", colors="b", labelsize=ticksize)
        if (count + 1) % ncols == 0:
            ax2.set_ylabel("% of oxic cells (-)", color="b", fontsize=axissize)
        else:
            ax2.set_yticklabels([])  # turning off secondary y axis yticklabels
        count = count + 1
    for ax in axes[:, 0]:
        ax.set_ylabel("DO (uM)", color="r", fontsize=axissize)
    for ax in axes[-1]:
        ax.set_xlabel("Y (cm)", fontsize=axissize)
    return fig


def biomassvscarryingcapacity(
    Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, cvars
):
    maxcap = 500
    yindex = list(range(51))
    cc = np.zeros([51, 1])
    df, massendtime, ma, bioconctime = biomasstimefunc(
        Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, cvars
    )
    c = [
        "firebrick",
        "darkblue",
        "darkgreen",
        "lightcoral",
        "cornflowerblue",
        "mediumseagreen",
    ]
    l = [
        "Active aerobes",
        "Active ammonia oxidizers",
        "Active nitrate reducers",
        "Inactive aerobes",
        "Inactive ammonia oxidizers",
        "Inactive nitrate reducers",
    ]
    bioconctimet = bioconctime[-1, :, :].transpose()
    cc[:, 0] = maxcap
    plt.figure()
    plt.stackplot(yindex, bioconctimet, colors=c, labels=l)
    plt.xlabel("Y (cm)")
    plt.ylabel("Average concentration (uM C)")
    #    plt.plot(yindex, cc, label = "Carrying capacity")
    plt.xlim(0, 50)
    #    plt.xlim(300,510)
    plt.title("Variance: " + str(Het) + " & Anisotropy: " + str(Anis))
    #    plt.gca().invert_yaxis()
    plt.legend(loc="lower right")

    return plt


def scatterbioflux(dataset1, dataset2, dataset3):
    dfall = pd.concat([dataset1, dataset2, dataset3], axis=0, ignore_index=True)
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    colseries = ["Reds", "Greens", "Blues"]
    ncols = len(Regimes)
    nrows = len(Chems)
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=[16, 10])
    plt.suptitle(
        "Change in removal of carbon and nitrogen with respect to change in biomass at steady state"
    )
    for i in Regimes:
        #    dfall0 = dfall[dfall['Time1']==0]
        df = dfall[dfall["Regime"] == i]
        col = 0
        for k in Chems:
            dfc = df[df["Chemb"] == k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.scatterplot(
                x="delbiomass%",
                y="del2massflux%",
                palette=colseries[Regimes.index(i)],
                data=dfc,
                ax=axes[colidx1][colidx2],
            )
            axes[colidx1][colidx2].set_xlabel("")
            axes[colidx1][colidx2].set_ylabel("")
            col = col + 1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:, 0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:, 2], Chems):
        ax.annotate(
            row,
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad + 270, 0),
            xycoords="axes fraction",
            textcoords="offset points",
            size="large",
            ha="left",
            va="center",
        )
    for ax, row in zip(axes[:, 0], Chems):
        if row == "Aerobes":
            ax.annotate(
                "DO",
                xy=(0, 0.5),
                xytext=(-ax.yaxis.labelpad - 70, 0),
                xycoords="axes fraction",
                textcoords="offset points",
                size="large",
                ha="right",
                va="center",
            )
        if row == "Ammonia oxidizers":
            ax.annotate(
                "Ammonium",
                xy=(0, 0.5),
                xytext=(-ax.yaxis.labelpad - 70, 0),
                xycoords="axes fraction",
                textcoords="offset points",
                size="large",
                ha="right",
                va="center",
            )
        if row == "Nitrate reducers":
            ax.annotate(
                "Nitrate",
                xy=(0, 0.5),
                xytext=(-ax.yaxis.labelpad - 70, 0),
                xycoords="axes fraction",
                textcoords="offset points",
                size="large",
                ha="right",
                va="center",
            )
    for ax in axes[2]:
        ax.set_xlabel("Relative difference in biomass (%)")
    fig.subplots_adjust(left=0.15, top=0.9)

    return fig


def boxVAflux(dataset1, dataset2, dataset3, chemseries):
    dfall = pd.concat([dataset1, dataset2, dataset3], axis=0, ignore_index=True)
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = chemseries
    colseries = ["Reds", "Greens", "Blues"]
    ncols = len(Regimes)
    nrows = len(Chems)
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=[16, 10])
    plt.suptitle("Change in removal of carbon and nitrogen at steady state")
    for i in Regimes:
        #    dfall0 = dfall[dfall['Time1']==0]
        df = dfall[dfall["Regime"] == i]
        col = 0
        for k in Chems:
            dfc = df[df["Chem"] == k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.boxplot(
                x="VA",
                y="del2massflux%",
                palette=colseries[Regimes.index(i)],
                data=dfc,
                ax=axes[colidx1][colidx2],
            )
            axes[colidx1][colidx2].set_xlabel("")
            axes[colidx1][colidx2].set_ylabel("")
            col = col + 1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:, 0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:, 2], Chems):
        ax.annotate(
            row,
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad + 270, 0),
            xycoords="axes fraction",
            textcoords="offset points",
            size="large",
            ha="left",
            va="center",
        )
    for ax in axes[2]:
        ax.set_xlabel("Variance x Anisotropy (-)")
    fig.subplots_adjust(left=0.15, top=0.9)

    return fig


def boxVAbio(dataset1, dataset2, dataset3):
    dfall = pd.concat([dataset1, dataset2, dataset3], axis=0, ignore_index=True)
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    colseries = ["Reds", "Greens", "Blues"]
    ncols = len(Regimes)
    nrows = len(Chems)
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=[16, 10])
    plt.suptitle("Change in biomass at steady state")
    for i in Regimes:
        #    dfall0 = dfall[dfall['Time1']==0]
        df = dfall[dfall["Regime"] == i]
        col = 0
        for k in Chems:
            dfc = df[df["Chemb"] == k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.boxplot(
                x="VA",
                y="delbiomass%",
                palette=colseries[Regimes.index(i)],
                data=dfc,
                ax=axes[colidx1][colidx2],
            )
            axes[colidx1][colidx2].set_xlabel("")
            axes[colidx1][colidx2].set_ylabel("")
            col = col + 1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:, 0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:, 2], Chems):
        ax.annotate(
            row,
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad + 270, 0),
            xycoords="axes fraction",
            textcoords="offset points",
            size="large",
            ha="left",
            va="center",
        )
    for ax in axes[2]:
        ax.set_xlabel("Variance x Anisotropy (-)")
    fig.subplots_adjust(left=0.15, top=0.9)

    return fig


def boxV_Aflux(dataset1, dataset2, dataset3, chemseries, imgsize):
    dfall = pd.concat([dataset1, dataset2, dataset3], axis=0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall["Variance"][i]) + ":" + str(dfall["Anisotropy"][i]))

    dfall["Xlabels"] = l
    dfall = dfall.sort_values(by=["Variance", "Anisotropy"])
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = chemseries
    colseries = ["Reds", "Greens", "Blues"]
    ncols = len(Regimes)
    nrows = len(Chems)
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=imgsize)
    plt.suptitle("Change in removal of carbon and nitrogen at steady state")
    for i in Regimes:
        #    dfall0 = dfall[dfall['Time1']==0]
        df = dfall[dfall["Regime"] == i]
        col = 0
        for k in Chems:
            dfc = df[df["Chem"] == k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.boxplot(
                x="Xlabels",
                y="fdelmassflux",
                palette=colseries[Regimes.index(i)],
                data=dfc,
                ax=axes[colidx1][colidx2],
            )
            axes[colidx1][colidx2].set_xlabel("")
            axes[colidx1][colidx2].set_ylabel("")
            col = col + 1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:, 0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:, 2], Chems):
        ax.annotate(
            row,
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad + 400, 0),
            xycoords="axes fraction",
            textcoords="offset points",
            size="large",
            ha="left",
            va="center",
        )
    for ax in axes[-1]:
        ax.set_xlabel("Variance : Anisotropy")
    fig.subplots_adjust(left=0.15, top=0.9)

    return fig


def boxV_Abio(dataset1, dataset2, dataset3, imgsize):
    dfall = pd.concat([dataset1, dataset2, dataset3], axis=0, ignore_index=True)
    #    dfall = pd.concat([equalbc, slowbc, fastbc], axis = 0, ignore_index=True)
    l = []
    for i in range(len(dfall)):
        l.append(str(dfall["Variance_b"][i]) + ":" + str(dfall["Anisotropy_b"][i]))

    dfall["Xlabels"] = l
    dfall = dfall.sort_values(by=["Variance_b", "Anisotropy_b"])
    Regimes = ["Slow", "Medium", "Fast"]
    Chems = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    colseries = ["Reds", "Greens", "Blues"]
    ncols = len(Regimes)
    nrows = len(Chems)
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=imgsize)
    plt.suptitle("Change in total biomass at steady state")
    for i in Regimes:
        #    dfall0 = dfall[dfall['Time1']==0]
        df = dfall[dfall["Regime"] == i]
        col = 0
        for k in Chems:
            dfc = df[df["Chemb"] == k]
            colidx1 = Chems.index(k)
            colidx2 = Regimes.index(i)
            dum = sns.boxplot(
                x="Xlabels",
                y="delbiomass%",
                palette=colseries[Regimes.index(i)],
                data=dfc,
                ax=axes[colidx1][colidx2],
            )
            axes[colidx1][colidx2].set_xlabel("")
            axes[colidx1][colidx2].set_ylabel("")
            col = col + 1
    for ax, col in zip(axes[0], Regimes):
        ax.set_title(col)
    for ax in axes[:, 0]:
        ax.set_ylabel("Relative difference (%)")
    for ax, row in zip(axes[:, 2], Chems):
        ax.annotate(
            row,
            xy=(0, 0.5),
            xytext=(-ax.yaxis.labelpad + 400, 0),
            xycoords="axes fraction",
            textcoords="offset points",
            size="large",
            ha="left",
            va="center",
        )
    for ax in axes[2]:
        ax.set_xlabel("Variance : Anisotropy")
    fig.subplots_adjust(left=0.15, top=0.9)

    return fig

def scatter_chem_regime (data, Xvariable, Xlabel, Yvariable, Ylabel, scaling):
    fig,ax = plt.subplots()
    for sp in list(data['Chem'].unique()):
        for r in list(data['Regime'].unique()):
            dfc = data[(data['Chem']==sp) & (data['Regime']==r)]
            m = markerseries[list(data['Chem'].unique()).index(sp)]
            plt.scatter(Xvariable, Yvariable, s = 100, c = defaultcolors[FlowRegimes.index(r)], linewidths = 1, 
                                                                         alpha = .7, edgecolor = 'k', data = dfc,
                                                                         marker = m, label = sp + " in " + r + " regime")
    if scaling == "Logx":
        plt.xscale("log")
    elif scaling == "Logy":
        plt.yscale("log")
    elif scaling == "True":
        plt.xscale("log")
        plt.yscale("log")
#        plt.axhline(10, linestyle = '--', color = 'gray')
        plt.axhline(100, linestyle = '--', color = 'gray')
    plt.ylabel (Ylabel)
    plt.xlabel (Xlabel)
    plt.legend()
    
    return plt