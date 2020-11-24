# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:57:47 2020

@author: khurana
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import matplotlib.patches as mpatches
import matplotlib as mpl
from analyses.saturated_transient import calcconcmasstime, calcsum_temp
from data_reader.data_processing import localmaxmin

Regimes = ["Slow", "Medium", "Fast"]
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
low_var = mpatches.Patch(color="silver", label="low")
mid_var = mpatches.Patch(color="darkgray", label="medium")
high_var = mpatches.Patch(color="black", label="high")
markerseries = ["d", "^", "o"]


def heatmapconcdist_temp(df, vars, Trial, gvarnames, d, fpre, title, timeidx):
    # Heatmaps
    from matplotlib.ticker import FuncFormatter

    fmt = lambda x, pos: "{:4.0f}".format(x)
    fmtvel = lambda x, pos: "{:1.1e}".format(x)
    titlesize = 20
    subtitlesize = 18
    vel = df[2, timeidx, :, :] * -1
    if "DOC" in gvarnames:
        fig, ax = plt.subplots(
            nrows=2, ncols=3, figsize=(8, 8), sharex=True, sharey=True
        )
        #        plt.tight_layout()
        plt.suptitle(title, ha="center", va="center", fontsize=titlesize)
        for i in range(len(vars)):
            name = sns.heatmap(
                df[vars[i] - 3, timeidx, :, :],
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
    elif "Aerobes" in gvarnames:
        fig, ax = plt.subplots(
            nrows=2, ncols=2, figsize=(7, 9), sharex=True, sharey=True
        )
        #        plt.tight_layout()
        plt.suptitle(title, ha="center", va="center", fontsize=titlesize)
        for i in range(len(vars)):
            name = sns.heatmap(
                df[vars[i] - 3, timeidx, :, :],
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


def plottimeseries(
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
    inttime,
    colorstime,
    nrows,
    ncols,
    figsize,
    chem,
    lylim,
    uylim,
):
    col = 0
    # figbig = plt.figure()
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
        for k in range(len(inttime)):
            ax[int(math.floor(col / ncols))][col % ncols].plot(
                conctime[inttime[k], 0:51, gvarnames.index(chem)],
                label=inttime[k] * 50,
                color=colorstime[k],
            )
            ax[int(math.floor(col / ncols))][col % ncols].set_title(
                Trial[Trial.index(i)]
            )
            ax[int(math.floor(col / ncols))][col % ncols].legend()
            ax[int(math.floor(col / ncols))][col % ncols].set_ylim(lylim, uylim)
        col = col + 1
    plt.legend()
    picname = (
        d
        + fpre
        + str(intsce[0])
        + "_"
        + str(intsce[-1])
        + "_"
        + chem
        + "_"
        + datetime.now().strftime("%d%m%Y%H%M")
        + ".png"
    )
    plt.savefig(picname)
    return plt


def plotoutflowtimeseries(
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
    inttime,
    colorstime,
    nrows,
    ncols,
    figsize,
    chem,
    lylim,
    uylim,
    start,
    end,
):
    col = 0
    # figbig = plt.figure()
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
        ax[int(math.floor(col / ncols))][col % ncols].plot(
            conctime[start:end, yout, gvarnames.index(chem)],
            label="Tracer",
            color="black",
        )
        ax[int(math.floor(col / ncols))][col % ncols].set_title(Trial[Trial.index(i)])
        #        ax[int(math.floor(col/ncols))][i%ncols].set_xlabel("Y")
        #        ax[int(math.floor(col/ncols))][i%ncols].set_ylabel("Concentration (uM)")
        ax[int(math.floor(col / ncols))][col % ncols].set_ylim(lylim, uylim)
        col = col + 1
    plt.legend()
    picname = (
        d + fpre + str(intsce[0]) + "_" + str(intsce[-1]) + "_" + chem + "_time.png"
    )
    plt.savefig(picname)
    return plt


def concprofilewithtime_temp(
    Trial,
    intsce,
    d,
    j,
    gvarnames,
    Het,
    Anis,
    gw,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    vars,
):
    titlesize = 20
    legendsize = 12
    ticksize = 15
    axissize = 15
    xindex = list(range(0, 1095 * 5, 5))
    if "DOC" in gvarnames:
        figbig, axes = plt.subplots(
            nrows=5, ncols=len(intsce), figsize=[30, 12], sharex=True
        )
        plt.suptitle("Concentration profile with time", fontsize=titlesize)
        for i in intsce:
            newd = d + j
            Titles = (
                "Variance:"
                + str(Het[Trial.index(i)])
                + " & Anisotropy:"
                + str(Anis[Trial.index(i)])
            )
            df, massendtime, masstime, conctime, Velocity, head = calcconcmasstime(
                Trial[Trial.index(i)],
                Het[Trial.index(i)],
                Anis[Trial.index(i)],
                gw,
                newd,
                fpre,
                fsuf,
                yin,
                yout,
                xleft,
                xright,
                vars,
                gvarnames,
            )
            axes[0][intsce.index(i)].plot(
                xindex, head, label="Velocity (m/d)", color="black"
            )
            axes[0][intsce.index(i)].set_title(Titles, fontsize=titlesize)
            axes[1][intsce.index(i)].plot(
                xindex, conctime[1:, yout, 0], label=gvarnames[0], color="black"
            )
            axes[2][intsce.index(i)].plot(
                xindex, conctime[1:, yout, 1], label=gvarnames[1], color="red"
            )
            axes[3][intsce.index(i)].plot(
                xindex, conctime[1:, yout, 2], label=gvarnames[2], color="blue"
            )
            axes[4][intsce.index(i)].plot(
                xindex, conctime[1:, yout, 3], label=gvarnames[3], color="green"
            )

            ax = axes[:, intsce.index(i)]
            if intsce.index(i) == 0:
                off = 75
                ax[0].annotate(
                    "Velocity",
                    xy=(0, 0.5),
                    xytext=(-ax[0].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[1].annotate(
                    str(gvarnames[0]),
                    xy=(0, 0.5),
                    xytext=(-ax[1].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[2].annotate(
                    str(gvarnames[1]),
                    xy=(0, 0.5),
                    xytext=(-ax[2].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[3].annotate(
                    str(gvarnames[2]),
                    xy=(0, 0.5),
                    xytext=(-ax[3].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[4].annotate(
                    str(gvarnames[3]),
                    xy=(0, 0.5),
                    xytext=(-ax[4].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
            else:
                ax[0].set_yticklabels([])
                ax[1].set_yticklabels([])
                ax[2].set_yticklabels([])
                ax[3].set_yticklabels([])
                ax[4].set_yticklabels([])

        for ax in axes[:, 0]:
            ax.tick_params(labelsize=ticksize)
        for ax in axes[4]:
            ax.set_xlabel("Time (days)", fontsize=axissize)
            ax.tick_params(labelsize=ticksize)
    elif "Aerobes" in gvarnames:
        figbig, axes = plt.subplots(
            nrows=4, ncols=len(intsce), figsize=[30, 12], sharex=True
        )
        plt.suptitle("Concentration profile with time", fontsize=titlesize)
        for i in intsce:
            newd = d + j
            Titles = (
                "Variance:"
                + str(Het[Trial.index(i)])
                + " & Anisotropy:"
                + str(Anis[Trial.index(i)])
            )
            df, massendtime, masstime, conctime, Velocity, head = calcconcmasstime(
                Trial[Trial.index(i)],
                Het[Trial.index(i)],
                Anis[Trial.index(i)],
                gw,
                newd,
                fpre,
                fsuf,
                yin,
                yout,
                xleft,
                xright,
                vars,
                gvarnames,
            )
            conctime = calcsum_temp(
                Trial[Trial.index(i)],
                Het[Trial.index(i)],
                Anis[Trial.index(i)],
                gw,
                newd,
                fpre,
                fsuf,
                yin,
                yout,
                xleft,
                xright,
                vars,
                gvarnames,
            )
            axes[0][intsce.index(i)].plot(
                xindex, head, label="Velocity (m/d)", color="black"
            )
            axes[0][intsce.index(i)].set_title(Titles, fontsize=titlesize)
            axes[1][intsce.index(i)].plot(
                xindex, conctime[:, 0], label=gvarnames[0], color="black"
            )
            axes[2][intsce.index(i)].plot(
                xindex, conctime[:, 1], label=gvarnames[1], color="red"
            )
            axes[3][intsce.index(i)].plot(
                xindex, conctime[:, 2], label=gvarnames[2], color="blue"
            )

            ax = axes[:, intsce.index(i)]
            if intsce.index(i) == 0:
                off = 75
                ax[0].annotate(
                    "Velocity",
                    xy=(0, 0.5),
                    xytext=(-ax[0].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[1].annotate(
                    str(gvarnames[0]),
                    xy=(0, 0.5),
                    xytext=(-ax[1].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[2].annotate(
                    str(gvarnames[1]),
                    xy=(0, 0.5),
                    xytext=(-ax[2].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[3].annotate(
                    str(gvarnames[2]),
                    xy=(0, 0.5),
                    xytext=(-ax[3].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
            else:
                ax[0].set_yticklabels([])
                ax[1].set_yticklabels([])
                ax[2].set_yticklabels([])
                ax[3].set_yticklabels([])

        for ax in axes[:, 0]:
            ax.tick_params(labelsize=ticksize)
        for ax in axes[3]:
            ax.set_xlabel("Time (days)", fontsize=axissize)
            ax.tick_params(labelsize=ticksize)
    return figbig


def shiftoxic_anoxic_temp(
    chem,
    Trial,
    intsce,
    d,
    j,
    gvarnames,
    Het,
    Anis,
    gw,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    vars,
):
    titlesize = 20
    legendsize = 12
    ticksize = 15
    axissize = 15
    yindex = list(range(51))
    nrows = 2
    ncols = 3
    figsize = [14, 8]
    colors = ["black", "red", "blue", "green", "orange", "grey"]
    times = [2 * 5, 180 * 5, 360 * 5, 500 * 5, 1000 * 5]
    # Plot all the scenarios excluding steady state scenarios
    #    intsce = [76, 73, 80, 84, 44, 63]
    col = 0
    figbig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    plt.suptitle(chem + " profile at different time points", fontsize=titlesize)
    for i in intsce:
        newd = d + j
        df, massendtime, masstime, conctime, Velocity, head = calcconcmasstime(
            Trial[Trial.index(i)],
            Het[Trial.index(i)],
            Anis[Trial.index(i)],
            gw,
            newd,
            fpre,
            fsuf,
            yin,
            yout,
            xleft,
            xright,
            vars,
            gvarnames,
        )
        for c in times:
            Titles = (
                "Variance:"
                + str(Het[Trial.index(i)])
                + " & Anisotropy:"
                + str(Anis[Trial.index(i)])
            )
            axes[int(math.floor(col / ncols))][col % ncols].plot(
                conctime[int(np.shape(conctime)[0] - c / 5), :, gvarnames.index(chem)],
                yindex,
                label=int(np.shape(conctime)[0] - c / 5) * 5,
                color=colors[times.index(c)],
            )
            axes[int(math.floor(col / ncols))][col % ncols].set_title(
                Titles, fontsize=axissize
            )
            axes[int(math.floor(col / ncols))][col % ncols].tick_params(
                labelsize=ticksize
            )
        axes[int(math.floor(col / ncols))][col % ncols].invert_yaxis()
        col = col + 1
    axes.flat[5].legend(
        (
            int(5 * (np.shape(conctime)[0] - times[0] / 5)),
            int(5 * (np.shape(conctime)[0] - times[1] / 5)),
            int(5 * (np.shape(conctime)[0] - times[2] / 5)),
            int(5 * (np.shape(conctime)[0] - times[3] / 5)),
            int(5 * (np.shape(conctime)[0] - times[4] / 5)),
        ),
        fontsize=legendsize,
    )
    for ax in axes[1]:
        ax.set_xlabel("Concentration of " + chem + " (uM)", size=axissize)
    for ax in axes[:, 0]:
        ax.set_ylabel("Y (cm)", size=axissize)

    return figbig


def generate_timeseries(
    Regimes,
    initseries,
    lastseries,
    Trial,
    Het,
    Anis,
    gw,
    d,
    fpre,
    vars,
    gvarnames,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    AFbiomassvars,
    AFbiomassgvarnames,
    biomassstate,
):
    for Reg in Regimes:
        titlesize = 35
        legendsize = 25
        ticksize = 30
        axissize = 25
        line = "--"
        dots = "dotted"
        biomasscenario = [
            "Active fixed",
            "Active mobile",
            "Inactive fixed",
            "Inactive mobile",
        ]
        prenum = int(
            (biomasscenario.index(biomassstate))
            * (len(AFbiomassgvarnames) / len(biomasscenario))
        )
        xindex = list(range(0, 1095 * 5, 5))
        Tforfpre = ["0/", "1/", "2/", "5/"]
        df0, massendtime0, masstime0, conctime0, Velocity0, head0 = calcconcmasstime(
            Trial,
            Het,
            Anis,
            gw,
            d + Reg + "AR_" + Tforfpre[0],
            fpre,
            fsuf,
            yin,
            yout,
            xleft,
            xright,
            vars,
            gvarnames,
        )
        sumalltime_total0 = calcsum_temp(
            Trial,
            Het,
            Anis,
            gw,
            d + Reg + "AR_" + Tforfpre[0],
            fpre,
            fsuf,
            yin,
            yout,
            xleft,
            xright,
            AFbiomassvars,
            AFbiomassgvarnames,
        )  # absolute numbers
        sumalltime0 = np.zeros(
            [np.shape(sumalltime_total0)[0], np.shape(sumalltime_total0)[1]]
        )
        for k in range(len(AFbiomassgvarnames)):
            sumalltime0[:, k] = sumalltime_total0[:, k] / np.sum(
                sumalltime_total0[:, :], axis=-1
            )
        h = []
        c = np.zeros([1095, np.shape(conctime0)[2]])
        s = np.zeros([1095, np.shape(sumalltime0)[1]])
        for idx in range(1095):
            h.append(head0[-1])
            s[idx, :] = sumalltime0[-1, :]
            c[idx, :] = conctime0[-1, yout, :]
        figbig, axes = plt.subplots(
            nrows=8,
            ncols=len(Tforfpre) - 1,
            figsize=[30, 20],
            sharey="row",
            sharex="col",
        )
        plt.suptitle(
            "Variance: " + str(Het) + " & Anisotropy: " + str(Anis), fontsize=titlesize
        )
        for j, init, last in zip(Tforfpre[1:], initseries, lastseries):
            newd = d + Reg + "AR_" + j
            Titles = "Time series: " + str(Tforfpre.index(j))
            df, massendtime, masstime, conctime, Velocity, head = calcconcmasstime(
                Trial,
                Het,
                Anis,
                gw,
                newd,
                fpre,
                fsuf,
                yin,
                yout,
                xleft,
                xright,
                vars,
                gvarnames,
            )
            sumalltime_total = calcsum_temp(
                Trial,
                Het,
                Anis,
                gw,
                newd,
                fpre,
                fsuf,
                yin,
                yout,
                xleft,
                xright,
                AFbiomassvars,
                AFbiomassgvarnames,
            )  # absolute numbers
            sumalltime = np.zeros(
                [np.shape(sumalltime_total)[0], np.shape(sumalltime_total)[1]]
            )
            concmeantime = np.zeros([1095, np.shape(conctime0)[2]])
            summeantime = np.zeros([1095, np.shape(sumalltime0)[1]])
            for k in range(len(AFbiomassgvarnames)):
                sumalltime[:, k] = sumalltime_total[:, k] / np.sum(
                    sumalltime_total[:, :], axis=-1
                )
            for idx in range(1095):
                concmeantime[idx, :] = np.mean(conctime[:, yout, :], axis=0)
                summeantime[idx, :] = np.mean(sumalltime[:, :], axis=0)
            axes[0][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                head[init:last],
                label="Groundwater head (m)",
                color="black",
            )
            axes[0][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                h[init:last],
                label="Groundwater head (m)",
                color="black",
                linestyle=line,
            )
            ymax, ymin, xposmax, xposmin = localmaxmin(
                head[init:last], h[init:last], init
            )
            axes[0][Tforfpre.index(j) - 1].plot(
                [xposmax, xposmax], [ymin, ymax], ":", c="grey"
            )
            axes[0][Tforfpre.index(j) - 1].plot(
                [xposmin, xposmin], [ymin, ymax], ":", c="grey"
            )
            axes[0][Tforfpre.index(j) - 1].annotate(
                xposmax,
                xy=(xposmax, ymax),
                xytext=(xposmax + 2, ymax - (ymax - ymin) / 2),
                size=legendsize,
            )
            axes[0][Tforfpre.index(j) - 1].annotate(
                xposmin,
                xy=(xposmin, ymin),
                xytext=(xposmin + 2, ymin + (ymax - ymin) / 2),
                size=legendsize,
            )
            axes[0][Tforfpre.index(j) - 1].set_title(Titles, fontsize=titlesize)
            axes[1][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                conctime[init + 1 : last + 1, yout, 0],
                label=gvarnames[0],
                color="black",
            )
            axes[1][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                concmeantime[init + 1 : last + 1, 0],
                label=gvarnames[0],
                color="black",
                linestyle=dots,
            )
            axes[1][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                c[init:last, 0],
                label=gvarnames[0],
                color="black",
                linestyle=line,
            )
            ymax, ymin, xposmax, xposmin = localmaxmin(
                conctime[init + 1 : last + 1, yout, 0], c[init:last, 0], init
            )
            axes[1][Tforfpre.index(j) - 1].plot(
                [xposmax, xposmax], [ymin, ymax], ":", c="grey"
            )
            axes[1][Tforfpre.index(j) - 1].plot(
                [xposmin, xposmin], [ymin, ymax], ":", c="grey"
            )
            axes[1][Tforfpre.index(j) - 1].annotate(
                xposmax,
                xy=(xposmax, ymax),
                xytext=(xposmax + 2, ymax - (ymax - ymin) / 2),
                size=legendsize,
            )
            axes[1][Tforfpre.index(j) - 1].annotate(
                xposmin,
                xy=(xposmin, ymin),
                xytext=(xposmin + 2, ymin + (ymax - ymin) / 2),
                size=legendsize,
            )
            axes[3][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                conctime[init + 1 : last + 1, yout, 1],
                label=gvarnames[1],
                color="red",
            )
            axes[3][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                c[init:last, 1],
                label=gvarnames[1],
                color="red",
                linestyle=line,
            )
            axes[3][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                concmeantime[init + 1 : last + 1, 1],
                label=gvarnames[1],
                color="red",
                linestyle=dots,
            )
            ymax, ymin, xposmax, xposmin = localmaxmin(
                conctime[init + 1 : last + 1, yout, 1], c[init:last, 1], init
            )
            axes[3][Tforfpre.index(j) - 1].plot(
                [xposmax, xposmax], [ymin, ymax], ":", c="grey"
            )
            axes[3][Tforfpre.index(j) - 1].plot(
                [xposmin, xposmin], [ymin, ymax], ":", c="grey"
            )
            axes[3][Tforfpre.index(j) - 1].annotate(
                xposmax,
                xy=(xposmax, ymax),
                xytext=(xposmax + 2, ymax - (ymax - ymin) / 2),
                size=legendsize,
            )
            axes[3][Tforfpre.index(j) - 1].annotate(
                xposmin,
                xy=(xposmin, ymin),
                xytext=(xposmin + 2, ymin + (ymax - ymin) / 2),
                size=legendsize,
            )
            axes[5][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                conctime[init + 1 : last + 1, yout, 2],
                label=gvarnames[2],
                color="blue",
            )
            axes[5][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                c[init:last, 2],
                label=gvarnames[2],
                color="blue",
                linestyle=line,
            )
            axes[5][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                concmeantime[init + 1 : last + 1, 2],
                label=gvarnames[2],
                color="blue",
                linestyle=dots,
            )
            ymax, ymin, xposmax, xposmin = localmaxmin(
                conctime[init + 1 : last + 1, yout, 2], c[init:last, 2], init
            )
            axes[5][Tforfpre.index(j) - 1].plot(
                [xposmax, xposmax], [ymin, ymax], ":", c="grey"
            )
            axes[5][Tforfpre.index(j) - 1].plot(
                [xposmin, xposmin], [ymin, ymax], ":", c="grey"
            )
            axes[5][Tforfpre.index(j) - 1].annotate(
                xposmax,
                xy=(xposmax, ymax),
                xytext=(xposmax + 2, ymax - (ymax - ymin) / 2),
                size=legendsize,
            )
            axes[5][Tforfpre.index(j) - 1].annotate(
                xposmin,
                xy=(xposmin, ymin),
                xytext=(xposmin + 2, ymin + (ymax - ymin) / 2),
                size=legendsize,
            )
            axes[7][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                conctime[init + 1 : last + 1 :, yout, 3],
                label=gvarnames[3],
                color="green",
            )
            axes[7][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                c[init:last, 3],
                label=gvarnames[3],
                color="green",
                linestyle=line,
            )
            axes[7][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                concmeantime[init + 1 : last + 1, 3],
                label=gvarnames[3],
                color="green",
                linestyle=dots,
            )
            ymax, ymin, xposmax, xposmin = localmaxmin(
                conctime[init + 1 : last + 1, yout, 3], c[init:last, 3], init
            )
            axes[7][Tforfpre.index(j) - 1].plot(
                [xposmax, xposmax], [ymin, ymax], ":", c="grey"
            )
            axes[7][Tforfpre.index(j) - 1].plot(
                [xposmin, xposmin], [ymin, ymax], ":", c="grey"
            )
            axes[7][Tforfpre.index(j) - 1].annotate(
                xposmax,
                xy=(xposmax, ymax),
                xytext=(xposmax + 2, ymax - (ymax - ymin) / 2),
                size=legendsize,
            )
            axes[7][Tforfpre.index(j) - 1].annotate(
                xposmin,
                xy=(xposmin, ymin),
                xytext=(xposmin + 2, ymin + (ymax - ymin) / 2),
                size=legendsize,
            )
            axes[2][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                sumalltime[init:last, prenum + 0],
                label=AFbiomassgvarnames[prenum + 0],
                color="red",
            )
            axes[2][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                s[init:last, prenum + 0],
                label=AFbiomassgvarnames[prenum + 0],
                color="red",
                linestyle=line,
            )
            axes[2][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                summeantime[init + 1 : last + 1, prenum + 0],
                label=AFbiomassgvarnames[prenum + 0],
                color="red",
                linestyle=dots,
            )
            ymax, ymin, xposmax, xposmin = localmaxmin(
                sumalltime[init:last, prenum + 0], s[init:last, prenum + 0], init
            )
            axes[2][Tforfpre.index(j) - 1].plot(
                [xposmax, xposmax], [ymin, ymax], ":", c="grey"
            )
            axes[2][Tforfpre.index(j) - 1].plot(
                [xposmin, xposmin], [ymin, ymax], ":", c="grey"
            )
            axes[2][Tforfpre.index(j) - 1].annotate(
                xposmax,
                xy=(xposmax, ymax),
                xytext=(xposmax + 2, ymax - (ymax - ymin) / 2),
                size=legendsize,
            )
            axes[2][Tforfpre.index(j) - 1].annotate(
                xposmin,
                xy=(xposmin, ymin),
                xytext=(xposmin + 2, ymin + (ymax - ymin) / 2),
                size=legendsize,
            )
            axes[4][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                sumalltime[init:last, prenum + 1],
                label=AFbiomassgvarnames[prenum + 1],
                color="blue",
            )
            axes[4][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                s[init:last, prenum + 1],
                label=AFbiomassgvarnames[prenum + 1],
                color="blue",
                linestyle=line,
            )
            axes[4][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                summeantime[init + 1 : last + 1, prenum + 1],
                label=AFbiomassgvarnames[prenum + 1],
                color="blue",
                linestyle=dots,
            )
            ymax, ymin, xposmax, xposmin = localmaxmin(
                sumalltime[init:last, 1], s[init:last, 1], init
            )
            axes[4][Tforfpre.index(j) - 1].plot(
                [xposmax, xposmax], [ymin, ymax], ":", c="grey"
            )
            axes[4][Tforfpre.index(j) - 1].plot(
                [xposmin, xposmin], [ymin, ymax], ":", c="grey"
            )
            axes[4][Tforfpre.index(j) - 1].annotate(
                xposmax,
                xy=(xposmax, ymax),
                xytext=(xposmax + 2, ymax - (ymax - ymin) / 2),
                size=legendsize,
            )
            axes[4][Tforfpre.index(j) - 1].annotate(
                xposmin,
                xy=(xposmin, ymin),
                xytext=(xposmin + 2, ymin + (ymax - ymin) / 2),
                size=legendsize,
            )
            axes[6][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                sumalltime[init:last, prenum + 2],
                label=AFbiomassgvarnames[prenum + 2],
                color="green",
            )
            axes[6][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                s[init:last, prenum + 2],
                label=AFbiomassgvarnames[prenum + 2],
                color="green",
                linestyle=line,
            )
            axes[6][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                summeantime[init + 1 : last + 1, prenum + 2],
                label=AFbiomassgvarnames[prenum + 2],
                color="green",
                linestyle=dots,
            )
            ymax, ymin, xposmax, xposmin = localmaxmin(
                sumalltime[init:last, prenum + 2], s[init:last, prenum + 2], init
            )
            axes[6][Tforfpre.index(j) - 1].plot(
                [xposmax, xposmax], [ymin, ymax], ":", c="grey"
            )
            axes[6][Tforfpre.index(j) - 1].plot(
                [xposmin, xposmin], [ymin, ymax], ":", c="grey"
            )
            axes[6][Tforfpre.index(j) - 1].annotate(
                xposmax,
                xy=(xposmax, ymax),
                xytext=(xposmax + 2, ymax - (ymax - ymin) / 2),
                size=legendsize,
            )
            axes[6][Tforfpre.index(j) - 1].annotate(
                xposmin,
                xy=(xposmin, ymin),
                xytext=(xposmin + 2, ymin + (ymax - ymin) / 2),
                size=legendsize,
            )

            ax = axes[:, Tforfpre.index(j) - 1]
            if Tforfpre.index(j) - 1 == 0:
                off = 150
                ax[0].annotate(
                    "Velocity",
                    xy=(0, 0.5),
                    xytext=(-ax[0].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[1].annotate(
                    str(gvarnames[0]),
                    xy=(0, 0.5),
                    xytext=(-ax[1].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[3].annotate(
                    str(gvarnames[1]),
                    xy=(0, 0.5),
                    xytext=(-ax[3].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[5].annotate(
                    str(gvarnames[2]),
                    xy=(0, 0.5),
                    xytext=(-ax[5].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[7].annotate(
                    str(gvarnames[3]),
                    xy=(0, 0.5),
                    xytext=(-ax[7].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[2].annotate(
                    str("Aerobes"),
                    xy=(0, 0.5),
                    xytext=(-ax[2].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[4].annotate(
                    "Ammonia\noxidizers",
                    xy=(0, 0.5),
                    xytext=(-ax[4].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[6].annotate(
                    "Nitrate\nreducers",
                    xy=(0, 0.5),
                    xytext=(-ax[6].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )

        for ax in axes[:, 0]:
            ax.tick_params(labelsize=ticksize)
        for ax in axes[7]:
            ax.set_xlabel("Time (days)", fontsize=axissize)
            ax.tick_params(labelsize=ticksize)
        picname = (
            "Z:/Saturated_flow/diffusion_transient/chemsand_"
            + biomassstate
            + "_biomasswithvel_temp_"
            + str(Tforfpre.index(j))
            + "_"
            + Reg
            + "_"
            + str(Trial)
            + "_ZOOMED.png"
        )
        plt.savefig(picname, dpi=300, bbox_inches="tight", pad_inches=0)

    return None


def generate_chem_timeseries(
    Regimes,
    initseries,
    lastseries,
    Trial,
    Het,
    Anis,
    gw,
    d,
    fpre,
    vars,
    gvarnames,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
):
    for Reg in Regimes:
        titlesize = 35
        legendsize = 25
        ticksize = 30
        axissize = 25
        line = "--"
        dots = "dotted"
        xindex = list(range(0, 1095 * 5, 5))
        Tforfpre = ["0/", "1/", "2/", "5/"]
        #        df0,massendtime0, masstime0, conctime0, Velocity0, head0 = calcconcmasstime (Trial, Het, Anis, gw, d+Reg + "AR_"+Tforfpre[0], fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
        #        h = []
        #        c = np.zeros([1095,np.shape(conctime0)[2]])
        #        for idx in range(1095):
        #            h.append(head0[-1])
        #            c[idx,:] = conctime0[-1,yout,:]
        figbig, axes = plt.subplots(
            nrows=5,
            ncols=len(Tforfpre) - 1,
            figsize=[30, 20],
            sharey="row",
            sharex="col",
        )
        plt.suptitle(
            "Variance: " + str(Het) + " & Anisotropy: " + str(Anis), fontsize=titlesize
        )
        for j, init, last in zip(Tforfpre[1:], initseries, lastseries):
            newd = d + Reg + "AR_" + j
            Titles = "Time series: " + str(Tforfpre.index(j))
            df, massendtime, masstime, conctime, Velocity, head = calcconcmasstime(
                Trial,
                Het,
                Anis,
                gw,
                newd,
                fpre,
                fsuf,
                yin,
                yout,
                xleft,
                xright,
                vars,
                gvarnames,
            )
            concmeantime = np.zeros([1095, np.shape(conctime)[2]])
            for idx in range(1095):
                concmeantime[idx, :] = np.mean(conctime[:, yout, :], axis=0)
            axes[0][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                head[init:last],
                label="Groundwater head (m)",
                color="black",
            )
            #            axes[0][Tforfpre.index(j)-1].plot(xindex[init:last], h[init:last], label = "Groundwater head (m)", color = "black", linestyle = line)
            #            ymax, ymin, xposmax, xposmin = localmaxmin(head[init:last],h[init:last], init)
            #            axes[0][Tforfpre.index(j)-1].plot([xposmax, xposmax], [ymin, ymax],':',  c ='grey')
            #            axes[0][Tforfpre.index(j)-1].plot([xposmin, xposmin], [ymin, ymax],':',  c ='grey')
            #            axes[0][Tforfpre.index(j)-1].annotate(xposmax, xy = (xposmax, ymax), xytext = (xposmax+2, ymax - (ymax-ymin)/2), size = legendsize)
            #            axes[0][Tforfpre.index(j)-1].annotate(xposmin, xy = (xposmin, ymin), xytext = (xposmin+2, ymin + (ymax-ymin)/2), size = legendsize)
            axes[0][Tforfpre.index(j) - 1].set_title(Titles, fontsize=titlesize)
            axes[1][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                conctime[init + 1 : last + 1, yout, 0],
                label=gvarnames[0],
                color="black",
            )
            axes[1][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                concmeantime[init + 1 : last + 1, 0],
                label=gvarnames[0],
                color="black",
                linestyle=dots,
            )
            #            axes[1][Tforfpre.index(j)-1].plot(xindex[init:last], c[init:last,0], label = gvarnames[0], color = "black", linestyle = line)
            #            ymax, ymin, xposmax, xposmin = localmaxmin(conctime[init+1:last+1,yout,0], c[init:last,0], init)
            #            axes[1][Tforfpre.index(j)-1].plot([xposmax, xposmax], [ymin, ymax],':',  c ='grey')
            #            axes[1][Tforfpre.index(j)-1].plot([xposmin, xposmin], [ymin, ymax],':',  c ='grey')
            #            axes[1][Tforfpre.index(j)-1].annotate(xposmax, xy = (xposmax, ymax), xytext = (xposmax+2, ymax-(ymax-ymin)/2), size = legendsize)
            #            axes[1][Tforfpre.index(j)-1].annotate(xposmin, xy = (xposmin, ymin), xytext = (xposmin+2, ymin+(ymax-ymin)/2), size = legendsize)
            axes[2][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                conctime[init + 1 : last + 1, yout, 1],
                label=gvarnames[1],
                color="red",
            )
            #            axes[2][Tforfpre.index(j)-1].plot(xindex[init:last], c[init:last,1], label = gvarnames[1], color = "red", linestyle = line)
            axes[2][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                concmeantime[init + 1 : last + 1, 1],
                label=gvarnames[1],
                color="red",
                linestyle=dots,
            )
            #            ymax, ymin, xposmax, xposmin = localmaxmin(conctime[init+1:last+1,yout,1], c[init:last,1], init)
            #            axes[2][Tforfpre.index(j)-1].plot([xposmax, xposmax], [ymin, ymax],':',  c ='grey')
            #            axes[2][Tforfpre.index(j)-1].plot([xposmin, xposmin], [ymin, ymax],':',  c ='grey')
            #            axes[2][Tforfpre.index(j)-1].annotate(xposmax, xy = (xposmax, ymax), xytext = (xposmax+2, ymax-(ymax-ymin)/2), size = legendsize)
            #            axes[2][Tforfpre.index(j)-1].annotate(xposmin, xy = (xposmin, ymin), xytext = (xposmin+2, ymin+(ymax-ymin)/2), size = legendsize)
            axes[3][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                conctime[init + 1 : last + 1, yout, 2],
                label=gvarnames[2],
                color="blue",
            )
            #            axes[3][Tforfpre.index(j)-1].plot(xindex[init:last], c[init:last,2], label = gvarnames[2], color = "blue", linestyle = line)
            axes[3][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                concmeantime[init + 1 : last + 1, 2],
                label=gvarnames[2],
                color="blue",
                linestyle=dots,
            )
            #            ymax, ymin, xposmax, xposmin = localmaxmin(conctime[init+1:last+1,yout,2], c[init:last,2], init)
            #            axes[3][Tforfpre.index(j)-1].plot([xposmax, xposmax], [ymin, ymax],':',  c ='grey')
            #            axes[3][Tforfpre.index(j)-1].plot([xposmin, xposmin], [ymin, ymax],':',  c ='grey')
            #            axes[3][Tforfpre.index(j)-1].annotate(xposmax, xy = (xposmax, ymax), xytext = (xposmax+2, ymax-(ymax-ymin)/2), size = legendsize)
            #            axes[3][Tforfpre.index(j)-1].annotate(xposmin, xy = (xposmin, ymin), xytext = (xposmin+2, ymin+(ymax-ymin)/2), size = legendsize)
            axes[4][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                conctime[init + 1 : last + 1 :, yout, 3],
                label=gvarnames[3],
                color="green",
            )
            #            axes[4][Tforfpre.index(j)-1].plot(xindex[init:last], c[init:last,3], label = gvarnames[3], color = "green", linestyle = line)
            axes[4][Tforfpre.index(j) - 1].plot(
                xindex[init:last],
                concmeantime[init + 1 : last + 1, 3],
                label=gvarnames[3],
                color="green",
                linestyle=dots,
            )
            #            ymax, ymin, xposmax, xposmin = localmaxmin(conctime[init+1:last+1,yout,3], c[init:last,3], init)
            #            axes[4][Tforfpre.index(j)-1].plot([xposmax, xposmax], [ymin, ymax],':',  c ='grey')
            #            axes[4][Tforfpre.index(j)-1].plot([xposmin, xposmin], [ymin, ymax],':',  c ='grey')
            #            axes[4][Tforfpre.index(j)-1].annotate(xposmax, xy = (xposmax, ymax), xytext = (xposmax+2, ymax-(ymax-ymin)/2), size = legendsize)
            #            axes[4][Tforfpre.index(j)-1].annotate(xposmin, xy = (xposmin, ymin), xytext = (xposmin+2, ymin+(ymax-ymin)/2), size = legendsize)

            ax = axes[:, Tforfpre.index(j) - 1]
            if Tforfpre.index(j) - 1 == 0:
                off = 150
                ax[0].annotate(
                    "Velocity",
                    xy=(0, 0.5),
                    xytext=(-ax[0].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[1].annotate(
                    str(gvarnames[0]),
                    xy=(0, 0.5),
                    xytext=(-ax[1].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[2].annotate(
                    str(gvarnames[1]),
                    xy=(0, 0.5),
                    xytext=(-ax[2].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[3].annotate(
                    str(gvarnames[2]),
                    xy=(0, 0.5),
                    xytext=(-ax[3].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )
                ax[4].annotate(
                    str(gvarnames[3]),
                    xy=(0, 0.5),
                    xytext=(-ax[4].yaxis.labelpad - off, 0),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    rotation="vertical",
                    size=axissize,
                )

        for ax in axes[:, 0]:
            ax.tick_params(labelsize=ticksize)
        for ax in axes[4]:
            ax.set_xlabel("Time (days)", fontsize=axissize)
            ax.tick_params(labelsize=ticksize)
        picname = (
            "Z:/Saturated_flow/diffusion_transient/chemswithvel_temp_"
            + str(Tforfpre.index(j))
            + "_"
            + Reg
            + "_"
            + str(Trial)
            + "_ZOOMED.png"
        )
        plt.savefig(picname, dpi=300, bbox_inches="tight", pad_inches=0)

    return None


def amplitude_biomass(data, Chemseries):
    Regimes = ["Slow", "Medium", "Fast"]
    #    colseries=["indianred","g","steelblue"]
    species = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    position = ["Immobile", "Mobile"]
    fig, axes = plt.subplots(
        nrows=2, ncols=3, figsize=[10, 6], sharey="row", sharex="col"
    )
    #    plt.suptitle("Change in biomass with respect to steady state conditions", fontsize = 20)
    for k in Chemseries:
        dfc = data[(data["Chem"] == k)]
        colidx1 = Chemseries.index(k)
        for i in Regimes:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp["Regime"] == i]
            l1 = axes.flat[colidx1].scatter(
                "fraction",
                "Min",  # c = "Time_series", cmap
                c=colseries[Regimes.index(i)],
                data=dfcr,
                label="Minimum",
                marker="^",
            )
            l2 = axes.flat[colidx1].scatter(
                "fraction",
                "Max",  # c = "Time_series", cmap
                c=colseries[Regimes.index(i)],
                data=dfcr,
                label="Maximum",
                marker="o",
            )
            axes.flat[colidx1].tick_params(axis="y", labelsize=15)
            axes.flat[colidx1].tick_params(axis="x", labelsize=15)
            axes.flat[colidx1].set_ylim((-4, 4))
            if (Chemseries.index(k) == len(Chemseries) - 1) & (
                Regimes.index(i) == len(Regimes) - 2
            ):
                handles1 = [l1, l2]
    #    fig.legend(handles1, ["Maximum decrease", "Maximum increase"], bbox_to_anchor=(0.9, 0.1),loc = 'center', title = "Value", fontsize = 12)
    #    fig.legend(handles=[red_patch, green_patch, blue_patch], bbox_to_anchor=(0.9, 0.5),loc = 'center', title = "Flow regime", fontsize = 12)
    #    fig.legend(handles=[low_var, mid_var, high_var], bbox_to_anchor=(0.7, 0.4),loc = 'center', title = "Variance in\nTime series", fontsize = 12)
    for ax, typsp in zip(axes[0, :], species):
        ax.set_title(typsp, fontsize=15)
    axes[0, -1].annotate(
        position[0],
        xy=(0, 0.5),
        xytext=(180, 0),
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
        xytext=(180, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )
    #    axes[0,0].set_ylabel(activity[0], fontsize = 15)
    #    axes[1,0].set_ylabel(activity[1], fontsize = 15)
    plt.annotate(
        "Normalized variation from steady state conditions",
        xy=(0, 2.2),
        xytext=(-450, -150),
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
    plt.tight_layout()

    return fig


def RMSamp_biomass(data, Chemseries, plotvar, yaxislabel, Regimestoplot):
    species = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    position = ["Immobile", "Mobile"]
    markerseries = ["d", "^", "o"]
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=[10, 6], sharex=True)
    for k in Chemseries:
        dfc = data[(data["Chem"] == k)]
        colidx1 = Chemseries.index(k)
        for i in Regimestoplot:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp["Regime"] == i]
            axes.flat[colidx1].scatter(
                "fraction",
                plotvar,
                c="Time_series",
                cmap=colseries[Regimes.index(i)],
                data=dfcr,
                marker=markerseries[Regimes.index(i)],
            )
            axes.flat[colidx1].tick_params(axis="y", labelsize=15)
            axes.flat[colidx1].tick_params(axis="x", labelsize=15)
    #    fig.legend(handles=[red_patch, green_patch, blue_patch], bbox_to_anchor=(0.9, 0.5),loc = 'center', title = "Flow regime", fontsize = 12)
    #    fig.legend(handles=[low_var, mid_var, high_var], bbox_to_anchor=(0.7, 0.4),loc = 'center', title = "Variance in\nTime series", fontsize = 12)
    for ax, typsp in zip(axes[0, :], species):
        ax.set_title(typsp, fontsize=15)
    axes[0, -1].annotate(
        position[0],
        xy=(0, 0.5),
        xytext=(180, 0),
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
        xytext=(180, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )
    plt.annotate(
        yaxislabel,
        xy=(0, 2.2),
        xytext=(-450, -150),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=18,
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
        fontsize=18,
    )
    #    plt.tight_layout()

    return fig


def amplitude_chem(data, Chemseries):
    Regimes = ["Slow", "Medium", "Fast"]
    fig, axes = plt.subplots(
        nrows=len(Chemseries), ncols=1, figsize=[4, 10], sharey="row", sharex="col"
    )
    plt.suptitle(
        "Change in flux averaged concentration with \nrespect to steady state conditions",
        fontsize=20,
    )
    for k in Chemseries:
        dfc = data[(data["Chem"] == k)]  # &(data['Min']>-0.9)&(data['Max']<1.5)]
        colidx1 = Chemseries.index(k)
        for i in Regimes:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp["Regime"] == i]
            l1 = axes.flat[colidx1].scatter(
                "fraction",
                "Min",
                c="Time_series",
                cmap=colseries[Regimes.index(i)],
                data=dfcr,
                label="Minimum",
                marker="^",
            )
            l2 = axes.flat[colidx1].scatter(
                "fraction",
                "Max",
                c="Time_series",
                cmap=colseries[Regimes.index(i)],
                data=dfcr,
                label="Maximum",
                marker="o",
            )
            axes.flat[colidx1].tick_params(axis="y", labelsize=15)
            axes.flat[colidx1].tick_params(axis="x", labelsize=15)
            if (Chemseries.index(k) == len(Chemseries) - 1) & (
                Regimes.index(i) == len(Regimes) - 2
            ):
                handles1 = [l1, l2]
    fig.legend(
        handles1,
        ["Maximum decrease", "Maximum increase"],
        bbox_to_anchor=(1.45, 0.1),
        loc="center",
        title="Value",
        fontsize=12,
    )
    fig.legend(
        handles=[red_patch, green_patch, blue_patch],
        bbox_to_anchor=(1.45, 0.25),
        loc="center",
        title="Flow regime",
        fontsize=12,
    )
    fig.legend(
        handles=[low_var, mid_var, high_var],
        bbox_to_anchor=(1.45, 0.4),
        loc="center",
        title="Variance in\nTime series",
    )
    for ax, typsp in zip(axes[:], Chemseries):
        ax.set_ylabel(typsp, fontsize=15)
    plt.annotate(
        "Amplitude as fraction of steady state conditions",
        xy=(0, 3.0),
        xytext=(-100, 0),
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
        xy=(0.1, -0.3),
        xytext=(100, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=15,
    )

    return fig


def correlation_delay_chem(data, Chemseries):
    Regimes = ["Slow", "Medium", "Fast"]
    fig1, axes1 = plt.subplots(nrows=len(Chemseries), ncols=1, figsize=[4, 10])
    plt.suptitle(
        "Change in flux averaged concentration with \nrespect to steady state conditions",
        fontsize=20,
    )
    for k in Chemseries:
        dfc = data[data["Chem"] == k]
        colidx1 = Chemseries.index(k)
        for i in Regimes:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp["Regime"] == i]
            paxis = axes1.flat[colidx1]
            l1 = axes1.flat[colidx1].scatter(
                "fraction",
                "Correlation",
                c="Time_series",
                cmap=colseries[Regimes.index(i)],
                data=dfcr,
                label="Correlation",
                marker="^",
            )
            paxis.tick_params(axis="y", labelsize=15)
            paxis.tick_params(axis="x", labelsize=15)
            if Chemseries.index(k) != len(Chemseries) - 1:
                paxis.set_xticklabels([])
                if Regimes.index(i) == len(Regimes) - 2:
                    handles1 = [l1]
    #    fig1.legend(handles1, ["Correlation"], bbox_to_anchor=(1.7, 0.1),loc = 'center', title = "Value")
    fig1.legend(
        handles=[red_patch, green_patch, blue_patch],
        bbox_to_anchor=(1.4, 0.15),
        loc="center",
        title="Flow regime",
        fontsize=12,
    )
    fig1.legend(
        handles=[low_var, mid_var, high_var],
        bbox_to_anchor=(1.4, 0.3),
        loc="center",
        title="Variance in\nTime series",
        fontsize=12,
    )
    for ax, typsp in zip(axes1[:], Chemseries):
        ax.set_ylabel(typsp, fontsize=15)
    plt.annotate(
        "Correlation",
        xy=(0, 2.0),
        xytext=(-100, 60),
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
        xy=(0.1, -0.3),
        xytext=(100, -20),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=15,
    )

    fig2, axes2 = plt.subplots(nrows=len(Chemseries), ncols=1, figsize=[4, 10])
    plt.suptitle(
        "Change in flux averaged concentration with \nrespect to steady state conditions",
        fontsize=20,
    )
    for k in Chemseries:
        dfc = data[data["Chem"] == k]
        colidx1 = Chemseries.index(k)
        for i in Regimes:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp["Regime"] == i]
            paxis = axes2.flat[colidx1]
            l1 = axes2.flat[colidx1].scatter(
                "fraction",
                "Delayfractionofbth",
                c="Time_series",
                cmap=colseries[Regimes.index(i)],
                data=dfcr,
                label="Correlation",
                marker="^",
            )
            paxis.tick_params(axis="y", labelsize=15)
            paxis.tick_params(axis="x", labelsize=15)
            paxis.set_ylim(
                np.min(dfc["Delayfractionofbth"]) - 0.1,
                np.max(dfc["Delayfractionofbth"]) + 0.1,
            )
            paxis.tick_params(axis="y", labelsize=15)
            if Chemseries.index(k) != len(Chemseries) - 1:
                paxis.set_xticklabels([])
                if Regimes.index(i) == len(Regimes) - 2:
                    handles1 = [l1]
    #    fig2.legend(handles1, ["Correlation"], bbox_to_anchor=(1.7, 0.1),loc = 'center', title = "Value")
    fig2.legend(
        handles=[red_patch, green_patch, blue_patch],
        bbox_to_anchor=(1.4, 0.15),
        loc="center",
        title="Flow regime",
        fontsize=12,
    )
    fig2.legend(
        handles=[low_var, mid_var, high_var],
        bbox_to_anchor=(1.4, 0.3),
        loc="center",
        title="Variance in\nTime series",
        fontsize=12,
    )
    for ax, typsp in zip(axes2[:], Chemseries):
        ax.set_ylabel(typsp, fontsize=15)
    plt.annotate(
        "Time lag with respect to fraction of breakthrough time",
        xy=(0, 1.5),
        xytext=(-100, 60),
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
        xy=(0.1, -0.3),
        xytext=(100, -20),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=15,
    )

    return fig1, fig2


def correlation_delay_biomass(data, Chemseries):
    Regimes = ["Slow", "Medium", "Fast"]
    species = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    position = ["Immobile", "Mobile"]
    #    activity = ["Active", "Inactive"]
    fig1, axes1 = plt.subplots(nrows=2, ncols=3, figsize=[10, 6], sharex="col")
    #    plt.suptitle("Change in biomass with respect to steady state conditions", fontsize = 20)
    for k in Chemseries:
        dfc = data[data["Chem"] == k]
        colidx1 = Chemseries.index(k)
        for i in Regimes:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp["Regime"] == i]
            paxis = axes1.flat[colidx1]
            l1 = paxis.scatter(
                "fraction",
                "Correlation",
                c="Time_series",
                cmap=colseries[Regimes.index(i)],
                data=dfcr,
                label="Correlation",
                marker="^",
            )
            paxis.tick_params(axis="y", labelsize=15)
            paxis.tick_params(axis="x", labelsize=15)
            if (Chemseries.index(k) == len(Chemseries) - 1) & (
                Regimes.index(i) == len(Regimes) - 2
            ):
                handles1 = [l1]
    fig1.legend(
        handles=[low_var, mid_var, high_var],
        bbox_to_anchor=(1.06, 0.4),
        loc="center",
        title="Variance in\nTime series",
        fontsize=12,
    )
    fig1.legend(
        handles=[red_patch, green_patch, blue_patch],
        bbox_to_anchor=(1.06, 0.2),
        loc="center",
        title="Flow regime",
        fontsize=12,
    )
    for ax, typsp in zip(axes1[0, :], species):
        ax.set_title(typsp, fontsize=15)
    axes1[0, -1].annotate(
        position[0],
        xy=(1, 0.5),
        xytext=(10, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )
    axes1[1, -1].annotate(
        position[1],
        xy=(1, 0.5),
        xytext=(10, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )
    fig1.tight_layout()
    plt.annotate(
        "Correlation",
        xy=(-0.3, 2.0),
        xytext=(-450, -150),
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
        xy=(0, 0),
        xytext=(-100, -50),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=15,
    )

    fig2, axes2 = plt.subplots(nrows=2, ncols=3, figsize=[10, 6], sharex="col")
    plt.suptitle(
        "Change in biomass with respect to steady state conditions", fontsize=20
    )
    for k in Chemseries:
        dfc = data[data["Chem"] == k]
        colidx1 = Chemseries.index(k)
        for i in Regimes:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp["Regime"] == i]
            paxis = axes2.flat[colidx1]
            l1 = paxis.scatter(
                "fraction",
                "Delayfractionofbth",
                c="Time_series",
                cmap=colseries[Regimes.index(i)],
                data=dfcr,
                label="Correlation",
                marker="^",
            )
            paxis.tick_params(axis="y", labelsize=15)
            paxis.tick_params(axis="x", labelsize=15)
            paxis.set_ylim(
                np.min(dfc["Delayfractionofbth"]) - 0.1,
                np.max(dfc["Delayfractionofbth"]) + 0.1,
            )
            paxis.tick_params(axis="y", labelsize=15)  # , labelcolor = 'b')
            if (Chemseries.index(k) == len(Chemseries) - 1) & (
                Regimes.index(i) == len(Regimes) - 2
            ):
                handles1 = [l1]
    fig2.legend(
        handles=[red_patch, green_patch, blue_patch],
        bbox_to_anchor=(0.5, 0.5),
        loc="center",
        title="Flow regime",
        fontsize=12,
    )
    fig2.legend(
        handles=[low_var, mid_var, high_var],
        bbox_to_anchor=(0.3, 0.5),
        loc="center",
        title="Variance in\nTime series",
        fontsize=12,
    )
    fig2.tight_layout()
    for ax, typsp in zip(axes2[0, :], species):
        ax.set_title(typsp, fontsize=15)
    axes2[0, -1].annotate(
        position[0],
        xy=(1, 0.5),
        xytext=(30, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )
    axes2[1, -1].annotate(
        position[1],
        xy=(1, 0.5),
        xytext=(30, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )
    plt.annotate(
        "Time lag with respect to fraction of breakthrough time",
        xy=(0, 2.2),
        xytext=(-480, -170),
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

    return fig1, fig2


def autocovariance_chem(data, Chemseries):
    Regimes = ["Slow", "Medium", "Fast"]
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=[11, 8], sharex=True)
    #    plt.suptitle("Autocovariance in flux averaged concentration", fontsize = 20)
    for k in Chemseries:
        dfc = data[(data["Chem"] == k)]
        colidx1 = Chemseries.index(k)
        for i in Regimes:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp["Regime"] == i]
            axes.flat[colidx1].scatter(
                "fraction",
                "Autocovariance",
                c="Time_series",
                cmap=colseries[Regimes.index(i)],
                data=dfcr,
                label="Time_series",
            )
            axes.flat[colidx1].tick_params(axis="y", labelsize=15)
            axes.flat[colidx1].tick_params(axis="x", labelsize=15)
            axes.flat[colidx1].set_title(k, fontsize=15)
            axes.flat[colidx1].set_yscale("log")
    fig.legend(
        handles=[red_patch, green_patch, blue_patch],
        bbox_to_anchor=(1.0, 0.25),
        loc="center",
        title="Flow regime",
        fontsize=12,
    )
    fig.legend(
        handles=[low_var, mid_var, high_var],
        bbox_to_anchor=(0.7, 0.4),
        loc="center",
        title="Variance in\nTime series",
        fontsize=12,
    )
    plt.annotate(
        "Autocovariance",
        xy=(-1.2, 1),
        xytext=(-100, 0),
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
        xy=(0, 0),
        xytext=(0, -50),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=15,
    )

    return fig


def amplitude_chem_2x2(data, Chemseries, Regimetoplot):
    Regimes = ["Slow"]  # ["Slow", "Medium", "Fast"]
    #    colseries = ["indianred","g","steelblue"]
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=[11, 8], sharex=True)
    #    plt.suptitle("Change in flux averaged concentration with \nrespect to steady state conditions", fontsize = 20)
    for k in Chemseries:
        dfc = data[(data["Chem"] == k)]  # &(data['Time_series']==3)]
        colidx1 = Chemseries.index(k)
        for i in Regimetoplot:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp["Regime"] == i]
            l1 = axes.flat[colidx1].scatter(
                "fraction",
                "Min",
                c="Time_series",
                cmap=colseries[Regimes.index(i)],
                data=dfcr,
                label="Minimum",
                marker="^",
            )
            l2 = axes.flat[colidx1].scatter(
                "fraction",
                "Max",
                c="Time_series",
                cmap=colseries[Regimes.index(i)],
                data=dfcr,
                label="Maximum",
                marker="o",
            )
            axes.flat[colidx1].tick_params(axis="y", labelsize=15)
            axes.flat[colidx1].tick_params(axis="x", labelsize=15)
            axes.flat[colidx1].set_title(k, fontsize=15)
            #            if((Chemseries.index(k) == len(Chemseries)-1) & (Regimes.index(i)==len(Regimes)-2)):
            handles1 = [l1, l2]
            if k == "DO":
                axes.flat[colidx1].set_ylim((-4.5, 4.5))
            else:
                axes.flat[colidx1].set_ylim((-0.55, 0.55))
    fig.legend(
        handles1,
        ["Maximum decrease", "Maximum increase"],
        bbox_to_anchor=(1.02, 0.1),
        loc="center",
        title="Value",
        fontsize=12,
    )
    fig.legend(
        handles=[red_patch, green_patch, blue_patch],
        bbox_to_anchor=(1.0, 0.25),
        loc="center",
        title="Flow regime",
        fontsize=12,
    )
    #    fig.legend(handles=[low_var, mid_var, high_var], bbox_to_anchor=(0.7, 0.4),loc = 'center', title = "Variance in\nTime series", fontsize = 12)
    plt.annotate(
        "Normalized variation from steady state conditions",
        xy=(-1.2, 1),
        xytext=(-100, 0),
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
        xy=(0, 0),
        xytext=(0, -50),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=15,
    )

    return fig


def RMS_chem_2x2(data, Chemseries, plotvar, yaxislabel, Regimetoplot):
    indices = list(Regimes.index(i) for i in Regimetoplot)
    patches = list(patchlist[i] for i in indices)
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=[11, 8], sharex=True)
    #    plt.suptitle("Change in flux averaged concentration with \nrespect to steady state conditions", fontsize = 20)
    for k in Chemseries:
        dfc = data[(data["Chem"] == k)]  # &(data['Time_series']==3)]
        colidx1 = Chemseries.index(k)
        for i in Regimetoplot:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp["Regime"] == i]
            axes.flat[colidx1].scatter(
                "fraction",
                plotvar,
                c="Time_series",
                cmap=colseries[Regimes.index(i)],
                data=dfcr,
                label="Minimum",
                marker=markerseries[Regimes.index(i)],
            )
            axes.flat[colidx1].tick_params(axis="y", labelsize=15)
            axes.flat[colidx1].tick_params(axis="x", labelsize=15)
            axes.flat[colidx1].set_title(k, fontsize=15)
    #            axes.flat[colidx1].set_ylim((0,60))
    #            if (k !="DO"):
    #                axes.flat[colidx1].set_ylim((-4.5,4.5))
    #                axes.flat[colidx1].set_yscale("log")
    #            else:
    #                axes.flat[colidx1].set_ylim((-0.55,0.55))
    fig.legend(
        handles=patches,
        bbox_to_anchor=(0.97, 0.25),
        loc="center",
        title="Flow regime",
        fontsize=12,
    )
    fig.legend(
        handles=[low_var, mid_var, high_var],
        bbox_to_anchor=(0.97, 0.5),
        loc="center",
        title="Variance in\nTime series",
        fontsize=12,
    )
    plt.annotate(
        yaxislabel,
        xy=(-1.2, 1.2),
        xytext=(-85, 5),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=18,
    )
    plt.annotate(
        "Fraction of breakthrough time in base case",
        xy=(0, 0),
        xytext=(0, -50),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=18,
    )

    return fig


def correlationdistributionchem(data, Chemseries):
    bins = [-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,
            0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1]
    centers = np.mean([bins[:-1], bins[1:]], axis=0)
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=[8, 8], sharex=True)
    for k in Chemseries:
        datac = data[data["Chem"] == k]
        ax1 = axes.flat[Chemseries.index(k)]
        ax1.bar(
            centers,
            np.histogram(datac[datac["Regime"] == "Slow"]["Crosscorrelation"], bins)[0],
            alpha=0.7,
            width=0.1,
            label="Slow flow",
            color="indianred",
        )  # , ax = ax1)
        ax1.bar(
            centers,
            np.histogram(datac[datac["Regime"] == "Medium"]["Crosscorrelation"], bins)[
                0
            ],
            alpha=0.7,
            width=0.1,
            label="Medium flow",
            color="g",
        )  # , ax = ax1)
        ax1.bar(
            centers,
            np.histogram(datac[datac["Regime"] == "Fast"]["Crosscorrelation"], bins)[0],
            alpha=0.7,
            width=0.1,
            label="Fast flow",
            color="steelblue",
        )  # , ax = ax1)
        ax1.set_title(k, fontsize=15)
        ax1.set_xlim((-1, 1))
        ax1.tick_params(axis="y", labelsize=15)
        ax1.tick_params(axis="x", labelsize=15)
        ax1.set_xlabel("")
    plt.legend()
    plt.annotate(
        "Number of scenarios",
        xy=(-1.2, 1),
        xytext=(-85, 5),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=18,
    )
    plt.annotate(
        "Peak cross-correlation between external forcing and\nresulting flux-averaged outlet concentration",
        xy=(0, 0),
        xytext=(0, -50),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=18,
    )

    return fig


def correlationdistributionbiomass(data, Chemseries):
    species = ["Aerobes", "Ammonia oxidizers", "Nitrate reducers"]
    position = ["Immobile", "Mobile"]
    bins = [
        -1.1,
        -1,
        -0.9,
        -0.8,
        -0.7,
        -0.6,
        -0.5,
        -0.4,
        -0.3,
        -0.2,
        -0.1,
        0,
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1,
        1.1,
    ]
    centers = np.mean([bins[:-1], bins[1:]], axis=0)
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=[10, 8], sharex="col")
    for k in Chemseries:
        datac = data[(data["Chem"] == k)]
        ax1 = axes.flat[Chemseries.index(k)]
        ax1.bar(
            centers,
            np.histogram(datac[datac["Regime"] == "Slow"]["Crosscorrelation"], bins)[0],
            alpha=0.7,
            width=0.1,
            label="Slow flow",
            color="indianred",
        )
        ax1.bar(
            centers,
            np.histogram(datac[datac["Regime"] == "Medium"]["Crosscorrelation"], bins)[
                0
            ],
            alpha=0.7,
            width=0.1,
            label="Medium flow",
            color="g",
        )
        ax1.bar(
            centers,
            np.histogram(datac[datac["Regime"] == "Fast"]["Crosscorrelation"], bins)[0],
            alpha=0.7,
            width=0.1,
            label="Fast flow",
            color="steelblue",
        )
        ax1.set_xlim((-1, 1))
        ax1.tick_params(axis="y", labelsize=15)
        ax1.tick_params(axis="x", labelsize=15)
        ax1.set_xlabel("")
    plt.legend()
    for ax, typsp in zip(axes[0, :], species):
        ax.set_title(typsp, fontsize=15)
    axes[0, -1].annotate(
        position[0],
        xy=(0, 0.5),
        xytext=(180, 0),
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
        xytext=(180, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=15,
    )
    plt.annotate(
        "Number of scenarios",
        xy=(0, 1.8),
        xytext=(-450, -150),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="left",
        va="center",
        rotation="vertical",
        fontsize=18,
    )
    plt.annotate(
        "Peak cross-correlation between external forcing and\nresulting biomass in the domain",
        xy=(-0.4, -0.3),
        xytext=(-50, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=18,
    )

    return fig


def amplitude_head(data):
    Regimes = ["Slow", "Medium", "Fast"]
    fig, axes = plt.subplots(ncols=2, figsize=[11, 8], sharex=True)
    dfc = data[(data["Chem"] == "Head") & (data["Regime"] == "Fast")]
    l1 = axes.flat[0].scatter(
        "fraction",
        "Min",
        c="Time_series",
        cmap=colseries[Regimes.index("Fast")],
        data=dfc,
        label="Minimum",
        marker="^",
    )
    l2 = axes.flat[1].scatter(
        "fraction",
        "Max",
        c="Time_series",
        cmap=colseries[Regimes.index("Fast")],
        data=dfc,
        label="Maximum",
        marker="o",
    )
    axes.flat[0].tick_params(axis="y", labelsize=15)
    axes.flat[0].tick_params(axis="x", labelsize=15)
    axes.flat[1].tick_params(axis="y", labelsize=15)
    axes.flat[1].set_title("Minimum", fontsize=15)
    axes.flat[1].set_title("Maximum", fontsize=15)
    fig.legend(
        handles=[low_var, mid_var, high_var],
        bbox_to_anchor=(0.7, 0.4),
        loc="center",
        title="Variance in\nTime series",
        fontsize=12,
    )
    plt.annotate(
        "Normalized variation from steady state conditions",
        xy=(-1.2, 1),
        xytext=(-100, 0),
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
        xy=(0, 0),
        xytext=(0, -50),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=15,
    )

    return fig


def RMS_head(data):
    Regimes = ["Slow", "Medium", "Fast"]
    for Reg in Regimes:
        dfc = data[(data["Regime"] == Reg)]
        plt.scatter(
            "fraction",
            "Max_amplitude",
            c="Time_series",
            cmap=colseries[Regimes.index(Reg)],
            data=dfc,
            label=Reg,
            marker="^",
        )
        plt.tick_params(axis="y", labelsize=15)
        plt.tick_params(axis="x", labelsize=15)
        plt.ylim((-0.001, 0.017))
    first_legend = plt.legend(
        handles=[red_patch, green_patch, blue_patch],
        bbox_to_anchor=(1.2, 0.2),
        loc="center",
        title="Flow regime",
        fontsize=12,
    )
    ax = plt.gca().add_artist(first_legend)
    plt.legend(
        handles=[low_var, mid_var, high_var],
        bbox_to_anchor=(1.2, 0.7),
        loc="center",
        title="Variance in\nTime series",
        fontsize=12,
    )
    plt.ylabel("RMS Amplitude", fontsize=15)
    plt.xlabel("Fraction of breakthrough time in base case", fontsize=15)

    return plt


def correlation_delay_chem_2x2(data, Chemseries):
    Regimes = ["Slow", "Medium", "Fast"]
    fig1, axes1 = plt.subplots(nrows=2, ncols=2, figsize=[11, 8], sharex=True)
    plt.suptitle(
        "Change in flux averaged concentration with \nrespect to steady state conditions",
        fontsize=20,
    )
    for k in Chemseries:
        dfc = data[data["Chem"] == k]
        colidx1 = Chemseries.index(k)
        for i in Regimes:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp["Regime"] == i]
            paxis = axes1.flat[colidx1]
            l1 = axes1.flat[colidx1].scatter(
                "fraction",
                "Correlation",
                c="Time_series",
                cmap=colseries[Regimes.index(i)],
                data=dfcr,
                label="Correlation",
                marker="^",
            )
            paxis.tick_params(axis="y", labelsize=15)
            paxis.tick_params(axis="x", labelsize=15)
            paxis.set_title(k, fontsize=15)
    #            if(Chemseries.index(k) != len(Chemseries)-1):
    #                paxis.set_xticklabels([])
    #                if (Regimes.index(i)==len(Regimes)-2):
    #                    handles1 = [l1]
    #    fig1.legend(handles1, ["Correlation"], bbox_to_anchor=(1.7, 0.1),loc = 'center', title = "Value")
    fig1.legend(
        handles=[red_patch, green_patch, blue_patch],
        bbox_to_anchor=(0.6, 0.6),
        loc="center",
        title="Flow regime",
        fontsize=12,
    )
    fig1.legend(
        handles=[low_var, mid_var, high_var],
        bbox_to_anchor=(0.55, 0.3),
        loc="center",
        title="Variance in\nTime series",
        fontsize=12,
    )
    #    for ax,typsp in zip(axes1[:], Chemseries):
    #        ax.set_ylabel(typsp, fontsize = 15)
    plt.annotate(
        "Correlation",
        xy=(-1.1, 1),
        xytext=(-80, 60),
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
        xy=(0, 0),
        xytext=(0, -50),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=15,
    )

    fig2, axes2 = plt.subplots(nrows=2, ncols=2, figsize=[11, 8], sharex=True)
    plt.suptitle(
        "Change in flux averaged concentration with \nrespect to steady state conditions",
        fontsize=20,
    )
    for k in Chemseries:
        dfc = data[data["Chem"] == k]
        colidx1 = Chemseries.index(k)
        for i in Regimes:
            dfctemp = dfc
            dfcr = dfctemp[dfctemp["Regime"] == i]
            paxis = axes2.flat[colidx1]
            l1 = axes2.flat[colidx1].scatter(
                "fraction",
                "Delayfractionofbth",
                c="Time_series",
                cmap=colseries[Regimes.index(i)],
                data=dfcr,
                label="Correlation",
                marker="^",
            )
            paxis.tick_params(axis="y", labelsize=15)
            paxis.tick_params(axis="x", labelsize=15)
            paxis.set_ylim(
                np.min(dfc["Delayfractionofbth"]) - 0.1,
                np.max(dfc["Delayfractionofbth"]) + 0.1,
            )
            paxis.tick_params(axis="y", labelsize=15)
            paxis.set_title(k, fontsize=15)
    #            if(Chemseries.index(k) != len(Chemseries)-1):
    #                paxis.set_xticklabels([])
    #                if (Regimes.index(i)==len(Regimes)-2):
    #                    handles1 = [l1]
    #    fig2.legend(handles1, ["Correlation"], bbox_to_anchor=(1.7, 0.1),loc = 'center', title = "Value")
    fig2.legend(
        handles=[red_patch, green_patch, blue_patch],
        bbox_to_anchor=(0.75, 0.75),
        loc="center",
        title="Flow regime",
        fontsize=12,
    )
    fig2.legend(
        handles=[low_var, mid_var, high_var],
        bbox_to_anchor=(0.75, 0.3),
        loc="center",
        title="Variance in\nTime series",
        fontsize=12,
    )
    #    for ax,typsp in zip(axes2[:], Chemseries):
    #        ax.set_ylabel(typsp, fontsize = 15)
    plt.annotate(
        "Time lag with respect to fraction of breakthrough time",
        xy=(-1.1, 0.8),
        xytext=(-80, 60),
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
        xy=(0.1, -0.3),
        xytext=(0, 0),
        xycoords="axes fraction",
        textcoords="offset points",
        size="large",
        ha="center",
        va="center",
        fontsize=15,
    )

    return fig1, fig2
