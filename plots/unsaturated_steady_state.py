# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 10:28:17 2020

@author: khurana
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from matplotlib.colors import LogNorm
from analyses.unsaturated_steady_state import calcoxiccells
from analyses.unsaturated_transient import (
    calcconcmasstime,
    biomasstimefunc,
    calcconcmasstimeX,
)
from data_reader.data_processing import tracerstudies

def plotdataindex(data, index, timindex, ylabel, legend):
    plt.plot(np.mean(data[index, timindex, :, :], axis=-1), label=legend)
    plt.ylabel(ylabel)
    plt.xlabel("Y (cm)")


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
):
    legendsize = 16
    axissize = 16
    ticksize = 14
    titlesize = 20
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
            cvars,
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

    elif "Aerobes" in gvarnames:
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


def plotoxiccellssdo(
    limit, Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars
):
    axissize = 14
    ticksize = 12
    titlesize = 14
    oxiccells = calcoxiccells(
        limit, Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars
    )
    nrows = 2
    ncols = 4
    figsize = [14, 8]
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
        )
        axes.flat[count].plot(conctime[-1, :, 1], "r-")
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