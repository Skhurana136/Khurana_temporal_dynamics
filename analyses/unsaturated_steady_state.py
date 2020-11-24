# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:56:44 2020

@author: khurana
"""

import numpy as np


def calcmassflux(
    Trial,
    Het,
    Anis,
    gw,
    directory,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    vars,
    gvarnames,
):
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    por = 0.2
    doc1 = 10 - gw
    Bmo1 = 9 - gw
    Bmn1 = 16 - gw
    Bms1 = 21 - gw
    Bma1 = 26 - gw
    Bimo1 = 14 - gw
    Bimn1 = 19 - gw
    Bims1 = 24 - gw
    Bima1 = 28 - gw
    POM1 = 30 - gw
    Amm1 = 12 - gw
    nitra1 = 17 - gw
    Nspecies = [Bmo1, Bmn1, Bms1, Bma1, Bimo1, Bimn1, Bims1, Bima1, POM1]
    Cspecies = [doc1, Bmo1, Bmn1, Bms1, Bma1, Bimo1, Bimn1, Bims1, Bima1, POM1]
    mf = np.zeros([len(Trial) * len(gvarnames), 9])
    for j in range(len(Trial)):
        di = directory + fpre + str(Trial[j]) + fsuf
        print(str(Trial[j]))
        #        di+fpre+str(Tforfpre[k])+str(Trial[j])+'_df'
        df = np.load(di + fpre + str(Trial[j]) + "_df.npy")
        veliredg = df[2, -1, yin, xright]
        veliledg = df[2, -1, yin, xleft]
        veloredg = df[2, -1, yout, xright]
        veloledg = df[2, -1, yout, xleft]
        veloelem = df[2, -1, yout, xleft + 1 : xright]
        velielem = df[2, -1, yin, xleft + 1 : xright]
        velelem = df[2, -1, yin + 1 : yout, xleft + 1 : xright]
        vellelem = df[2, -1, yin + 1 : yout, xleft]
        velrelem = df[2, -1, yin + 1 : yout, xright]
        satielem = df[4, -1, yin, xleft + 1 : xright]
        satoelem = df[4, -1, yout, xleft + 1 : xright]
        satiredg = df[4, -1, yin, xright]
        satiledg = df[4, -1, yin, xright]
        satoledg = df[4, -1, yout, xleft]
        satoredg = df[4, -1, yout, xright]
        satlelem = df[4, -1, yin + 1 : yout, xleft]
        satrelem = df[4, -1, yin + 1 : yout, xright]
        satelem = df[4, -1, yin + 1 : yout, xleft + 1 : xright]
        for i in range(len(gvarnames)):
            idx = j * len(gvarnames) + i
            mf[idx, 0] = j
            mf[idx, 1] = Het[j]
            mf[idx, 2] = Anis[j]
            mf[idx, 3] = i
            if gvarnames[i] == "Nitrogen":
                ninlet = 0
                noutlet = 0
                for n in Nspecies:
                    ninlet = (
                        ninlet
                        + (
                            df[n - 3, -1, yin, xleft] * satiledg * veliledg * vedge
                            + df[n - 3, -1, yin, xright] * satiredg * veliredg * vedge
                            + sum(
                                df[n - 3, -1, yin, xleft + 1 : xright]
                                * satielem
                                * velielem
                                * velem
                            )
                        )
                        / por
                    )
                    noutlet = (
                        noutlet
                        + (
                            df[n - 3, -1, yout, xleft] * satoledg * veloledg * vedge
                            + df[n - 3, -1, yout, xright] * satoredg * veloredg * vedge
                            + sum(
                                df[n - 3, -1, yout, xleft + 1 : xright]
                                * satoelem
                                * veloelem
                                * velem
                            )
                        )
                        / por
                    )
                ninlet = (
                    ninlet / 10
                    + (
                        df[Amm1 - 3, -1, yin, xleft] * satiledg * veliledg * vedge
                        + df[Amm1 - 3, -1, yin, xright] * satiredg * veliredg * vedge
                        + sum(
                            df[Amm1 - 3, -1, yin, xleft + 1 : xright]
                            * satielem
                            * velielem
                            * velem
                        )
                        + df[nitra1 - 3, -1, yin, xleft] * satiledg * veliledg * vedge
                        + df[nitra1 - 3, -1, yin, xright] * satiredg * veliredg * vedge
                        + sum(
                            df[nitra1 - 3, -1, yin, xleft + 1 : xright]
                            * satielem
                            * velielem
                            * velem
                        )
                    )
                    / por
                )
                noutlet = (
                    noutlet / 10
                    + (
                        df[Amm1 - 3, -1, yout, xleft] * satoledg * veloledg * vedge
                        + df[Amm1 - 3, -1, yout, xright] * satoredg * veloredg * vedge
                        + sum(
                            df[Amm1 - 3, -1, yout, xleft + 1 : xright]
                            * satoelem
                            * veloelem
                            * velem
                        )
                        + df[nitra1 - 3, -1, yout, xleft] * satoledg * veloledg * vedge
                        + df[nitra1 - 3, -1, yout, xright] * satoredg * veloredg * vedge
                        + sum(
                            df[nitra1 - 3, -1, yout, xleft + 1 : xright]
                            * satoelem
                            * veloelem
                            * velem
                        )
                    )
                    / por
                )
                mf[idx, 4] = abs(
                    ninlet
                )  # /(sum(velielem)*velem + (veliledg+veliredg)*vedge)
                mf[idx, 5] = abs(
                    noutlet
                )  # /(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
            elif gvarnames[i] == "TOC":
                cinlet = 0
                coutlet = 0
                for c in Cspecies:
                    cinlet = (
                        cinlet
                        + (
                            df[c - 3, -1, yin, xleft] * satiledg * veliledg * vedge
                            + df[c - 3, -1, yin, xright] * satiredg * veliredg * vedge
                            + sum(
                                df[c - 3, -1, yin, xleft + 1 : xright]
                                * satielem
                                * velielem
                                * velem
                            )
                        )
                        / por
                    )
                    coutlet = (
                        coutlet
                        + (
                            df[c - 3, -1, yout, xleft] * satoledg * veloledg * vedge
                            + df[c - 3, -1, yout, xright] * satoredg * veloredg * vedge
                            + sum(
                                df[c - 3, -1, yout, xleft + 1 : xright]
                                * satoelem
                                * veloelem
                                * velem
                            )
                        )
                        / por
                    )
                mf[idx, 4] = abs(
                    cinlet
                )  # /(sum(velielem)*velem + (veliledg+veliredg)*vedge)
                mf[idx, 5] = abs(
                    coutlet
                )  # /(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
            else:
                mf[idx, 4] = abs(
                    (
                        (
                            df[vars[i] - 3, -1, yin, xleft]
                            * satiledg
                            * veliledg
                            * vedge
                            + df[vars[i] - 3, -1, yin, xright]
                            * satiredg
                            * veliredg
                            * vedge
                            + sum(
                                df[vars[i] - 3, -1, yin, xleft + 1 : xright]
                                * satielem
                                * velielem
                                * velem
                            )
                        )
                    )
                    / por
                )  # /(sum(velielem)*velem + (veliledg+veliredg)*vedge)
                mf[idx, 5] = abs(
                    (
                        (
                            df[vars[i] - 3, -1, yout, xleft]
                            * satoledg
                            * veloledg
                            * vedge
                            + df[vars[i] - 3, -1, yout, xright]
                            * satoredg
                            * veloredg
                            * vedge
                            + sum(
                                df[vars[i] - 3, -1, yout, xleft + 1 : xright]
                                * satoelem
                                * veloelem
                                * velem
                            )
                        )
                    )
                    / por
                )  # /(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
            # calculating removal in mass
            mf[idx, 6] = mf[idx, 4] - mf[idx, 5]
            # normalizing removal in mass with incoming mass
            mf[idx, 7] = (mf[idx, 4] - mf[idx, 5]) / mf[idx, 4]
            # comparing removal with homogeneous case
            mf[idx, 8] = mf[idx, 6] / mf[i, 6]
    return mf


def calcsum(
    Trial,
    Het,
    Anis,
    gw,
    directory,
    fpre,
    fsuf,
    yin,
    yout,
    xleft,
    xright,
    vars,
    gvarnames,
):
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    por = 0.2
    sumall = np.zeros([len(Trial) * len(vars), 9])
    for j in range(len(Trial)):
        di = directory + fpre + str(Trial[j]) + fsuf
        print(str(Trial[j]))
        df = np.load(di + fpre + str(Trial[j]) + "_df.npy")
        satielem = df[4, np.shape(df)[1] - 1, yin, xleft + 1 : xright]
        satoelem = df[4, np.shape(df)[1] - 1, yout, xleft + 1 : xright]
        satiredg = df[4, np.shape(df)[1] - 1, yin, xright]
        satiledg = df[4, np.shape(df)[1] - 1, yin, xright]
        satoledg = df[4, np.shape(df)[1] - 1, yout, xleft]
        satoredg = df[4, np.shape(df)[1] - 1, yout, xright]
        satlelem = df[4, np.shape(df)[1] - 1, yin + 1 : yout, xleft]
        satrelem = df[4, np.shape(df)[1] - 1, yin + 1 : yout, xright]
        satelem = df[4, np.shape(df)[1] - 1, yin + 1 : yout, xleft + 1 : xright]
        for i in range(len(gvarnames)):
            idx = j * len(gvarnames) + i
            sumall[idx, 0] = Trial[1] + j
            sumall[idx, 1] = Het[j]
            sumall[idx, 2] = Anis[j]
            sumall[idx, 3] = i
            sumall[idx, 4] = (
                (
                    (
                        df[vars[i] - 3, np.shape(df)[1] - 1, yin, xleft] * satiledg
                        + df[vars[i] - 3, np.shape(df)[1] - 1, yout, xleft] * satoledg
                        + df[vars[i] - 3, np.shape(df)[1] - 1, yin, xright] * satiredg
                        + df[vars[i] - 3, np.shape(df)[1] - 1, yout, xright] * satoredg
                    )
                    * vedge
                    * vedge
                    + (
                        sum(
                            df[
                                vars[i] - 3,
                                np.shape(df)[1] - 1,
                                yin,
                                xleft + 1 : xright,
                            ]
                            * satielem
                        )
                        + sum(
                            df[
                                vars[i] - 3,
                                np.shape(df)[1] - 1,
                                yout,
                                xleft + 1 : xright,
                            ]
                            * satoelem
                        )
                        + sum(
                            df[vars[i] - 3, np.shape(df)[1] - 1, yin + 1 : yout, xleft]
                            * satlelem
                        )
                        + sum(
                            df[vars[i] - 3, np.shape(df)[1] - 1, yin + 1 : yout, xright]
                            * satrelem
                        )
                    )
                    * vedge
                    * velem
                )
                + sum(
                    sum(
                        df[
                            vars[i] - 3,
                            np.shape(df)[1] - 1,
                            yin + 1 : yout,
                            xleft + 1 : xright,
                        ]
                        * satelem
                    )
                )
                * velem
                * velem
            ) * por
            sumall[idx, 5] = (sumall[idx, 4]) / sumall[i, 4]
        total = np.sum(sumall[idx - 11 : idx + 1, 4])
        for k in range(len(gvarnames)):
            idxk = j * len(gvarnames) + k
            sumall[idxk, 6] = sumall[idxk, 4] / total
            # comparing the fraction of species with homogeneous case in the uniform flow rate scenario
            sumall[idxk, 7] = sumall[idxk, 6] / sumall[k, 6]
            # comparing the fraction of species with homogeneous case in the varying flow rate scenario
            sumall[idxk, 8] = sumall[idxk, 6] / sumall[j * len(gvarnames) + k, 6]
    return sumall


def calcoxiccells(
    limit, Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars
):
    oxiccells = np.zeros([len(Trial), 51, 1])
    for j in range(len(Trial)):
        df, massendtime, masstime, conctime, Velocity, head = calcconcmasstime(
            Trial[j], Het[j], Anis[j], gw, d, fpre, fsuf, yin, yout, xleft, xright, vars
        )
        c = []
        for k in range(51):
            c.append(
                np.count_nonzero(df[vars[1] - 3, np.shape(df)[1] - 1, k, :] > limit)
            )
        #        print (Trial[j], c)
        #        print ("Calculating number of oxic cells: ", c)
        oxiccells[j, :, 0] = c

    return oxiccells
