# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:57:00 2020

@author: khurana
"""
import numpy as np

def calcconcmasstime(
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
    di = directory + fpre + str(Trial) + fsuf
    print(str(Trial))
    df = np.load(di + fpre + str(Trial) + "_df.npy")
    conctime = np.zeros([np.shape(df)[1], 51, len(gvarnames)])
    mftime = np.zeros([np.shape(df)[1], 51, len(gvarnames)])
    veliredg = df[2, 1:, yin, xright]
    veliledg = df[2, 1:, yin, xleft]
    veloredg = df[2, 1:, yout, xright]
    veloledg = df[2, 1:, yout, xleft]
    veloelem = df[2, 1:, yout, xleft + 1 : xright]
    velielem = df[2, 1:, yin, xleft + 1 : xright]
    velelem = df[2, 1:, yin + 1 : yout, xleft + 1 : xright]
    vellelem = df[2, 1:, yin + 1 : yout, xleft]
    velrelem = df[2, 1:, yin + 1 : yout, xright]
    satielem = df[4, 1:, yin, xleft + 1 : xright]
    satoelem = df[4, 1:, yout, xleft + 1 : xright]
    satlelem = df[4, 1:, yin + 1 : yout, xleft]
    satrelem = df[4, 1:, yin + 1 : yout, xright]
    satiredg = df[4, 1:, yin, xright]
    satiledg = df[4, 1:, yin, xleft]
    satoledg = df[4, 1:, yout, xleft]
    satoredg = df[4, 1:, yout, xright]
    satelem = df[4, 1:, yin + 1 : yout, xleft + 1 : xright]
    for i in range(len(gvarnames)):
        if gvarnames[i] == "Nitrogen":
            ninlet = 0
            noutlet = 0
            nminlet = 0
            nmoutlet = 0
            for n in Nspecies:
                ninlet = ninlet + (
                    df[n - 3, 1:, yin, xleft] * satiledg * veliledg * vedge
                    + df[n - 3, 1:, yin, xright] * satiredg * veliredg * vedge
                    + np.sum(
                        df[n - 3, 1:, yin, xleft + 1 : xright]
                        * satielem
                        * velielem
                        * velem,
                        axis=-1,
                    )
                ) / (
                    vedge * (satiredg * veliredg + satiledg * veliledg)
                    + np.sum(velem * satielem * velielem, axis=-1)
                )
                noutlet = noutlet + (
                    df[n - 3, 1:, yout, xleft] * satoledg * veloledg * vedge
                    + df[n - 3, 1:, yout, xright] * satoredg * veloredg * vedge
                    + np.sum(
                        df[n - 3, 1:, yout, xleft + 1 : xright]
                        * satoelem
                        * veloelem
                        * velem,
                        axis=-1,
                    )
                ) / (
                    vedge * (satoredg * veloredg + satoledg * veloledg)
                    + np.sum(velem * satoelem * veloelem, axis=-1)
                )
                nminlet = nminlet + (
                    df[n - 3, 1:, yin, xleft] * satiledg * veliledg * vedge
                    + df[n - 3, 1:, yin, xright] * satiredg * veliredg * vedge
                    + np.sum(
                        df[n - 3, 1:, yin, xleft + 1 : xright]
                        * satielem
                        * velielem
                        * velem,
                        axis=-1,
                    )
                )  # /(vedge*(satiredg*veliredg+satiledg*veliledg) + np.sum(velem*satielem*velielem, axis = -1))
                nmoutlet = nmoutlet + (
                    df[n - 3, 1:, yout, xleft] * satoledg * veloledg * vedge
                    + df[n - 3, 1:, yout, xright] * satoredg * veloredg * vedge
                    + np.sum(
                        df[n - 3, 1:, yout, xleft + 1 : xright]
                        * satoelem
                        * veloelem
                        * velem,
                        axis=-1,
                    )
                )  # /(vedge*(satoredg*veloredg+satoledg*veloledg) + np.sum(velem*satoelem*veloelem, axis = -1))
            conctime[1:, yin, i] = ninlet / 10 + (
                df[Amm1 - 3, 1:, yin, xleft] * satiledg * veliledg * vedge
                + df[Amm1 - 3, 1:, yin, xright] * satiredg * veliredg * vedge
                + np.sum(
                    df[Amm1 - 3, 1:, yin, xleft + 1 : xright]
                    * satielem
                    * velielem
                    * velem,
                    axis=-1,
                )
                + df[nitra1 - 3, 1:, yin, xleft] * satiledg * veliledg * vedge
                + df[nitra1 - 3, 1:, yin, xright] * satiredg * veliredg * vedge
                + np.sum(
                    df[nitra1 - 3, 1:, yin, xleft + 1 : xright]
                    * satielem
                    * velielem
                    * velem,
                    axis=-1,
                )
            ) / (
                vedge * (satiredg * veliredg + satiledg * veliledg)
                + np.sum(velem * satielem * velielem, axis=-1)
            )
            conctime[1:, yout, i] = noutlet / 10 + (
                df[Amm1 - 3, 1:, yout, xleft] * satoledg * veloledg * vedge
                + df[Amm1 - 3, 1:, yout, xright] * satoredg * veloredg * vedge
                + np.sum(
                    df[Amm1 - 3, 1:, yout, xleft + 1 : xright]
                    * satoelem
                    * veloelem
                    * velem,
                    axis=-1,
                )
                + df[nitra1 - 3, 1:, yout, xleft] * satoledg * veloledg * vedge
                + df[nitra1 - 3, 1:, yout, xright] * satoredg * veloredg * vedge
                + np.sum(
                    df[nitra1 - 3, 1:, yout, xleft + 1 : xright]
                    * satoelem
                    * veloelem
                    * velem,
                    axis=-1,
                )
            ) / (
                vedge * (satoredg * veloredg + satoledg * veloledg)
                + np.sum(velem * satoelem * veloelem, axis=-1)
            )
            mftime[1:, yin, i] = np.abs(
                ninlet / 10
                + (
                    df[Amm1 - 3, 1:, yin, xleft] * satiledg * veliledg * vedge
                    + df[Amm1 - 3, 1:, yin, xright] * satiredg * veliredg * vedge
                    + np.sum(
                        df[Amm1 - 3, 1:, yin, xleft + 1 : xright]
                        * satielem
                        * velielem
                        * velem,
                        axis=-1,
                    )
                    + df[nitra1 - 3, 1:, yin, xleft] * satiledg * veliledg * vedge
                    + df[nitra1 - 3, 1:, yin, xright] * satiredg * veliredg * vedge
                    + np.sum(
                        df[nitra1 - 3, 1:, yin, xleft + 1 : xright]
                        * satielem
                        * velielem
                        * velem,
                        axis=-1,
                    )
                )
            )  # /(vedge*(satiredg*veliredg+satiledg*veliledg) + np.sum(velem*satielem*velielem, axis = -1))
            mftime[1:, yout, i] = np.abs(
                noutlet / 10
                + (
                    df[Amm1 - 3, 1:, yout, xleft] * satoledg * veloledg * vedge
                    + df[Amm1 - 3, 1:, yout, xright] * satoredg * veloredg * vedge
                    + np.sum(
                        df[Amm1 - 3, 1:, yout, xleft + 1 : xright]
                        * satoelem
                        * veloelem
                        * velem,
                        axis=-1,
                    )
                    + df[nitra1 - 3, 1:, yout, xleft] * satoledg * veloledg * vedge
                    + df[nitra1 - 3, 1:, yout, xright] * satoredg * veloredg * vedge
                    + np.sum(
                        df[nitra1 - 3, 1:, yout, xleft + 1 : xright]
                        * satoelem
                        * veloelem
                        * velem,
                        axis=-1,
                    )
                )
            )  # /(vedge*(satoredg*veloredg+satoledg*veloledg) + np.sum(velem*satoelem*veloelem, axis = -1))
        elif gvarnames[i] == "TOC":
            cinlet = 0
            coutlet = 0
            cminlet = 0
            cmoutlet = 0
            for c in Cspecies:
                cinlet = cinlet + (
                    df[c - 3, 1:, yin, xleft] * satiledg * veliledg * vedge
                    + df[c - 3, 1:, yin, xright] * satiredg * veliredg * vedge
                    + np.sum(
                        df[c - 3, 1:, yin, xleft + 1 : xright]
                        * satielem
                        * velielem
                        * velem,
                        axis=-1,
                    )
                ) / (
                    vedge * (satiredg * veliredg + satiledg * veliledg)
                    + np.sum(velem * satielem * velielem, axis=-1)
                )
                coutlet = coutlet + (
                    df[c - 3, 1:, yout, xleft] * satoledg * veloledg * vedge
                    + df[c - 3, 1:, yout, xright] * satoredg * veloredg * vedge
                    + np.sum(
                        df[c - 3, 1:, yout, xleft + 1 : xright]
                        * satoelem
                        * veloelem
                        * velem,
                        axis=-1,
                    )
                ) / (
                    vedge * (satoredg * veloredg + satoledg * veloledg)
                    + np.sum(satoelem * velem * veloelem, axis=-1)
                )
                cminlet = cminlet + (
                    df[c - 3, 1:, yin, xleft] * satiledg * veliledg * vedge
                    + df[c - 3, 1:, yin, xright] * satiredg * veliredg * vedge
                    + np.sum(
                        df[c - 3, 1:, yin, xleft + 1 : xright]
                        * satielem
                        * velielem
                        * velem,
                        axis=-1,
                    )
                )  # /(vedge*(satiredg*veliredg+satiledg*veliledg) + np.sum(velem*satielem*velielem, axis = -1))
                cmoutlet = cmoutlet + (
                    df[c - 3, 1:, yout, xleft] * satoledg * veloledg * vedge
                    + df[c - 3, 1:, yout, xright] * satoredg * veloredg * vedge
                    + np.sum(
                        df[c - 3, 1:, yout, xleft + 1 : xright]
                        * satoelem
                        * veloelem
                        * velem,
                        axis=-1,
                    )
                )  # /(vedge*(satoredg*veloredg+satoledg*veloledg) + np.sum(satoelem*velem*veloelem, axis = -1))
            conctime[
                1:, yin, i
            ] = cinlet  # /(sum(velielem)*velem + (veliledg+veliredg)*vedge)
            conctime[
                1:, yout, i
            ] = coutlet  # /(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
            mftime[1:, yin, i] = np.abs(cminlet)
            mftime[1:, yout, i] = np.abs(coutlet)
        else:
            conctime[1:, yin, i] = (
                (
                    df[vars[i] - 3, 1:, yin, xleft] * satiledg * veliledg
                    + df[vars[i] - 3, 1:, yin, xright] * satiredg * veliredg
                )
                * (vedge)
                + (
                    np.sum(
                        df[vars[i] - 3, 1:, yin, xleft + 1 : xright]
                        * satielem
                        * velielem,
                        axis=-1,
                    )
                )
                * velem
            ) / (
                vedge * satiledg * (veliredg + veliledg)
                + np.sum(satielem * velem * velielem, axis=-1)
            )
            conctime[1:, yout, i] = (
                (
                    df[vars[i] - 3, 1:, yout, xleft] * satoledg * veloledg
                    + df[vars[i] - 3, 1:, yout, xright] * satoredg * veloredg
                )
                * (vedge)
                + (
                    np.sum(
                        df[vars[i] - 3, 1:, yout, xleft + 1 : xright]
                        * satoelem
                        * veloelem,
                        axis=-1,
                    )
                )
                * velem
            ) / (
                vedge * satoledg * (veloredg + veloledg)
                + np.sum(satoelem * velem * veloelem, axis=-1)
            )
            conctime[1:, yin + 1 : yout, i] = (
                np.sum(
                    df[vars[i] - 3, 1:, yin + 1 : yout, xleft + 1 : xright]
                    * satelem
                    * velem
                    * velelem,
                    axis=-1,
                )
                + (
                    df[vars[i] - 3, 1:, yin + 1 : yout, xleft] * satlelem * vellelem
                    + df[vars[i] - 3, 1:, yin + 1 : yout, xright] * satrelem * velrelem
                )
                * vedge
            ) / (
                vedge * (vellelem * satlelem + velrelem * satrelem)
                + np.sum(velem * velelem * satelem, axis=-1)
            )
            mftime[1:, yin, i] = np.abs(
                (
                    df[vars[i] - 3, 1:, yin, xleft] * satiledg * veliledg
                    + df[vars[i] - 3, 1:, yin, xright] * satiredg * veliredg
                )
                * (vedge)
                + (
                    np.sum(
                        df[vars[i] - 3, 1:, yin, xleft + 1 : xright]
                        * satielem
                        * velielem,
                        axis=-1,
                    )
                )
                * velem
            )  # /(vedge*satiledg*(veliredg+veliledg) + np.sum(satielem*velem*velielem, axis = -1))
            mftime[1:, yout, i] = np.abs(
                (
                    df[vars[i] - 3, 1:, yout, xleft] * satoledg * veloledg
                    + df[vars[i] - 3, 1:, yout, xright] * satoredg * veloredg
                )
                * (vedge)
                + (
                    np.sum(
                        df[vars[i] - 3, 1:, yout, xleft + 1 : xright]
                        * satoelem
                        * veloelem,
                        axis=-1,
                    )
                )
                * velem
            )  # /(vedge*satoledg*(veloredg+veloledg) + np.sum(satoelem*velem*veloelem, axis = -1))
            mftime[1:, yin + 1 : yout, i] = np.abs(
                np.sum(
                    df[vars[i] - 3, 1:, yin + 1 : yout, xleft + 1 : xright]
                    * satelem
                    * velem
                    * velelem,
                    axis=-1,
                )
                + (
                    df[vars[i] - 3, 1:, yin + 1 : yout, xleft] * satlelem * vellelem
                    + df[vars[i] - 3, 1:, yin + 1 : yout, xright] * satrelem * velrelem
                )
                * vedge
            )  # /(vedge*(vellelem*satlelem+velrelem*satrelem) + np.sum(velem*velelem*satelem, axis = -1))
    TotalFlow = (veliledg + veloledg + veliredg + veloredg) * vedge + (
        np.sum(vellelem)
        + np.sum(velrelem)
        + np.sum(velelem)
        + np.sum(velielem)
        + np.sum(veloelem)
    ) * velem
    #    Velocity = np.mean ([InVelocity, OutVelocity, MidVelocity])
    Velocity = df[2, np.shape(df)[1] - 1, :, :]
    Headinlettime = np.mean(df[2, 1:, yin, :], axis=-1) * -1
    return df, conctime, mftime, np.mean(Velocity), Headinlettime


def biomasstimefunc(
    Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, biomassvars
):
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    biomassendtime = np.zeros([len(biomassvars)])
    biomassendtimey = np.zeros([51, len(biomassvars)])
    biomassendtimey[:, 0] = range(51)
    di = d + fpre + str(Trial) + fsuf
    print(str(Trial))
    df = np.load(di + fpre + str(Trial) + "_df.npy")
    biomasstime = np.zeros([np.shape(df)[1] - 1, 51, len(biomassvars)])
    bioconctime = np.zeros([np.shape(df)[1] - 1, 51, len(biomassvars)])
    satielem = df[4, 1:, yin, xleft + 1 : xright]
    satoelem = df[4, 1:, yout, xleft + 1 : xright]
    satlelem = df[4, 1:, yin + 1 : yout, xleft]
    satrelem = df[4, 1:, yin + 1 : yout, xright]
    satiredg = df[4, 1:, yin, xright]
    satiledg = df[4, 1:, yin, xleft]
    satoledg = df[4, 1:, yout, xleft]
    satoredg = df[4, 1:, yout, xright]
    satelem = df[4, 1:, yin + 1 : yout, xleft + 1 : xright]
    for i in range(len(biomassvars)):
        biomassendtime[i] = (
            (
                df[biomassvars[i] - 3, np.shape(df)[1] - 1, yin, xleft] * satiledg
                + df[biomassvars[i] - 3, np.shape(df)[1] - 1, yin, xright] * satiredg
                + df[biomassvars[i] - 3, np.shape(df)[1] - 1, yout, xleft] * satoledg
                + df[biomassvars[i] - 3, np.shape(df)[1] - 1, yout, xright] * satoredg
            )
            * (vedge ** 2)
            + sum(
                sum(
                    df[
                        biomassvars[i] - 3,
                        np.shape(df)[1] - 1,
                        yin + 1 : yout - 1,
                        xleft + 1 : xright - 1,
                    ]
                    * satelem
                    * (velem ** 2)
                )
            )
            + (
                sum(
                    df[
                        biomassvars[i] - 3,
                        np.shape(df)[1] - 1,
                        yin,
                        xleft + 1 : xright - 1,
                    ]
                    * satielem
                )
                + sum(
                    df[
                        biomassvars[i] - 3,
                        np.shape(df)[1] - 1,
                        yout,
                        xleft + 1 : xright - 1,
                    ]
                    * satoelem
                )
                + sum(
                    df[
                        biomassvars[i] - 3,
                        np.shape(df)[1] - 1,
                        yin + 1 : yout - 1,
                        xleft,
                    ]
                    * satlelem
                )
                + sum(
                    df[
                        biomassvars[i] - 3,
                        np.shape(df)[1] - 1,
                        yin + 1 : yout - 1,
                        xright,
                    ]
                    * satrelem
                )
            )
            * velem
            * vedge
        )
        biomassendtimey[yin, i] = (
            (
                df[biomassvars[i] - 3, np.shape(df)[1] - 1, yin, xleft] * satiledg
                + df[biomassvars[i] - 3, np.shape(df)[1] - 1, yin, xright] * satiredg
            )
            * (vedge ** 2)
            + (
                sum(
                    df[biomassvars[i] - 3, np.shape(df)[1] - 1, yin, xleft + 1 : xright]
                    * satielem
                )
            )
            * velem
            * vedge
        )
        biomassendtimey[yout, i] = (
            (
                df[biomassvars[i] - 3, np.shape(df)[1] - 1, yout, xleft] * satoledg
                + df[biomassvars[i] - 3, np.shape(df)[1] - 1, yout, xright] * satoredg
            )
            * (vedge ** 2)
            + (
                sum(
                    df[
                        biomassvars[i] - 3,
                        np.shape(df)[1] - 1,
                        yout,
                        xleft + 1 : xright,
                    ]
                    * satoelem
                )
            )
            * velem
            * vedge
        )
        biomassendtimey[yin + 1 : yout, i] = (
            sum(
                sum(
                    df[
                        biomassvars[i] - 3,
                        np.shape(df)[1] - 1,
                        yin + 1 : yout,
                        xleft + 1 : xright,
                    ]
                    * satelem
                    * (velem ** 2)
                )
            )
            + (
                sum(
                    df[biomassvars[i] - 3, np.shape(df)[1] - 1, yin + 1 : yout, xleft]
                    * satlelem
                )
                + sum(
                    df[biomassvars[i] - 3, np.shape(df)[1] - 1, yin + 1 : yout, xright]
                    * satrelem
                )
            )
            * velem
            * vedge
        )
        biomasstime[:, yin, i] = (
            (
                df[biomassvars[i] - 3, 1:, yin, xleft] * satiledg
                + df[biomassvars[i] - 3, 1:, yin, xright] * satiredg
            )
            * (vedge ** 2)
            + (
                np.sum(
                    df[biomassvars[i] - 3, 1:, yin, xleft + 1 : xright] * satielem,
                    axis=-1,
                )
            )
            * velem
            * vedge
        )
        biomasstime[:, yout, i] = (
            (
                df[biomassvars[i] - 3, 1:, yout, xleft] * satoledg
                + df[biomassvars[i] - 3, 1:, yout, xright] * satoredg
            )
            * (vedge ** 2)
            + (
                np.sum(
                    df[biomassvars[i] - 3, 1:, yout, xleft + 1 : xright] * satoelem,
                    axis=-1,
                )
            )
            * velem
            * vedge
        )
        biomasstime[:, yin + 1 : yout, i] = (
            np.sum(
                df[biomassvars[i] - 3, 1:, yin + 1 : yout, xleft + 1 : xright]
                * satelem
                * (velem ** 2),
                axis=-1,
            )
            + (
                (
                    df[biomassvars[i] - 3, 1:, yin + 1 : yout, xleft] * satlelem
                    + df[biomassvars[i] - 3, 1:, yin + 1 : yout, xright]
                )
                * satrelem
            )
            * velem
            * vedge
        )
        bioconctime[:, yin, i] = (
            (
                df[biomassvars[i] - 3, 1:, yin, xleft] * satiledg
                + df[biomassvars[i] - 3, 1:, yin, xright] * satiredg
            )
            * (vedge ** 2)
            + (
                np.sum(
                    df[biomassvars[i] - 3, 1:, yin, xleft + 1 : xright - 1] * satielem,
                    axis=-1,
                )
            )
            * velem
            * vedge
        ) / (vbc * vedge)
        bioconctime[:, yout, i] = (
            (
                df[biomassvars[i] - 3, 1:, yout, xleft] * satoledg
                + df[biomassvars[i] - 3, 1:, yout, xright] * satoredg
            )
            * (vedge ** 2)
            + (
                np.sum(
                    df[biomassvars[i] - 3, 1:, yout, xleft + 1 : xright] * satoelem,
                    axis=-1,
                )
            )
            * velem
            * vedge
        ) / (vbc * vedge)
        bioconctime[:, yin + 1 : yout, i] = (
            np.sum(
                df[biomassvars[i] - 3, 1:, yin + 1 : yout, xleft + 1 : xright]
                * satelem
                * (velem ** 2),
                axis=-1,
            )
            + (
                (
                    df[biomassvars[i] - 3, 1:, yin + 1 : yout, xleft] * satlelem
                    + df[biomassvars[i] - 3, 1:, yin + 1 : yout, xright]
                )
                * satrelem
            )
            * velem
            * vedge
        ) / (vbc * velem)

    for i in range(len(biomassvars)):
        biomassendtime[i] = (
            (
                df[biomassvars[i] - 3, np.shape(df)[1] - 2, yin, xleft]
                * satiledg[int(np.shape(df)[1]) - 2]
                + df[biomassvars[i] - 3, np.shape(df)[1] - 2, yin, xright]
                * satiredg[int(np.shape(df)[1]) - 2]
                + df[biomassvars[i] - 3, np.shape(df)[1] - 2, yout, xleft]
                * satoledg[int(np.shape(df)[1]) - 2]
                + df[biomassvars[i] - 3, np.shape(df)[1] - 2, yout, xright]
                * satoredg[int(np.shape(df)[1]) - 2]
            )
            * (vedge ** 2)
            + sum(
                sum(
                    df[
                        biomassvars[i] - 3,
                        np.shape(df)[1] - 2,
                        yin + 1 : yout,
                        xleft + 1 : xright,
                    ]
                    * satelem[int(np.shape(df)[1]) - 2, :, :]
                    * (velem ** 2)
                )
            )
            + (
                sum(
                    df[biomassvars[i] - 3, np.shape(df)[1] - 2, yin, xleft + 1 : xright]
                    * satielem[np.shape(df)[1] - 2, :]
                )
                + sum(
                    df[
                        biomassvars[i] - 3,
                        np.shape(df)[1] - 2,
                        yout,
                        xleft + 1 : xright,
                    ]
                    * satoelem[np.shape(df)[1] - 2, :]
                )
                + sum(
                    df[biomassvars[i] - 3, np.shape(df)[1] - 2, yin + 1 : yout, xleft]
                    * satlelem[np.shape(df)[1] - 2, :]
                )
                + sum(
                    df[biomassvars[i] - 3, np.shape(df)[1] - 2, yin + 1 : yout, xright]
                    * satrelem[np.shape(df)[1] - 2, :]
                )
            )
            * velem
            * vedge
        )
        biomassendtimey[yin, i] = (
            (
                df[biomassvars[i] - 3, np.shape(df)[1] - 2, yin, xleft]
                * satiledg[np.shape(df)[1] - 2]
                + df[biomassvars[i] - 3, np.shape(df)[1] - 2, yin, xright]
                * satiredg[np.shape(df)[1] - 2]
            )
            * (vedge ** 2)
            + (
                sum(
                    df[biomassvars[i] - 3, np.shape(df)[1] - 2, yin, xleft + 1 : xright]
                    * satielem[np.shape(df)[1] - 2, :]
                )
            )
            * velem
            * vedge
        )
        biomassendtimey[yout, i] = (
            (
                df[biomassvars[i] - 3, np.shape(df)[1] - 2, yout, xleft]
                * satoledg[np.shape(df)[1] - 2]
                + df[biomassvars[i] - 3, np.shape(df)[1] - 2, yout, xright]
                * satoredg[np.shape(df)[1] - 2]
            )
            * (vedge ** 2)
            + (
                sum(
                    df[
                        biomassvars[i] - 3,
                        np.shape(df)[1] - 2,
                        yout,
                        xleft + 1 : xright,
                    ]
                    * satoelem[np.shape(df)[1] - 2, :]
                )
            )
            * velem
            * vedge
        )
        biomassendtimey[yin + 1 : yout, i] = (
            sum(
                sum(
                    df[
                        biomassvars[i] - 3,
                        np.shape(df)[1] - 2,
                        yin + 1 : yout,
                        xleft + 1 : xright,
                    ]
                    * satelem[np.shape(df)[1] - 2, :, :]
                    * (velem ** 2)
                )
            )
            + (
                sum(
                    df[biomassvars[i] - 3, np.shape(df)[1] - 2, yin + 1 : yout, xleft]
                    * satlelem[np.shape(df)[1] - 2, :]
                )
                + sum(
                    df[biomassvars[i] - 3, np.shape(df)[1] - 2, yin + 1 : yout, xright]
                    * satrelem[np.shape(df)[1] - 2, :]
                )
            )
            * velem
            * vedge
        )
        biomasstime[:, yin, i] = (
            (
                df[biomassvars[i] - 3, 1:, yin, xleft] * satiledg
                + df[vars[i] - 3, 1:, yin, xright] * satiredg
            )
            * (vedge ** 2)
            + (
                np.sum(
                    df[biomassvars[i] - 3, 1:, yin, xleft + 1 : xright] * satielem,
                    axis=-1,
                )
            )
            * velem
            * vedge
        )
        biomasstime[:, yout, i] = (
            (
                df[biomassvars[i] - 3, 1:, yout, xleft] * satoledg
                + df[vars[i] - 3, 1:, yout, xright] * satoredg
            )
            * (vedge ** 2)
            + (
                np.sum(
                    df[biomassvars[i] - 3, 1:, yout, xleft + 1 : xright] * satoelem,
                    axis=-1,
                )
            )
            * velem
            * vedge
        )
        biomasstime[:, yin + 1 : yout, i] = (
            np.sum(
                df[biomassvars[i] - 3, 1:, yin + 1 : yout, xleft + 1 : xright]
                * satelem
                * (velem ** 2),
                axis=-1,
            )
            + (
                df[biomassvars[i] - 3, 1:, yin + 1 : yout, xleft] * satlelem
                + df[biomassvars[i] - 3, 1:, yin + 1 : yout, xright] * satrelem
            )
            * velem
            * vedge
        )
        bioconctime[:, yin, i] = (
            (
                df[biomassvars[i] - 3, 1:, yin, xleft] * satiledg
                + df[vars[i] - 3, 1:, yin, xright] * satiredg
            )
            * (vedge ** 2)
            + np.sum(
                df[biomassvars[i] - 3, 1:, yin, xleft + 1 : xright] * satielem, axis=-1
            )
            * velem
            * vedge
        ) / (vbc * vedge)
        bioconctime[:, yout, i] = (
            (
                df[biomassvars[i] - 3, 1:, yout, xleft] * satoledg
                + df[vars[i] - 3, 1:, yout, xright] * satoredg
            )
            * (vedge ** 2)
            + (
                np.sum(
                    df[biomassvars[i] - 3, 1:, yout, xleft + 1 : xright] * satoelem,
                    axis=-1,
                )
            )
            * velem
            * vedge
        ) / (vbc * vedge)
        bioconctime[:, yin + 1 : yout, i] = (
            np.sum(
                df[biomassvars[i] - 3, 1:, yin + 1 : yout, xleft + 1 : xright]
                * satelem
                * (velem ** 2),
                axis=-1,
            )
            + (
                df[biomassvars[i] - 3, 1:, yin + 1 : yout, xleft] * satlelem
                + df[biomassvars[i] - 3, 1:, yin + 1 : yout, xright] * satrelem
            )
            * velem
            * vedge
        ) / (vbc * velem)
    return df, biomassendtime, biomasstime, bioconctime


def calcconcmasstimeX(
    Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars
):
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    massendtime = np.zeros([len(vars)])
    massendtimey = np.zeros([51, len(vars) + 1])
    massendtimey[:, 0] = range(51)
    di = d + fpre + str(Trial) + fsuf
    print(str(Trial))
    df = np.load(di + fpre + str(Trial) + "_df.npy")
    masstime = np.zeros([np.shape(df)[1], 31, len(vars) + 1])
    conctime = np.zeros([np.shape(df)[1], 31, len(vars) + 1])
    veliredg = df[2, 1:, yin, xright]
    veliledg = df[2, 1:, yin, xleft]
    veloredg = df[2, 1:, yout, xright]
    veloledg = df[2, 1:, yout, xleft]
    veloelem = df[2, 1:, yout, xleft + 1 : xright]
    velielem = df[2, 1:, yin, xleft + 1 : xright]
    velelem = df[2, 1:, yin + 1 : yout, xleft + 1 : xright]
    vellelem = df[2, 1:, yin + 1 : yout, xleft]
    velrelem = df[2, 1:, yin + 1 : yout, xright]
    satielem = df[4, 1:, yin, xleft + 1 : xright]
    satoelem = df[4, 1:, yout, xleft + 1 : xright]
    satlelem = df[4, 1:, yin + 1 : yout, xleft]
    satrelem = df[4, 1:, yin + 1 : yout, xright]
    satiredg = df[4, 1:, yin, xright]
    satiledg = df[4, 1:, yin, xleft]
    satoledg = df[4, 1:, yout, xleft]
    satoredg = df[4, 1:, yout, xright]
    satelem = df[4, 1:, yin + 1 : yout, xleft + 1 : xright]
    for i in range(len(vars)):
        #         massendtime[i] = (df[vars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[vars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg + df[vars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[vars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*(vedge**2) + sum(sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xleft+1:xright-1]*satelem*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright-1]*satielem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright-1]*satoelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xleft]*satlelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xright]*satrelem))*velem*vedge
        #         massendtimey[yin,i+1] = (df[vars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[vars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg)*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright-1]*satielem))*velem*vedge
        #         massendtimey[yout,i+1] = (df[vars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[vars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright-1]*satoelem))*velem*vedge
        #         massendtimey[yin+1:yout, i+1] = sum(sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft+1:xright-1]*satelem*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft]*satlelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xright]*satrelem))*velem*vedge
        masstime[1:, xleft, i + 1] = (
            (
                df[vars[i] - 3, 1:, yin, xleft] * satiledg
                + df[vars[i] - 3, 1:, yout, xleft] * satoledg
            )
            * (vedge ** 2)
            + (np.sum(df[vars[i] - 3, 1:, yin + 1 : yout, xleft] * satielem, axis=-1))
            * velem
            * vedge
        )
        masstime[1:, xright, i + 1] = (
            (
                df[vars[i] - 3, 1:, yin, xright] * satiredg
                + df[vars[i] - 3, 1:, yout, xright] * satoredg
            )
            * (vedge ** 2)
            + (np.sum(df[vars[i] - 3, 1:, yin + 1 : yout, xright] * satoelem, axis=-1))
            * velem
            * vedge
        )
        masstime[1:, xleft + 1 : xright, i + 1] = (
            np.sum(
                df[vars[i] - 3, 1:, yin + 1 : yout, xleft + 1 : xright]
                * satelem
                * (velem ** 2),
                axis=-2,
            )
            + (
                (
                    df[vars[i] - 3, 1:, yin, xleft + 1 : xright] * satlelem
                    + df[vars[i] - 3, 1:, yout, xleft + 1 : xright]
                )
                * satrelem
            )
            * velem
            * vedge
        )
        conctime[1:, xleft, i + 1] = (
            (
                df[vars[i] - 3, 1:, yin, xleft] * satiledg * veliledg
                + df[vars[i] - 3, 1:, yout, xleft] * satoledg * veloledg
            )
            * (vedge)
            + (
                np.sum(
                    df[vars[i] - 3, 1:, yin + 1 : yout, xleft] * satlelem * vellelem,
                    axis=-1,
                )
            )
            * velem
        ) / (vedge * (veloledg + veliledg) + np.sum(velem * vellelem, axis=-1))
        conctime[1:, xright, i + 1] = (
            (
                df[vars[i] - 3, 1:, yin, xright] * satiredg * veliredg
                + df[vars[i] - 3, 1:, yout, xright] * satoredg * veloredg
            )
            * (vedge)
            + (
                np.sum(
                    df[vars[i] - 3, 1:, yin + 1 : yout, xright] * satrelem * velrelem,
                    axis=-1,
                )
            )
            * velem
        ) / (vedge * (veloredg + veliredg) + np.sum(velem * velrelem, axis=-1))
        conctime[1:, xleft + 1 : xright, i + 1] = (
            np.sum(
                df[vars[i] - 3, 1:, yin + 1 : yout, xleft + 1 : xright]
                * satelem
                * velem
                * velelem,
                axis=-2,
            )
            + (
                df[vars[i] - 3, 1:, yin, xleft + 1 : xright] * satielem * velielem
                + df[vars[i] - 3, 1:, yout, xleft + 1 : xright] * satoelem * veloelem
            )
            * vedge
        ) / (vedge * (veloelem + velielem) + np.sum(velem * velelem, axis=-2))

    for i in range(len(vars)):
        massendtime[i] = (
            (
                df[vars[i] - 3, np.shape(df)[1] - 2, yin, xleft]
                * satiledg[int(np.shape(df)[1]) - 2]
                + df[vars[i] - 3, np.shape(df)[1] - 2, yin, xright]
                * satiredg[int(np.shape(df)[1]) - 2]
                + df[vars[i] - 3, np.shape(df)[1] - 2, yout, xleft]
                * satoledg[int(np.shape(df)[1]) - 2]
                + df[vars[i] - 3, np.shape(df)[1] - 2, yout, xright]
                * satoredg[int(np.shape(df)[1]) - 2]
            )
            * (vedge ** 2)
            + sum(
                sum(
                    df[
                        vars[i] - 3,
                        np.shape(df)[1] - 2,
                        yin + 1 : yout,
                        xleft + 1 : xright,
                    ]
                    * satelem[int(np.shape(df)[1]) - 2, :, :]
                    * (velem ** 2)
                )
            )
            + (
                sum(
                    df[vars[i] - 3, np.shape(df)[1] - 2, yin, xleft + 1 : xright]
                    * satielem[np.shape(df)[1] - 2, :]
                )
                + sum(
                    df[vars[i] - 3, np.shape(df)[1] - 2, yout, xleft + 1 : xright]
                    * satoelem[np.shape(df)[1] - 2, :]
                )
                + sum(
                    df[vars[i] - 3, np.shape(df)[1] - 2, yin + 1 : yout, xleft]
                    * satlelem[np.shape(df)[1] - 2, :]
                )
                + sum(
                    df[vars[i] - 3, np.shape(df)[1] - 2, yin + 1 : yout, xright]
                    * satrelem[np.shape(df)[1] - 2, :]
                )
            )
            * velem
            * vedge
        )
        massendtimey[yin, i + 1] = (
            (
                df[vars[i] - 3, np.shape(df)[1] - 2, yin, xleft]
                * satiledg[np.shape(df)[1] - 2]
                + df[vars[i] - 3, np.shape(df)[1] - 2, yin, xright]
                * satiredg[np.shape(df)[1] - 2]
            )
            * (vedge ** 2)
            + (
                sum(
                    df[vars[i] - 3, np.shape(df)[1] - 2, yin, xleft + 1 : xright]
                    * satielem[np.shape(df)[1] - 2, :]
                )
            )
            * velem
            * vedge
        )
        massendtimey[yout, i + 1] = (
            (
                df[vars[i] - 3, np.shape(df)[1] - 2, yout, xleft]
                * satoledg[np.shape(df)[1] - 2]
                + df[vars[i] - 3, np.shape(df)[1] - 2, yout, xright]
                * satoredg[np.shape(df)[1] - 2]
            )
            * (vedge ** 2)
            + (
                sum(
                    df[vars[i] - 3, np.shape(df)[1] - 2, yout, xleft + 1 : xright]
                    * satoelem[np.shape(df)[1] - 2, :]
                )
            )
            * velem
            * vedge
        )
        massendtimey[yin + 1 : yout, i + 1] = (
            sum(
                sum(
                    df[
                        vars[i] - 3,
                        np.shape(df)[1] - 2,
                        yin + 1 : yout,
                        xleft + 1 : xright,
                    ]
                    * satelem[np.shape(df)[1] - 2, :, :]
                    * (velem ** 2)
                )
            )
            + (
                sum(
                    df[vars[i] - 3, np.shape(df)[1] - 2, yin + 1 : yout, xleft]
                    * satlelem[np.shape(df)[1] - 2, :]
                )
                + sum(
                    df[vars[i] - 3, np.shape(df)[1] - 2, yin + 1 : yout, xright]
                    * satrelem[np.shape(df)[1] - 2, :]
                )
            )
            * velem
            * vedge
        )
        masstime[1:, yin, i + 1] = (
            (
                df[vars[i] - 3, 1:, yin, xleft] * satiledg
                + df[vars[i] - 3, 1:, yin, xright] * satiredg
            )
            * (vedge ** 2)
            + (np.sum(df[vars[i] - 3, 1:, yin, xleft + 1 : xright] * satielem, axis=-1))
            * velem
            * vedge
        )
        masstime[1:, yout, i + 1] = (
            (
                df[vars[i] - 3, 1:, yout, xleft] * satoledg
                + df[vars[i] - 3, 1:, yout, xright] * satoredg
            )
            * (vedge ** 2)
            + (
                np.sum(
                    df[vars[i] - 3, 1:, yout, xleft + 1 : xright] * satoelem, axis=-1
                )
            )
            * velem
            * vedge
        )
        masstime[1:, yin + 1 : yout, i + 1] = (
            np.sum(
                df[vars[i] - 3, 1:, yin + 1 : yout, xleft + 1 : xright]
                * satelem
                * (velem ** 2),
                axis=-1,
            )
            + (
                df[vars[i] - 3, 1:, yin + 1 : yout, xleft] * satlelem
                + df[vars[i] - 3, 1:, yin + 1 : yout, xright] * satrelem
            )
            * velem
            * vedge
        )
        conctime[1:, yin, i + 1] = (
            (
                df[vars[i] - 3, 1:, yin, xleft] * satiledg
                + df[vars[i] - 3, 1:, yin, xright] * satiredg
            )
            * (vedge ** 2)
            + np.sum(df[vars[i] - 3, 1:, yin, xleft + 1 : xright] * satielem, axis=-1)
            * velem
            * vedge
        ) / (vbc * vedge)
        conctime[1:, yout, i + 1] = (
            (
                df[vars[i] - 3, 1:, yout, xleft] * satoledg
                + df[vars[i] - 3, 1:, yout, xright] * satoredg
            )
            * (vedge ** 2)
            + (
                np.sum(
                    df[vars[i] - 3, 1:, yout, xleft + 1 : xright] * satoelem, axis=-1
                )
            )
            * velem
            * vedge
        ) / (vbc * vedge)
        conctime[1:, yin + 1 : yout, i + 1] = (
            np.sum(
                df[vars[i] - 3, 1:, yin + 1 : yout, xleft + 1 : xright]
                * satelem
                * (velem ** 2),
                axis=-1,
            )
            + (
                df[vars[i] - 3, 1:, yin + 1 : yout, xleft] * satlelem
                + df[vars[i] - 3, 1:, yin + 1 : yout, xright] * satrelem
            )
            * velem
            * vedge
        ) / (vbc * velem)

    TotalFlow = (veliledg + veloledg + veliredg + veloredg) * vedge + (
        np.sum(vellelem)
        + np.sum(velrelem)
        + np.sum(velelem)
        + np.sum(velielem)
        + np.sum(veloelem)
    ) * velem
    #    Velocity = np.mean ([InVelocity, OutVelocity, MidVelocity])
    Velocity = df[2, np.shape(df)[1] - 1, :, :]
    return df, massendtime, masstime, conctime, np.mean(Velocity)
