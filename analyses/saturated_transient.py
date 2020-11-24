# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:56:16 2020

@author: khurana
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def calcconcmasstimenew (numpyarray,yin,yout,xleft,xright, nodesinydirection, gvarnames,flowregime):
    import data_reader.data_processing as proc
    species = proc.speciesdict(flowregime)
    mobilespecies = list(t for t in species.keys() if (species[t]['Location'] == "Mobile") and (species[t]['State'] != "Dissolved"))

    vedge = 0.005
    velem = 0.01

    df = numpyarray
    
    conctime = np.zeros([np.shape(df)[1], nodesinydirection, len(gvarnames)])    
    veliredg = df[2, 1:, yin, xright]
    veliledg = df[2, 1:, yin, xleft]
    veloredg = df[2, 1:, yout, xright]
    veloledg = df[2, 1:, yout, xleft]
    veloelem = df[2, 1:, yout, xleft + 1 : xright]
    velielem = df[2, 1:, yin, xleft + 1 : xright]
    velelem = df[2, 1:, yin + 1 : yout, xleft + 1 : xright]
    vellelem = df[2, 1:, yin + 1 : yout, xleft]
    velrelem = df[2, 1:, yin + 1 : yout, xright]
    if flowregime == "Saturated":
        satielem = 1
        satoelem = 1
        satlelem = 1
        satrelem = 1
        satiredg = 1
        satiledg = 1
        satoledg = 1
        satoredg = 1
        satelem = 1
    elif flowregime == "Unsaturated":
        satiredg = df[4, 1:, yin, xright]
        satiledg = df[4, 1:, yin, xleft]
        satoredg = df[4, 1:, yout, xright]
        satoledg = df[4, 1:, yout, xleft]
        satoelem = df[4, 1:, yout, xleft + 1 : xright]
        satielem = df[4, 1:, yin, xleft + 1 : xright]
        satlelem = df[4, 1:, yin + 1 : yout, xleft]
        satrelem = df[4, 1:, yin + 1 : yout, xright]
        satelem = df[4, 1:, yin + 1 : yout, xleft + 1 : xright]
    for i in gvarnames:
        idx = gvarnames.index(i)
        if i == "Nitrogen":
           Nspecies = mobilespecies
           ninlet = 0
           noutlet = 0
           for n in Nspecies:
               ninlet = ninlet + (
                   df[species[n]['TecIndex'], 1:, yin, xleft] * satiledg * veliledg * vedge
                   + df[species[n]['TecIndex'], 1:, yin, xright] * satiredg * veliredg * vedge
                   + np.sum(
                       df[species[n]['TecIndex'], 1:, yin, xleft + 1 : xright]
                       * satielem
                       * velielem
                       * velem,
                       axis=-1,
                   )
               ) / (
                   vedge * (veliredg + veliledg)
                   + np.sum(velem * velielem, axis=-1)
               )
               noutlet = noutlet + (
                   df[species[n]['TecIndex'], 1:, yout, xleft] * satoledg * veloledg * vedge
                   + df[species[n]['TecIndex'], 1:, yout, xright] * satoredg * veloredg * vedge
                   + np.sum(
                       df[species[n]['TecIndex'], 1:, yout, xleft + 1 : xright]
                       * satoelem
                       * veloelem
                       * velem,
                       axis=-1,
                   )
                   ) / (
                   vedge * (veloredg + veloledg)
                   + np.sum(velem * veloelem, axis=-1)
               )
           sumin = 0
           sumout = 0
           for r in ["Ammonium", "Nitrate"]:
               sumin = sumin + (
                   df[species[r]['TecIndex'], 1:, yin, xleft] * satiledg * veliledg * vedge
                   + df[species[r]['TecIndex'], 1:, yin, xright] * satiredg * veliredg * vedge
                   + np.sum(
                       df[species[r]['TecIndex'], 1:, yin, xleft + 1 : xright]
                       * satielem
                       * velielem
                       * velem,
                       axis=-1,
                   )
               ) / (
                   vedge * (veliredg + veliledg)
                   + np.sum(velem * velielem, axis=-1)
               )
               sumout = sumout + (
                   df[species[r]['TecIndex'], 1:, yout, xleft] * satoledg * veloledg * vedge
                   + df[species[r]['TecIndex'], 1:, yout, xright] * satoredg * veloredg * vedge
                   + np.sum(
                       df[species[r]['TecIndex'], 1:, yout, xleft + 1 : xright]
                       * satoelem
                       * veloelem
                       * velem,
                       axis=-1,
                   )
                   ) / (
                   vedge * (veloredg + veloledg)
                   + np.sum(velem * veloelem, axis=-1))
            
           conctime[1:, yin, idx] = ninlet / 10 + sumin
           conctime[1:, yout, idx] = noutlet / 10 + sumout  
        elif i == "TOC":
           cinlet = 0
           coutlet = 0
           for c in list(mobilespecies + ["DOC"]):
               cinlet = cinlet + (
                       df[species[c]['TecIndex'], 1:, yin, xleft] * satiledg * veliledg * vedge
                       + df[species[c]['TecIndex'], 1:, yin, xright] * satiredg * veliredg * vedge
                       + np.sum(
                        df[species[c]['TecIndex'], 1:, yin, xleft + 1 : xright]
                        * satielem
                        * velielem
                        * velem,
                        axis=-1)) / (vedge * (veliredg + veliledg)+ np.sum(velem * velielem, axis=-1))
               coutlet = coutlet + (
                    df[species[c]['TecIndex'], 1:, yout, xleft] * satoledg * veloledg * vedge
                    + df[species[c]['TecIndex'], 1:, yout, xright] * satoredg * veloredg * vedge
                    + np.sum(
                        df[species[c]['TecIndex'], 1:, yout, xleft + 1 : xright]
                        * satoelem
                        * veloelem
                        * velem,
                        axis=-1,
                    )
                ) / (
                    vedge * (veloredg + veloledg)
                    + np.sum(velem * veloelem, axis=-1)
                )
           conctime[1:, yin, idx] = cinlet
           conctime[1:, yout, idx] = coutlet
        else:
           conctime[1:, yin, idx] = (
                    (
                        df[species[i]['TecIndex'], 1:, yin, xleft] * satiledg * veliledg
                        + df[species[i]['TecIndex'], 1:, yin, xright] * satiredg * veliredg
                    )
                    * (vedge)
                    + (np.sum(df[species[i]['TecIndex'], 1:, yin, xleft + 1 : xright]
                            * satielem
                            * velielem,
                            axis=-1))* velem) / (vedge * (veliredg + veliledg) + velem*np.sum(velielem, axis=-1))
           conctime[1:, yout, idx] = (
                    (
                        df[species[i]['TecIndex'], 1:, yout, xleft] * satoledg * veloledg
                        + df[species[i]['TecIndex'], 1:, yout, xright] * satoredg * veloredg
                    )
                    * (vedge)
                    + (
                        np.sum(
                            df[species[i]['TecIndex'], 1:, yout, xleft + 1 : xright]
                            * satoelem
                            * veloelem,
                            axis=-1,
                        )
                    )
                    * velem
                ) / (vedge * (veloredg + veloledg) + np.sum(velem * veloelem, axis=-1))
           conctime[1:, yin + 1 : yout, idx] = (
                    np.sum(
                        df[species[i]['TecIndex'], 1:, yin + 1 : yout, xleft + 1 : xright]
                        * satelem
                        * velem
                        * velelem,
                        axis=-1,
                    )
                    + (
                        df[species[i]['TecIndex'], 1:, yin + 1 : yout, xleft] * satlelem * vellelem
                        + df[species[i]['TecIndex'], 1:, yin + 1 : yout, xright]
                        * satrelem
                        * velrelem
                    )
                    * vedge
                ) / (vedge * (vellelem + velrelem) + np.sum(velem * velelem, axis=-1))

    TotalFlow = (veliledg + veloledg + veliredg + veloredg) * vedge + (
        np.sum(vellelem)
        + np.sum(velrelem)
        + np.sum(velelem)
        + np.sum(velielem)
        + np.sum(veloelem)
    ) * velem
    
    Headinlettime = np.mean(df[2, 1:, yin, :], axis=-1) * -1
    
    return conctime, TotalFlow, Headinlettime

def calcsum_temp(data, yin,yout,xleft,xright,gvarnames,flowregime):
    import data_reader.data_processing as proc
    species = proc.speciesdict(flowregime)
    vedge = 0.005
    velem = 0.01
    por = 0.2
    vbc = 0.5 * 0.3
    df = data
    sumalltime = np.zeros([np.shape(df)[1] - 1, len(gvarnames)])
    if flowregime == "Saturated":
        satielem = 1
        satoelem = 1
        satlelem = 1
        satrelem = 1
        satiredg = 1
        satiledg = 1
        satoledg = 1
        satoredg = 1
        satelem = 1
    else:
        satiredg = df[4, 1:, yin, xright]
        satiledg = df[4, 1:, yin, xleft]
        satoredg = df[4, 1:, yout, xright]
        satoledg = df[4, 1:, yout, xleft]
        satoelem = df[4, 1:, yout, xleft + 1 : xright]
        satielem = df[4, 1:, yin, xleft + 1 : xright]
        satlelem = df[4, 1:, yin + 1 : yout, xleft]
        satrelem = df[4, 1:, yin + 1 : yout, xright]
        satelem = df[4, 1:, yin + 1 : yout, xleft + 1 : xright]
    for g in gvarnames:
        b = gvarnames.index(g)
        sumalltime[:, b] = ((df[species[g]['TecIndex'], 1:, yin, xleft] * satiledg
                    + df[species[g]['TecIndex'], 1:, yout, xleft] * satoledg
                    + df[species[g]['TecIndex'], 1:, yin, xright] * satiredg
                    + df[species[g]['TecIndex'], 1:, yout, xright] * satoredg)* vedge**2
                + (np.sum(df[species[g]['TecIndex'], 1:, yin, xleft + 1 : xright] * satielem,axis=-1)
                    + np.sum(df[species[g]['TecIndex'], 1:, yout, xleft + 1 : xright] * satoelem,axis=-1)
                    + np.sum(df[species[g]['TecIndex'], 1:, yin + 1 : yout, xleft] * satlelem,axis=-1)
                    + np.sum(df[species[g]['TecIndex'], 1:, yin + 1 : yout, xright] * satrelem,axis=-1)
                    )* vedge * velem
            + np.sum(np.sum(df[species[g]['TecIndex'], 1:, yin + 1 : yout, xleft + 1 : xright]* satelem, axis=-1),
                     axis=-1)* velem**2)*por/vbc

    return sumalltime


def conc_norm_amplitude(data, benchmark, yin, yout, xleft, xright, nodesinydirection, gvarnames, flowregime):
    
    massflux_amplitude = np. zeros([len(gvarnames)])
    massflux_amplitude_basecase = np. zeros([len(gvarnames)])
    maxwhere = np. zeros([len(gvarnames)])
    basemaxwhere = np. zeros([len(gvarnames)])
    normavgconcout = np.zeros([np.shape(data)[1], len(gvarnames)])
    baseavgconcout = np.zeros([np.shape(data)[1], len(gvarnames)])
    
    conctime, TotalFlow, Headinlettime = calcconcmasstimenew (data,yin,yout,xleft,xright, nodesinydirection, gvarnames,flowregime)
    baseconctime, baseTotalFlow, baseHeadinlettime = calcconcmasstimenew (benchmark,yin,yout,xleft,xright, nodesinydirection, gvarnames,flowregime)
    
    for g in gvarnames:
        normavgconcout[:, gvarnames.index(g)] = conctime[:, yout, gvarnames.index(g)] / np.mean(conctime[:, yout, gvarnames.index(g)])
        baseavgconcout[:, gvarnames.index(g)] = conctime[:, yout, gvarnames.index(g)] / np.mean(baseconctime[-1, yout, gvarnames.index(g)])
    
    for g in gvarnames:
        f, Pxx_spec = signal.periodogram(normavgconcout[:, gvarnames.index(g)], scaling="spectrum")
        massflux_amplitude[gvarnames.index(g)] = np.sqrt(Pxx_spec.max())
        maxwhere[gvarnames.index(g)] = np.argmax(Pxx_spec)
        basef, basePxx_spec = signal.periodogram(baseavgconcout[:, gvarnames.index(g)], scaling="spectrum")
        massflux_amplitude_basecase[gvarnames.index(g)] = np.sqrt(basePxx_spec.max())
        basemaxwhere[gvarnames.index(g)] = np.argmax(basePxx_spec)

    return massflux_amplitude, maxwhere, massflux_amplitude_basecase, basemaxwhere

def mass_norm_amplitude(data, benchmark, yin, yout, xleft, xright, nodesinydirection, gvarnames, flowregime):
    
    mass_amplitude = np. zeros([len(gvarnames)])
    basemass_amplitude = np. zeros([len(gvarnames)])
    maxwhere = np. zeros([len(gvarnames)])
    basemaxwhere = np. zeros([len(gvarnames)])
    normavgmass = np.zeros([np.shape(data)[1]-1, len(gvarnames)])
    baseavgmass = np.zeros([np.shape(data)[1]-1, len(gvarnames)])
    
    masstime = calcsum_temp(data, 0, -1, 0, -1, gvarnames, "Saturated")
    basemasstime = calcsum_temp(benchmark, 0, -1, 0, -1, gvarnames, "Saturated")
    
    for g in gvarnames:
        normavgmass[:, gvarnames.index(g)] = masstime[:, gvarnames.index(g)] / np.mean(masstime[:, gvarnames.index(g)])
        baseavgmass[:, gvarnames.index(g)] = masstime[:, gvarnames.index(g)] / basemasstime[-1, gvarnames.index(g)]
    
    for g in gvarnames:
        f, Pxx_spec = signal.periodogram(normavgmass[:, gvarnames.index(g)], scaling="spectrum")
        mass_amplitude[gvarnames.index(g)] = np.sqrt(Pxx_spec.max())
        maxwhere[gvarnames.index(g)] = np.argmax(Pxx_spec)
        bf, basePxx_spec = signal.periodogram(baseavgmass[:, gvarnames.index(g)], scaling="spectrum")
        basemass_amplitude[gvarnames.index(g)] = np.sqrt(basePxx_spec.max())
        basemaxwhere[gvarnames.index(g)] = np.argmax(basePxx_spec)

    return mass_amplitude, maxwhere, basemass_amplitude, basemaxwhere

def correlation(numpyarray,yin,yout,xleft,xright, nodesinydirection, gvarnames,flowregime):
    df = numpyarray
    conctime0, TotalFlow, Headinlettime0 = calcconcmasstimenew (df,yin,yout,xleft,xright, nodesinydirection, gvarnames,flowregime)
    corrchem = np.zeros([2 * np.shape(Headinlettime0)[0], (len(gvarnames))])
    normavgconcout = np.zeros([np.shape(df)[1], len(gvarnames)])
    
    for g in gvarnames:
        k = gvarnames.index(g)
        normavgconcout[:, k] = conctime0[:, yout, k] - np.mean(conctime0[:, yout, k])
        normheadin = Headinlettime0 - np.mean(Headinlettime0)
        for k in range(len(gvarnames)):
            corrchem[:, k] = np.correlate(normavgconcout[:, k], normheadin, "full") / (
                (np.std(conctime0[:, yout, k]))
                * (np.std(Headinlettime0) * np.shape(Headinlettime0)[0])
            )

    return corrchem, Headinlettime0

def mass_correlation(numpyarray,yin,yout,xleft,xright, nodesinydirection, gvarnames,flowregime):
    data = numpyarray
    masstime = calcsum_temp (data,yin,yout,xleft,xright, gvarnames,flowregime)
    conctime, TotalFlow, Headinlettime0 = calcconcmasstimenew (data,yin,yout,xleft,xright, nodesinydirection, gvarnames,flowregime)
    corrchem = np.zeros([2 * np.shape(Headinlettime0)[0]-1, (len(gvarnames))])
    normavgmass = np.zeros([np.shape(data)[1]-1, len(gvarnames)])
    
    for g in gvarnames:
        k = gvarnames.index(g)
        normavgmass[:, k] = masstime[:, k] - np.mean(masstime[:, k])
        normheadin = Headinlettime0 - np.mean(Headinlettime0)
        for k in range(len(gvarnames)):
            corrchem[:, k] = np.correlate(normavgmass[:, k], normheadin, "full") / (
                (np.std(masstime[:, k]))
                * (np.std(Headinlettime0) * np.shape(Headinlettime0)[0])
            )

    return corrchem, Headinlettime0

#def calcmft_temp(Tforfpre, Trial, gvarnames, d, yin, yout, xleft, xright, fpre, fsuf, Het, Anis):
def calcmft_temp(numpyarray, yin, yout, xleft, xright, gvarnames, flowregime):            
    import data_reader.data_processing as proc
    species = proc.speciesdict(flowregime)
    mobilespecies = list(t for t in species.keys() if (species[t]['Location'] == "Mobile") and (species[t]['State'] != "Dissolved"))

    vedge = 0.005
    velem = 0.01
    por = 0.2
    
    massfluxin = np.zeros([len(gvarnames)])
    massfluxout = np.zeros([len(gvarnames)])
    df = numpyarray
    veliredg = df[2, 1:, yin, xright]
    veliledg = df[2, 1:, yin, xleft]
    veloredg = df[2, 1:, yout, xright]
    veloledg = df[2, 1:, yout, xleft]
    veloelem = df[2, 1:, yout, xleft + 1 : xright]
    velielem = df[2, 1:, yin, xleft + 1 : xright]
#    velelem = df[2, 1:, yin + 1 : yout, xleft + 1 : xright]
#    vellelem = df[2, 1:, yin + 1 : yout, xleft]
#    velrelem = df[2, 1:, yin + 1 : yout, xright]
    if flowregime == "Saturated":
        satielem = 1
        satoelem = 1
        satlelem = 1
        satrelem = 1
        satiredg = 1
        satiledg = 1
        satoledg = 1
        satoredg = 1
        satelem = 1
    elif flowregime == "Unsaturated":
        satiredg = df[4, 1:, yin, xright]
        satiledg = df[4, 1:, yin, xleft]
        satoredg = df[4, 1:, yout, xright]
        satoledg = df[4, 1:, yout, xleft]
        satoelem = df[4, 1:, yout, xleft + 1 : xright]
        satielem = df[4, 1:, yin, xleft + 1 : xright]
#        satlelem = df[4, 1:, yin + 1 : yout, xleft]
#        satrelem = df[4, 1:, yin + 1 : yout, xright]
#        satelem = df[4, 1:, yin + 1 : yout, xleft + 1 : xright]
    for i in gvarnames:
        idx = gvarnames.index(i)
        if i == "Nitrogen":
            ninlet = 0
            noutlet = 0
            for n in mobilespecies:
                ninlet = (
                        ninlet
                        + (sum(df[species[n]['TecIndex'], 1:, yin, xleft]*satiledg*veliledg*vedge)
                            + sum(df[species[n]['TecIndex'], 1:, yin, xright]*satiredg*veliredg*vedge)
                            + sum(np.sum(df[species[n]['TecIndex'], 1:, yin, xleft + 1 : xright]*satielem*velielem*velem,axis=-1)))/por)
                noutlet = (
                            noutlet
                        + (sum(df[species[n]['TecIndex'], 1:, yout, xleft]*satoledg*veloledg*vedge)
                            + sum(df[species[n]['TecIndex'], 1:, yout, xright]*satoredg*veloredg*vedge)
                            + sum(np.sum(df[species[n]['TecIndex'], 1:, yout, xleft + 1 : xright]*satoelem*veloelem*velem,axis=-1)))/por)
            sumin = 0
            sumout = 0
            for r in ["Ammonium", "Nitrate"]:                            
                rin = (
                        sum(df[species[r]['TecIndex'], 1:, yin, xleft]*satiledg*veliledg*vedge)
                + sum(df[species[r]['TecIndex'], 1:, yin, xright]*satiredg*veliredg*vedge)
                + sum(np.sum(df[species[r]['TecIndex'], 1:, yin, xleft + 1 : xright]*satielem*velielem*velem,axis=-1)))/por
                rout = (
                        sum(df[species[r]['TecIndex'], 1:, yout, xleft]*satoledg*veloledg*vedge)
                + sum(df[species[r]['TecIndex'], 1:, yout, xright]*satoredg*veloredg*vedge)
                + sum(                            np.sum(
                        df[species[r]['TecIndex'], 1:, yout, xleft + 1 : xright]
                        * satoelem
                        * veloelem
                        * velem,
                        axis=-1,
                        )
                )) / por
                sumin = sumin + rin
                sumout = sumout + rout
            massfluxin[idx] = ninlet/10 + sumin
            massfluxout[idx] = noutlet/10 + sumout
        elif i == "TOC":
            cinlet = 0
            coutlet = 0
            for c in list(mobilespecies + ["DOC"]):
                cinlet = (
                        cinlet
                        + (
                        sum(                            df[species[c]['TecIndex'], 1:, yin, xleft]
                            * satiledg
                            * veliledg
                            * vedge
                        )
                        + sum(                                df[species[c]['TecIndex'], 1:, yin, xright]
                                * satiredg
                                * veliredg
                                * vedge
                            )
                        + sum(                                np.sum(
                                    df[species[c]['TecIndex'], 1:, yin, xleft + 1 : xright]
                                    * satielem
                                    * velielem
                                    * velem,
                                    axis=-1,
                                )
                            )
                        )
                        / por
                    )
                coutlet = (
                    coutlet
                    + (
                        sum(                            df[species[c]['TecIndex'], 1:, yout, xleft]
                            * satoledg
                            * veloledg
                            * vedge
                            )
                        + sum(                                df[species[c]['TecIndex'], 1:, yout, xright]
                                * satoredg
                                * veloredg
                                * vedge
                            )
                        + sum(                                np.sum(
                                    df[species[c]['TecIndex'], 1:, yout, xleft + 1 : xright]
                                    * satoelem
                                    * veloelem
                                    * velem,
                                    axis=-1,
                                )
                            )
                        )
                        / por
                    )
            massfluxin[idx] = cinlet
            massfluxout[idx] = coutlet
        else:
            massfluxin[idx] = (
                sum(                        df[species[i]['TecIndex'], 1:, yin, xleft]
                        * satiledg
                        * veliledg
                        * vedge
                    )
                + sum(                        df[species[i]['TecIndex'], 1:, yin, xright]
                        * satiredg
                        * veliredg
                        * vedge
                    )
                + sum(                        np.sum(
                            df[species[i]['TecIndex'], 1:, yin, xleft + 1 : xright]
                            * satielem
                            * velielem
                            * velem,
                            axis=-1,
                        )
                    )
                    ) / por
            massfluxout[idx] = (                    sum(                        df[species[i]['TecIndex'], 1:, yout, xleft]
                        * satoledg
                        * veloledg
                        * vedge
                    )
                    + sum(                        df[species[i]['TecIndex'], 1:, yout, xright]
                        * satoredg
                        * veloredg
                        * vedge
                    )
                    + sum(                        np.sum(                            df[species[i]['TecIndex'], 1:, yout, xleft + 1 : xright]
                            * satoelem
                            * veloelem
                            * velem,
                                axis=-1,
                        )
                    )
                ) / por
    return massfluxin, massfluxout

def calcmassfluxtime(
    Trialhead, Trial, Het, Anis, gw, d, fpre, fsuf, yin, yout, xleft, xright, vars
):
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    mft = np.ones([len(Trial) * len(vars), 8])
    dire = d + fpre + str(Trialhead)
    for j in range(len(Trial)):
        print(str(Trial[j]))
        di = dire + str(Trial[j]) + fsuf
        df = np.load(di + fpre + str(Trialhead) + str(Trial[j]) + "_df.npy")
        veliredg = df[2, 1:, yin, xright]
        veliledg = df[2, 1:, yin, xleft]
        veloredg = df[2, 1:, yout, xright]
        veloledg = df[2, 1:, yout, xleft]
        veloelem = df[2, 1:, yout, xleft + 1 : xright]
        velielem = df[2, 1:, yin, xleft + 1 : xright]
        velelem = df[2, 1:, yin + 1 : yout, xleft + 1 : xright]
        vellelem = df[2, 1:, yin + 1 : yout, xleft]
        velrelem = df[2, 1:, yin + 1 : yout, xright]
        if gw == 1:
            satielem = 1
            satoelem = 1
            satlelem = 1
            satrelem = 1
            satiredg = 1
            satiledg = 1
            satoledg = 1
            satoredg = 1
            satelem = 1
        else:
            satiredg = df[4, 1:, yin, xright]
            satiledg = df[4, 1:, yin, xleft]
            satoredg = df[4, 1:, yout, xright]
            satoledg = df[4, 1:, yout, xleft]
            satoelem = df[4, 1:, yout, xleft + 1 : xright]
            satielem = df[4, 1:, yin, xleft + 1 : xright]
            satlelem = df[4, 1:, yin + 1 : yout, xleft]
            satrelem = df[4, 1:, yin + 1 : yout, xright]
            satelem = df[4, 1:, yin + 1 : yout, xleft + 1 : xright]
        for i in range(len(vars)):
            idx = j * len(vars) + i
            if type(Trial[1]) == str:
                mft[idx, 0] = j
            else:
                mft[idx, 0] = Trial[1] + j
            mft[idx, 1] = Het[j]
            mft[idx, 2] = Anis[j]
            mft[idx, 3] = i
            mft[idx, 4] = (
                sum(df[vars[i] - 3, 1:, yin, xleft] * satiledg * veliledg * vedge)
                + sum(df[vars[i] - 3, 1:, yin, xright] * satiredg * veliredg * vedge)
                + sum(
                    np.sum(
                        df[vars[i] - 3, 1:, yin, xleft + 1 : xright]
                        * satielem
                        * velielem
                        * velem,
                        axis=-1,
                    )
                )
            ) / sum(np.sum(velielem * velem, axis=-1) + (veliledg + veliredg) * vedge)
            mft[idx, 5] = (
                sum(df[vars[i] - 3, 1:, yout, xleft] * satoledg * veloledg * vedge)
                + sum(df[vars[i] - 3, 1:, yout, xright] * satoredg * veloredg * vedge)
                + sum(
                    np.sum(
                        df[vars[i] - 3, 1:, yout, xleft + 1 : xright]
                        * satoelem
                        * veloelem
                        * velem,
                        axis=-1,
                    )
                )
            ) / sum((np.sum(veloelem * velem, axis=-1) + (veloledg + veloredg) * vedge))
            mft[idx, 6] = (mft[idx, 4] - mft[idx, 5]) / mft[idx, 4]
            mft[idx, 7] = (mft[idx, 6] - mft[i, 6]) / mft[i, 6]
    return mft


def calcaggrestime(mf):
    H = sorted(np.unique(mf[:, 2, ...]))
    A = sorted(np.unique(mf[..., 3]))
    C = sorted(np.unique(mf[..., 6]))
    Ch = sorted(np.unique(mf[..., 11]))
    meanresults = []
    stdresults = []
    covresults = []
    for i in H:
        remH = mf[np.where(mf[..., 2] == i)]
        for j in A:
            remHA = remH[np.where(remH[..., 3] == j)]
            for k in C:
                remHAC = remHA[np.where(remHA[..., 6] == k)]
                for l in Ch:
                    meanresults.append(
                        [
                            np.mean(
                                remHAC[
                                    np.where(remHAC[..., np.shape(mf[1])[0] - 1] == l)
                                ][..., np.shape(mf[1])[0] - 3]
                            ),
                            np.mean(
                                remHAC[
                                    np.where(remHAC[..., np.shape(mf[1])[0] - 1] == l)
                                ][..., np.shape(mf[1])[0] - 2]
                            ),
                            i,
                            j,
                            k,
                            l,
                        ]
                    )
                    stdresults.append(
                        [
                            np.std(
                                remHAC[
                                    np.where(remHAC[..., np.shape(mf[1])[0] - 1] == l)
                                ][..., np.shape(mf[1])[0] - 3]
                            ),
                            np.std(
                                remHAC[
                                    np.where(remHAC[..., np.shape(mf[1])[0] - 1] == l)
                                ][..., np.shape(mf[1])[0] - 2]
                            ),
                            i,
                            j,
                            k,
                            l,
                        ]
                    )
                    covresults.append(
                        [
                            np.cov(
                                remHAC[
                                    np.where(remHAC[..., np.shape(mf[1])[0] - 1] == l)
                                ][..., np.shape(mf[1])[0] - 3]
                            ),
                            np.cov(
                                remHAC[
                                    np.where(remHAC[..., np.shape(mf[1])[0] - 1] == l)
                                ][..., np.shape(mf[1])[0] - 2]
                            ),
                            i,
                            j,
                            k,
                            l,
                        ]
                    )
    cleanmresults = [x for x in meanresults if str(x[0]) != "nan"]
    cleansresults = [x for x in stdresults if str(x[0]) != "nan"]
    cleancresults = [x for x in covresults if str(x[0]) != "nan"]
    meanresultsarr = np.array(cleanmresults)
    stdresultsarr = np.array(cleansresults)
    covresultsarr = np.array(cleancresults)
    return meanresultsarr, stdresultsarr, covresultsarr


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
    massendtime = np.zeros([len(gvarnames)])
    massendtimey = np.zeros([yout + 1, len(gvarnames)])
    massendtimey[:, 0] = range(yout + 1)
    di = directory + fpre + str(Trial) + fsuf
    print(str(Trial))
    df = np.load(di + fpre + str(Trial) + "_df.npy")
    masstime = np.zeros([np.shape(df)[1], yout + 1, len(gvarnames)])
    conctime = np.zeros([np.shape(df)[1], yout + 1, len(gvarnames)])
    veliredg = df[2, 1:, yin, xright]
    veliledg = df[2, 1:, yin, xleft]
    veloredg = df[2, 1:, yout, xright]
    veloledg = df[2, 1:, yout, xleft]
    veloelem = df[2, 1:, yout, xleft + 1 : xright]
    velielem = df[2, 1:, yin, xleft + 1 : xright]
    velelem = df[2, 1:, yin + 1 : yout, xleft + 1 : xright]
    vellelem = df[2, 1:, yin + 1 : yout, xleft]
    velrelem = df[2, 1:, yin + 1 : yout, xright]
    if gw == 1:
        satielem = 1
        satoelem = 1
        satlelem = 1
        satrelem = 1
        satiredg = 1
        satiledg = 1
        satoledg = 1
        satoredg = 1
        satelem = 1
        for i in range(len(gvarnames)):
            if gvarnames[i] == "Nitrogen":
                ninlet = 0
                noutlet = 0
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
                        vedge * (veliredg + veliledg)
                        + np.sum(velem * velielem, axis=-1)
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
                        vedge * (veloredg + veloledg)
                        + np.sum(velem * veloelem, axis=-1)
                    )
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
                ) / (vedge * (veliredg + veliledg) + np.sum(velem * velielem, axis=-1))
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
                ) / (vedge * (veloredg + veloledg) + np.sum(velem * veloelem, axis=-1))
            elif gvarnames[i] == "TOC":
                cinlet = 0
                coutlet = 0
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
                        vedge * (veliredg + veliledg)
                        + np.sum(velem * velielem, axis=-1)
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
                        vedge * (veloredg + veloledg)
                        + np.sum(velem * veloelem, axis=-1)
                    )
                conctime[
                    1:, yin, i
                ] = cinlet  # /(sum(velielem)*velem + (veliledg+veliredg)*vedge)
                conctime[
                    1:, yout, i
                ] = coutlet  # /(sum(veloelem)*velem + (veloledg+veloredg)*vedge)
            else:
                #                massendtime[i] = (df[vars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[vars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg + df[vars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[vars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*(vedge**2) + sum(sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xleft+1:xright-1]*satelem*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright-1]*satielem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright-1]*satoelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xleft]*satlelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xright]*satrelem))*velem*vedge
                #                massendtimey[yin,i+1] = (df[vars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[vars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg)*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright-1]*satielem))*velem*vedge
                #                massendtimey[yout,i+1] = (df[vars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[vars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright-1]*satoelem))*velem*vedge
                #                massendtimey[yin+1:yout, i+1] = sum(sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft+1:xright-1]*satelem*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft]*satlelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xright]*satrelem))*velem*vedge
                #                masstime[1:,yin,i] = (df[vars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yin,xleft+1:xright]*satielem, axis=-1))*velem*vedge
                #                masstime[1:,yout,i] = (df[vars[i]-3,1:,yout,xleft]*satoledg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1))*velem*vedge
                #                masstime[1:,yin+1:yout, i] = np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + ((df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[vars[i]-3,1:,yin+1:yout,xright])*satrelem)*velem*vedge
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
                ) / (vedge * (veliredg + veliledg) + np.sum(velem * velielem, axis=-1))
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
                ) / (vedge * (veloredg + veloledg) + np.sum(velem * veloelem, axis=-1))
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
                        + df[vars[i] - 3, 1:, yin + 1 : yout, xright]
                        * satrelem
                        * velrelem
                    )
                    * vedge
                ) / (vedge * (vellelem + velrelem) + np.sum(velem * velelem, axis=-1))
    else:
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
            #           massendtime[i] = (df[vars[i]-3,np.shape(df)[1]-2,yin,xleft]*satiledg[int(np.shape(df)[1])-2] + df[vars[i]-3,np.shape(df)[1]-2,yin,xright]*satiredg[int(np.shape(df)[1])-2] + df[vars[i]-3,np.shape(df)[1]-2,yout,xleft]*satoledg[int(np.shape(df)[1])-2] + df[vars[i]-3,np.shape(df)[1]-2,yout,xright]*satoredg[int(np.shape(df)[1])-2])*(vedge**2) + sum(sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft+1:xright]*satelem[int(np.shape(df)[1])-2,:,:]*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yin,xleft+1:xright]*satielem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yout,xleft+1:xright]*satoelem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft]*satlelem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xright]*satrelem[np.shape(df)[1]-2,:]))*velem*vedge
            #           massendtimey[yin,i] = (df[vars[i]-3,np.shape(df)[1]-2,yin,xleft]*satiledg[np.shape(df)[1]-2] + df[vars[i]-3,np.shape(df)[1]-2,yin,xright]*satiredg[np.shape(df)[1]-2])*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yin,xleft+1:xright]*satielem[np.shape(df)[1]-2,:]))*velem*vedge
            #           massendtimey[yout,i] = (df[vars[i]-3,np.shape(df)[1]-2,yout,xleft]*satoledg[np.shape(df)[1]-2] + df[vars[i]-3,np.shape(df)[1]-2,yout,xright]*satoredg[np.shape(df)[1]-2])*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yout,xleft+1:xright]*satoelem[np.shape(df)[1]-2,:]))*velem*vedge
            #           massendtimey[yin+1:yout, i] = sum(sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft+1:xright]*satelem[np.shape(df)[1]-2,:,:]*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xleft]*satlelem[np.shape(df)[1]-2,:]) + sum(df[vars[i]-3,np.shape(df)[1]-2,yin+1:yout,xright]*satrelem[np.shape(df)[1]-2,:]))*velem*vedge
            #           masstime[1:,yin,i] = (df[vars[i]-3,1:,yin,xleft]*satiledg + df[vars[i]-3,1:,yin,xright]*satiredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yin,xleft+1:xright]*satielem, axis=-1))*velem*vedge
            #           masstime[1:,yout,i] = (df[vars[i]-3,1:,yout,xleft]*satoledg + df[vars[i]-3,1:,yout,xright]*satoredg)*(vedge**2) + (np.sum(df[vars[i]-3,1:,yout,xleft+1:xright]*satoelem, axis = -1))*velem*vedge
            #           masstime[1:,yin+1:yout, i] = np.sum(df[vars[i]-3,1:,yin+1:yout,xleft+1:xright]*satelem*(velem**2),axis=-1) + (df[vars[i]-3,1:,yin+1:yout,xleft]*satlelem + df[vars[i]-3,1:,yin+1:yout,xright]*satrelem)*velem*vedge
            conctime[1:, yin, i] = (
                (
                    df[vars[i] - 3, 1:, yin, xleft] * satiledg
                    + df[vars[i] - 3, 1:, yin, xright] * satiredg
                )
                * (vedge ** 2)
                + np.sum(
                    df[vars[i] - 3, 1:, yin, xleft + 1 : xright] * satielem, axis=-1
                )
                * velem
                * vedge
            ) / (vbc * vedge)
            conctime[1:, yout, i] = (
                (
                    df[vars[i] - 3, 1:, yout, xleft] * satoledg
                    + df[vars[i] - 3, 1:, yout, xright] * satoredg
                )
                * (vedge ** 2)
                + (
                    np.sum(
                        df[vars[i] - 3, 1:, yout, xleft + 1 : xright] * satoelem,
                        axis=-1,
                    )
                )
                * velem
                * vedge
            ) / (vbc * vedge)
            conctime[1:, yin + 1 : yout, i] = (
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
    Velocity = df[2, - 1, :, :]
    Headinlettime = np.mean(df[2, 1:, yin, :], axis=-1) * -1
    return df, massendtime, masstime, conctime, TotalFlow, Headinlettime


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
    if gw == 1:
        satielem = 1
        satoelem = 1
        satlelem = 1
        satrelem = 1
        satiredg = 1
        satiledg = 1
        satoledg = 1
        satoredg = 1
        satelem = 1
        for i in range(len(vars)):
            #            massendtime[i] = (df[vars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[vars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg + df[vars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[vars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*(vedge**2) + sum(sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xleft+1:xright-1]*satelem*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright-1]*satielem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright-1]*satoelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xleft]*satlelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout-1,xright]*satrelem))*velem*vedge
            #            massendtimey[yin,i+1] = (df[vars[i]-3,np.shape(df)[1]-1,yin,xleft]*satiledg + df[vars[i]-3,np.shape(df)[1]-1,yin,xright]*satiredg)*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin,xleft+1:xright-1]*satielem))*velem*vedge
            #            massendtimey[yout,i+1] = (df[vars[i]-3,np.shape(df)[1]-1,yout,xleft]*satoledg + df[vars[i]-3,np.shape(df)[1]-1,yout,xright]*satoredg)*(vedge**2) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yout,xleft+1:xright-1]*satoelem))*velem*vedge
            #            massendtimey[yin+1:yout, i+1] = sum(sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft+1:xright-1]*satelem*(velem**2))) + (sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xleft]*satlelem) + sum(df[vars[i]-3,np.shape(df)[1]-1,yin+1:yout,xright]*satrelem))*velem*vedge
            masstime[1:, xleft, i + 1] = (
                (
                    df[vars[i] - 3, 1:, yin, xleft] * satiledg
                    + df[vars[i] - 3, 1:, yout, xleft] * satoledg
                )
                * (vedge ** 2)
                + (
                    np.sum(
                        df[vars[i] - 3, 1:, yin + 1 : yout, xleft] * satielem, axis=-1
                    )
                )
                * velem
                * vedge
            )
            masstime[1:, xright, i + 1] = (
                (
                    df[vars[i] - 3, 1:, yin, xright] * satiredg
                    + df[vars[i] - 3, 1:, yout, xright] * satoredg
                )
                * (vedge ** 2)
                + (
                    np.sum(
                        df[vars[i] - 3, 1:, yin + 1 : yout, xright] * satoelem, axis=-1
                    )
                )
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
                        df[vars[i] - 3, 1:, yin + 1 : yout, xleft]
                        * satlelem
                        * vellelem,
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
                        df[vars[i] - 3, 1:, yin + 1 : yout, xright]
                        * satrelem
                        * velrelem,
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
                    + df[vars[i] - 3, 1:, yout, xleft + 1 : xright]
                    * satoelem
                    * veloelem
                )
                * vedge
            ) / (vedge * (veloelem + velielem) + np.sum(velem * velelem, axis=-2))
    else:
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
                + (
                    np.sum(
                        df[vars[i] - 3, 1:, yin, xleft + 1 : xright] * satielem, axis=-1
                    )
                )
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
                        df[vars[i] - 3, 1:, yout, xleft + 1 : xright] * satoelem,
                        axis=-1,
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
                + np.sum(
                    df[vars[i] - 3, 1:, yin, xleft + 1 : xright] * satielem, axis=-1
                )
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
                        df[vars[i] - 3, 1:, yout, xleft + 1 : xright] * satoelem,
                        axis=-1,
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


def biomasstimefunc(numpyarray, yin, yout, xleft, xright, nodesinydirection, gvarnames, flowregime):
    import data_reader.data_processing as proc
    vedge = 0.005
    velem = 0.01
    vbc = 0.3
    df = numpyarray
    biomasstime = np.zeros([np.shape(df)[1]-1, nodesinydirection,  len(gvarnames)])
    bioconctime = np.zeros([np.shape(df)[1]-1, nodesinydirection, len(gvarnames)])
    species = proc.speciesdict(flowregime)
    
    if flowregime == "Saturated":
        satielem = 1
        satoelem = 1
        satlelem = 1
        satrelem = 1
        satiredg = 1
        satiledg = 1
        satoledg = 1
        satoredg = 1
        satelem = 1
    else:
        satielem = df[4, 1:, yin, xleft + 1 : xright]
        satoelem = df[4, 1:, yout, xleft + 1 : xright]
        satlelem = df[4, 1:, yin + 1 : yout, xleft]
        satrelem = df[4, 1:, yin + 1 : yout, xright]
        satiredg = df[4, 1:, yin, xright]
        satiledg = df[4, 1:, yin, xleft]
        satoledg = df[4, 1:, yout, xleft]
        satoredg = df[4, 1:, yout, xright]
        satelem = df[4, 1:, yin + 1 : yout, xleft + 1 : xright]
    for b in gvarnames:
        i = gvarnames.index(b)
        biomasstime[:, yin, i] = (
                (
                    df[species[b]['TecIndex'], 1:, yin, xleft] * satiledg
                    + df[species[b]['TecIndex'], 1:, yin, xright] * satiredg
                )
                * (vedge ** 2)
                + (
                    np.sum(
                        df[species[b]['TecIndex'], 1:, yin, xleft + 1 : xright] * satielem,
                        axis=-1,
                    )
                )
                * velem
                * vedge
            )
        biomasstime[:, yout, i] = (
                (
                    df[species[b]['TecIndex'], 1:, yout, xleft] * satoledg
                    + df[species[b]['TecIndex'], 1:, yout, xright] * satoredg
                )
                * (vedge ** 2)
                + (
                    np.sum(
                        df[species[b]['TecIndex'], 1:, yout, xleft + 1 : xright] * satoelem,
                        axis=-1,
                    )
                )
                * velem
                * vedge
            )
        biomasstime[:, yin + 1 : yout, i] = (
                np.sum(
                    df[species[b]['TecIndex'], 1:, yin + 1 : yout, xleft + 1 : xright]
                    * satelem
                    * (velem ** 2),
                    axis=-1,
                )
                + (
                    df[species[b]['TecIndex'], 1:, yin + 1 : yout, xleft] * satlelem
                    + df[species[b]['TecIndex'], 1:, yin + 1 : yout, xright] * satrelem
                )
                * velem
                * vedge
            )
        bioconctime[:, yin, i] = (
                (
                    df[species[b]['TecIndex'], 1:, yin, xleft] * satiledg
                    + df[species[b]['TecIndex'], 1:, yin, xright] * satiredg
                )
                * (vedge ** 2)
                + np.sum(
                    df[species[b]['TecIndex'], 1:, yin, xleft + 1 : xright] * satielem,
                    axis=-1,
                )
                * velem
                * vedge
            ) / (vbc * vedge)
        bioconctime[:, yout, i] = (
                (
                    df[species[b]['TecIndex'], 1:, yout, xleft] * satoledg
                    + df[species[b]['TecIndex'], 1:, yout, xright] * satoredg
                )
                * (vedge ** 2)
                + (
                    np.sum(
                        df[species[b]['TecIndex'], 1:, yout, xleft + 1 : xright] * satoelem,
                        axis=-1,
                    )
                )
                * velem
                * vedge
            ) / (vbc * vedge)
        bioconctime[:, yin + 1 : yout, i] = (
                np.sum(
                    df[species[b]['TecIndex'], 1:, yin + 1 : yout, xleft + 1 : xright]
                    * satelem
                    * (velem ** 2),
                    axis=-1,
                )
                + (
                    df[species[b]['TecIndex'], 1:, yin + 1 : yout, xleft] * satlelem
                    + df[species[b]['TecIndex'], 1:, yin + 1 : yout, xright] * satrelem
                )
                * velem
                * vedge
            ) / (vbc * velem)
    return biomasstime, bioconctime

def correlationanalysis(
    directory,
    numberofTempscenarios,
    Trial,
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
    gvarnames,
    AFbiomassvars,
    AFbiomassgvarnames,
):
    Regimes = ["Slow", "Equal", "Fast"]
    amplitudec = np.zeros(
        [len(Trial) * (numberofTempscenarios) * len(Regimes) * (len(gvarnames)), 8]
    )
    amplitudeb = np.zeros(
        [
            len(Trial)
            * (numberofTempscenarios)
            * len(Regimes)
            * (len(AFbiomassgvarnames)),
            8,
        ]
    )
    conccorr = np.zeros(
        [len(Trial) * (numberofTempscenarios) * len(Regimes) * (len(gvarnames)), 8]
    )
    conccorrb = np.zeros(
        [
            len(Trial)
            * (numberofTempscenarios)
            * len(Regimes)
            * (len(AFbiomassgvarnames)),
            8,
        ]
    )
    plotg = ["DOC", "DO", "TOC", "Nitrogen"]
    for Reg in Regimes:
        print(Reg)
        Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
        for p in Tforfpre[1:]:
            print(p)
            for t in Trial:
                print(t)
                #                df = np.load(directory+p+fpre+str(t)+fpre+str(t)+"_df.npy")
                df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = calcconcmasstime(
                    t,
                    Het[Trial.index(t)],
                    Anis[Trial.index(t)],
                    gw,
                    directory + p,
                    fpre,
                    fsuf,
                    yin,
                    yout,
                    xleft,
                    xright,
                    vars,
                    gvarnames,
                )
                sumall0 = calcsum_temp(
                    t,
                    Het[Trial.index(t)],
                    Anis[Trial.index(t)],
                    gw,
                    directory + p,
                    fpre,
                    fsuf,
                    yin,
                    yout,
                    xleft,
                    xright,
                    AFbiomassvars,
                    AFbiomassgvarnames,
                )
                normavgconcout = np.zeros([np.shape(df0)[1], len(gvarnames)])
                normsum = np.zeros([np.shape(df0)[1] - 1, len(AFbiomassgvarnames)])
                spratio = np.zeros([np.shape(df0)[1] - 1, len(AFbiomassgvarnames)])
                normspratio = np.zeros([np.shape(df0)[1] - 1, len(AFbiomassgvarnames)])
                for k in range(len(gvarnames)):
                    normavgconcout[:, k] = conctime0[:, yout, k] - np.mean(
                        conctime0[:, yout, k]
                    )
                normheadin = Headinlettime0 - np.mean(Headinlettime0)
                fig, axes = plt.subplots(nrows=2, ncols=2, figsize=[11, 8], sharex=True)
                plt.suptitle(Reg + p + str(t))
                for k in range(len(gvarnames)):
                    idxc = (
                        Regimes.index(Reg)
                        * (len(Tforfpre) - 1)
                        * len(Trial)
                        * (len(gvarnames))
                        + (Tforfpre.index(p) - 1) * len(Trial) * (len(gvarnames))
                        + Trial.index(t) * (len(gvarnames))
                        + k
                    )
                    if t == "H":
                        conccorr[idxc, 0] = Trial.index(t)
                    else:
                        conccorr[idxc, 0] = t
                    conccorr[idxc, 1] = Het[Trial.index(t)]
                    conccorr[idxc, 2] = Anis[Trial.index(t)]
                    conccorr[idxc, 3] = Tforfpre.index(p)
                    conccorr[idxc, 4] = Regimes.index(Reg)
                    conccorr[idxc, 5] = k
                    conccorrarray = np.correlate(
                        normavgconcout[:, k], normheadin, "full"
                    ) / (
                        (np.std(conctime0[:, yout, k]))
                        * (np.std(Headinlettime0) * np.shape(Headinlettime0)[0])
                    )
                    conccorr[idxc, 6] = conccorrarray[
                        np.shape(Headinlettime0)[0]
                        - 1
                        + np.argmax(
                            np.abs(conccorrarray[np.shape(Headinlettime0)[0] - 1 :])
                        )
                    ]
                    conccorr[idxc, 7] = (
                        np.argmax(
                            np.abs(conccorrarray[np.shape(Headinlettime0)[0] - 1 :])
                        )
                    ) * 5
                    amplitudec[idxc, :6] = conccorr[idxc, :6]
                    amplitudec[idxc, 6] = normavgconcout[
                        np.argmin(normavgconcout[1:, k]), k
                    ] / np.mean(conctime0[:, yout, k])
                    amplitudec[idxc, 7] = normavgconcout[
                        np.argmax(normavgconcout[1:, k]), k
                    ] / np.mean(conctime0[:, yout, k])
                    if gvarnames[k] in plotg:
                        axes.flat[plotg.index(gvarnames[k])].stem(
                            conccorrarray[np.shape(Headinlettime0)[0] - 1 :]
                        )
                for k in range(len(AFbiomassgvarnames)):
                    normsum[:, k] = sumall0[:, k] - np.mean(sumall0[:, k])
                    spratio[:, k] = sumall0[:, k] / np.sum(sumall0[:, :], axis=-1)
                for k in range(len(AFbiomassgvarnames)):
                    normspratio[:, k] = spratio[:, k] - np.mean(spratio[:, k])
                for k in range(len(AFbiomassgvarnames)):
                    idxb = (
                        Regimes.index(Reg)
                        * (len(Tforfpre) - 1)
                        * len(Trial)
                        * (len(AFbiomassgvarnames))
                        + (Tforfpre.index(p) - 1)
                        * len(Trial)
                        * (len(AFbiomassgvarnames))
                        + Trial.index(t) * (len(AFbiomassgvarnames))
                        + k
                    )
                    if t == "H":
                        conccorrb[idxb, 0] = Trial.index(t)
                    else:
                        conccorrb[idxb, 0] = t
                    conccorrb[idxb, 1] = Het[Trial.index(t)]
                    conccorrb[idxb, 2] = Anis[Trial.index(t)]
                    conccorrb[idxb, 3] = Tforfpre.index(p)
                    conccorrb[idxb, 4] = Regimes.index(Reg)
                    conccorrb[idxb, 5] = k
                    conccorrarrayb = np.correlate(
                        normspratio[:, k], normheadin, "full"
                    ) / (
                        (np.std(spratio[:, k]))
                        * (np.std(Headinlettime0) * np.shape(Headinlettime0)[0])
                    )
                    conccorrb[idxb, 6] = conccorrarrayb[
                        np.shape(Headinlettime0)[0]
                        - 1
                        + np.argmax(
                            np.abs(conccorrarrayb[np.shape(Headinlettime0)[0] - 1 :])
                        )
                    ]
                    conccorrb[idxb, 7] = (
                        np.argmax(
                            np.abs(conccorrarrayb[np.shape(Headinlettime0)[0] - 1 :])
                        )
                    ) * 5
                    amplitudeb[idxb, :6] = conccorrb[idxb, :6]
                    amplitudeb[idxb, 6] = normsum[
                        np.argmin(normsum[1:, k]), k
                    ] / np.mean(sumall0[:, k])
                    amplitudeb[idxb, 7] = normsum[
                        np.argmax(normsum[1:, k]), k
                    ] / np.mean(sumall0[:, k])

    return amplitudec, conccorr, amplitudeb, conccorrb


def norm_amplitude(
    directory,
    numberofTempscenarios,
    Trial,
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
    gvarnames,
    AFbiomassvars,
    AFbiomassgvarnames,
):
    Regimes = ["Slow", "Equal", "Fast"]
    Ampchem = np.zeros(
        [len(Regimes) * len(Trial) * numberofTempscenarios * len(gvarnames), 8]
    )
    Ampbio = np.zeros(
        [len(Regimes) * len(Trial) * numberofTempscenarios * len(AFbiomassgvarnames), 8]
    )
    for Reg in Regimes:
        print(Reg)
        Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
        for p in Tforfpre[1:]:
            print(p)
            for t in Trial:
                print(t)
                path = directory + p + fpre + str(t) + fsuf + fpre + str(t)
                df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = calcconcmasstime(
                    t,
                    Het[Trial.index(t)],
                    Anis[Trial.index(t)],
                    gw,
                    directory + p,
                    fpre,
                    fsuf,
                    yin,
                    yout,
                    xleft,
                    xright,
                    vars,
                    gvarnames,
                )
                sumall0 = calcsum_temp(
                    t,
                    Het[Trial.index(t)],
                    Anis[Trial.index(t)],
                    gw,
                    directory + p,
                    fpre,
                    fsuf,
                    yin,
                    yout,
                    xleft,
                    xright,
                    AFbiomassvars,
                    AFbiomassgvarnames,
                )
                normsum = np.zeros([np.shape(df0)[1] - 1, len(AFbiomassgvarnames)])
                spratio = np.zeros([np.shape(df0)[1] - 1, len(AFbiomassgvarnames)])
                normspratio = np.zeros([np.shape(df0)[1] - 1, len(AFbiomassgvarnames)])
                normavgconcout = np.zeros([np.shape(df0)[1], len(gvarnames)])
                normheadin = Headinlettime0 - np.mean(Headinlettime0)
                for k in range(len(gvarnames)):
                    normavgconcout[:, k] = conctime0[:, yout, k] / np.mean(
                        conctime0[:, yout, k]
                    )
                for k in range(len(gvarnames)):
                    f, Pxx_spec = signal.periodogram(
                        normavgconcout[:, k], scaling="spectrum"
                    )
                    idxc = (
                        Regimes.index(Reg)
                        * (len(Tforfpre) - 1)
                        * len(Trial)
                        * (len(gvarnames))
                        + (Tforfpre.index(p) - 1) * len(Trial) * (len(gvarnames))
                        + Trial.index(t) * (len(gvarnames))
                        + k
                    )
                    if t == "H":
                        Ampchem[idxc, 0] = Trial.index(t)
                    else:
                        Ampchem[idxc, 0] = t
                    Ampchem[idxc, 1] = Het[Trial.index(t)]
                    Ampchem[idxc, 2] = Anis[Trial.index(t)]
                    Ampchem[idxc, 3] = Tforfpre.index(p)
                    Ampchem[idxc, 4] = Regimes.index(Reg)
                    Ampchem[idxc, 5] = k
                    Ampchem[idxc, 6] = np.sqrt(Pxx_spec.max())
                    Ampchem[idxc, 7] = np.sqrt(Pxx_spec.mean())
                for k in range(len(AFbiomassgvarnames)):
                    normsum[:, k] = sumall0[:, k] / np.mean(sumall0[:, k])
                    spratio[:, k] = sumall0[:, k] / np.sum(sumall0[:, :], axis=-1)
                for k in range(len(AFbiomassgvarnames)):
                    normspratio[:, k] = spratio[:, k] / np.mean(spratio[:, k])
                for k in range(len(AFbiomassgvarnames)):
                    f, Pxx_spec = signal.periodogram(
                        normspratio[:, k], scaling="spectrum"
                    )
                    idxb = (
                        Regimes.index(Reg)
                        * (len(Tforfpre) - 1)
                        * len(Trial)
                        * (len(AFbiomassgvarnames))
                        + (Tforfpre.index(p) - 1)
                        * len(Trial)
                        * (len(AFbiomassgvarnames))
                        + Trial.index(t) * (len(AFbiomassgvarnames))
                        + k
                    )
                    if t == "H":
                        Ampbio[idxb, 0] = Trial.index(t)
                    else:
                        Ampbio[idxb, 0] = t
                    Ampbio[idxb, 1] = Het[Trial.index(t)]
                    Ampbio[idxb, 2] = Anis[Trial.index(t)]
                    Ampbio[idxb, 3] = Tforfpre.index(p)
                    Ampbio[idxb, 4] = Regimes.index(Reg)
                    Ampbio[idxb, 5] = k
                    Ampbio[idxb, 6] = np.sqrt(Pxx_spec.max())
                    Ampbio[idxb, 7] = np.sqrt(Pxx_spec.mean())

    return Ampchem, Ampbio


def head_amplitude(
    directory,
    numberofTempscenarios,
    Trial,
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
    gvarnames,
):
    Regimes = ["Slow", "Equal", "Fast"]
    amplitude = np.zeros([len(Trial) * (numberofTempscenarios) * len(Regimes), 8])
    idxc = 0
    for Reg in Regimes:
        print(Reg)
        Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
        for p in Tforfpre[1:]:
            print(p)
            for t in Trial:
                print(t)
                df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = calcconcmasstime(
                    t,
                    Het[Trial.index(t)],
                    Anis[Trial.index(t)],
                    gw,
                    directory + p,
                    fpre,
                    fsuf,
                    yin,
                    yout,
                    xleft,
                    xright,
                    vars,
                    gvarnames,
                )
                normheadin = Headinlettime0 - np.mean(Headinlettime0)
                if t == "H":
                    amplitude[idxc, 0] = Trial.index(t)
                else:
                    amplitude[idxc, 0] = t
                amplitude[idxc, 1] = Het[Trial.index(t)]
                amplitude[idxc, 2] = Anis[Trial.index(t)]
                amplitude[idxc, 3] = Tforfpre.index(p)
                amplitude[idxc, 4] = Regimes.index(Reg)
                amplitude[idxc, 5] = 0
                amplitude[idxc, 6] = normheadin[np.argmin(normheadin[1:])] / np.mean(
                    Headinlettime0
                )
                amplitude[idxc, 7] = normheadin[np.argmax(normheadin[1:])] / np.mean(
                    Headinlettime0
                )
                idxc = idxc + 1

    return amplitude
