# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 21:18:49 2020

@author: khurana
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 12:26:03 2019

@author: khurana
"""
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
Trial = list(t for t,values in scdict.items())
Het = list(values['Het'] for t,values in scdict.items())
Anis = list(values['Anis'] for t,values in scdict.items())


# Variations:
#notlist = [43,54]
#Trial = list(t for t,values in scdict.items() if t not in notlist)
#Het = list(values['Het'] for t,values in scdict.items() if t not in notlist)
#Anis = list(values['Anis'] for t,values in scdict.items() if t not in notlist)

Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
nScenarios = 3

# Scan for files
for j in Tforfpre:
    print(j)
    for k in Trial:
        file_name = fpre[1:] + str(k) + "_df.npy"
        cur_dir = (
            directory + j + fpre + str(k) + fsuf
        )  # Dir from where search starts can be replaced with any path
        file_list = os.listdir(cur_dir)
        if file_name in file_list:
            #                print (file_name, " File Exists in: ", cur_dir)
            continue
        else:
            print(file_name, "File not found: ", cur_dir)

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

# Spectral density/spectrogram function plot
Regimes = ["Slow", "Equal", "Fast"]
plotg = ["DOC", "DO", "TOC", "Nitrogen"]
fsample = 73
from scipy import signal
from scipy import fft

fig, axes = plt.subplots(nrows=3, ncols=3, figsize=[20, 15], sharey=True, sharex=True)
plt.suptitle("Spectogram: Head at inlet")
count = 0
for Reg in Regimes:
    print(Reg)
    Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
    for p in Tforfpre[1:]:
        print(p)
        for t in Trial[:1]:
            print(t)
            path = directory + p + fpre + str(t) + fsuf + fpre + str(t)
            #            df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = sta.calcconcmasstime(t, Het[Trial.index(t)], Anis[Trial.index(t)], gw, directory+p, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
            df = np.load(path + "_df.npy")
            Headinlettime0 = np.mean(df[2, 1:, yin, :], axis=-1) * -1
            freqs, times, spect = signal.spectrogram(
                Headinlettime0, fsample, scaling="spectrum"
            )
            print(np.shape(spect))
            #            axes.flat[count].imshow(spect, aspect = 'auto', cmap = 'hot_r', origin = 'lower')
            axes.flat[count].pcolormesh(times, freqs, spect)
            axes.flat[count].set_title(p)
            axes.flat[count].set_ylim([0, 2])
            axes.flat[count].set_ylabel("Frequency (per year)")
            axes.flat[count].set_xlabel("Time (years)")
            count = count + 1
plt.savefig(
    directory + "spectrogram__head_09062020.png",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0,
)

# fastfourier transform - power spectrum - head at inlet
N = 1095
dt = 5 / 365
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=[20, 15], sharex=True, sharey="col")
plt.suptitle("FFT: Head at inlet")
count = 0
for Reg in Regimes:
    print(Reg)
    Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
    for p in Tforfpre[1:]:
        print(p)
        for t in Trial[:1]:
            print(t)
            path = directory + p + fpre + str(t) + fsuf + fpre + str(t)
            #            df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = sta.calcconcmasstime(t, Het[Trial.index(t)], Anis[Trial.index(t)], gw, directory+p, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
            df = np.load(path + "_df.npy")
            Headinlettime0 = np.mean(df[2, 1:, yin, :], axis=-1) * -1
            normheadin = (Headinlettime0 - np.mean(Headinlettime0)) / np.std(
                Headinlettime0
            )
            fhat = np.fft.fft(normheadin, N)
            PSD = fhat * np.conj(fhat) / N
            freq = (1 / (N * dt)) * np.arange(N)
            L = np.arange(1, np.floor(N / 2), dtype="int")
            print(np.shape(yf))
            axes.flat[count].semilogy(
                freq[L], PSD[L], color="r", LineWidth=2, label="Noisy"
            )
            axes.flat[count].set_xlim((0, 5))
            axes.flat[count].grid()
            axes.flat[count].set_title(p)
            axes.flat[count].set_ylabel("Power")
            axes.flat[count].set_xlabel("Frequency (per year)")
            count = count + 1
plt.savefig(
    directory + "fft__head_trimmed.png", dpi=300, bbox_inches="tight", pad_inches=0
)

# fastfourier transform - power spectrum - concentration at outlet inlet
N = 1095
dt = 5 / 365
for Reg in Regimes:
    print(Reg)
    Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
    for p in Tforfpre[1:]:
        print(p)
        for k in plotg:
            fig, axes = plt.subplots(
                nrows=4, ncols=3, figsize=[20, 15], sharey=True, sharex=True
            )
            plt.suptitle(p + k)
            for t in Trial:
                print(t)
                path = directory + p + fpre + str(t) + fsuf + fpre + str(t)
                df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = sta.calcconcmasstime(
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
                normavgconcout = np.zeros([np.shape(df0)[1], len(gvarnames)])
                normavgconcout[:, gvarnames.index(k)] = (
                    conctime0[:, yout, gvarnames.index(k)]
                    - np.mean(conctime0[:, yout, gvarnames.index(k)])
                ) / np.std(conctime0[:, yout, gvarnames.index(k)])
                fhat = np.fft.fft(normavgconcout[:, gvarnames.index(k)], N)
                PSD = fhat * np.conj(fhat) / N
                freq = (1 / (N * dt)) * np.arange(N)
                L = np.arange(1, np.floor(N / 2), dtype="int")
                axes.flat[Trial.index(t)].semilogy(
                    freq[L], PSD[L], color="r", LineWidth=2, label="Noisy"
                )
                axes.flat[Trial.index(t)].set_xlim((0, 3))
                axes.flat[Trial.index(t)].grid()
                axes.flat[Trial.index(t)].set_ylim(bottom=0.1)
                axes.flat[Trial.index(t)].set_ylabel("Power")
                axes.flat[Trial.index(t)].set_xlabel("Frequency (per year)")
                axes.flat[Trial.index(t)].set_title(t)
            plt.savefig(
                directory + p + "_med_var_fft_trimmed_" + k + "_chem.png",
                dpi=300,
                bbox_inches="tight",
                pad_inches=0.01,
            )

for Reg in Regimes:
    print(Reg)
    Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
    for p in Tforfpre[1:]:
        print(p)
        for k in plotg:
            fig, axes = plt.subplots(
                nrows=8, ncols=6, figsize=[20, 15], sharey=True, sharex=True
            )
            plt.suptitle(p + k)
            for t in Trial:
                print(t)
                path = directory + p + fpre + str(t) + fsuf + fpre + str(t)
                df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = sta.calcconcmasstime(
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
                freqs, times, spect = signal.spectrogram(
                    conctime0[:, -1, gvarnames.index(k)], fsample
                )
                #                axes.flat[Trial.index(t)].imshow(spect, aspect = 'auto')#, cmap = 'hot_r', origin = 'lower')
                axes.flat[Trial.index(t)].pcolormesh(times, freqs, spect)
                axes.flat[Trial.index(t)].set_title(t)
                axes.flat[Trial.index(t)].set_ylim([0, 2])
            plt.ylabel("Frequency (per year)")
            plt.xlabel("Time (years)")
            plt.savefig(
                directory + p + "_spectrogram_" + k + "_chem.png",
                dpi=300,
                bbox_inches="tight",
                pad_inches=0.01,
            )

# Welch
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=[20, 15], sharey=True, sharex=True)
plt.suptitle("Spectogram: Head at inlet")
count = 0
for Reg in Regimes:
    print(Reg)
    Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
    for p in Tforfpre[1:]:
        print(p)
        for t in Trial[:1]:
            print(t)
            path = directory + p + fpre + str(t) + fsuf + fpre + str(t)
            df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = sta.calcconcmasstime(
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
            f, Pxx_den = signal.welch(Headinlettime0, fsample)
            axes.flat[count].semilogy(f, Pxx_den)
            #            axes.flat[count].set_ylim([0.5e-6, 1])
            axes.flat[count].set_xlabel("frequency [Hz]")
            axes.flat[count].set_ylabel("PSD [V**2/Hz]")
            count = count + 1
plt.savefig(
    directory + "spectrogram_welch__head.png",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0,
)

for Reg in Regimes:
    print(Reg)
    Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
    for p in Tforfpre[1:]:
        print(p)
        for k in plotg:
            fig, axes = plt.subplots(
                nrows=8, ncols=6, figsize=[20, 15], sharey=True, sharex=True
            )
            plt.suptitle(p + k)
            for t in Trial:
                print(t)
                path = directory + p + fpre + str(t) + fsuf + fpre + str(t)
                df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = sta.calcconcmasstime(
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
                f, Pxx_den = signal.welch(conctime0[:, -1, gvarnames.index(k)], fsample)
                axes.flat[Trial.index(t)].semilogy(f, Pxx_den)
                axes.flat[Trial.index(t)].set_title(t)
            plt.savefig(
                directory + p + "_spectrogram_welch_" + k + "_chem.png",
                dpi=300,
                bbox_inches="tight",
                pad_inches=0,
            )

# Cross spectral density
for Reg in Regimes:
    print(Reg)
    Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
    for p in Tforfpre[1:]:
        print(p)
        for k in plotg:
            fig, axes = plt.subplots(
                nrows=8, ncols=6, figsize=[20, 15], sharey=True, sharex=True
            )
            plt.suptitle(p + k + "_CSD")
            for t in Trial:
                print(t)
                path = directory + p + fpre + str(t) + fsuf + fpre + str(t)
                df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = sta.calcconcmasstime(
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
                normavgconcout = np.zeros([np.shape(df0)[1], len(gvarnames)])
                normavgconcout[:, gvarnames.index(k)] = conctime0[
                    :, yout, gvarnames.index(k)
                ] - np.mean(conctime0[:, yout, gvarnames.index(k)])
                normheadin = Headinlettime0 - np.mean(Headinlettime0)
                f, Pxy = signal.csd(
                    normavgconcout[:, gvarnames.index(k)], normheadin, fsample
                )
                axes.flat[Trial.index(t)].semilogy(f, np.abs(Pxy))
                axes.flat[Trial.index(t)].set_title(t)
            plt.savefig(
                directory + p + "_csd_welch_" + k + "_chem.png",
                dpi=300,
                bbox_inches="tight",
                pad_inches=0,
            )

fig, axes = plt.subplots(nrows=3, ncols=3, figsize=[20, 15], sharey=True, sharex=True)
plt.suptitle("Cross spectral density: Head at inlet")
count = 0
for Reg in Regimes:
    print(Reg)
    Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
    for p in Tforfpre[1:]:
        print(p)
        for t in Trial[:1]:
            print(t)
            path = directory + p + fpre + str(t) + fsuf + fpre + str(t)
            df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = sta.calcconcmasstime(
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
            f, Pxx_den = signal.csd(normheadin, normheadin, fsample)
            axes.flat[count].semilogy(f, Pxx_den)
            #            axes.flat[count].set_ylim([0.5e-6, 1])
            axes.flat[count].set_xlabel("frequency [Hz]")
            axes.flat[count].set_ylabel("PSD [V**2/Hz]")
            count = count + 1
plt.savefig(
    directory + "csd_welch__head.png", dpi=300, bbox_inches="tight", pad_inches=0
)

# Normalized periodogram and RMS amplitude
RMSampchem, RMSampbio = sta.norm_amplitude(
    directory,
    3,
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
)
ampdf = pd.DataFrame(
    RMSampchem,
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Time_series",
        "Regime",
        "Chem",
        "Max_amplitude",
        "Average_amplitude",
    ],
)
ampdf = proc.processdataframe(ampdf, gvarnames)
discard, bth = proc.tracerstudies()
ampbth = pd.merge(
    ampdf, bth[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"]
)
ampbth.to_csv(
    "Z:/Saturated_flow/diffusion_transient/Normalized_RMSamplitude_chem.csv", sep="\t"
)
ampbth = pd.read_csv(
    "Z:/Saturated_flow/diffusion_transient/Normalized_RMSamplitude_chem.csv", sep="\t"
)
ampbth["Max_amplitude%"] = ampbth["Max_amplitude"] * 100
ampbth["Average_amplitude%"] = ampbth["Average_amplitude"] * 100
plotg = ["DOC", "DO", "TOC", "Nitrogen"]
dummy = stp.RMS_chem_2x2(
    ampbth[ampbth["Trial"] != "52"],
    plotg,
    "Max_amplitude%",
    "Sensitivity",
    ["Medium", "Fast"],
)
dummy.savefig(
    directory + "Normalized_Maximum%_RMS_amp.png",
    dpi=300,
    bbox_inches="tight",
    pad=0.01,
)
dummy = stp.RMS_chem_2x2(
    ampbth[ampbth["Trial"] != "52"],
    plotg,
    "Average_amplitude%",
    "Average RMS amplitude",
)
dummy.savefig(
    directory + "Normalized_Average%_RMS_amp.png",
    dpi=300,
    bbox_inches="tight",
    pad=0.01,
)

ampdf = pd.DataFrame(
    RMSampbio,
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Time_series",
        "Regime",
        "Chem",
        "Max_amplitude",
        "Average_amplitude",
    ],
)
ampdf = proc.processdataframe(ampdf, AFbiomassgvarnames)
discard, bth = proc.tracerstudies()
ampbth = pd.merge(
    ampdf, bth[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"]
)
ampbth.to_csv(
    "Z:/Saturated_flow/diffusion_transient/Normalized_RMSamplitude_biomass.csv",
    sep="\t",
)
Chemseries = [
    "Active fixed Aerobes",
    "Active fixed Ammonia oxidizers",
    "Active fixed Nitrate reducers",
    "Active mobile Aerobes",
    "Active mobile Ammonia oxidizers",
    "Active mobile Nitrate reducers",
]
ampbth = pd.read_csv(
    "Z:/Saturated_flow/diffusion_transient/Normalized_RMSamplitude_biomass.csv",
    sep="\t",
)
ampbth["Max_amplitude%"] = ampbth["Max_amplitude"] * 100
ampbth["Average_amplitude%"] = ampbth["Average_amplitude"] * 100
amplitude = stp.RMSamp_biomass(
    ampbth[ampbth["Trial"] != "52"],
    Chemseries,
    "Max_amplitude%",
    "Sensitivity",
    ["Slow", "Medium", "Fast"],
)
amplitude.savefig(
    directory + "Normalized_Maximum%_RMS_amp_biomass.png",
    dpi=300,
    bbox_inches="tight",
    pad=0.01,
)

plotg = ["Head"]
RMSamplitude = np.zeros(
    [len(Regimes) * len(Trial) * (len(Tforfpre) - 1) * len(plotg), 8]
)
for Reg in Regimes:
    print(Reg)
    Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
    for p in Tforfpre[1:]:
        print(p)
        for t in Trial:
            print(t)
            path = directory + p + fpre + str(t) + fsuf + fpre + str(t)
            df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = sta.calcconcmasstime(
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
            normavgconcout = np.zeros([np.shape(df0)[1], len(gvarnames)])
            normheadin = Headinlettime0 - np.mean(Headinlettime0)
            for k in range(len(plotg)):
                f, Pxx_spec = signal.periodogram(normheadin, scaling="spectrum")
                idxc = (
                    Regimes.index(Reg) * (len(Tforfpre) - 1) * len(Trial) * (len(plotg))
                    + (Tforfpre.index(p) - 1) * len(Trial) * (len(plotg))
                    + Trial.index(t) * (len(plotg))
                    + k
                )
                if t == "H":
                    RMSamplitude[idxc, 0] = Trial.index(t)
                else:
                    RMSamplitude[idxc, 0] = t
                RMSamplitude[idxc, 1] = Het[Trial.index(t)]
                RMSamplitude[idxc, 2] = Anis[Trial.index(t)]
                RMSamplitude[idxc, 3] = Tforfpre.index(p)
                RMSamplitude[idxc, 4] = Regimes.index(Reg)
                RMSamplitude[idxc, 5] = k
                RMSamplitude[idxc, 6] = np.sqrt(Pxx_spec.max())
                RMSamplitude[idxc, 7] = np.sqrt(Pxx_spec.mean())

Headampdf = pd.DataFrame(
    RMSamplitude,
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Time_series",
        "Regime",
        "Chem",
        "Max_amplitude",
        "Average_amplitude",
    ],
)
Headampdf = proc.processdataframe(Headampdf, ["Head"])
discard, bth = proc.tracerstudies()
Headampbth = pd.merge(
    Headampdf, bth[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"]
)
Headampbth.to_csv(
    "Z:/Saturated_flow/diffusion_transient/RMSamplitude_Head.csv", sep="\t"
)

# Correlation function plot
Regimes = ["Slow", "Equal", "Fast"]
plotg = ["DOC", "DO", "TOC", "Nitrogen"]
for Reg in Regimes:
    print(Reg)
    Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
    for p in Tforfpre[1:]:
        print(p)
        for k in plotg:
            fig, axes = plt.subplots(nrows=8, ncols=6, figsize=[20, 15], sharey=True)
            plt.suptitle(p + k)
            for t in Trial:
                print(t)
                path = directory + p + fpre + str(t) + fsuf + fpre + str(t)
                acfchem, acfbiomass, Headinlettime = sta.autocorrelation(
                    directory,
                    p,
                    t,
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
                )
                y = acfchem[
                    np.shape(Headinlettime)[0]
                    - 1 : np.shape(Headinlettime)[0]
                    - 1
                    + 300,
                    gvarnames.index(k),
                ]
                my_color = np.where(((y >= 0.5) | (y <= -0.5)), "orange", "skyblue")
                axes.flat[Trial.index(t)].set_title(str(t))
                axes.flat[Trial.index(t)].vlines(
                    range(300), ymin=0, ymax=y, color=my_color, alpha=0.4
                )
                axes.flat[Trial.index(t)].scatter(
                    range(300), y, color=my_color, s=1, alpha=1
                )
                axes.flat[Trial.index(t)].set_title(t)
                axes.flat[Trial.index(t)].set_ylim((-1, 1))
            plt.savefig(
                directory + p + "_" + k + "_chem_till300.png",
                dpi=300,
                bbox_inches="tight",
                pad_inches=0,
            )

# Time analysis for correlation
plotg = ["Head"]
Regimes = ["Slow", "Equal", "Fast"]
Headautocorr = np.zeros(
    [len(Regimes) * len(Trial) * (len(Tforfpre) - 1) * len(plotg), 8]
)
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=[10, 7], sharey=True, sharex=True)
plt.suptitle("Autocorrelation function of external forcing")
count = 0
for Reg in Regimes:
    print(Reg)
    Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
    for p in Tforfpre[1:]:
        print(p)
        for k in plotg:
            for t in Trial:
                print(t)
                path = directory + p + fpre + str(t) + fsuf + fpre + str(t)
                #                df0, massendtime0, masstime0, conctime0, Velocity0, Headinlettime0 = sta.calcconcmasstime(t, Het[Trial.index(t)], Anis[Trial.index(t)], gw, directory+p, fpre, fsuf, yin, yout, xleft, xright, vars, gvarnames)
                df = np.load(path + "_df.npy")
                Headinlettime0 = np.mean(df[2, 1:, yin, :], axis=-1) * -1
                normheadin = Headinlettime0 - np.mean(Headinlettime0)
                autocorr = np.correlate(normheadin, normheadin, "full") / (
                    np.std(Headinlettime0)
                    * np.std(Headinlettime0)
                    * np.shape(Headinlettime0)[0]
                )
                idxc = (
                    Regimes.index(Reg) * (len(Tforfpre) - 1) * len(Trial) * (len(plotg))
                    + (Tforfpre.index(p) - 1) * len(Trial) * (len(plotg))
                    + Trial.index(t) * (len(plotg))
                )
                if t == "H":
                    Headautocorr[idxc, 0] = Trial.index(t)
                else:
                    Headautocorr[idxc, 0] = t
                Headautocorr[idxc, 1] = Het[Trial.index(t)]
                Headautocorr[idxc, 2] = Anis[Trial.index(t)]
                Headautocorr[idxc, 3] = Tforfpre.index(p)
                Headautocorr[idxc, 4] = Regimes.index(Reg)
                Headautocorr[idxc, 5] = plotg.index(k)
                Headautocorr[idxc, 6] = autocorr[np.argmax(np.abs(autocorr))]
                Headautocorr[idxc, 7] = np.mean(autocorr)
                y = autocorr[int((len(autocorr) + 1) / 2) - 1 :]
                my_color = np.where(((y >= 0.7) | (y <= -0.7)), "orange", "skyblue")
                axes.flat[count].vlines(
                    range(int((len(autocorr) + 1) / 2)),
                    ymin=0,
                    ymax=y,
                    color=my_color,
                    alpha=0.4,
                )
                axes.flat[count].scatter(
                    range(int((len(autocorr) + 1) / 2)), y, color=my_color, s=1, alpha=1
                )
                axes.flat[count].set_ylim((-1, 1))
                if count < 3:
                    axes.flat[count].set_title("Time series: " + str(count + 1))
                if count % 3 == 0:
                    axes.flat[count].set_ylabel(Regimes[int(count / 3)])
                if count > 5:
                    axes.flat[count].set_xlabel("Lag in time (days x 5)")
                count = count + 1
plt.savefig(
    directory + p + "_" + k + "_4thattempt_full.png",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0,
)

Ampchem = np.zeros([len(Regimes) * len(Trial) * 3 * len(gvarnames), 9])
Ampbio = np.zeros([len(Regimes) * len(Trial) * 3 * len(AFbiomassgvarnames), 9])
for Reg in Regimes:
    print(Reg)
    Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
    for p in Tforfpre[1:]:
        print(p)
        for t in Trial:
            print(t)
            path = directory + p + fpre + str(t) + fsuf + fpre + str(t)
            acfchem, acfbiomass, Headinlettime = sta.autocorrelation(
                directory,
                p,
                t,
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
            )
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
                maxchem = np.argmax(
                    np.abs(acfchem[np.shape(Headinlettime)[0] - 1 :, k])
                )
                ychem = acfchem[np.shape(Headinlettime)[0] - 1 + maxchem :, k]
                memorychem = np.where((ychem <= 0.7))
                if t == "H":
                    Ampchem[idxc, 0] = Trial.index(t)
                else:
                    Ampchem[idxc, 0] = t
                Ampchem[idxc, 1] = Het[Trial.index(t)]
                Ampchem[idxc, 2] = Anis[Trial.index(t)]
                Ampchem[idxc, 3] = Tforfpre.index(p)
                Ampchem[idxc, 4] = Regimes.index(Reg)
                Ampchem[idxc, 5] = k
                Ampchem[idxc, 6] = maxchem * 5
                Ampchem[idxc, 7] = memorychem[0][0] * 5
                Ampchem[idxc, 8] = acfchem[np.shape(Headinlettime)[0] - 1 + maxchem, k]
            for k in range(len(AFbiomassgvarnames)):
                idxc = (
                    Regimes.index(Reg)
                    * (len(Tforfpre) - 1)
                    * len(Trial)
                    * (len(AFbiomassgvarnames))
                    + (Tforfpre.index(p) - 1) * len(Trial) * (len(AFbiomassgvarnames))
                    + Trial.index(t) * (len(AFbiomassgvarnames))
                    + k
                )
                maxbio = np.argmax(
                    np.abs(acfbiomass[np.shape(Headinlettime)[0] - 1 :, k])
                )
                ybio = acfbiomass[np.shape(Headinlettime)[0] - 1 + maxbio :, k]
                memorybio = np.where((ybio <= 0.7))
                if t == "H":
                    Ampbio[idxc, 0] = Trial.index(t)
                else:
                    Ampbio[idxc, 0] = t
                Ampbio[idxc, 1] = Het[Trial.index(t)]
                Ampbio[idxc, 2] = Anis[Trial.index(t)]
                Ampbio[idxc, 3] = Tforfpre.index(p)
                Ampbio[idxc, 4] = Regimes.index(Reg)
                Ampbio[idxc, 5] = k
                Ampbio[idxc, 6] = maxbio * 5
                Ampbio[idxc, 7] = memorybio[0][0] * 5
                Ampbio[idxc, 8] = acfbiomass[np.shape(Headinlettime)[0] - 1 + maxbio, k]

Ampchem = pd.DataFrame(
    Ampchem,
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Time_series",
        "Regime",
        "Chem",
        "Delay",
        "Memory",
        "Crosscorrelation",
    ],
)
Ampchem = proc.processdataframe(Ampchem, gvarnames)
bth1, bth = proc.tracerstudies()
dfall2 = pd.merge(
    Ampchem, bth[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"]
)
dfall2.to_csv(
    "Z:/Saturated_flow/diffusion_transient/crosschor_memory_chem.csv", sep="\t"
)
dfall2 = pd.read_csv(
    "Z:/Saturated_flow/diffusion_transient/crosschor_memory_chem.csv", sep="\t"
)
dummy = stp.correlationdistributionchem(dfall2[dfall2["Trial"] != "52"], plotg)
dummy.savefig(
    directory + "cross_correlation_distribution_chem.png",
    dpi=300,
    bbox_inches="tight",
    pad=0.01,
)
dfall2["Norm_delay"] = dfall2["Delay"] / (dfall2["Time"] / dfall2["fraction"])
dfall2["Norm_memory"] = dfall2["Memory"] / (dfall2["Time"] / dfall2["fraction"])
dummy = stp.RMS_chem_2x2(
    dfall2[dfall2["Trial"] != "52"], plotg, "Norm_memory", "Normalized Memory", ["Fast"]
)
dummy.savefig(
    directory + "Norm_memory_chem.png", dpi=300, bbox_inches="tight", pad=0.01
)
dummy = stp.RMS_chem_2x2(
    dfall2[dfall2["Trial"] != "52"],
    plotg,
    "Memory",
    "Memory (days)",
    ["Medium", "Fast"],
)
dummy.savefig(directory + "Memory_chem.png", dpi=300, bbox_inches="tight", pad=0.01)
dummy = stp.RMS_chem_2x2(
    dfall2[dfall2["Trial"] != "52"], plotg, "Norm_delay", "Normalized delay", ["Fast"]
)
dummy.savefig(directory + "Norm_delay_chem.png", dpi=300, bbox_inches="tight", pad=0.01)
dummy = stp.RMS_chem_2x2(
    dfall2[dfall2["Trial"] != "52"], plotg, "Delay", "Delay", ["Fast"]
)
dummy.savefig(directory + "Delay_chem.png", dpi=300, bbox_inches="tight", pad=0.01)

dfall = dfall2[dfall2["Chem"] == "Nitrogen"]
plt.hist2d(dfall["fraction"], dfall["Memory"], bins=20, cmap=plt.cm.BuGn_r)

Ampbio = pd.DataFrame(
    Ampbio,
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Time_series",
        "Regime",
        "Chem",
        "Delay",
        "Memory",
        "Crosscorrelation",
    ],
)
Ampbio = proc.processdataframe(Ampbio, AFbiomassgvarnames)
bth1, bth = proc.tracerstudies()
dfallb = pd.merge(
    Ampbio, bth[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"]
)
dfallb.to_csv(
    "Z:/Saturated_flow/diffusion_transient/crosschor_memory_biomass.csv", sep="\t"
)
dfallb = pd.read_csv(
    "Z:/Saturated_flow/diffusion_transient/crosschor_memory_biomass.csv", sep="\t"
)
dummy = stp.correlationdistributionbiomass(dfallb[dfallb["Trial"] != "52"], Chemseries)
dummy.savefig(
    directory + "cross_correlation_distribution_biomass.png",
    dpi=300,
    bbox_inches="tight",
    pad=0.01,
)
dfallb = pd.read_csv(
    "Z:/Saturated_flow/diffusion_transient/crosschor_memory_biomass.csv", sep="\t"
)
dfallb["Norm_delay"] = dfallb["Delay"] / (dfallb["Time"] / dfallb["fraction"])
dfallb["Norm_memory"] = dfallb["Memory"] / (dfallb["Time"] / dfallb["fraction"])
Chemseries = [
    "Active fixed Aerobes",
    "Active fixed Ammonia oxidizers",
    "Active fixed Nitrate reducers",
    "Active mobile Aerobes",
    "Active mobile Ammonia oxidizers",
    "Active mobile Nitrate reducers",
]
# dfall = dfallb[dfallb["Regime"]=="Medium"]
dummy = stp.RMSamp_biomass(
    dfallb[dfallb["Trial"] != "52"],
    Chemseries,
    "Norm_memory",
    "Normalized Memory",
    ["Slow"],
)
dummy.savefig(
    directory + "Norm_memory_biomass.png", dpi=300, bbox_inches="tight", pad=0.01
)
dummy = stp.RMSamp_biomass(
    dfallb[dfallb["Trial"] != "52"],
    Chemseries,
    "Memory",
    "Memory (days)",
    ["Slow", "Medium", "Fast"],
)
dummy.savefig(directory + "Memory_biomass.png", dpi=300, bbox_inches="tight", pad=0.01)
dummy = stp.RMSamp_biomass(
    dfallb[dfallb["Trial"] != "52"], Chemseries, "Norm_delay", "Normalized delay"
)
dummy.savefig(
    directory + "Norm_delay_biomass.png", dpi=300, bbox_inches="tight", pad=0.01
)
dummy = stp.RMSamp_biomass(
    dfallb[dfallb["Trial"] != "52"], Chemseries, "Delay", "Delay"
)
dummy.savefig(directory + "Delay_biomass.png", dpi=300, bbox_inches="tight", pad=0.01)

Regimes = ["Slow", "Medium", "Fast"]
plotg = ["DOC", "DO", "TOC", "Nitrogen"]

for Reg in Regimes:
    print(Reg)
    Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
    for p in Tforfpre[1:]:
        print(p)
        for t in Trial:
            print(t)
            path = directory + p + fpre + str(t) + fsuf + fpre + str(t)
            acfchem, acfbiomass, Headinlettime = sta.autocorrelation(
                directory,
                p,
                t,
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
            )
            fig, axes = plt.subplots(nrows=2, ncols=2, figsize=[11, 8], sharex=True)
            plt.suptitle(
                Reg + p + str(t) + str(Het[Trial.index(t)]) + str(Anis[Trial.index(t)])
            )
            for k in plotg:
                axes.flat[plotg.index(k)].stem(
                    acfchem[np.shape(Headinlettime)[0] - 1 :, gvarnames.index(k)],
                    use_line_collection=True,
                )
                axes.flat[plotg.index(k)].set_title(k)
            plt.savefig(directory + p + fsuf + fpre + str(t) + "_chem.png")

import statsmodels.tsa.stattools as smtools

Regimes = ["Slow", "Equal", "Fast"]
numberofTempscenarios = 3
autocov = np.zeros(
    [len(Trial) * (numberofTempscenarios) * len(Regimes) * (len(gvarnames)), 7]
)
crosscov = np.zeros(
    [len(Trial) * (numberofTempscenarios) * len(Regimes) * (len(gvarnames)), 7]
)
for Reg in Regimes:
    d = r"Z:/Saturated_flow/diffusion_transient/"
    Tforfpre = [Reg + "AR_0", Reg + "AR_1/", Reg + "AR_2/", Reg + "AR_5/"]
    for p in Tforfpre:
        newd = d + p
        for t in Trial:
            df, massendtime, masstime, conctime, Velocity, head = sta.calcconcmasstime(
                t,
                Het[Trial.index(t)],
                Anis[Trial.index(t)],
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
            for i in gvarnames:
                idxc = (
                    Regimes.index(Reg)
                    * (len(Tforfpre) - 1)
                    * len(Trial)
                    * (len(gvarnames))
                    + (Tforfpre.index(p) - 1) * len(Trial) * (len(gvarnames))
                    + Trial.index(t) * (len(gvarnames))
                    + gvarnames.index(i)
                )
                x = conctime[:, yout, gvarnames.index(i)]
                y = conctime[:, yout, gvarnames.index("Head")]
                acovcal = np.max(smtools.acovf(x, demean=True))
                ccovcal = np.max(smtools.ccovf(x, y, demean=True))
                if t == "H":
                    autocov[idxc, 0] = Trial.index(t)
                else:
                    autocov[idxc, 0] = t
                autocov[idxc, 1] = Het[Trial.index(t)]
                autocov[idxc, 2] = Anis[Trial.index(t)]
                autocov[idxc, 3] = Tforfpre.index(p)
                autocov[idxc, 4] = Regimes.index(Reg)
                autocov[idxc, 5] = gvarnames.index(i)
                autocov[idxc, 6] = acovcal
                crosscov[idxc, :6] = autocov[idxc, :6]
                crosscov[idxc, 6] = ccovcal

autocovdf = pd.DataFrame(
    autocov,
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Time_series",
        "Regime",
        "Chem",
        "Autocovariance",
    ],
)
autocovdf = proc.processdataframe(autocovdf, gvarnames)
bth1, bth = proc.tracerstudies()
dfall2 = pd.merge(
    autocovdf, bth[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"]
)
dfall2.to_csv("Z:/Saturated_flow/diffusion_transient/autocovariance_chem.csv", sep="\t")
data = pd.read_csv(
    "Z:/Saturated_flow/diffusion_transient/autocovariance_chem.csv", sep="\t"
)
data = data[data["Trial"] != "52"]
plotg = ["DOC", "TOC", "DO", "Nitrogen"]
dummy = stp.autocovariance_chem(data, plotg)
colseries = ["Reds", "Greens", "Blues"]
Regimes = ["Slow", "Medium", "Fast"]
xall = dfall2[dfall2["Chem"] == "Head"]
for k in plotg:
    dfall3 = dfall2[dfall2["Chem"] == k]
    #    fig, host = plt.subplots()
    #    fig.subplots_adjust(right=0.75)
    #    par1 = host.twinx()
    plt.figure()
    for Reg in ["Slow", "Medium", "Fast"]:
        y = dfall3[dfall3["Regime"] == Reg]
        x = xall[xall["Regime"] == Reg]
        plt.scatter(
            y["fraction"],
            y["Autocovariance"],
            c=y["Time_series"],
            cmap=colseries[Regimes.index(Reg)],
            marker="o",
        )
        #        par1.scatter(y["fraction"],x["Autocovariance"], c = y["Time_series"], cmap = colseries[Regimes.index(Reg)], marker = '^')
        #        par1.set_yscale("log")
        plt.yscale("log")
        plt.title(k)

ccovdf = pd.DataFrame(
    crosscov,
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Time_series",
        "Regime",
        "Chem",
        "Crosscovariance",
    ],
)
ccovdf = proc.processdataframe(ccovdf, gvarnames)
dfall2 = pd.merge(
    ccovdf, bth[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"]
)
dfall2.to_csv(
    "Z:/Saturated_flow/diffusion_transient/crosscovariance_chem.csv", sep="\t"
)
plotg = ["DOC", "TOC", "DO", "Nitrogen", "Head"]
colseries = ["Reds", "Greens", "Blues"]
Regimes = ["Slow", "Medium", "Fast"]
xall = dfall2[dfall2["Chem"] == "Head"]
for k in plotg:
    dfall3 = dfall2[dfall2["Chem"] == k]
    #    fig, host = plt.subplots()
    #    fig.subplots_adjust(right=0.75)
    #    par1 = host.twinx()
    plt.figure()
    for Reg in ["Slow", "Medium", "Fast"]:
        y = dfall3[dfall3["Regime"] == Reg]
        x = xall[xall["Regime"] == Reg]
        plt.scatter(
            y["fraction"],
            y["Crosscovariance"],
            c=y["Time_series"],
            cmap=colseries[Regimes.index(Reg)],
            marker="o",
        )
        #        par1.scatter(y["fraction"],x["Autocovariance"], c = y["Time_series"], cmap = colseries[Regimes.index(Reg)], marker = '^')
        #        par1.set_yscale("log")
        plt.yscale("log")
        plt.title(k)

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
    #    for p in Tforfpre[1:]:
    #        for t in Trial:
    #            sumall0 = sta.calcsum_temp(t, Het[Trial.index(t)], Anis[Trial.index(t)], gw, d+p, fpre, fsuf, yin, yout, xleft, xright, AFbiomassvars, AFbiomassgvarnames)
    #            spratio = np.zeros([1095,len(AFbiomassgvarnames)])
    #            for k in range(len(AFbiomassgvarnames)):
    #                spratio[:,k] = sumall0[:,k]/np.sum(sumall0[:,:], axis = -1)
    # Delete following file writing lines if df.to_csv seems to be working fine.
    # Writing the results in a csv file
    #    fnamemf = "X:/diffusion_"+Reg+"_Transient_"+str(Trial[0])+"_"+str(Trial[-1])+"_"+"MF_"+str(gw) + ".csv"
    #    csvfile= open(fnamemf, "w")
    #    writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
    #    writer.writerow(["Sno","Trial", "Variance", "Anisotropy", "Chem", "Inlet_total_mass_flux", "Outlet_mass_flux", "Removal", "Removal_bc_Ratio", "RemovalRatio_homo_Time","RemovalRatio_hetero_Time","Time", "Regime"])
    #    idx = 0
    #    for l in range(len(Tforfpre)):
    #        for j in range(len(Trial)):
    #            for i in range(len(gvarnames)):
    #                idx = (l*len(Trial)+j)*len(gvarnames)+i
    #                if (Reg == "Equal"):
    #                    writer.writerow([idx, Trial[j], Het[j],Anis[j], gvarnames[i], mft[idx,4], mft[idx,5], mft[idx,6], mft[idx,7],mft[idx,10],mft[idx,8],mft[idx,9],"Medium"])
    #                else:
    #                    writer.writerow([idx, Trial[j], Het[j],Anis[j], gvarnames[i], mft[idx,4], mft[idx,5], mft[idx,6], mft[idx,7],mft[idx,10],mft[idx,8],mft[idx,9],Reg])
    #    csvfile.close()
    #
    #    fnamebiosum = "X:/diffusion_"+Reg+"_Transient"+str(Trial[0])+"_"+str(Trial[-1])+"_"+"biomass_"+str(gw) +".csv"
    #    csvfile= open(fnamebiosum, "w")
    #    writer = csv.writer(csvfile, delimiter='\t', quotechar='\t', quoting=csv.QUOTE_MINIMAL, lineterminator = '\n')
    #    writer.writerow(["Sno","Trial", "Variance", "Anisotropy", "Chem", "Total_biomass_umoles", "TotalRatio_bc", "TotalRatio_hetero_Time",  "TotalRatio_homo_Time", "Species_fraction","SpeciesRatio_bc", "SpeciesRatio_hetero_Time", "SpeciesRatio_homo_Time","Time", "Regime"])
    #    idxb = 0
    #    for l in range(len(Tforfpre)):
    #        for j in range(len(Trial)):
    #            for i in AFbiomassgvarnames:
    #                idx = (l*len(Trial) + j)*len(AFbiomassgvarnames)+AFbiomassgvarnames.index(i)
    #                if (Reg == "Equal"):
    #                    writer.writerow([idx+1, Trial[j], Het[j], Anis[j], i, calcsum[idx,4], calcsum[idx,5], calcsum[idx,6],calcsum[idx,7],calcsum[idx,9],calcsum[idx,10],calcsum[idx,11],calcsum[idx,12],calcsum[idx,8], "Medium"])
    #                else:
    #                    writer.writerow([idx+1, Trial[j], Het[j], Anis[j], i, calcsum[idx,4], calcsum[idx,5], calcsum[idx,6],calcsum[idx,7],calcsum[idx,9],calcsum[idx,10],calcsum[idx,11],calcsum[idx,12],calcsum[idx,8], Reg])
    #    csvfile.close()
    # ----------
    print("Files written, processing as dataframes ...")
    if Reg == "Fast":
        #        Reg = "Fast"
        #        fnamemf = "Fast_Transient_H_59_MF_1_time.csv"
        #        fnamebiosum = "Fast_Transient_AnaerobicH_59_biomass_1_time.csv"
        fastb = proc.processchembiomassfiles(fnamebiomasssum, Reg)
        #        fastmbc = sk.processchembiomassfiles(fnamebiosummob, Reg)
        fastc = proc.processchemfiles(fnamemf, Reg)
    elif Reg == "Slow":
        #        Reg = "Slow"
        #        fnamemf = "X:/Saturated_flow/Anaerobic_Slow_Transient_H_59_MF_1.csv"
        #        fnamebiosum = "X:/Saturated_flow/Anaerobic_Slow_Transient_AnaerobicH_59_biomass_1.csv"
        #        fnamebiosummob = "X:/Saturated_flow/Anaerobic_Slow_Transient_AnaerobicH_59__mobile_biomass_1.csv"
        slowb = proc.processchembiomassfiles(fnamebiomasssum, Reg)
        #        slowmbc = sk.processchembiomassfiles(fnamebiosummob, Reg)
        slowc = proc.processchemfiles(fnamemf, Reg)
    elif Reg == "Equal":
        #        Reg = "Equal"
        #        fnamemf = "Equal_Transient_H_59_MF_1_time.csv"
        #        fnamebiosum = "Equal_Transient_AnaerobicH_59_biomass_1_time.csv"
        equalb = proc.processchembiomassfiles(fnamebiomasssum, "Medium")
        #        equalmbc = sk.processchembiomassfiles(fnamebiosummob, "Medium")
        equalc = proc.processchemfiles(fnamemf, "Medium")
    print("Dataframes processed, move to next series")
mastermf2 = proc.processdataframe(mastermf, gvarnames)
masterbiomass2 = proc.processdataframe(masterbiomass, AFbiomassgvarnames)
mastermf.to_csv(
    "Z:/Saturated_flow/diffusion_transient/mass_flux_temporal_impact.csv", sep="\t"
)
masterbiomass.to_csv(
    "Z:/Saturated_flow/diffusion_transient/biomass_temporal_impact.csv", sep="\t"
)
# Aggregated results
# Mass flux and biomass with average residence time
# dummy = sssp.scatterrestime_flux_temp(fastc, equalc, slowc, ["DOC", "DO", "Nitrogen", "TOC"], "Removal_bc_Ratio", "Ratio of removal with base case in uniform flow rate")
# dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_flux_bc_uf.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
# dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_flux_bc_uf.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
# dummy = sssp.scatterrestime_flux_temp(fastc, equalc, slowc, ["DOC", "DO", "Nitrogen", "TOC"], "RemovalRatio_homo_Time", "Ratio of removal with base case in varying flow rate")
# dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_flux_bc_vf.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
# dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_flux_bc_vf.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
# dummy = stp.scatterrestime_flux_temp(fastc, equalc, slowc, ["DOC", "DO", "Nitrogen", "TOC"], "RemovalRatio_hetero_Time", "Ratio of removal with heterogeneous case in varying flow rate")
# dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_flux_hc_vf.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
# dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_flux_hc_vf.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
# dummy = stp.scatterrestime_flux_temp_singleaxis(fastc, equalc, slowc, ["DOC", "DO", "Nitrogen", "TOC"], "RemovalRatio_homo_Time", "Ratio of removal with base case in varying flow rate")
# dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_flux_bv_vf_1col.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
# dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_flux_bc_vf_1col.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
# dummy = stp.scatterrestime_biomass_temp(fastb, slowb, equalb, AFbiomassgvarnames, "SpeciesRatio_hetero_Time", "Ratio of species fraction with heterogeneous case in varying flow rate")
# dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_biomass_temp_speciesfraction_hc_vf.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
# dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_biomass_temp_speciesfraction_hc_vf.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
# dummy = stp.scatterrestime_biomass_temp(fastb, slowb, equalb, AFbiomassgvarnames, "SpeciesRatio_homo_Time", "Ratio of species fraction with homogeneous case in varying flow rate")
# dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_biomass_temp_speciesfraction_bc_vf.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
# dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_biomass_temp_speciesfraction_bc_vf.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
# dummy = stp.scatterrestime_biomass_temp(fastb, slowb, equalb, AFbiomassgvarnames, "SpeciesRatio_bc", "Ratio of species fraction with homogeneous case in uniform flow rate")
# dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_biomass_temp_speciesfraction_bc_uf.png", dpi = 300, bbox_inches = 'tight', pad_inches = 0)
# dummy.savefig("Z:/Saturated_flow/diffusion_transient/scatterrestime_biomass_temp_speciesfraction_bc_uf.pdf", dpi = 300, bbox_inches = 'tight', pad_inches = 0)

# Correlation analysis
# loading tracer study results for comparison with steady state results later on
bth1, bth = proc.tracerstudies()
amplitudev = sta.head_amplitude(
    directory,
    nScenarios,
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
)
corrdfv = pd.DataFrame(
    amplitudev,
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Time_series",
        "Regime",
        "Chem",
        "Min",
        "Max",
    ],
)
corrdfv = proc.processdataframe(corrdfv, ["Head"])

amplitudec, conccorr, amplitudeb, conccorrb = sta.correlationanalysis(
    directory,
    nScenarios,
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
)
corrdf = pd.DataFrame(
    amplitudec,
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Time_series",
        "Regime",
        "Chem",
        "Min",
        "Max",
    ],
)
corrdf = proc.processdataframe(corrdf, gvarnames)
corrd = pd.concat([corrdf, corrdfv], axis=0)
dfall2 = pd.merge(
    corrd, bth[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"]
)
dfall2.to_csv("Z:/Saturated_flow/diffusion_transient/amplitude%_chem.csv", sep="\t")

corrdf = pd.DataFrame(
    amplitudeb,
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Time_series",
        "Regime",
        "Chem",
        "Min",
        "Max",
    ],
)
corrdf = proc.processdataframe(corrdf, AFbiomassgvarnames)
corrd = pd.concat([corrdf, corrdfv], axis=0)
dfall2b = pd.merge(
    corrd, bth[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"]
)
dfall2b.to_csv("Z:/Saturated_flow/diffusion_transient/amplitude%_biomass.csv", sep="\t")

corrdf = pd.DataFrame(
    conccorr,
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Time_series",
        "Regime",
        "Chem",
        "Correlation",
        "Delay",
    ],
)
corrdf = proc.processdataframe(corrdf, gvarnames)
dfall2 = pd.merge(
    corrdf, bth[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"]
)
dfall2["Delayfractionofbth"] = dfall2["Delay"] / dfall2["Time"]
dfall2.to_csv("Z:/Saturated_flow/diffusion_transient/correlation_chem.csv", sep="\t")

corrdf = pd.DataFrame(
    conccorrb,
    columns=[
        "Trial",
        "Variance",
        "Anisotropy",
        "Time_series",
        "Regime",
        "Chem",
        "Correlation",
        "Delay",
    ],
)
corrdf = proc.processdataframe(corrdf, AFbiomassgvarnames)
dfall2b = pd.merge(
    corrdf, bth[["Trial", "Regime", "fraction", "Time"]], on=["Trial", "Regime"]
)
dfall2b["Delayfractionofbth"] = dfall2b["Delay"] / dfall2b["Time"]
dfall2b.to_csv(
    "Z:/Saturated_flow/diffusion_transient/correlation_biomass.csv", sep="\t"
)

Chemseries = [
    "Active fixed Aerobes",
    "Active fixed Ammonia oxidizers",
    "Active fixed Nitrate reducers",
    "Active mobile Aerobes",
    "Active mobile Ammonia oxidizers",
    "Active mobile Nitrate reducers",
]
data_biomass = pd.read_csv(
    "Z:/Saturated_flow/diffusion_transient/amplitude%_biomass.csv",
    sep="\t",
    header="infer",
)
amplitude = stp.amplitude_biomass(data_biomass, Chemseries)
amplitude.savefig(
    "Z:/Saturated_flow/diffusion_transient/amplitude_biomass_singletimeseries.png",
    dpi=300,
    bbox_inches="tight",
)  # , pad_inches = 0)
amplitude.savefig(
    "Z:/Saturated_flow/diffusion_transient/amplitude_biomass_singletimeseries.pdf",
    dpi=300,
    bbox_inches="tight",
)  # , pad_inches = 0.1)

data_biomass = pd.read_csv(
    "Z:/Saturated_flow/diffusion_transient/correlation_biomass.csv",
    sep="\t",
    header="infer",
)
corrplot, delayplot = stp.correlation_delay_biomass(data_biomass, Chemseries)
corrplot.savefig(
    "Z:/Saturated_flow/diffusion_transient/corr_biomass.png",
    dpi=300,
    bbox_inches="tight",
)  # , pad_inches = 0)
corrplot.savefig(
    "Z:/Saturated_flow/diffusion_transient/corr_biomass.pdf",
    dpi=300,
    bbox_inches="tight",
)  # , pad_inches = 0.1)
delayplot.savefig(
    "Z:/Saturated_flow/diffusion_transient/delay_biomass.png",
    dpi=300,
    bbox_inches="tight",
)  # , pad_inches = 0)
delayplot.savefig(
    "Z:/Saturated_flow/diffusion_transient/delay_biomass.pdf",
    dpi=300,
    bbox_inches="tight",
)  # , pad_inches = 0.1)

Chemseries = ["DOC", "DO", "Nitrogen", "TOC", "Head"]
data_chem = pd.read_csv(
    "Z:/Saturated_flow/diffusion_transient/amplitude%_chem.csv",
    sep="\t",
    header="infer",
)
amplitude = stp.amplitude_chem(data_chem, Chemseries)
amplitude.savefig(
    "Z:/Saturated_flow/diffusion_transient/amplitude_chem.png",
    dpi=300,
    bbox_inches="tight",
)  # , pad_inches = 0)
amplitude.savefig(
    "Z:/Saturated_flow/diffusion_transient/amplitude_chem.pdf",
    dpi=300,
    bbox_inches="tight",
)  # , pad_inches = 0.1)

Chemseries = ["DOC", "Head", "Nitrogen", "TOC"]
amplitude = stp.amplitude_chem_2x2(data_chem, Chemseries)
amplitude.savefig(
    "Z:/Saturated_flow/diffusion_transient/amplitude_chem_square_singletimeseries.png",
    dpi=300,
    bbox_inches="tight",
)  # , pad_inches = 0)
amplitude.savefig(
    "Z:/Saturated_flow/diffusion_transient/amplitude_chem_square_singletimeseries.pdf",
    dpi=300,
    bbox_inches="tight",
)  # , pad_inches = 0.1)

amplitudehead = stp.amplitude_head(data_chem)
data_chem = pd.read_csv(
    "Z:/Saturated_flow/diffusion_transient/correlation_chem.csv",
    sep="\t",
    header="infer",
)
corrplot, delayplot = stp.correlation_delay_chem(data_chem, Chemseries)
corrplot.savefig(
    "Z:/Saturated_flow/diffusion_transient/corr_chem.png", dpi=300, bbox_inches="tight"
)  # , pad_inches = 0)
corrplot.savefig(
    "Z:/Saturated_flow/diffusion_transient/corr_chem.pdf", dpi=300, bbox_inches="tight"
)  # , pad_inches = 0)
delayplot.savefig(
    "Z:/Saturated_flow/diffusion_transient/delay_chem.pdf", dpi=300, bbox_inches="tight"
)  # , pad_inches = 0.1)
delayplot.savefig(
    "Z:/Saturated_flow/diffusion_transient/delay_chem.png", dpi=300, bbox_inches="tight"
)  # , pad_inches = 0.1)

corrplot, delayplot = stp.correlation_delay_chem_2x2(data_chem, Chemseries)
corrplot.savefig(
    "Z:/Saturated_flow/diffusion_transient/corr_chem_square.png",
    dpi=300,
    bbox_inches="tight",
)  # , pad_inches = 0)
corrplot.savefig(
    "Z:/Saturated_flow/diffusion_transient/corr_chem_square.pdf",
    dpi=300,
    bbox_inches="tight",
)  # , pad_inches = 0)
delayplot.savefig(
    "Z:/Saturated_flow/diffusion_transient/delay_chem_square.pdf",
    dpi=300,
    bbox_inches="tight",
)  # , pad_inches = 0.1)
delayplot.savefig(
    "Z:/Saturated_flow/diffusion_transient/delay_chem_square.png",
    dpi=300,
    bbox_inches="tight",
)  # , pad_inches = 0.1)


Regimes = ["Slow", "Equal", "Fast"]
intsce = ["H", 44, 76, 73, 80, 84, 63]
intsce = [76, 73, 80, 84, 63]
for i in intsce:
    initseries = [500, 430, 600]
    lastseries = [700, 630, 800]
    stp.generate_chem_timeseries(
        Regimes,
        initseries,
        lastseries,
        Trial[Trial.index(i)],
        Het[Trial.index(i)],
        Anis[Trial.index(i)],
        gw,
        directory,
        fpre,
        vars,
        gvarnames,
        fsuf,
        yin,
        yout,
        xleft,
        xright,
    )
    stp.generate_timeseries(
        Regimes,
        initseries,
        lastseries,
        Trial[Trial.index(i)],
        Het[Trial.index(i)],
        Anis[Trial.index(i)],
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
        "Active fixed",
    )
    stp.generate_timeseries(
        Regimes,
        initseries,
        lastseries,
        Trial[Trial.index(i)],
        Het[Trial.index(i)],
        Anis[Trial.index(i)],
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
        "Active mobile",
    )
    stp.generate_timeseries(
        Regimes,
        initseries,
        lastseries,
        Trial[Trial.index(i)],
        Het[Trial.index(i)],
        Anis[Trial.index(i)],
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
        "Inactive fixed",
    )
    stp.generate_timeseries(
        Regimes,
        initseries,
        lastseries,
        Trial[Trial.index(i)],
        Het[Trial.index(i)],
        Anis[Trial.index(i)],
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
        "Inactive mobile",
    )

fig, ax = plt.subplots(nrows=len(Regimes), ncols=1, figsize=[8, 10])
for a, Reg in zip(ax, Regimes):
    f = "chemsamndbiomasswithvel_temp_3_" + Reg + "_H_ZOOMED.png"
    im = plt.imread(d + f)
    a.imshow(im)
    if Reg == "Equal":
        a.annotate(
            "Medium flow",
            xy=(0, 0.5),
            xytext=(-ax[1].yaxis.labelpad - 10, 0),
            xycoords="axes fraction",
            textcoords="offset points",
            ha="left",
            va="center",
            rotation="vertical",
            size=15,
        )
    else:
        a.annotate(
            Reg + " flow",
            xy=(0, 0.5),
            xytext=(-ax[1].yaxis.labelpad - 10, 0),
            xycoords="axes fraction",
            textcoords="offset points",
            ha="left",
            va="center",
            rotation="vertical",
            size=15,
        )
    #    a.set_title (Reg+" flow", fontsize = 15)
    a.axis("off")
plt.savefig(
    "//tsclient/D/Saturated_flow/diffusion/homogeneous_profiles_chem_time.pdf",
    dpi=300,
    bbox_inches="tight",
)
plt.savefig(
    "//tsclient/D/Saturated_flow/diffusion/homogeneous_profiles_chem_time.png",
    dpi=300,
    bbox_inches="tight",
)

intsce = [80, 84, 73, 63]
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=[15, 10])
axlinear = ax.flatten()
for a, i in zip(axlinear, intsce):
    f = "chemsamndbiomasswithvel_temp_3_Equal_" + str(i) + "_ZOOMED.png"
    im = plt.imread(d + f)
    a.imshow(im)
    a.axis("off")
plt.savefig(
    "//tsclient/D/Saturated_flow/diffusion/Medium_heterogeneous_profiles_chem_time.pdf",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0,
)
plt.savefig(
    "//tsclient/D/Saturated_flow/diffusion/Medium_Lheterogeneous_profiles_chem_time.png",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0,
)

for Reg in Regimes:
    newd = d  # + Reg + "AR_"
    for j in Tforfpre[1:]:
        # Concentration profiles
        intsce = ["H", 76, 73, 80, 84, 63]
        for chem in ["DO", "DOC"]:
            dum = sk.shiftoxic_anoxic_temp(
                chem,
                Trial,
                intsce,
                newd,
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
            )
            picname = (
                "Z:/Saturated_flow/diffusion_transient/chemprofile_temp_"
                + chem
                + str(Tforfpre.index(j))
                + "_"
                + Reg
                + "_.png"
            )
            dum.savefig(picname, dpi=300, bbox_inches="tight", pad_inches=0)

        dummy = sk.concprofilewithtime_temp(
            Trial,
            intsce,
            newd,
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
        )
        picname = (
            "Z:/Saturated_flow/diffusion_transient/chemswithvel_temp_"
            + str(Tforfpre.index(j))
            + "_"
            + Reg
            + "_.png"
        )
        dummy.savefig(picname, dpi=300, bbox_inches="tight", pad_inches=0)

fig, ax = plt.subplots(ncols=len(Regimes), nrows=1, figsize=[14, 8])
for a, Reg in zip(ax, Regimes):
    f = "chemprofile_temp_DO3_" + Reg + "_.png"
    im = plt.imread(d + f)
    a.imshow(im)
    if Reg == "Equal":
        #        a.annotate("Medium flow", xy = (0,0.5), xytext = (-ax[1].yaxis.labelpad-10,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = 15)
        a.set_title("Medium flow", fontsize=10)
    else:
        #        a.annotate(Reg+" flow", xy = (0,0.5), xytext = (-ax[1].yaxis.labelpad-10,0), xycoords = 'axes fraction', textcoords='offset points', ha = 'left', va = 'center', rotation = 'vertical', size = 15)
        a.set_title(Reg + " flow", fontsize=10)
    #    a.set_title (Reg+" flow", fontsize = 15)
    a.axis("off")
plt.savefig(
    "//tsclient/D/Saturated_flow/diffusion/heterogeneous_profiles_DO_time.pdf",
    dpi=300,
    bbox_inches="tight",
)
plt.savefig(
    "//tsclient/D/Saturated_flow/diffusion/heterogeneous_profiles_DO_time.png",
    dpi=300,
    bbox_inches="tight",
)

# Heatmaps at select time intervals
intsce = [63, 84]
timepoints = [0, 180, 360, 1000]
hmgvarnames = ["DOC", "DO", "Ammonium", "Nitrate", "Sulphate"]
hmbgvarnames = ["Aerobes", "Ammonium oxidizers", "Nitrate reducers"]
indices = list(gvarnames.index(i) for i in hmgvarnames)
hmvars = list(vars[i] for i in indices)
hmbvars = AFbiomassvars[0:3]
Regimes = ["Slow", "Equal"]
for Reg in Regimes:
    Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
    for p in Tforfpre[1:]:
        d = r"Z:/Saturated_flow/diffusion_transient/" + p
    for k in intsce:
        df = np.load(d + fpre + str(k) + fpre + str(k) + "_df.npy")
        for time in timepoints:
            title = (
                "Variance "
                + str(Het[Trial.index(k)])
                + " : Anisotropy "
                + str(Anis[Trial.index(k)])
                + " at Time: "
                + str(time * 5)
                + " days"
            )
            heatmapc = stp.heatmapconcdist_temp(
                df, hmvars, k, hmgvarnames, d, fpre, title, time
            )
            picname = (
                "Z:/Saturated_flow/diffusion_transient/"
                + Reg
                + str(k)
                + "_"
                + "heatmap"
                + "_conc_"
                + str(time * 5)
                + ".png"
            )
            #            pictosave = heatmapc.get_figure()
            heatmapc.savefig(picname, dpi=300)
            heatmapb = stp.heatmapconcdist_temp(
                df, hmbvars, k, hmbgvarnames, d, fpre, title, time
            )
            picname = (
                "Z:/Saturated_flow/diffusion_transient/"
                + Reg
                + str(k)
                + "_"
                + "heatmap"
                + "_biomass_"
                + str(time * 5)
                + ".png"
            )
            #            pictosave = heatmapb.get_figure()
            heatmapb.savefig(picname, dpi=300)

# calculating average velocities
for Reg in ["Slow", "Equal"]:
    print(Reg)
    Tforfpre = [Reg + "AR_0", Reg + "AR_1", Reg + "AR_2", Reg + "AR_5"]
    for j in Tforfpre:
        print(j)
        if Tforfpre.index(j) == 0:
            pass
        else:
            intsce = [43, 44, 45, 52, 53, 54, 61, 62, 63, 73, 74]
            for i in intsce:
                newd = d + j
                df, massendtime, masstime, conctime, Velocity, head = sk.calcconcmasstime(
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
                )
                Velocity = df[2, 1:, :, :]
                print(i, ": ", np.mean(Velocity))

for Reg in ["SlowAR_5"]:
    d = r"Z:/Saturated_flow/diffusion_transient/" + Reg + "/"
    for j in range(len(Trial)):
        di = d + fpre + str(Trial[j]) + fsuf
        print(di)
        fwithd = di + filename
        print("Reading tec file....")
        size, steps, Headers, D = rdr.readTecfile(fwithd)
        print("Converting to array....")
        df = rdr.Converttomarr(D)
        print("Saving numpy array....")
        np.save(di + fpre + str(Trial[j]) + "_df", df)
        # Test for correct orientation of the data
        for i in range(np.shape(df)[0]):
            print(Headers[i + 3], np.mean(df[i, steps - 1, 0, :]))
