# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 14:07:23 2021

@author: khurana
"""

import os
#import plots.saturated_steady_state as sssp
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib as mpl  
from data_reader import data_processing as proc
import matplotlib.gridspec as gridspec
import cv2
from PIL import Image

#File paths
DIR = "E:/Zenodo_temporal_dynamics"
results_dir = os.path.join(DIR,"images")
raw_directory = DIR
os.chdir(raw_directory)
#Standard color and font options
legendkw = {'fontsize' : 14}
labelkw = {'labelsize' : 14}
secondlabelkw = {'labelsize' : 16}
suptitlekw = {'fontsize' : 18}
titlekw = {'fontsize' : 16}
mpl.rc('font',family='Arial')

Regimes = ["Slow", "Equal", "Fast"]
trialist = proc.masterscenarios()
Trial = ["50", "73", "63"]
species = proc.speciesdict("Saturated")
gvarnames = ["DO", "DOC", "Ammonium", "Nitrate"]
velindex = 2
colorscheme = 'YlGnBu'
columntitles = ["Velocity\ndistribution pattern", "Slow\nflow", "Medium\nflow", "Fast\nflow"]

def plotfigure(timidx, timseries):
    spmin = np.zeros((len(gvarnames),3))
    spmax = np.zeros((len(gvarnames),3))
    velmin = np.zeros((len(Trial),1))
    velmax = np.zeros((len(Trial),1))
    for t in Trial:
        for r in Regimes:
            file = os.path.join(raw_directory, r + "AR_"+timseries+"_NS-A"+str(t)+"_df.npy")
            data = np.load(file)
            for g in gvarnames:
                spmin[gvarnames.index(g), Regimes.index(r)] = np.min(data[species[g]['TecIndex'],1:,:,:])
                spmax[gvarnames.index(g), Regimes.index(r)] = np.max(data[species[g]['TecIndex'],1:,:,:])
            
    for t in Trial:
        file = os.path.join(raw_directory, "EqualAR_"+timseries+"_NS-A"+str(t)+"_df.npy")
        data = np.load(file)
        velmin[Trial.index(t)] = np.min(abs(data[velindex,1:,:,:]))
        velmax[Trial.index(t)] = np.max(abs(data[velindex,1:,:,:]))
    
    fig = plt.figure(figsize=(14, 14))
    outer = gridspec.GridSpec(3, 4, wspace=0.2, hspace=0.2)
    pad = 200
    for t in Trial:
        file = os.path.join(raw_directory, "EqualAR_"+timseries+"_NS-A"+str(t)+"_df.npy")
        data = np.load(file)
        left = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer[4*Trial.index(t)], wspace=0.3, hspace=0.1)
        axe = plt.Subplot(fig, left[0])
        velocity = abs(data[velindex, timidx, :, :])
        sns.heatmap(velocity, cmap = colorscheme, ax = axe, cbar = False,
                    vmax = velmax[Trial.index(t)], vmin = velmin[Trial.index(t)])
        axe.set_ylabel ("Variance: " + str(trialist[t]["Het"])+ " &\nAnisotropy: " + str(trialist[t]["Anis"]),
                        rotation = "vertical", ha = "center", **titlekw)
        axe.set_xticks([])
        axe.set_yticks([])
        fig.add_subplot(axe)
        for r in Regimes:
            i = Trial.index(t)*len(Regimes) + Regimes.index(r) + Trial.index(t) + 1
            if i%4 != 0:
                inner = gridspec.GridSpecFromSubplotSpec(2, 2,
                                                     subplot_spec=outer[i], wspace=0.15, hspace=0.17)
                file = os.path.join(raw_directory, r + "AR_"+timseries+"_NS-A"+str(t)+"_df.npy")
                data = np.load(file)
                for g in gvarnames:
                    axe = plt.Subplot(fig, inner[gvarnames.index(g)])
                    sns.heatmap (data[species[g]["TecIndex"], timidx, :, :], cmap = colorscheme, ax= axe,
                                 cbar = False, #fig.add_axes([0.25*Regimes.index(r), 0.1*gvarnames.index(g), 0.03, 0.4]),
                                 vmin = spmin[gvarnames.index(g), Regimes.index(r)],
                                 vmax = spmax[gvarnames.index(g),Regimes.index(r)])
                    axe.set_title(g, ha = "center", **legendkw)
                    axe.set_xticks([])
                    axe.set_yticks([])
                    fig.add_subplot(axe)
    for a in range(4):
        plt.annotate(columntitles[a], xy=(0.15, 0.92), xytext=(0.0 + pad*a, 0),
                     xycoords='figure fraction', textcoords='offset points',
                     size='large', ha='center', va='baseline',
                     **titlekw)
    plt.annotate(str(timidx*5) + " days", xy=(0.9, 0.97), xytext=(0.0, 0),
                     xycoords='figure fraction', textcoords='offset points',
                     size='large', ha='center', va='baseline',
                     **titlekw)
    fig.show()
    
    return fig
    
timseries = np.asarray(list(range(1097)))
timeseries = ["0","1","2","5"]
series = "1"

for time in timseries[835:]:
    plotfigure(time, series)
    picname = os.path.join(results_dir,"tim1","anime_"+str(time)+".png")
    plt.savefig(picname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)

firstimg = cv2.imread(os.path.join(results_dir,"tim1","anime_1005.png"))
height,width,layers= firstimg.shape

# Resizing of the images to give
# them same width and height 
#for file in os.listdir('.'):
#    if file.endswith(".jpg") or file.endswith(".jpeg") or file.endswith("png"):
for j in timseries[:1006]:
    imgfile = os.path.join(results_dir,"tim1","anime_"+str(j)+".png")
    # opening image using PIL Image
    im = Image.open(imgfile) 
    # resizing 
    imResize = im.resize((width, height), Image.ANTIALIAS) 
    imResize.save(imgfile, 'png', dpi = (300,300)) # setting quality
        
fourcc = cv2.VideoWriter_fourcc(*'H264')#*'mp4v')
videofile = os.path.join(results_dir, "tim1",'video_till1005.mp4')
video=cv2.VideoWriter(videofile, fourcc, 12,(width,height))

for j in timseries[:1006]:
    imgfile = os.path.join(results_dir,"tim1","anime_"+str(j)+".png")
    img = cv2.imread(imgfile)
    video.write(img)

cv2.destroyAllWindows()
video.release()