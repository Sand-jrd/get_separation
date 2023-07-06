#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 10:26:06 2023

@author: sand-jrd
"""
from vip_hci.fits import write_fits, open_fits
import matplotlib.pyplot as plt
import numpy as np
from mustard.utils import circle
from vip_hci.preproc import frame_pad, frame_shift, frame_crop
import os

import json
from vip_hci.stats.distances import cube_distance
from matplotlib.colors import LogNorm

# %% Folders 
star = "HD100546"
path = "../"+star+"/"

est_name = "GreeDS_full_H2_no_ref_30_1_1disks"

figsize=(9.5,8.5)

# %% Arguments 

plsca = 0.01226

perX = 99.8
perL = 99
perD = 99


shape = 501

rnd_nb = 2 # Precision
# %% Some functions

def get_point():
    
    pos = []
    
    pos = plt.ginput(n=1)
    plt.pause(1)
    pos = pos[0]

    siftx = (pos[0] - shape/2)*plsca
    sifty = (pos[1] - shape/2)*plsca
    
    print("Point coordinate in pixels: " + str(pos))
    print("Point coordinate in arcsec: " + str((siftx,sifty)))
    print("Raduis : " + str(np.sqrt(siftx**2 + sifty**2)) +"\n")

    tkt_sift = 10
    plt.text(pos[0]+np.sign(siftx)*tkt_sift,pos[1]+np.sign(sifty)*tkt_sift, str((round(siftx, rnd_nb), round(sifty, rnd_nb))), color="white", size=10)
    
    plt.plot([pos[0],shape/2], [pos[1],shape/2], color="white", lw=2 )
    plt.text((pos[0]+ shape/2)/2 +tkt_sift,(pos[1] + shape/2)/2 +tkt_sift, str(round(np.sqrt(siftx**2 + sifty**2), rnd_nb)), color="white", size=10, rotation=np.rad2deg(np.arctan(sifty/siftx)))
    

def draw_circle():
    plt.pause(1)

    pos = plt.ginput(n=1)
    plt.pause(1)

    pos = pos[0]
    siftx = (pos[0] - shape/2)*plsca
    sifty = (pos[1] - shape/2)*plsca
    
    print("Point coordinate in pixels: " + str(pos))
    print("Point coordinate in arcsec: " + str((siftx,sifty)))
    print("Raduis : " + str(np.sqrt(siftx**2 + sifty**2)) +"\n")
    
    tkt_sift = 15
    
    r = np.sqrt((pos[0] - shape/2)**2 + (pos[1] - shape/2) **2)
    
    ang = np.linspace(0,2*np.pi,50)
    
    cx = (r*np.cos(ang)) +shape/2
    cy = (r*np.sin(ang)) +shape/2
    
    plt.plot(cx,cy,color="white", lw=2 )
    plt.text(pos[0]+np.sign(siftx)*tkt_sift,pos[1]+np.sign(sifty)*tkt_sift, str(round(np.sqrt(siftx**2 + sifty**2), rnd_nb)), color="white", size=10, rotation=np.rad2deg(np.arctan(sifty/siftx))+np.pi)

def circle(shape: tuple, r: float, offset=(0.5, 0.5)):
    """ Create circle of 1 in a 2D matrix of zeros"
       
       Parameters
       ----------

       shape : tuple
           shape x,y of the matrix
       
       r : float
           radius of the circle
       offset : (optional) float
           offset from the center
       
       Returns
       -------
       M : ndarray
           Zeros matrix with a circle filled with ones
       
    """
    assert(len(shape) == 2 or len(shape) == 3)
    if isinstance(offset, (int, float)): offset = (offset, offset)

    nb_f  = shape[0]  if len(shape) == 3 else 0
    shape = shape[1:] if len(shape) == 3 else shape

    M = np.zeros(shape)
    w, l = shape
    for x in range(0, w):
        for y in range(0, l):
            if pow(x - (w // 2) + offset[0], 2) + pow(y - (l // 2) + offset[1], 2) < pow(r, 2):
                M[x, y] = 1

    if nb_f: M = np.tile(M, (nb_f, 1, 1))

    return M

def genere_args(Tm, M, Dtype):
    """Tm is the thing you want to plot.
    If multiplte plot, M is the ref to scale colorbar
    Dtype is to choose between disk, psf or residualplot arguments"""
    if Dtype == "X": 
        arg = {"cmap":"magma"}
        per = perX
        vmin = np.percentile(M[np.where(M>0)], 10)
        Tm[np.where(Tm<=0)] = vmin
        M[np.where(M<=0)] = vmin
        Tm[np.where(Tm<=0)] = vmin
        arg["norm"]= LogNorm(vmin=vmin, vmax=np.percentile(M, per))
    elif Dtype == "L": 
        arg = {"cmap":"jet"}
        per = perL
        arg["vmax"]=np.percentile(Tm, per)
    else :
        arg = {"cmap":"seismic"}
        per = perD
        arg["vmax"]=np.percentile(M, per)
    arg["X"]=Tm
    return arg

# %% Open stuff , reshape

X =  open_fits(path+est_name)[-1]

if X.shape[0] > shape : X = frame_crop(X,shape)

pup = circle((shape,shape),shape//2) -circle((shape,shape),10)
X =  pup*X

#%%` Plots

xticks_lab = np.array([-2,-1.5,-1,-0.5, 0, 0.5, 1, 1.5, 2])
yticks_lab = xticks_lab
xticks = (xticks_lab/plsca) + shape/2
yticks = (yticks_lab/plsca) + shape/2

fig = plt.figure(star, figsize=figsize)

heatmap = plt.imshow(**genere_args(X, X, "X"))
plt.gca().invert_yaxis()
plt.xticks(xticks, labels = xticks_lab)
plt.yticks(yticks, labels = yticks_lab)
plt.colorbar()


# %% Do something

get_point()
get_point()
draw_circle()

plt.savefig(star+".png")


