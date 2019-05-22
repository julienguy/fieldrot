#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys


# pixsize = 9um, platescale ~ 70 um/arcsec -> pixsize = 0.13 arcsec
# for seeing = 1 arcsec, sigma = 1./2.35/0.13 ~ 3 pix

def xysig(stamp,first_guess_sigma=3.) :
    n0=stamp.shape[0]
    n1=stamp.shape[1]
    
    exclusion_radius = min(n0,n1)/4.
    x = np.tile(np.arange(n1),(n0,1))
    y = np.tile(np.arange(n0),(n1,1)).T
    
    # first guess of center
    cx = n1/2.
    cy = n0/2.
    
    # saved values
    scx=0.
    scy=0.
    
    for i in range(40) :
        r = np.sqrt((x-cx)**2+(y-cy)**2)
        bkg = np.median(stamp[r>exclusion_radius])
        nstamp = stamp-bkg
        nstamp = nstamp/float(np.sum(nstamp))
        cx = np.sum(x*nstamp)
        cy = np.sum(y*nstamp)
        diff = max(np.abs(cx-scx),np.abs(cy-scy))
        scx = cx+0
        scy = cy+0
        print(i,cx,cy,diff)
        if diff < 0.005 : break
    
    pstamp = nstamp*(nstamp>0)
    pstamp /= float(np.sum(pstamp))
    
    

    # first guess
    sx = first_guess_sigma + 0.
    sy = first_guess_sigma + 0.
    
    # saved values
    ssx=0.
    ssy=0.

    # iterative measure, weighted by a Gaussian 
    # to limit impact of noise
    for i in range(40) :
        weight = np.exp(-(x-cx)**2/sx**2/2-(y-cy)**2/sy**2/2)
        norme = np.sum(pstamp*weight)
        # sqrt(2) because gaussian**2 here
        cx = np.sum(x*pstamp*weight)/norme
        cy = np.sum(y*pstamp*weight)/norme
        sx = np.sqrt(2.)*np.sqrt(np.sum((x-cx)**2*pstamp*weight)/norme)
        sy = np.sqrt(2.)*np.sqrt(np.sum((y-cy)**2*pstamp*weight)/norme)
        diff = max(np.abs(sx-ssx),np.abs(sy-ssy))
        ssx = sx+0
        ssy = sy+0
        print(i,cx,cy,sx,sy,diff)
        if diff < 0.005 : break
    
    #sx = np.sqrt(np.sum((x-cx)**2*pstamp))
    #sy = np.sqrt(np.sum((y-cy)**2*pstamp))
    #print(cx,cy,sx,sy)
    return cx,cy,sx,sy




def main() :
    # test 
    
    n0=21
    n1=31
    stamp=np.zeros((n0,n1))
    xc = 21
    yc = 10.
    x = np.tile(np.arange(n1),(n0,1))-xc
    y = np.tile(np.arange(n0),(n1,1)).T-yc
    sx = 3
    sy = 4
    stamp = np.exp(-x**2/sx**2/2-y**2/sy**2/2)
    cx,cy,sx,sy = xysig(stamp)
    print(cx,cy,sx,sy)
    
    noise = 0.1*np.random.normal(size=stamp.shape)
    cx,cy,sx,sy = xysig(stamp+noise)
    print(cx,cy,sx,sy)
    
    plt.subplot(2,1,1)
    plt.imshow(stamp,origin=0)
    plt.subplot(2,1,2)
    plt.imshow(stamp+noise,origin=0)
    plt.show()


if __name__ == '__main__' :

    main()


    
    
