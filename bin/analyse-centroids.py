#!/usr/bin/env python

import os, sys
import json
import argparse

import numpy as np
from astropy.table import Table
import fitsio
import matplotlib.pyplot as plt
from fieldrot.rotlib import *


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", type=str, nargs = "*", required = True, help="input centroids file(s)")
parser.add_argument("-o", "--outfile", type=str, required = False, help="output file")
parser.add_argument("--plot",action="store_true")


args = parser.parse_args()

# collect results
header_keys=['NIGHT','EXPID','MOUNTHA','TARGTRA','TARGTDEC','AIRMASS']
res={}
for k in header_keys :
    res[k]=[]
res["ROTANGLE"]=[]
res["ROTANGLERMS"]=[]
res["NFRAMES"]=[]
res["NSTARS"]=[]
res["NCAMS"]=[]
res["RMSX"]=[]
res["RMSY"]=[]
res["RMSXB"]=[]
res["RMSYB"]=[]

for filename in args.infile :
    
     #- Load centroids
    cx = json.load(open(filename))

    #- ra, dec aren't in centroids file, so see if there is a guide file
    expid = cx['expid']
    dirname = os.path.dirname(filename)
    guidefile = os.path.join(dirname, 'guide-{:08d}.fits.fz'.format(expid))
    if os.path.exists(guidefile):
        hdr = fitsio.read_header(guidefile, 0)
    else:
        hdr = None
        print("no guide file",guidefile)
        continue
    
    pmfile = os.path.join(dirname, 'pm-{:08d}.fits'.format(expid))
    if not os.path.exists(pmfile):
        print("no PM file",pmfile)
        continue 
    try :
        PMGSTARS = fitsio.read(pmfile,"PMGSTARS")
    except OSError :
        print("error reading",pmfile)
        continue

    #- Convert json into Table; convert pixel offsets to arcsec
    rows = list()
    nframes = len(cx['frames'])
    for i, frame in enumerate(cx['frames'].values()):
        #- Use CIX and ROI=1 for the combined guider offset values
        ci = 'CIX'
        roi = 1
        
        x1=[]
        y1=[]
        x2=[]
        y2=[]
        
        #- Now loop over whatever cameras were included in this file
        for key, value in frame.items():
            if not key.startswith('CI'):
                continue
            ci = key[0:3]
            roi = int(key.split('_')[1])
            dico = dict(frame=i, ci=ci, roi=roi)
            for k in ['x_expected','y_expected','x_centroid','y_centroid'] :
                dico[k]=value[k]
            rows.append(dico)

    if len(rows)==0 : continue

    data = Table(rows,
             names=['frame','ci','roi','x_expected','y_expected','x_centroid','y_centroid'],
             dtype=(int, str, int, float, float, float, float))
    
    cameras=np.unique(data['ci'])
    ncam=np.sum(cameras!="CIC")
    if ncam<2 :
        print("ignore data with only",cameras)
        continue

    frames=np.unique(data['frame'])
    
    # central RA Dec
    cra=hdr["TARGTRA"]
    cdec=hdr["TARGTDEC"]
    # central HA
    cha=hdr["MOUNTHA"]
    lst=cha+cra
    
    nentries=len(data['x_expected'])
    data['tx_expected']=np.zeros(nentries,dtype=float)
    data['ty_expected']=np.zeros(nentries,dtype=float)
    data['tx_centroid']=np.zeros(nentries,dtype=float)
    data['ty_centroid']=np.zeros(nentries,dtype=float)
    

    #print(PMGSTARS.dtype.names)
    for cam in cameras :
        ii=(data['ci']==cam)
        stars=np.unique(data['roi'][ii])
        print("{} nstars={}".format(cam,len(stars)))
        
        

        # match to pm catalogPMGSTARS["GFA_LOC"]
        #print(np.unique(PMGSTARS["GFA_LOC"]))
        ii=np.where(PMGSTARS["GFA_LOC"]==bytes(cam, 'utf-8'))[0]
        print("{} n pm stars={}".format(cam,np.sum(ii)))
        # pixel coordinates ?
        py=PMGSTARS["COL"][ii]
        px=PMGSTARS["ROW"][ii]
        
        pra=PMGSTARS["RA"][ii]
        pdec=PMGSTARS["DEC"][ii]
        # convert to ha
        pha=lst-pra
        # convert to tangent plane
        ptx,pty = hadec2xy(pha,pdec,cha,cdec)

        # fit transfo from pixel to tangent plane
        w = np.ones(px.size)
        A=np.zeros((3,3))
        B=np.zeros((3))
        A[0,0] = np.sum(w)
        A[1,0] = A[0,1] = np.sum(w*px)
        A[1,1] = np.sum(w*px**2)
        A[1,2] = A[2,1] = np.sum(w*px*py)
        A[2,0] = A[0,2] = np.sum(w*py)
        A[2,2] = np.sum(w*py**2)
        Ai     = np.linalg.inv(A)
        B[0]   = np.sum(w*ptx)
        B[1]   = np.sum(w*px*ptx)
        B[2]   = np.sum(w*py*ptx)
        
        Tx = Ai.dot(B)
        B[0]   = np.sum(w*pty)
        B[1]   = np.sum(w*px*pty)
        B[2]   = np.sum(w*py*pty)
        Ty = Ai.dot(B)
        

        #plt.figure("all")
        #plt.plot(ptx,pty,"o",alpha=0.2)
        #tx=Tx[0]+Tx[1]*px+Tx[2]*py
        #ty=Ty[0]+Ty[1]*px+Ty[2]*py
        #plt.plot(tx,ty,"x",c="gray",alpha=0.2)
        
        

       
            
        ok=(data['x_centroid']>=-10000)&(data['ci']==cam)
        
        px=data['x_expected'][ok]
        py=data['y_expected'][ok]
        data['tx_expected'][ok]=Tx[0]+Tx[1]*px+Tx[2]*py
        data['ty_expected'][ok]=Ty[0]+Ty[1]*px+Ty[2]*py
        
        px=data['x_centroid'][ok]
        py=data['y_centroid'][ok]
        data['tx_centroid'][ok]=Tx[0]+Tx[1]*px+Tx[2]*py
        data['ty_centroid'][ok]=Ty[0]+Ty[1]*px+Ty[2]*py
        
        #plt.plot(tx,ty,"+")
     
    if np.sum( data['tx_expected']==0 )>0 :
        print("problem with tx_expected=0")
        continue
    
    
    #print("I can now measure the rotation angle")
    data['rotangle']=np.zeros(nentries,dtype=float)
    
    frames=np.unique(data['frame'])
    rotangle_of_frames=np.zeros(frames.size)
    rotangle_rms_of_frames=np.zeros(frames.size)
    mx_of_frames=np.zeros(frames.size)
    my_of_frames=np.zeros(frames.size)
    rmsx_of_frames=np.zeros(frames.size)
    rmsy_of_frames=np.zeros(frames.size)
    rmsxb_of_frames=np.zeros(frames.size)
    rmsyb_of_frames=np.zeros(frames.size)
    for j,frame in enumerate(frames) :
        ii=(data['frame']==frame)&(data['tx_expected']!=0)&(data['ci']!="CIC")
        if np.sum(ii)==0 : continue
        #print(np.unique(data['ci'][ii]))
        x1=data['tx_expected'][ii]
        y1=data['ty_expected'][ii]
        x2=data['tx_centroid'][ii]
        y2=data['ty_centroid'][ii]
        dx = np.mean(x2-x1)
        dy = np.mean(y2-y1)
        x2 -= dx
        y2 -= dy
        angles = np.arcsin((x1*y2-x2*y1)/np.sqrt((x1**2+y1**2)*(x2**2+y2**2)))/D2R
        angle  = np.mean(angles)
        rotangle_of_frames[j] = angle
        n=len(angles)
        if n>1 :
            rotangle_rms_of_frames[j] = np.std(angles)/np.sqrt(n-1.)
        data['rotangle'][ii] = angle
    if np.sum(rotangle_of_frames==0)>0 :
        print("skip this one")
        continue
    mangle = np.median(rotangle_of_frames[rotangle_of_frames!=0])
    
    
    for j,frame in enumerate(frames) :
        mdx2=[]
        mdy2=[]
        mdx3=[]
        mdy3=[]
        bad=False
        for cam in cameras :
            if cam=="CIC" : continue
            ii=(data['ci']==cam)&(data['frame']==frame)&(data['tx_expected']!=0)
            if np.sum(ii)<2 : # at least 2 stars
                bad=True
                break
            x1=data['tx_expected'][ii]
            y1=data['ty_expected'][ii]
            x2=data['tx_centroid'][ii]-dx
            y2=data['ty_centroid'][ii]-dy
            x3=x2*np.cos(mangle*D2R)+y2*np.sin(mangle*D2R)
            y3=x2*np.sin(-mangle*D2R)+y2*np.cos(mangle*D2R)
            angle12  = np.mean(np.arcsin((x1*y2-x2*y1)/np.sqrt((x1**2+y1**2)*(x2**2+y2**2)))/D2R)
            angle13  = np.mean(np.arcsin((x1*y3-x3*y1)/np.sqrt((x1**2+y1**2)*(x3**2+y3**2)))/D2R)
            #print(cam,angle*3600,angle12*3600,"->",angle13*3600)
            mdx2.append(np.mean(x2-x1))
            mdy2.append(np.mean(y2-y1))
            mdx3.append(np.mean(x3-x1))
            mdy3.append(np.mean(y3-y1))
            
        if bad : 
            continue
        
        if len(mdx2)>1 :
            mx_of_frames[j] = np.mean(mdx2)
            my_of_frames[j] = np.mean(mdy2)
            rmsx_of_frames[j] = np.std(mdx2)
            rmsy_of_frames[j] = np.std(mdy2)
            rmsxb_of_frames[j] = np.std(mdx3)
            rmsyb_of_frames[j] = np.std(mdy3)

        plt.figure("mdx")
        nn=len(mdy2)
        for i in range(nn) :
            plt.scatter(j,mdy2[i],color=["b","g","r","k"][i])
    
        
        
    if args.plot :
        figurename = os.path.basename(filename).split(".")[0]
        plt.figure(figurename)
        plt.errorbar(frames,rotangle_of_frames*3600.,rotangle_rms_of_frames*3600.,fmt="o")
        plt.xlabel("Frame number")
        plt.ylabel("Rotation angle (centroid-PM) (arcsec)")
        plt.grid()
        #plt.plot(frames,rotangle_of_frames,"o")

    mangle = np.median(rotangle_of_frames)
    rms = np.std(rotangle_of_frames)
    for k in header_keys :
        res[k].append(hdr[k])
    res["ROTANGLE"].append(mangle)
    res["ROTANGLERMS"].append(rms)
    res["NFRAMES"].append(nframes)
    
    nstars=0
    for cam in cameras :
        if cam=="CIC" : continue
        nstars += np.unique(data['roi'][data['ci']==cam]).size
    res["NSTARS"].append(nstars)
    res["NCAMS"].append(ncam)
    
    res["RMSX"].append(np.mean(rmsx_of_frames))
    res["RMSY"].append(np.mean(rmsy_of_frames))
    res["RMSXB"].append(np.mean(rmsxb_of_frames))
    res["RMSYB"].append(np.mean(rmsyb_of_frames))
    print("angle= {} deg = {} arcsec, nstars={}".format(mangle,mangle*3600,nstars))

if args.outfile :
    line=""
    tmp=[]
    for i,k in enumerate(res.keys()) :
        line=line+" "+k
        tmp.append(res[k])
    np.savetxt(args.outfile,np.array(tmp).T,header=line)
    print("wrote",args.outfile)



if args.plot :
    plt.show()

