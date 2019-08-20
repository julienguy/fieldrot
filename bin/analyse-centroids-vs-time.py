#!/usr/bin/env python

import os, sys
import json
import argparse

import numpy as np
from astropy.table import Table
import fitsio
import matplotlib.pyplot as plt
from fieldrot.rotlib import *
import astropy.io.fits as pyfits

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", type=str, nargs = "*", required = True, help="input centroids file(s)")
parser.add_argument("-o", "--outfile", type=str, required = False, help="output file")
parser.add_argument("--plot-per-exposure",action="store_true")
parser.add_argument("--title-suffix",type=str, default=None)



args = parser.parse_args()

# collect results
header_keys=['NIGHT','EXPID','MOUNTHA','TARGTRA','TARGTDEC','AIRMASS']
res={}
for k in header_keys :
    res[k]=[]
res["DX"]=[]
res["DY"]=[]


# correlation function 
nc=100  
s_cdra=np.zeros(nc)
s_cddec=np.zeros(nc)
s_entries=np.zeros(nc)

pid_params=None
guidtime=None
guidecam=None
timelapse=11.
night=None
nights=None

for filename in args.infile :
    
    print(filename)
    
     #- Load centroids
    cx = json.load(open(filename))

    #- ra, dec aren't in centroids file, so see if there is a guide file
    expid = cx['expid']
    dirname = os.path.dirname(filename)
    guidefile = os.path.join(dirname, 'guide-rois-{:08d}.fits.fz'.format(expid))
    if not os.path.exists(guidefile):
        guidefile = os.path.join(dirname, 'guide-{:08d}.fits.fz'.format(expid))
    print(guidefile)
    hdr = None
    if os.path.exists(guidefile):
        try : 
            hdr = fitsio.read_header(guidefile, 0)
        except Exception as e:
            print("failed to read guide file",guidefile)
    else :
        print("no guide file",guidefile)
        
    if hdr :
        print("TCSPIDEC=",hdr['TCSPIDEC'])
        current_pid_params = "{} {} {}".format(hdr['TCSGRA'],hdr['TCSMFRA'],hdr['TCSPIDEC'])
        if pid_params is not None :
            if current_pid_params !=  pid_params :
                print("DIFFERENT PID params:", current_pid_params, pid_params)
                print("skipping",filename)
                continue

        if not 'GUIDTIME' in hdr :
            print ("no GUIDTIME?")
            hdr['GUIDTIME']=5

        if guidtime is not None :
            if hdr['GUIDTIME'] != guidtime :
                print("DIFFERENT GUIDTIME:",hdr['GUIDTIME'],guidtime)
                print("skipping",filename)
                continue

        if guidecam is not None :
            if hdr['GUIDECAM'] != guidecam :
                print("DIFFERENT GUIDECAM:",hdr['GUIDECAM'],guidecam)
                #print("skipping",filename)
                #continue

        if night is None :
            night  = hdr['NIGHT']
            nights = "{}".format(night)
        if night is not None :
            if hdr['NIGHT'] != night :
                night = hdr['NIGHT']
                nights += ",{}".format(night)

        if False and guidtime is None :
            print("look for the time sequence of one guiding image to get an accurate measure of the time lapse")

            h = pyfits.open(guidefile)
            for hdu in ["CINT","CIWT","CIST","CIET","CICT"] :
                if hdu in h :
                    table=h[hdu].data
                    #print(table)
                    tmp=table["MJD-OBS"]
                    tmp-=tmp[0]
                    tmp *= 24.*3600.
                    #print(tmp)
                    timelapse = np.median(np.gradient(tmp))
                    print("timelapse=",timelapse,"sec")
                    #times=h[hdu].data["TIME-OBS"]
                    #print(times)
                    #sys.exit(12)




        pid_params = current_pid_params
        guidtime = hdr['GUIDTIME']
        guidecam = hdr['GUIDECAM']
    else :
        guidtime = 5 # guess
        guidecam = "unknown"
        tcspidec = "unknown"
        timelapse = 11 # guess



    print("number of frames x cameras =",len(cx['frames'].values()))
    if len(cx['frames'].values())<3 :
        print("skip this one")
        continue

    #- Convert json into Table; convert pixel offsets to arcsec
    rows = list()
    
    
    #- platescale values taken from Echo22 design in
    #- desimodel/data/focalplane/platescale.txt
    center_platescale = 67.5 # um/arcsec
    radial_platescale = 76.3 # um/arcsec at r=400 mm
    az_platescale = 70.3 # um/arcsec at r=400m
    pixsize = 9 # um

    for i, frame in enumerate(cx['frames'].values()):
        for key, value in frame.items():
            if not key.startswith('CI'):
                continue
            ci = key[0:3]
            star = int(key.split('_')[1])+1000*(ci=="CIN")+2000*(ci=="CIW")+3000*(ci=="CIS")+4000*(ci=="CIE")
            x_error = value['x_error']
            y_error = value['y_error']
            if ci == 'CIC':
                ra_err = pixsize * x_error / center_platescale
                dec_err = -pixsize * y_error / center_platescale
            elif ci == 'CIS':
                ra_err = -pixsize * x_error / az_platescale
                dec_err = pixsize * y_error / radial_platescale
            elif ci == 'CIN':
                ra_err = pixsize * x_error / az_platescale
                dec_err = -pixsize * y_error / radial_platescale
            elif ci == 'CIE':
                dec_err = -pixsize * x_error / az_platescale
                ra_err = -pixsize * y_error / radial_platescale
            elif ci == 'CIW':
                dec_err = pixsize * x_error / az_platescale
            rows.append(dict(
                frame=i, ci=ci, star=star,
                dx=x_error, dy=y_error,
                dra=ra_err, ddec=dec_err,
                tcs_dra=frame["tcs_correction_ra"],
                tcs_ddec=frame["tcs_correction_dec"])
                    )
                
    if len(rows)==0 : continue
    
    data = Table(rows,
                 names=['frame','ci','star','dx','dy','dra','ddec','tcs_dra','tcs_ddec'],
                 dtype=(int, str, int, float, float, float, float, float, float))
    
    frames=np.unique(data['frame'])
    nf=len(frames)
    print("number of frames=",nf)
    if nf<3 :
        continue
    #print(data)
    dra=np.zeros(nf)
    ddec=np.zeros(nf)
    tcs_dra=np.zeros(nf)
    tcs_ddec=np.zeros(nf)
    era=np.zeros(nf)
    edec=np.zeros(nf)
    for f,frame in enumerate(frames) :
        ii=(data['frame']==frame)
        dra[f]  = np.mean(data['dra'][ii])
        ddec[f] = np.mean(data['ddec'][ii])
        tcs_dra[f]  = np.mean(data['tcs_dra'][ii])
        tcs_ddec[f] = np.mean(data['tcs_ddec'][ii])
        if ii.size>1 :
            era[f] = np.std(data['dra'][ii])/np.sqrt(ii.size-1)
            edec[f] = np.std(data['ddec'][ii])/np.sqrt(ii.size-1)

    
    print("rms along ra = ",np.std(dra),"arcsec")
    print("rms along dec = ",np.std(ddec),"arcsec")
    print("median stat. error along ra = ",np.median(era),"arcsec")
    print("median stat. error along dec = ",np.median(edec),"arcsec")
    
    title=os.path.basename(filename).split(".")[0]
    if args.plot_per_exposure :
        plt.figure(title+"-ra")
        plt.errorbar(frames*timelapse,dra,era,fmt="o-",color="k")
        plt.plot(frames*timelapse,-tcs_dra,"--",color="r")
        plt.grid()
        plt.xlabel("time (sec) ({:2.1f} sec between frames)".format(timelapse))
        plt.ylabel("dx (arcsec) (along ra)")
        print("rms dx(dra)=",np.std(dra))
    
        plt.figure(title+"-dec")
        plt.errorbar(frames*timelapse,ddec,edec,fmt="o-",color="k")
        plt.plot(frames*timelapse,-tcs_ddec,"--",color="r")
        plt.xlabel("time (sec) ({:2.1f} sec between frames)".format(timelapse))
        plt.ylabel("dy (arcsec) (along dec)")
        plt.grid()
        print("rms dy(ddec)=",np.std(ddec))

        if False :
            for ci in np.unique(data['ci']) :
                ii=(data['ci']==ci)
                print(ci)
                plt.figure(title+"-ra")
                plt.plot(data['frame'][ii],data['dra'][ii],"o-",alpha=0.4)
                plt.figure(title+"-dec")
                plt.plot(data['frame'][ii],data['ddec'][ii],"o-",alpha=0.4)

    if np.std(dra)<0.00001 or np.std(dra)>12 or np.std(ddec)>12 :
        print("don't add {} because odd rms of residuals".format(title))
        continue
    
    dra -= np.mean(dra)
    ddec -= np.mean(ddec)
    for i in range(nc) :
        if nf>i :
            s_cdra[i] += np.mean(dra[0:nf-i]*dra[i:nf])
            s_cddec[i] += np.mean(ddec[0:nf-i]*ddec[i:nf])
            s_entries[i] += 1
    
    
if np.sum(s_entries>0) == 0 :
    print("empty")
    sys.exit(12)

cdra=np.zeros(s_cdra.shape)
cddec=np.zeros(s_cddec.shape)
cdra[s_entries>0] = s_cdra[s_entries>0]/s_entries[s_entries>0]
cddec[s_entries>0] = s_cddec[s_entries>0]/s_entries[s_entries>0]
    
title="guiding-errors-correlation"
if args.title_suffix is not None :
    title += "-"+args.title_suffix
plt.figure(title)


times=np.arange(cdra.size)*timelapse
a=plt.subplot(1,1,1)
a.plot(times,cdra,label="along RA")
a.plot(times,cddec,label="along Dec")
#a.legend(loc="upper right",title="{} with PID={} guiding with {}".format(nights,pid_params,guidecam))
a.legend(loc="upper right",title="PID={}".format(pid_params))
a.grid()
a.set_xlabel("time (sec) ({:2.1f} sec between frames)".format(timelapse))
a.set_ylabel("correlation (arcsec**2)")

if 0 :
    a=plt.subplot(2,1,2)
    a.plot(np.sqrt(np.abs(cdra)))
    a.plot(np.sqrt(np.abs(cddec)))
    a.grid()
    a.set_xlabel("frame")
    a.set_ylabel("sqrt(|correlation|) (arcsec)")

plt.show()
    
