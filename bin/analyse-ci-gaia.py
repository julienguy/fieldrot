#!/usr/bin/env python


"""
precession done with ICRS -> FK5 coord transform from astropy ...
aberration done with ICRS -> GCRS coord transform from astropy ...
need to compare with Mike's code
precessiom, etc: https://desi.lbl.gov/DocDB/cgi-bin/private/ShowDocument?docid=4957
ADC corr:  https://desi.lbl.gov/DocDB/cgi-bin/private/ShowDocument?docid=5095

"""







import os, sys
import json
import argparse
import datetime
import numpy as np
from astropy.table import Table
import fitsio
import matplotlib.pyplot as plt
from fieldrot.rotlib import *


from astropy.coordinates import SkyCoord,FK5,GCRS
import astropy.time
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", type=str, nargs= "*", required = True, help="input centroids gaia file(s) ci-????????_gaia-a.fits")
parser.add_argument("-o", "--outfile", type=str, required = False, help="output file")
parser.add_argument("--plot",action="store_true")
parser.add_argument("-r", "--reffile", type=str, required = False, default="v0001/20190403/ci-00003852/ci-00003852_gaia-a.fits", help="input centroids gaia file to define transfo")
parser.add_argument("--precessed",action="store_true")

D2R=np.pi/180.

args = parser.parse_args()

# collect results,
header_keys=['NIGHT','EXPID','EXPTIME','MJD-OBS','TARGTRA','TARGTDEC','TARGTEL','TARGTAZ','SKYRA','SKYDEC','MOUNTHA','MOUNTAZ','MOUNTEL','MOUNTDEC','MNTOFFD','MNTOFFR','PARALLAC','AIRMASS','CI-T1','CI-T2','CI-T3','CI-T4','CI-T5','TAIRITMP','TAIROTMP','TAIRTEMP','TTRUETBT','HUMIDITY','PRESSURE','OUTTEMP','ADCCORR','HEXPOS','TCSST','ADC1PHI','ADC2PHI']
res={}
for k in header_keys :
    if k == 'HEXPOS' :
        res['HEXPOS1']=[]
        res['HEXPOS2']=[]
        res['HEXPOS3']=[]
        res['HEXPOS4']=[]
        res['HEXPOS5']=[]
        res['HEXPOS6']=[]
    else :
        res[k]=[]
res["THETA"]=[]
res["THETARMS"]=[]
res["CONTRAST"]=[]
res["NSTARS"]=[]
res["NSELECTED"]=[]
res["NCAMS"]=[]
res["TARGTHA"]=[]
res["POINTINGDX"]=[]
res["POINTINGDY"]=[]
res["SIGSEEING"]=[]

Tx=None
Ty=None

filenames = np.append(args.reffile,args.infile)

for filename in filenames :
    print("filename=",filename)
    table = fitsio.read(filename)
    night=os.path.basename(os.path.dirname(os.path.dirname(filename)))
    expid=os.path.basename(filename).replace("_gaia-a.fits","").replace("ci-","")
    #print("night=",night,"expid=",expid)
    
    ofilename="/project/projectdirs/desi/spectro/data/{}/{}/ci-{}.fits.fz".format(night,expid,expid)
    header = fitsio.read_header(ofilename,1)
    #print(header)
    if not "TARGTRA" in header :
        print ("no TARGTRA")
        continue
    
    cra=header["TARGTRA"]
    cdec=header["TARGTDEC"]
    ra=table["RA_EXTERN"]
    dec=table["DEC_EXTERN"]
    sigma=np.sqrt((table["SIG_MAJOR_PIX"]**2+table["SIG_MINOR_PIX"]**2)/2.)
    

    if args.precessed :

        gcrs=GCRS(obstime=astropy.time.Time(header["MJD-OBS"],format="mjd")) # precession? + aberration
        fk5=FK5(equinox=astropy.time.Time(header["MJD-OBS"],format="mjd")) # precession

        c = SkyCoord(cra,cdec, frame='icrs', unit='deg')
        c_ab = c.transform_to(gcrs)
        dra=c_ab.ra.value-c.ra.value
        ddec=c_ab.dec.value-c.dec.value
        c = SkyCoord(cra+dra,cdec+ddec, frame='icrs', unit='deg').transform_to(fk5)
        cra = c.ra.value
        cdec = c.dec.value
        
        c = SkyCoord(ra,dec, frame='icrs', unit='deg')
        c_ab = c.transform_to(gcrs)
        dra=c_ab.ra.value-c.ra.value
        ddec=c_ab.dec.value-c.dec.value
        c = SkyCoord(ra+dra,dec+ddec, frame='icrs', unit='deg').transform_to(fk5)
        ra = c.ra.value
        dec = c.dec.value
        
    cha=header["MOUNTHA"]
    #print("cha=",cha)
    #tmp=header["ST"].split(":")
    #lst1=(float(tmp[0])+float(tmp[1])/60.+float(tmp[2])/3600.)*180./24. # deg
    #cha=lst-cra
    lst=cha+cra
    #print(lst1,lst)
    #print(table.dtype.names)
    #sys.exit(12)
    x=table["XCENTROID"]
    y=table["YCENTROID"]
    cam=table["CAMERA"]
    cin=table["CI_NUMBER"]
    cont=table["CONTRAST"]
    ha=lst-ra

    # redefine cha cdec such that they are at the center of the
    # four outer CI cameras
    dx=[]
    dy=[]
    for c in np.unique(cam) :
        if c=="CIC": continue
        cdx=np.mean(ra[cam==c])
        cdy=np.mean(dec[cam==c])
        dx.append(cdx)
        dy.append(cdy)
    cra=np.mean(dx)
    cdec=np.mean(dy)
        
    # convert to tangent plane
    tx,ty = hadec2xy(ha,dec,cha,cdec)
    
    # fit transfo from pixel to tangent plane
    if Tx is None :
        Tx=dict()
        Ty=dict()
        for c in np.unique(cam) :
            A=np.zeros((3,3))
            B=np.zeros((3))
            ii=(cam==c) # stars in that camera
            w = np.ones(x[ii].size)
            A[0,0] = np.sum(w)
            A[1,0] = A[0,1] = np.sum(w*x[ii])
            A[1,1] = np.sum(w*x[ii]**2)
            A[1,2] = A[2,1] = np.sum(w*x[ii]*y[ii])
            A[2,0] = A[0,2] = np.sum(w*y[ii])
            A[2,2] = np.sum(w*y[ii]**2)
            Ai     = np.linalg.inv(A)
            B[0]   = np.sum(w*tx[ii])
            B[1]   = np.sum(w*x[ii]*tx[ii])
            B[2]   = np.sum(w*y[ii]*tx[ii])
            Tx[c] = Ai.dot(B)
            B[0]   = np.sum(w*ty[ii])
            B[1]   = np.sum(w*x[ii]*ty[ii])
            B[2]   = np.sum(w*y[ii]*ty[ii])
            Ty[c] = Ai.dot(B)


    # exclude the central camera
    #print(cin)
    #print(table.dtype.names)
    #print(table["CAMERA"][cin==3])
    #print(np.unique(cam))
    #sys.exit(12)
    cic=b'CIC'
    ok=cam!=cic
    x=x[ok]
    y=y[ok]
    tx=tx[ok]
    ty=ty[ok]
    ha=ha[ok]
    dec=dec[ok]
    cam=cam[ok]
    cont=cont[ok]
    ncam=np.unique(cam).size
    nstars=x.size
    mx=np.zeros(nstars) # "measured" x
    my=np.zeros(nstars) # "measured" y
    
    for c in np.unique(cam) :
        ii=(cam==c) # stars in that camera
        mx[ii]=Tx[c][0]+Tx[c][1]*x[ii]+Tx[c][2]*y[ii]
        my[ii]=Ty[c][0]+Ty[c][1]*x[ii]+Ty[c][2]*y[ii]

    # centering of field, with an actual rotation
    pointing_dx=0
    pointing_dy=0
    
    for loop in range(3) :
        # first average per camera to avoid artificially weighting more one side
        dx=[]
        dy=[]
        for c in np.unique(cam) :
            cdx=np.mean(mx[cam==c]-tx[cam==c])
            cdy=np.mean(my[cam==c]-ty[cam==c])
            #print(c,cdx,cdy)
            dx.append(cdx)
            dy.append(cdy)
        dx=np.mean(dx)
        dy=np.mean(dy)
        if loop == 0 :
            pointing_dx=dx
            pointing_dy=dy
        # print("{} dx= {} dy= {}".format(loop,dx,dy))
        cha2,cdec2=xy2hadec(-dx,-dy,cha,cdec)
        cha=cha2
        cdec=cdec2
        # convert to tangent plane
        tx,ty = hadec2xy(ha,dec,cha,cdec)
    
    if args.plot : # a plot for fun
        fig = plt.figure(os.path.basename(filename).split(".")[0])
        plt.subplot(1,2,1)
        plt.plot(tx,ty,"o")
        plt.plot(mx,my,"x")
        #plt.show()

    if 0 : # test
        plt.plot(tx,ty,"o")
        i=np.argmax(ha)
        plt.plot(tx[i],ty[i],"o",label="max HA")
        i=np.where(cam==b"CIW")[0]
        plt.plot(tx[i],ty[i],"o",label="CIW")
        i=np.argmax(dec)
        plt.plot(tx[i],ty[i],"o",label="max Dec")
        plt.xlabel("x (increasing HA, sky coordinates projected on instrument??)")
        plt.ylabel("y (increasing Dec)")
        a=20./180.*np.pi
        ca=np.cos(a)
        sa=np.sin(a)
        rx=ca*tx-sa*ty
        ry=sa*tx+ca*ty
        plt.plot(rx,ry,".",label="theta=20deg")
        plt.plot([0,1.5*ca],[0,1.5*sa],c="k")
        plt.plot([0,1.5],[0,0],"--",c="k")
        plt.text(0.5,0.06,"theta")
        plt.text(0.,-.5,"theta: rotation angle, counter-clockwise,",
                 horizontalalignment='center',verticalalignment='center')
        plt.text(0.,-.7,"of the star positions in the CI images",
                 horizontalalignment='center',verticalalignment='center')
        ab=np.arcsin((tx*ry-rx*ty)/np.sqrt((tx**2+ty**2)*(rx**2+ry**2)))/D2R # degree
        print("angle=",np.mean(ab))
        plt.legend()
        plt.show()
        sys.exit(12)
    
        
    
    angles = np.arcsin((tx*my-mx*ty)/np.sqrt((tx**2+ty**2)*(mx**2+my**2)))/D2R # degree
    theta  = np.median(angles)
    if args.plot :
        plt.subplot(1,2,2)
        plt.hist(angles[np.abs(angles)<0.2],bins=30)
        plt.axvline(theta,color="k")
        plt.xlabel("Theta (deg)")
        print("theta=",theta)
    plt.show()
    
    theta  = np.median(angles)
    rms  = 1.48*np.median(np.abs(angles-theta))
    ok   = np.abs(angles-theta)<3.*rms
    nselected = np.sum(ok)

    if nselected == 0 :
        print("no stars left...")
        continue
    
    theta  = np.mean(angles[ok])
    rms  = np.std(angles[ok])
    # check all keys are here
    bad=False
    for k in header_keys :
        if not k in  header :
            print("missing",k)
            if k == "ADC1PHI" or k == "ADC2PHI" :
                print("continue anyway")
            else :
                bad=True
                break
    if bad :
        continue
    
    for k in header_keys :   
        if k == "HEXPOS" :
            for i,v in enumerate(header[k]) :
                res["HEXPOS{}".format(i+1)].append(float(v))
        elif k == "TCSST" :
            vals=header[k].split(":")
            lst=(float(vals[0])+float(vals[1])/60.+float(vals[2])/3600.)*180/12.
            res[k].append(lst)

        elif k == "ADCCORR" :
            if header[k] is False :
                res[k].append(0)
            else :
                res[k].append(1)
        else :
            if k in header :
                res[k].append(header[k])
            else :
                res[k].append(-1)
    
    ha=res["TCSST"][-1]-res["TARGTRA"][-1]
    ha=ha%360.
    if ha>180 : ha -= 360.
    res["TARGTHA"].append(ha)
    res["THETA"].append(theta)
    res["THETARMS"].append(rms)
    res["NSTARS"].append(nstars)
    res["NSELECTED"].append(nselected)
    res["NCAMS"].append(ncam)
    if len(cont)>0 :
        res["CONTRAST"].append(np.min(cont))
    else :
        res["CONTRAST"].append(0.)
    res["POINTINGDX"].append(pointing_dx)
    res["POINTINGDY"].append(pointing_dy)
    res["SIGSEEING"].append(np.median(sigma))
    
    print(theta)


if args.outfile :
    if 0 :
        line=""
        tmp=[]
        for k in res.keys() :
            line=line+" "+k
            tmp.append(res[k])
            np.savetxt(args.outfile,np.array(tmp).T,header=line)
    else :
        file=open(args.outfile,"w")
        line="#"
        for k in res.keys() :
            line=line+" "+k
        file.write(line+"\n")
        file.write("# {}\n".format(str(datetime.datetime.now())))
        file.write("# From the files /project/projectdirs/desi/users/ameisner/CI/reduced/v0001/YYYYMMDD/ci-EXPID/ci-CEXPID_gaia-a.fits from Aaron Meisner\n")
        file.write("# Theta (in deg) measures the apparent rotation of stars in the CI images\n")
        file.write("# Reference exposure defining Theta=0 : {}\n".format(args.reffile))
        for i in range(len(res["EXPID"])) :
            line=""
            for  k in res.keys() :
                if k in ["NIGHT","EXPID","NSTARS","ADCCORR","NSELECTED","NCAMS"] :
                    line += "{:d} ".format(int(res[k][i]))
                else :
                    try :
                        line += "{:5.4f} ".format(res[k][i])
                    except ValueError :
                        line += "0 "
                        print(sys.exc_info())
            file.write(line+"\n")
        file.close()
    print("wrote",args.outfile)

