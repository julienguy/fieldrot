#!/usr/bin/env python

import fitsio
import numpy as np
import matplotlib.pyplot as plt
import glob
import os,sys
import argparse

from fieldrot.rotlib import *


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="fit field rotation")

parser.add_argument('-i','--infile', type = str, default = None, required = True, nargs="*", help = 'path to many ci-xxxxx_catalog.fits from Aaron')
parser.add_argument('-o','--outfile',type = str, default = None, required = True, help="output ascii file")

args = parser.parse_args()

if  "DESI_SPECTRO_DATA" not in os.environ :
    print("need DESI_SPECTRO_DATA env. variable, like /project/projectdirs/desi/spectro/data at NERSC")
    sys.exit(12)


"""
('xcen_init',
 'ycen_init',
 'detmap_xmin',
 'detmap_xmax',
 'detmap_ymin',
 'detmap_ymax',
 'detmap_peak',
 'no_centroid_refinement',
 'centroid_refinement_fail',
 'xcentroid',
 'ycentroid',
 'sig_major_pix',
 'sig_minor_pix',
 'ellipticity',
 'pos_angle',
 'dxcentroid',
 'dycentroid',
 'min_edge_dist_pix',
 'dq_flags',
 'aper_sum_bkgsub_0',
 'aper_bkg_0',
 'aperture_sum_err_0',
 'aper_sum_bkgsub_1',
 'aper_bkg_1',
 'aperture_sum_err_1',
 'aper_sum_bkgsub_2',
 'aper_bkg_2',
 'aperture_sum_err_2',
 'aper_sum_bkgsub_3',
 'aper_bkg_3',
 'aperture_sum_err_3',
 'aper_sum_bkgsub_4',
 'aper_bkg_4',
 'aperture_sum_err_4',
 'aper_sum_bkgsub_5',
 'aper_bkg_5',
 'aperture_sum_err_5',
 'aper_sum_bkgsub_6',
 'aper_bkg_6',
 'aperture_sum_err_6',
 'aper_sum_bkgsub_7',
 'aper_bkg_7',
 'aperture_sum_err_7',
 'aper_ell_sum_bkgsub_0',
 'aper_ell_bkg_0',
 'aperture_ell_sum_err_0',
 'aper_ell_sum_bkgsub_1',
 'aper_ell_bkg_1',
 'aperture_ell_sum_err_1',
 'aper_ell_sum_bkgsub_2',
 'aper_ell_bkg_2',
 'aperture_ell_sum_err_2',
 'aper_ell_sum_bkgsub_3',
 'aper_ell_bkg_3',
 'aperture_ell_sum_err_3',
 'aper_ell_sum_bkgsub_4',
 'aper_ell_bkg_4',
 'aperture_ell_sum_err_4',
 'aper_ell_sum_bkgsub_5',
 'aper_ell_bkg_5',
 'aperture_ell_sum_err_5',
 'aper_ell_sum_bkgsub_6',
 'aper_ell_bkg_6',
 'aperture_ell_sum_err_6',
 'aper_ell_sum_bkgsub_7',
 'aper_ell_bkg_7',
 'aperture_ell_sum_err_7',
 'aper_sum_bkgsub_fib',
 'aper_bkg_fib',
 'aperture_sum_err_fib',
 'sky_annulus_area_pix',
 'sky_annulus_median',
 'ra',
 'dec',
 'MJD_OBS',
 'det_id',
 'camera',
 'ci_number')
"""

#filenames=sorted(glob.glob("ci-*/ci-*_catalog.fits"))
filenames=args.infile
print(filenames)
#filenames=[filenames[0],filenames[-1]]

cameras=None

# tangent plane coordinates for first exposure
ctx0=dict()
cty0=dict()

# transformation coeff to tangent plane
T=dict()

keys="EXPID NIGHT TARGTRA TARGTDEC SKYRA SKYDEC MOUNTHA MOUNTAZ MOUNTEL MOUNTDEC PARALLAC".split()

results=[]
for filename in filenames :
    print(filename)
    cat=fitsio.read(filename)
    reduced_filename=filename.replace("catalog","reduced")
    head=fitsio.read_header(reduced_filename)
    night=head["NIGHT"]
    expid=head["EXPID"]
    # need the raw file to get full header
    raw_filename="{}/{}/{:08d}/ci-{:08d}.fits.fz".format(os.environ["DESI_SPECTRO_DATA"],night,expid,expid)
    print(raw_filename)
    head=fitsio.read_header(raw_filename,1)
    if cameras is None :
        cameras=np.unique(cat['camera'])
        print("cameras:",cameras)
        x=cat['xcentroid']
        y=cat['ycentroid']
        xe=cat['dxcentroid']
        ye=cat['dycentroid']
        
        # use only first image ra,dec to define a fixed tangent plane transfo per camera 
        # (x,y) -> tx,ty = a*x+b*y+c	
        ra=cat['ra']
        ## HA = LST-RA !!!!
        ha = np.mean(ra)-ra # sign flip is crutial here, offset has no effect
        
        dec=cat['dec']
        cha=[]
        cdec=[]
        for cam in cameras :
            ii=cat['camera']==cam
            cha.append(np.mean(ha[ii]))
            cdec.append(np.mean(dec[ii]))
        cha=np.mean(cha)
        cdec=np.mean(cdec)
        tx,ty = hadec2xy(ha,dec,cha,cdec)
        
        for cam in cameras :
            ii=(cat['camera']==cam)&(cat['dxcentroid']>0)&(cat['dycentroid']>0)&(cat['sig_major_pix']>2)
            x=cat['xcentroid'][ii]
            y=cat['ycentroid'][ii]
            xe=cat['dxcentroid'][ii]
            ye=cat['dycentroid'][ii]
            ttx=tx[ii]
            tty=ty[ii]

            print(cam,"tx=",np.mean(ttx),"ty=",np.mean(tty))
            
            w = 1./(xe**2+0.02**2)
            A=np.zeros((3,3))
            B=np.zeros((3))
            A[0,0] = np.sum(w)
            A[1,0] = A[0,1] = np.sum(w*x)
            A[1,1] = np.sum(w*x**2)
            A[1,2] = A[2,1] = np.sum(w*x*y)
            A[2,0] = A[0,2] = np.sum(w*y)
            A[2,2] = np.sum(w*y**2)
            Ai     = np.linalg.inv(A)
            B[0]   = np.sum(w*ttx)
            B[1]   = np.sum(w*x*ttx)
            B[2]   = np.sum(w*y*ttx)
            T[cam]=dict()
            T[cam]["x"] = Ai.dot(B)
            B[0]   = np.sum(w*tty)
            B[1]   = np.sum(w*x*tty)
            B[2]   = np.sum(w*y*tty)
            T[cam]["y"] = Ai.dot(B)

            print(cam,"tx={}+{}*x+{}*y".format(T[cam]["x"][0],T[cam]["x"][1],T[cam]["x"][2]))
            print(cam,"ty={}+{}*x+{}*y".format(T[cam]["y"][0],T[cam]["y"][1],T[cam]["y"][2]))
            


            if 0 :
                p = T[cam]["x"]
                plt.figure()
                plt.plot(ttx,p[0]+p[1]*x+p[2]*y,"o")
                plt.plot(ttx,ttx,"-")
                p = T[cam]["y"]
                plt.figure()
                plt.plot(tty,p[0]+p[1]*x+p[2]*y,"o")
                plt.plot(tty,tty,"-")
                plt.show()
            
            p = T[cam]["x"]
            ctx0[cam] =  p[0]+p[1]*x+p[2]*y
            p = T[cam]["y"]
            cty0[cam] =  p[0]+p[1]*x+p[2]*y
        #sys.exit(12) # debug
        continue

    # compute tangent plane coords for this exposure
        
    angles_of_cams = []
    errors_of_cams = []

    for cam in cameras :
        ii=(cat['camera']==cam)&(cat['dxcentroid']>0)&(cat['dycentroid']>0)&(cat['dxcentroid']<0.5)&(cat['dycentroid']<0.5)&(cat['sig_major_pix']>2)
        
        if np.sum(ii)==0 : continue
        
        x=cat['xcentroid'][ii]
        y=cat['ycentroid'][ii]
        xe=cat['dxcentroid'][ii]
        ye=cat['dycentroid'][ii]

        # select stars 
        #ms = np.median(cat['sig_major_pix'][ii])
        #rms = 1.48*np.median(np.abs(cat['sig_major_pix'][ii]-ms))
        #ok = 
        #plt.figure()
        #plt.plot(cat['sig_major_pix'][ii],"o")
        
        p = T[cam]["x"]
        ctx =  p[0]+p[1]*x+p[2]*y
        p = T[cam]["y"]
        cty =  p[0]+p[1]*x+p[2]*y

        n0=ctx0[cam].size
        nc=ctx.size

        dist2 = (np.tile(ctx,(n0,1))-np.tile(ctx0[cam],(nc,1)).T)**2 \
                + (np.tile(cty,(n0,1))-np.tile(cty0[cam],(nc,1)).T)**2
        #print(ctx.shape)
        #print(ctx0[cam].shape)
        
        j = np.argmin(dist2,axis=0)
        
        # matched coordinates
        mctx0=ctx0[cam][j]
        mcty0=cty0[cam][j]
        dist = np.sqrt((ctx-mctx0)**2+(cty-mcty0)**2)
        ok = np.where(dist<5./3600.)[0]
        r0=np.sqrt(mctx0[ok]**2+mcty0[ok]**2)
        rc=np.sqrt(ctx[ok]**2+cty[ok]**2)
        angles = np.arcsin((cty[ok]*mctx0[ok] - ctx[ok]*mcty0[ok])/r0/rc)/D2R # deg
        mangle = np.median(angles)
        rms    = 1.48*np.median(np.abs(angles-mangle))
        for loop in range(5) :
            ii     = np.abs(angles-mangle)<2.5*rms
            if loop<1 :
                mangle = np.median(angles[ii])
                rms    = 1.48*np.median(np.abs(angles[ii]-mangle))
            else :
                mangle = np.mean(angles[ii])
                rms    = np.std(angles[ii])
            #print(loop,mangle,rms,np.sum(ii))
        
        err    = rms/np.sqrt(np.sum(ii)-1)
        print("{} angle = {} +- {} deg (rms={})".format(cam,mangle,err,rms))
        angles_of_cams.append(mangle)
        errors_of_cams.append(err)
        if False :
            plt.figure()
            plt.plot(mctx0[ok],mcty0[ok],"o")
            plt.plot(ctx[ok],cty[ok],"X")
            print("mean x y ",np.mean(ctx[ok]),np.mean(cty[ok]))
            print("mean angle=",mangle)
            print("mean dx (recent - first) = ",np.mean(ctx[ok]-mctx0[ok]))
            print("mean dy (recent - first) = ",np.mean(cty[ok]-mcty0[ok]))
            plt.show()
            

        if False :
            plt.figure()
            plt.plot(r0,angles,"o")
            plt.plot(r0[ii],angles[ii],"o")
            plt.plot(r0,np.ones(r0.size)*mangle)
        
    #print(head)
    res=[]
    for k in keys :
        res.append(head[k])
    for c in range(4) :
        res.append(angles_of_cams[c])
        res.append(errors_of_cams[c])
    res.append(np.mean(angles_of_cams))
    res.append(np.mean(errors_of_cams)/2.)
    results.append(res)


results=np.array(results)
np.savetxt(args.outfile,results,header="EXPID NIGHT TARGTRA TARGTDEC SKYRA SKYDEC MOUNTHA MOUNTAZ MOUNTEL MOUNTDEC PARALLAC THETA1 ERR1  THETA2 ERR2 THETA3 ERR3 THETA4 ERR4 THETA_DEG_EOFN_RELATIVE THETA_ERROR")


    #plt.show()
    #sys.exit(12)

#ii=cat['camera']==cameras[0]
#xx.append(cat['xcentroid'])
#yy.append(cat['ycentroid'])


plt.show()

