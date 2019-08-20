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
parser.add_argument('-r','--reffile', type = str, default = None, required = False, help = 'reference frame from this ci-xxxxx_catalog.fits, default is the first of the input series')
parser.add_argument('-o','--outfile',type = str, default = None, required = True, help="output ascii file")
parser.add_argument('--maxdist',type=float,default=1,required=False,help="max allowed separation between stars in match")
parser.add_argument('--maxerr',type=float,default=0.3,required=False,help="max allowed uncertainty on position")
parser.add_argument('--anchor-to-north',action = 'store_true')

args = parser.parse_args()

if args.reffile is None :
    args.reffile = args.infile[0]

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


# tangent plane coordinates for first exposure
ctx0=dict()
cty0=dict()
# pixel coordinates of first exposure (to discard stars on edges)
cx0=dict()
cy0=dict()
flux0=dict()

# transformation coeff to tangent plane
T=dict()

keys="EXPID NIGHT TARGTRA TARGTDEC SKYRA SKYDEC MOUNTHA MOUNTAZ MOUNTEL MOUNTDEC PARALLAC".split()


if True :
    filename = args.reffile
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
    reference_lst = np.mean(ra) # this is just a reference to get consistent ha, the actual value does not matter
    
    ha = reference_lst-ra # sign flip is crutial here, offset has no effect
    
    dec=cat['dec']
    cha=[]
    cdec=[]
    for cam in cameras :
        ii=cat['camera']==cam
        cha.append(np.mean(ha[ii]))
        cdec.append(np.mean(dec[ii]))
    reference_ha=np.mean(cha)
    reference_dec=np.mean(cdec)
    tx,ty = hadec2xy(ha,dec,reference_ha,reference_dec)
    
    for cam in cameras :
        ii=(cat['camera']==cam)&(cat['dxcentroid']>0)&(cat['dycentroid']>0)&(cat['dxcentroid']<0.5)&(cat['dycentroid']<0.5)&(cat['sig_major_pix']>2)
        x=cat['xcentroid'][ii]
        y=cat['ycentroid'][ii]
        xe=cat['dxcentroid'][ii]
        ye=cat['dycentroid'][ii]
        ttx=tx[ii]
        tty=ty[ii]
        
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

        cx0[cam] = x
        cy0[cam] = y
        flux0[cam] =cat['aper_sum_bkgsub_3'][ii]
        




results=[]
for filename in args.infile :

    if filename == args.reffile :
        continue
    
    cat=fitsio.read(filename)
    reduced_filename=filename.replace("catalog","reduced")
    head=fitsio.read_header(reduced_filename)
    night=head["NIGHT"]
    expid=head["EXPID"]
    # need the raw file to get full header
    raw_filename="{}/{}/{:08d}/ci-{:08d}.fits.fz".format(os.environ["DESI_SPECTRO_DATA"],night,expid,expid)
    print(raw_filename)
    head=fitsio.read_header(raw_filename,1)
    

    # compute tangent plane coords for this exposure
        
    angles_of_cams = []
    errors_of_cams = []
    x0_of_cams = []
    y0_of_cams = []
    xc_of_cams = []
    yc_of_cams = []
    dx_of_cams = []
    dy_of_cams = []
    cam_of_cams = []
    
    for cam in cameras :
        
        ii=(cat['camera']==cam)&(cat['dxcentroid']>0)&(cat['dycentroid']>0)&(cat['dxcentroid']<args.maxerr)&(cat['dycentroid']<args.maxerr)&(cat['sig_major_pix']>2)&(cat['sig_minor_pix']>1)
        
        if np.sum(ii)==0 : 
            print(cam,'no stars!')
            continue
        
        x=cat['xcentroid'][ii]
        y=cat['ycentroid'][ii]
        xe=cat['dxcentroid'][ii]
        ye=cat['dycentroid'][ii]

        p = T[cam]["x"]
        ctx  =  p[0]+p[1]*x+p[2]*y
        ctxe =  np.sqrt((p[1]*xe)**2+(p[2]*ye)**2)
        p = T[cam]["y"]
        cty  =  p[0]+p[1]*x+p[2]*y
        ctye =  np.sqrt((p[1]*xe)**2+(p[2]*ye)**2)
        n0=ctx0[cam].size
        nc=ctx.size

        

        if 1 : # use astrometry of this exposure for match
            tx,ty = hadec2xy(reference_lst - cat['ra'][ii],cat['dec'][ii],reference_ha,reference_dec)

            k=np.argsort(flux0[cam])[::-1]
            n0=nc+3
            ktx0=ctx0[cam][k][:n0]
            kty0=cty0[cam][k][:n0]
            n0=ktx0.size

            dx=0.
            dy=0.
            for loop in range(3) :
                dist2 = (np.tile(tx,(n0,1))-np.tile(ktx0-dx,(nc,1)).T)**2 \
                        + (np.tile(ty,(n0,1))-np.tile(kty0-dy,(nc,1)).T)**2
                j = np.argmin(dist2,axis=0)
                dx = np.median(ktx0[j]-tx)
                dy = np.median(kty0[j]-ty)
                
                mctx0=ktx0[j]
                mcty0=kty0[j]
                        
            dist = np.sqrt((tx-(mctx0-dx))**2+(ty-(mcty0-dy))**2)
            mx0=cx0[cam][j]
            my0=cy0[cam][j]
            # selection on distance of match and avoid edges of CCD
            margin=1
            ok = np.where((dist<args.maxdist/3600.)&(mx0>margin)&(mx0<(3072-margin))&(my0>margin)&(my0<(2048-margin)))[0]
            print(cam,"number of matched stars=",ok.size)
            
            if False :
                plt.plot(ctx0[cam],cty0[cam],"o",c="b",alpha=0.2)
                plt.plot(ctx0[cam][k][:n0],cty0[cam][k][:n0],"o",c="b")
                plt.plot(tx,ty,"o",c="orange")
                for i in ok :
                    plt.plot([mctx0[i],tx[i]],[mcty0[i],ty[i]],"-",c="black")
                    plt.plot([mctx0[i],ctx[i]],[mcty0[i],cty[i]],"--",c="red")
                    
                plt.show()

        else : # use tranfo from reference exposure for the match
            dist2 = (np.tile(ctx,(n0,1))-np.tile(ctx0[cam],(nc,1)).T)**2 \
                + (np.tile(cty,(n0,1))-np.tile(cty0[cam],(nc,1)).T)**2
            j = np.argmin(dist2,axis=0)
            # matched coordinates
            mctx0=ctx0[cam][j]
            mcty0=cty0[cam][j]
            dist = np.sqrt((ctx-mctx0)**2+(cty-mcty0)**2)
            mx0=x0[cam][j]
            my0=y0[cam][j]
            # selection on distance of match and avoid edges of CCD
            ok = np.where((dist<args.maxdist/3600.)&(mx0>10)&(mx0<(3072-10))&(my0>10)&(my0<(2048-10)))[0]
        
        if ok.size<3 :
            print(cam,"not enough matched star:",ok.size)
            continue
        x0_of_cams.append(mctx0[ok])
        y0_of_cams.append(mcty0[ok])
        xc_of_cams.append(ctx[ok])
        yc_of_cams.append(cty[ok])
        dx_of_cams.append(np.median(ctx[ok]-mctx0[ok]))
        dy_of_cams.append(np.median(cty[ok]-mcty0[ok]))
        cam_of_cams.append(cam)

    if len(dx_of_cams)<2 :
        print("skip this exposure because only {} cams with enough stars".format(len(dx_of_cams)))
        continue

    # will translate before measuring angles
    dx=np.mean(dx_of_cams)
    dy=np.mean(dy_of_cams)
    print("mean dx,dy of cams = ",dx,dy)
    
    for c,cam in enumerate(cam_of_cams) :
        x0 = x0_of_cams[c]
        y0 = y0_of_cams[c]
        xc = xc_of_cams[c]-dx
        yc = yc_of_cams[c]-dy
        
        r0=np.sqrt(x0**2+y0**2)
        rc=np.sqrt(xc**2+yc**2)
        angles = np.arcsin((yc*x0 - xc*y0)/r0/rc)/D2R # deg
        mangle = np.median(angles)
        rms    = 1.48*np.median(np.abs(angles-mangle))
        err    = 1.25*rms/np.sqrt(len(x0_of_cams)-1)
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
            plt.errorbar(r0,angles,errors,fmt="o")
            plt.figure()
            plt.plot(r0,angles,"o")
            #plt.errorbar(r0,angles,errors,fmt="o")
            plt.plot(r0[kk],angles[kk],"o")
            plt.plot(r0,np.ones(r0.size)*mangle)
            plt.figure()
            plt.plot(x[ok],y[ok],"o")
            plt.plot(x[ok][kk],y[ok][kk],"o")
            plt.show()


    angles_of_cams = np.array(angles_of_cams)
    errors_of_cams = np.array(errors_of_cams)    
 
    if np.std(angles_of_cams) == 0 or np.mean(angles_of_cams) == 0 :
        print("nothing interesting for this exposure")
        continue
    
    cam_of_cams = np.array(cam_of_cams)
    #print(cam_of_cams)
    cc=np.where(cam_of_cams!=b'CIC')[0]
    #print(cc)
    #print("n outer cams=",cc.size)
    print("outer cams=",cam_of_cams[cc])
    
    if cc.size<2 :
        print("skip this exposure")
        continue
        
    res=[]
    for k in keys :
        res.append(head[k])
    
    # I don't want to do a weighted mean because I want to weight
    # equally the cameras.
    mtheta=np.mean(angles_of_cams[cc])
    
    # approximate error
    chi2=np.sum((angles_of_cams[cc]-mtheta)**2/errors_of_cams[cc]**2)
    rchi2=chi2/3.
    mthetaerr=np.sqrt(np.sum(errors_of_cams[cc]**2))/cc.size
    if rchi2>1 :
        mthetaerr *= np.sqrt(rchi2)
    res.append(mtheta)
    res.append(mthetaerr)
    print("angle= {:6.5f} +- {:6.5f}".format(mtheta,mthetaerr))
    results.append(res)


results=np.array(results)
np.savetxt(args.outfile,results,header="EXPID NIGHT TARGTRA TARGTDEC SKYRA SKYDEC MOUNTHA MOUNTAZ MOUNTEL MOUNTDEC PARALLAC THETA_DEG_EOFN_RELATIVE THETA_ERROR")

plt.show()

