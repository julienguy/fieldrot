#!/usr/bin/env python

import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import argparse

from fieldrot.rotlib import *

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
description="Fit field rotation data")
parser.add_argument('-i','--infile', type = str, default = None, required = True, 
                    help = 'path to ASCII file with sky rotation data (like sky-ci-rotation.txt)')
parser.add_argument('--opposite', action = 'store_true',
                    help = 'interpret rotation angle as rot of instrument wrt sky instead of sky wrt instrument, needed for Aarons angles')
parser.add_argument('--scan', action = 'store_true',
                    help = 'scan polar misalignment')
parser.add_argument('--fit-torque', action = 'store_true',
                    help = 'fit torque on top of polar misalignement')
parser.add_argument('--fix-b', type = float, help = 'fix b')
parser.add_argument('--fix-a', type = float, help = 'fix a')
parser.add_argument('--fixed-torque', action = 'store_true',
                    help = 'use fixed torque model (hardcoded parameters)')
parser.add_argument('--fit-offset', action = 'store_true',
                    help = 'fit an offset')
parser.add_argument('--with-boresight-offset', action = 'store_true',
                    help = 'include boresight offset')

args = parser.parse_args()


tmp    = np.linspace(0,2*np.pi,10)
x2d  = np.cos(tmp)
y2d  = np.sin(tmp)


filename=args.infile

print("loading",filename)

tmp=np.loadtxt(filename).T
table={}
for line in open(filename).readlines() :
    if line.find("# EXPID")>=0 :
        keys=(line.replace("#","").strip()).split()
        break
#print(keys)
i=0
for k in keys :
    if len(k)>0 :
        table[k]=tmp[i]
        i+=1
print(table.keys())


#for i in range(table["HEXPOS1"].size) :
#    print(table["NIGHT"][i],table["HEXPOS1"][i],table["HEXPOS2"][i],table["HEXPOS3"][i],table["HEXPOS4"][i],table["HEXPOS5"][i],table["HEXPOS6"][i])

configs = []
if "HEXPOS1" in table :
    configs.append( (np.abs(table["HEXPOS1"]+2000)<2.) & (np.abs(table["HEXPOS2"]+4000)<2.)  )
    #configs.append( (np.abs(table["HEXPOS1"]-2000)<2.) & (np.abs(table["HEXPOS2"]+4000)<2.)  )
    #configs.append( (np.abs(table["HEXPOS1"]-5740)<2.) & (np.abs(table["HEXPOS2"]+4000)<2.)  )
    #configs.append( (np.abs(table["HEXPOS1"]+3790)<2.) & (np.abs(table["HEXPOS2"]+4000)<2.)  )
    if 0 :
        configs.append( (np.abs(table["HEXPOS1"]-5740)<2.) & (np.abs(table["HEXPOS2"]+4000)<2.)  )
        configs.append( (np.abs(table["HEXPOS1"]+2000)<2.) & (np.abs(table["HEXPOS2"]-6700)<2.)  )
        configs.append( (np.abs(table["HEXPOS1"]+1617.0)<2.) & (np.abs(table["HEXPOS2"]+3815.5)<2.)  )
        configs.append( (np.abs(table["HEXPOS1"]-2000)<2.) & (np.abs(table["HEXPOS2"]+4000)<2.)  )
        configs.append( (np.abs(table["HEXPOS1"]+3790)<2.) & (np.abs(table["HEXPOS2"]+4000)<2.)  )
else :
    configs = [np.ones(table["EXPID"].size)]
nconfig = len(configs)


if "THETA_DEG_EOFN" in table :
    # recalibration over nights
    #plt.figure()
    #plt.plot(table["NIGHT"],table["THETA_DEG_EOFN_RELATIVE"]-table["THETA_DEG_EOFN"],"o")
    nights=np.unique(table["NIGHT"])
    mdiff=np.mean(table["THETA_DEG_EOFN_RELATIVE"]-table["THETA_DEG_EOFN"])
    for night in nights :
        ii=np.where(table["NIGHT"]==night)[0]
        night_ddiff = np.mean(table["THETA_DEG_EOFN_RELATIVE"][ii]-table["THETA_DEG_EOFN"][ii])-mdiff
        table["THETA_DEG_EOFN_RELATIVE"][ii] -= night_ddiff
    #plt.plot(table["NIGHT"],table["THETA_DEG_EOFN_RELATIVE"]-table["THETA_DEG_EOFN"],"o")


# latitude
# 31:57:50.5 DESI-3182-v4
lat = 31.96403 #deg
north_axis = unit_vector(0,0)
dec_axis   = unit_vector(np.pi/2,-np.pi/2.)  

    
for sample,ii in enumerate(configs ) :

    if np.sum(ii)<2 : continue
    ii = np.where(ii)[0]

    if "HEXPOS1" in table :
        HEXPOS1 = int(np.floor(np.mean(table["HEXPOS1"][ii])))
        HEXPOS2 = int(np.floor(np.mean(table["HEXPOS2"][ii])))
    else :
        HEXPOS1 = HEXPOS2 = 0.
    
    ca=np.cos(lat*D2R)
    sa=np.sin(lat*D2R)
    cp=np.cos(table["MOUNTHA"][ii]*D2R)
    sp=np.sin(table["MOUNTHA"][ii]*D2R)
    ct=np.cos((90.-table["MOUNTDEC"][ii])*D2R)
    st=np.sin((90.-table["MOUNTDEC"][ii])*D2R)
    x=ct*cp*ca-st*sa
    y=sp*ca
    if "ROT" in table :
        z=table["ROT"][ii]
        if args.opposite :
            print("with angle = ROT , no opposite")
            sys.exit(12)
        ze=table["ERR"][ii]
    else :
        z=table["THETA_DEG_EOFN_RELATIVE"][ii]
        if args.opposite :
            z *= -1
            print("!!!!!!!!!!!!! OPPOSITE SIGN FOR THETA_DEG_EOFN_RELATIVE !!!!!!!!!!")
        ze=table["THETA_ERROR"][ii]
            
    w=1./(ze**2+0.003**2) #np.ones(ii.size)
    with_offset = True
    if with_offset :
        A=np.zeros((3,3))
        B=np.zeros(3)
        A[0,0]=np.sum(w)
        A[0,1]=A[1,0]=np.sum(w*x)
        A[1,1]=np.sum(w*x**2)
        A[0,2]=A[2,0]=np.sum(w*y)
        A[1,2]=A[2,1]=np.sum(w*x*y)
        A[2,2]=np.sum(w*y**2)
    else :
        A=np.zeros((2,2))
        B=np.zeros(2)
        A[0,0]=np.sum(w*x**2)
        A[0,1]=A[1,0]=np.sum(w*x*y)
        A[1,1]=np.sum(w*y**2)

    if args.fix_a is not None :
        print("add prior fixing a={}".format(args.fix_a))
        A[1,1] += 1e12 
    if args.fix_b is not None :
        print("add prior fixing b={}".format(args.fix_b))
        A[2,2] += 1e12 
    Ai = np.linalg.inv(A)

    # refraction + polar misalignment, best for of aaron's angle data
    # telescope_polar_axis_dha = 343  ;  telescope_polar_axis_ddec = 0.036 # degrees
    # because of yet another sign flip
    telescope_polar_axis_dha = 10.  ;  telescope_polar_axis_ddec = 0.032 # degrees
    # best fit without b term
    #telescope_polar_axis_dha = 240  ;  telescope_polar_axis_ddec = 0.02 # degrees
    
    
    fit_offset = True
    
    if 1 :
        # use actual measurements of the pole
        # from Dick Joyce
        # "in the northern hemisphere, positive MA means that the pole of the mounting is to the right of due north"
        # "in the northern hemisphere, positive ME means that the pole of the mounting is below the true (unrefracted) north"
        MA=25/3600.
        ME=-103/3600. 
        ha,dec=altaz2hadec(alt=lat-ME,az=MA,lat=lat)
        telescope_polar_axis_dha = ha
        telescope_polar_axis_ddec = 90 - dec
        print("telescope_polar_axis_dha  = ",telescope_polar_axis_dha )
        print("telescope_polar_axis_ddec = ",telescope_polar_axis_ddec )

    # for actual polar-misalignment
    boresight_offset_dx = -0.021 # deg , real one is indeed -0.013
    boresight_offset_dy = 0. # deg
    
    # best is
    if 0 :
        telescope_polar_axis_dha = 10.  ;  telescope_polar_axis_ddec = 0.032 # degrees
        boresight_offset_dx = -0.013 # deg
        boresight_offset_dy = 0.0 # deg
    
    scale=1
    
    tha_a=[]
    tdec_a=[]
    scale_a=[]
    rms_a=[]
    chi2_a=[]
    dx_a=[]
    dy_a=[]
    ##################################
    
    if args.scan :
        thas=[telescope_polar_axis_dha]
        tdecs=[telescope_polar_axis_ddec]
        scales=[scale]
        bsdxs=[boresight_offset_dx]
        bsdys=[boresight_offset_dy]


        #bsdxs=np.linspace(-0.08,0.,81)
        #bsdys=np.linspace(-0.04,0.04,41)        
        #thas=np.linspace(0,360,37)
        #tdecs=np.linspace(0.02,0.04,21)
        scales=np.linspace(0,1.5,16)
        
    else :
        thas=[telescope_polar_axis_dha]
        tdecs=[telescope_polar_axis_ddec]
        scales=[scale]
        bsdxs=[boresight_offset_dx]
        bsdys=[boresight_offset_dy]

    for bsdx in bsdxs :
        for bsdy in bsdys :
            for rscale in scales :
                for tha  in thas :
                    for tdec in tdecs :
                        # see rotlib tests for verification
                        pol_ha  = tha
                        pol_dec = 90-tdec
                        t,p = hadec2thetaphi(pol_ha,pol_dec)
                        cross = np.cross(unit_vector(t,p),unit_vector(0,0))
                        norme = np.sqrt(np.sum(cross**2))
                        if norme > 0 :
                            drot = rotation_matrix(cross/norme,np.arcsin(norme))
                        else :
                            drot = np.eye(3)

                        ptheta1=np.zeros(ii.size)
                        ptheta=np.zeros(ii.size)

                        for j,i in enumerate(ii) :

                            # atmospheric refraction = 79. arcsec at a zenith angle of 60 deg, DESI-3182
                            el = table["MOUNTEL"][i]
                            paral = table["PARALLAC"][i]
                            hacen  = table["MOUNTHA"][i]
                            deccen = table["MOUNTDEC"][i]
                            if bsdx !=0 or bsdy !=0 :
                                # shift field center
                                hacen,deccen=xy2hadec(bsdx,bsdy,hacen,deccen)    
                            ha1,dec1   = xy2hadec(x2d,y2d,hacen,deccen)
                            ha2,dec2   = rotation_hadec(ha1,dec1,drot)
                            alt,az     = hadec2altaz(ha2,dec2,lat)

                            R = 79./3600.*np.tan(30.*D2R)/np.tan(alt*D2R)  # deg , refraction per point in field
                            R *= rscale
                            alt += R # apply true direct refraction where it applies

                            ha2,dec2   = altaz2hadec(alt,az,lat)

                            if False : # adjust pointing, guide on one star (it's not needed, does not change rotation)
                                ddec=dec2[0]-dec1[0]
                                dha=ha2[0]-ha1[0]
                                dha=ha2[0]-ha1[0]
                                if dha>180: dha  -= 360.
                                if dha<-180: dha  += 360.
                                print("step",0,dha,ddec)
                                for loop in range(1) :
                                    mha=ha2[0]
                                    adjustment_rot = rotation_matrix(north_axis,(mha-dha)*D2R).dot(rotation_matrix(dec_axis,-ddec*D2R).dot(rotation_matrix(north_axis,-mha*D2R)))
                                    ha2,dec2       = rotation_hadec(ha2,dec2,adjustment_rot)
                                    ddec=dec2[0]-dec1[0]
                                    dha=ha2[0]-ha1[0]
                                    if dha>180: dha  -= 360.
                                    if dha<-180: dha  += 360.
                                    #dha=np.mean(dha)
                                    print("step",loop+1,dha,ddec)

                            # rotation angle, counterclockwise, due to polar misalignment, + boresight misalignment
                            ptheta1[j]  = np.mean(measure_field_rotation(ha1,dec1,ha2,dec2,0,0)) ####

                            if args.with_boresight_offset :
                                ptheta[j]   = np.mean(measure_field_rotation(ha1,dec1,ha2,dec2,bsdx,bsdy)) ####
                            else :
                                ptheta[j]   = ptheta1[j]
                                
                        if args.fit_torque :
                            if with_offset :
                                B[0]=np.sum(w*(z-ptheta))
                                B[1]=np.sum(w*x*(z-ptheta))
                                B[2]=np.sum(w*y*(z-ptheta))
                            else :
                                B[0]=np.sum(w*x*(z-ptheta))
                                B[1]=np.sum(w*y*(z-ptheta))
                            if args.fix_a :
                                B[1]+=1e12*args.fix_a
                            if args.fix_b :
                                B[2]+=1e12*args.fix_b
                                
                            params = Ai.dot(B)
                            
                            if with_offset : 
                                c = params[0]
                                a = params[1]
                                b = params[2]
                                model=a*x+b*y+c
                            else :
                                a = params[0]
                                b = params[1]
                                model=a*x+b*y
                        else :
                            a=b=c=0
                            model = 0*z

                        if args.fixed_torque :
                            # after changing sign of aaron's !
                            a= 0.01677551402438771
                            b= -0.016771011680446243
                            c= -0.02275433847462704
                            #a= -0.018394127123866604
                            #b= 0.07552874818502238
                            #c= 0.0206462726669795
                            model=a*x+b*y+c
                            #if args.fit_offset :
                            #    model += np.mean(z-ptheta-model)

                        rms0 = np.std(z-ptheta)
                        if fit_offset :
                            delta = np.sum(w*(z-ptheta-model))/np.sum(w)
                        rms  = np.std(z-ptheta-model)
                        chi2  = np.sum(w*(z-ptheta-model-delta)**2)/w.size
                        print(tha,tdec,"rms(ptheta):",np.std(ptheta),"rms(res):",rms0,"->",rms)
                        tha_a.append(tha)
                        tdec_a.append(tdec)
                        rms_a.append(rms)
                        chi2_a.append(chi2)
                        scale_a.append(rscale)
                        dx_a.append(bsdx)
                        dy_a.append(bsdy)

                
    i=np.argmin(chi2_a)
    print("HEXPOS=({},{},) a={:4.3f} b={:4.3f} rms={:4.3f}".format(HEXPOS1,HEXPOS2,a,b,rms_a[i]))

    if sample > 0 :
        print("skip plot")
        continue
    
    print("minimum tha = {:6.5f} deg; tdec = {:6.5f} deg = {:3.1f} arcsec; dx = {:6.5f} deg; dy = {:6.5f} deg; rms = {:6.5f} deg = {:3.1f} arcsec ; chi2={:3.2f}".format(tha_a[i],tdec_a[i],tdec_a[i]*3600.,dx_a[i],dy_a[i],rms_a[i],rms_a[i]*3600.,chi2_a[i]))
    alt,az=hadec2altaz(tha_a[i],90.-tdec_a[i],lat)
    if az>180. : az -=360.
    print("corresponding to alt = lat + {:5.4f} arcsec; az = {:5.4f} arcsec".format((alt-lat)*3600.,az*3600.))
    ma = az*3600.
    me = (lat-alt)*3600.
    
    print("DESI-3182 31:59:30.5 00:00:27 => alt = lat + 100 arcsec; az = 27 arcsec")
    alt2=lat+100./3600.
    az2=27./3600.

    v1 = unit_vector((90.-alt)*D2R,az*D2R)
    v2 = unit_vector((90.-alt2)*D2R,az2*D2R)
    c12 = np.cross(v1,v2)
    angle = np.arcsin(np.sqrt(np.sum(c12**2)))/D2R
    print("difference = {} deg = {} arcsec".format(angle,angle*3600.))
    
    print("scale of refraction = ",scale_a[i])
    
    print("a=",a)
    print("b=",b)
    if with_offset: print("c=",c)
    print("rms=",rms," deg = ",rms*3600.," arcsec")

    # remove mean of boresight offset
    ptheta -= np.mean(ptheta-ptheta1)
    
    
    if 0 :
        if args.fit_torque or args.fixed_torque :
            plt.figure()
            x=model+ptheta
            plt.errorbar(x,z,ze,fmt="o")
            plt.xlabel("rotation from polar misalignment + refraction + FP torque (deg)")
            plt.ylabel("measured rotation angle on sky + offset (deg)")
            i=np.argsort(x); plt.plot(x,x,"-",c="gray")
            plt.grid()

            plt.figure()
            x=model
            plt.errorbar(x,z-ptheta,ze,fmt="o")
            plt.xlabel("rotation from FP torque (deg)")
            plt.ylabel("measured rotation angle on sky + offset - polar misalignment + refraction (deg)")
            i=np.argsort(x); plt.plot(x,x,"-",c="gray")
            plt.grid()

            plt.figure()
            x=ptheta
            plt.errorbar(x,z-model,ze,fmt="o")
            plt.xlabel("rotation from polar misalignment + refraction (deg)")
            plt.ylabel("measured rotation angle on sky - FP torque (deg)")
            i=np.argsort(x); plt.plot(x,x,"-",c="gray")
            plt.grid()
        else :
            plt.figure()
            x=ptheta
            
            plt.errorbar(x,z,ze,fmt="o",label="20190419 expid 8240-8274")
            if args.with_boresight_offset :
                plt.xlabel("rotation from polar misalignment + refraction + boresight offset (deg)")
            else :
                plt.xlabel("rotation from polar misalignment + refraction (deg)")
            
            plt.ylabel("measured rotation angle on sky + offset (deg)")
            
            offset = np.mean(z-x)
            i=np.argsort(x); plt.plot(x,x+offset,"-",c="gray",label="slope=1, arb. offset")
    
            plt.grid()
            plt.legend(title="Polar axis MA={:3.1f} ME={:3.1f} arcsec".format(ma,me))

    if 0 :
        plt.figure()
        plt.scatter(table["MOUNTEL"][ii],z-model-ptheta,c=table["MOUNTDEC"][ii])
        plt.colorbar()
        plt.xlabel("elevation (deg)")
        plt.ylabel("measured rotation angle on sky - model")
        plt.grid()

    mstyle="-"
    conversion = 1.
    conversion = 3600.
    if conversion == 1: unit="deg"
    elif conversion == 3600: unit="arcsec"
    
    night=int(np.median(table["NIGHT"][ii]))
    emin=int(np.min(table["EXPID"][ii]))
    emax=int(np.max(table["EXPID"][ii]))    
    label  = "{} expid {}-{}".format(night,emin,emax)
    if args.fit_offset :
        offset = np.mean(z-ptheta)
        print("offset=",offset)
        label  += " + offset"
    else :
        offset = 0.
    if 1 : # vs HA
        title=filename.split(".")[0]+"-vs-ha"
        fig=plt.figure(title)

        
        
        plt.scatter(table["MOUNTHA"][ii],(z-offset)*conversion,label=label,alpha=0.8)
        plt.plot(table["MOUNTHA"][ii],ptheta1*conversion,mstyle,c="red",label="polar misalignment + refraction")
        if args.with_boresight_offset :
            plt.plot(table["MOUNTHA"][ii],ptheta*conversion,"v",c="green",label="polar misalignment + refraction + boresight offset")
        if args.fit_torque or args.fixed_torque :
            #plt.plot(table["MOUNTHA"][ii],model,"-",label="FP torque")
            label="polar misalignment + refraction + FP torque"
            label+=" (a={:4.3f} b={:4.3f} deg)".format(a,b)
            plt.plot(table["MOUNTHA"][ii],(ptheta+model-offset)*conversion,mstyle,c='orange',label=label)
        legend_title='with ME={:d}" MA={:d}"'.format(int(me),int(ma))
        #if args.fit_torque or args.fixed_torque : legend_title+=" a={:4.3f} b={:4.3f} deg".format(a,b)
        plt.legend(title=legend_title)
        plt.xlabel("HA (deg)")
        if offset :
            plt.ylabel("measured rotation angle on sky + offset ({})".format(unit))
        else :
            plt.ylabel("measured rotation angle on sky ({})".format(unit))
        plt.grid()
        fig.savefig(title+".png")
        print("wrote",title+".png")


    if 1 : # vs EXPID
        title=filename.split(".")[0]+"-vs-expid"
        fig=plt.figure(title)

        
        
        plt.scatter(table["EXPID"][ii],(z-offset)*conversion,label=label,alpha=0.8)
        plt.plot(table["EXPID"][ii],ptheta1*conversion,mstyle,c="red",label="polar misalignment + refraction")
        if args.with_boresight_offset :
            plt.plot(table["EXPID"][ii],ptheta*conversion,"v",c="green",label="polar misalignment + refraction + boresight offset")
        if args.fit_torque or args.fixed_torque :
            #plt.plot(table["EXPID"][ii],model,"-",label="FP torque")
            label="polar misalignment + refraction + FP torque"
            label+=" (a={:4.3f} b={:4.3f} deg)".format(a,b)
            plt.plot(table["EXPID"][ii],(ptheta+model-offset)*conversion,mstyle,c='orange',label=label)
        legend_title='with ME={:d}" MA={:d}"'.format(int(me),int(ma))
        #if args.fit_torque or args.fixed_torque : legend_title+=" a={:4.3f} b={:4.3f} deg".format(a,b)
        plt.legend(title=legend_title)
        plt.xlabel("EXPID")
        if offset :
            plt.ylabel("measured rotation angle on sky + offset ({})".format(unit))
        else :
            plt.ylabel("measured rotation angle on sky ({})".format(unit))
        plt.grid()
        fig.savefig(title+".png")
        print("wrote",title+".png")

    if 0 : # vs dec
        title=filename.split(".")[0]+"-vs-dec"
        fig=plt.figure(title)
        plt.scatter(table["MOUNTDEC"][ii],z-offset,label=label)
        plt.plot(table["MOUNTDEC"][ii],ptheta,mstyle,c="red",label="polar misalignment + refraction")
        if args.fit_torque or args.fixed_torque :
            plt.plot(table["MOUNTDEC"][ii],ptheta+model,mstyle,c='orange',label="polar misalignment + refraction + FP torque")
        plt.legend()
        plt.xlabel("DEC (deg)")
        if offset :
            plt.ylabel("measured rotation angle on sky + offset ({})".format(unit))
        else :
            plt.ylabel("measured rotation angle on sky ({})".format(unit))
        plt.grid()
        fig.savefig(title+".png")
        print("wrote",title+".png")
        
    if 0 : # vs alt
        title=filename.split(".")[0]+"-vs-el"
        fig=plt.figure(title)
        plt.scatter(table["MOUNTEL"][ii],z-offset,label=label)
        plt.plot(table["MOUNTEL"][ii],ptheta,mstyle,c="red",label="polar misalignment + refraction")
        if args.fit_torque or args.fixed_torque :
            plt.plot(table["MOUNTEL"][ii],ptheta+model,mstyle,c='orange',label="polar misalignment + refraction + FP torque")
        plt.legend()
        plt.xlabel("ELEVATION (deg)")
        if offset :
            plt.ylabel("measured rotation angle on sky + offset ({})".format(unit))
        else :
            plt.ylabel("measured rotation angle on sky ({})".format(unit))
        plt.grid()
        fig.savefig(title+".png")
        print("wrote",title+".png")
        
plt.show()
