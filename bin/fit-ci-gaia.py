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
parser.add_argument('-i','--infile', type = str, nargs = "*" , default = None, required = True, 
                    help = 'path to ASCII files with sky rotation data (like sky-ci-rotation.txt)')
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
parser.add_argument('--fit-offset', action = 'store_true', default =True, help = 'fit an offset')
parser.add_argument('--with-boresight-offset', action = 'store_true', default =True, 
                    help = 'include boresight offset')

args = parser.parse_args()


#tmp    = np.linspace(0,2*np.pi,10)
tmp    = [0,np.pi/2.,np.pi,3*np.pi/2.]
x2d  = np.cos(tmp)
y2d  = np.sin(tmp)


table=None

for j,filename in enumerate(args.infile) :

    print("loading",filename)

    tmp=np.loadtxt(filename).T
    for line in open(filename).readlines() :
        if line.find("#")>=0 :
            keys=(line.replace("#","").strip()).split()
            break
    if table is None :
        table={}
        i=0
        for k in keys :
            if len(k)>0 :
                table[k]=tmp[i]
                i+=1
        table["fileid"]=j*np.ones(tmp.shape[1])
    else :
        i=0
        for k in keys :
            if len(k)>0 :
                table[k]=np.append(table[k],tmp[i])
                i+=1

        table["fileid"]=np.append(table["fileid"],j*np.ones(tmp.shape[1]))


if "TARGTHA" in table :
    RA="TARGTRA"
    HA="TARGTHA"
    DEC="TARGTDEC"
else :
    RA="TARGTRA"
    HA="MOUNTHA"
    DEC="TARGTDEC"

#plt.plot(table["MOUNTHA"],table["TCSST"]-table["TARGTRA"]-table["MOUNTHA"],"o")
#tmp=(table["TCSST"]-table["TARGTRA"]-table["TARGTHA"])%360.
#tmp[tmp>180] -= 360
#plt.plot(table["MOUNTHA"],tmp,"o")
#plt.show()

#plt.plot(table["NIGHT"],table["HEXPOS2"],"o")
#plt.show()


#table["SIGSEEING"] *= (0.12*2.35) # convert to arcsec FWHM
#plt.plot(table["EXPID"],table["SIGSEEING"],"o")
#plt.show()

adcphi=(table["ADC1PHI"]-table["ADC2PHI"])%360
adcphi[adcphi>180] -= 360

#plt.plot(table["NIGHT"],adcphi,"o")
#plt.hist(adcphi,bins=1000)
#plt.show()

print(table.keys())
fileids = np.unique(table["fileid"])


#for i in range(table["HEXPOS1"].size) :
#    print(table["NIGHT"][i],table["HEXPOS1"][i],table["HEXPOS2"][i],table["HEXPOS3"][i],table["HEXPOS4"][i],table["HEXPOS5"][i],table["HEXPOS6"][i])

configs = []
if "HEXPOS1" in table :
    
    configs.append( (table["NIGHT"]<=20190505) & (np.abs(table["HEXPOS6"])<10.)  )
    #configs.append( (table["NIGHT"]>20190505) & (np.abs(table["HEXPOS6"])<10.)  )

    #configs.append( (np.abs(table["HEXPOS2"])<1000) & (np.abs(table["HEXPOS6"])<10.) & (table["NIGHT"]>20190505)  )
    #configs.append( (np.abs(table["HEXPOS2"]+4000)<500) )
    #configs.append( (np.abs(table["HEXPOS2"])<1000)  & (table["NIGHT"]>20190505) )
    #configs.append( (np.abs(table["HEXPOS6"])<10.) )
    
else :
    configs = [np.ones(table["EXPID"].size)]


for config in configs : # more selection
    t=table
    
    err=t["THETARMS"]/np.sqrt(t["NSELECTED"])
    config &= (err<0.05)&(t["NCAMS"]==4)&(t["ADCCORR"]==0)  # 4 outer cameras have matched stars
    config &= (t["NSTARS"]>10) # at least nstars (total of 4 outer cameras)
    config &= (t["CONTRAST"]>2) # quality of match to gain
    #config &= (t["SIGSEEING"]>0.5)&(t["SIGSEEING"]<3.) # quality of match to gain
    config &= (t["AIRMASS"]<1.8) # I need to do a better job with airmass
    config &= (np.abs(adcphi)<2.)
    decmax=75
    config &=(t[DEC]<decmax)
    config &=(t["EXPID"]!=10329)&(t["EXPID"]!=10311) # bad PSF
    
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
    cp=np.cos(table[HA][ii]*D2R)
    sp=np.sin(table[HA][ii]*D2R)
    ct=np.cos((90.-table[DEC][ii])*D2R)
    st=np.sin((90.-table[DEC][ii])*D2R)
    x=ct*cp*ca-st*sa
    y=sp*ca
    if "ROT" in table :
        z=table["ROT"][ii]
        if args.opposite :
            print("with angle = ROT , no opposite")
            sys.exit(12)
        ze=table["ERR"][ii]
    elif "THETA" in table :
        z=table["THETA"][ii]
        ze=table["THETARMS"][ii]/np.sqrt(table["NSELECTED"][ii])
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

    
    fit_offset = True
    boresight_offset_dx = 0.
    boresight_offset_dy = 0.
    scale = 1.
    
    # "in the northern hemisphere, positive MA means that the pole of the mounting is to the right of due north"
    # "in the northern hemisphere, positive ME means that the pole of the mounting is below the true (unrefracted) north"
    # not clear any of this is accurate
    MA=25/3600. # polar axis azimuth offset (positive being to the right of true north in the Northern Hemisphere)
    ME=-103/3600. # polar axis elevation offset (positive being below the true pole)
    
    if 0 :
        # use MA ME this  
        ha,dec=altaz2hadec(alt=lat-ME,az=MA,lat=lat)
        telescope_polar_axis_dha = ha
        telescope_polar_axis_ddec = 90 - dec
        print("telescope_polar_axis_dha  = ",telescope_polar_axis_dha )
        print("telescope_polar_axis_ddec = ",telescope_polar_axis_ddec )

    
    
    # best is for external rotation
    if 1 :
        telescope_polar_axis_dha = 1.  ;  telescope_polar_axis_ddec = 0.0335 ; boresight_offset_dx = -0.017 # deg, best fit for dec<75 (has no outliers) HEXPOS2=-4000
        #telescope_polar_axis_dha = 1.  ;  telescope_polar_axis_ddec = 0.03 ; boresight_offset_dx = -0.0185 # deg# degrees , best fit for dec<75 (has no outliers) HEXPOS2=0, night>20190505
            
    
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


        bsdxs=np.linspace(-0.02,0.01,31)
        #bsdys=np.linspace(-0.04,0.04,41)        
        #thas=np.linspace(-10,10,21)
        #tdecs=np.linspace(0.03,0.04,21)
        #scales=np.linspace(0,1.5,16)
        
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
                            hacen  = table[HA][i]
                            deccen = table[DEC][i]
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
                            delta = np.zeros(z.size)
                            for fid in fileids :
                                jj=np.where(table["fileid"][ii]==fid)[0]
                                d = np.sum(w[jj]*(z[jj]-ptheta[jj]-model[jj]))/np.sum(w[jj])
                                delta[jj]=d
                                #print(jj.size,"delta=",delta[jj])
                        rms  = np.std(z-ptheta-model)
                        chi2  = np.sum(w*(z-ptheta-model-delta)**2)/w.size
                        print('tha=',tha,'tdec=',tdec,'bsdx=',bsdx,'bsdy=',bsdy,"rms(ptheta):",np.std(ptheta),"rms(res):",rms0,"->",rms)
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


    mstyle="."
    conversion = 1.
    conversion = 3600.
    if conversion == 1: unit="deg"
    elif conversion == 3600: unit="arcsec"
    
    for x in [HA,DEC] :
        title="field-rotation-vs-"+x
        fig=plt.figure(title)
        first=True
        offset = np.zeros(z.size)
        nightmin=int(np.min(table["NIGHT"][ii]))
        nightmax=int(np.max(table["NIGHT"][ii]))
        emin=int(np.min(table["EXPID"][ii]))
        emax=int(np.max(table["EXPID"][ii]))    
        label  = "measurements: nights {}-{}".format(nightmin,nightmax)
        if args.fit_offset :
            offset = np.mean(z-ptheta)
            print("file #{} offset={}".format(fid,np.mean(offset)))
            label  += " + offset"
            
        plt.scatter(table[x][ii],(z-offset)*conversion,label=label,alpha=0.8)
        mlabel="model: polar misalignment + refraction"
        plt.plot(table[x][ii],ptheta*conversion,mstyle,c="k",label=mlabel)
        legend_title='with ME={:d}" MA={:d}"'.format(int(me),int(ma))
        plt.legend(title=legend_title)
        plt.xlabel(x)
        if args.fit_offset :
            plt.ylabel("measured rotation angle on sky + offset ({})".format(unit))
        else :
            plt.ylabel("measured rotation angle on sky ({})".format(unit))
        plt.grid()
        fig.savefig(title+".png")
        print("wrote",title+".png")
    if True :
        #for x in ["NIGHT","HEXPOS1","HEXPOS2","HEXPOS3","HEXPOS4","HEXPOS5","HEXPOS6"] :
        #for x in ["NIGHT","NSELECTED","CONTRAST","SIGSEEING"]:
        for x in ["AIRMASS","EXPID"]:
        #for x in ["CONTRAST","POINTINGDX","POINTINGDY","NSELECTED","OUTTEMP","TAIRITMP","TAIRTEMP","TTRUETBT","PRESSURE","CI-T1"] :
            title="field-rotation-residuals-vs-"+x
            fig=plt.figure(title)
            first=True
            offset = np.zeros(z.size)
            nightmin=int(np.min(table["NIGHT"][ii]))
            nightmax=int(np.max(table["NIGHT"][ii]))
            emin=int(np.min(table["EXPID"][ii]))
            emax=int(np.max(table["EXPID"][ii]))    
            label  = "measurements: nights {}-{}".format(nightmin,nightmax)
            if args.fit_offset :
                offset = np.mean(z-ptheta)
                print("file #{} offset={}".format(fid,np.mean(offset)))
                label  += " + offset"
            plt.scatter(table[x][ii],(z-offset-ptheta)*conversion,label=label,alpha=0.8)
            legend_title='with ME={:d}" MA={:d}"'.format(int(me),int(ma))
            plt.legend(title=legend_title)
            plt.xlabel(x)
            plt.ylabel("measured rotation angle residuals")
            plt.grid()
            fig.savefig(title+".png")
            print("wrote",title+".png")
            if x == "HEXPOS6" :
                for j in np.where(np.abs(table[x][ii])>20)[0] :
                    print("{}/{} {}={} Theta={} arcsec".format(int(table["NIGHT"][ii][j]),int(table["EXPID"][ii][j]),x,table[x][ii][j],(z-offset-ptheta)[j]*3600.))
            if x == "NIGHT" :
                for j in np.where(table[x][ii]==20190430)[0] :
                    print("{}/{} {}={} Theta res={} arcsec".format(int(table["NIGHT"][ii][j]),int(table["EXPID"][ii][j]),x,table[x][ii][j],(z-offset-ptheta)[j]*3600.))
            if x == DEC :
                for j in np.where(table[x][ii]==32.)[0] :
                    print("{}/{} {}={} Theta res={} arcsec".format(int(table["NIGHT"][ii][j]),int(table["EXPID"][ii][j]),x,table[x][ii][j],(z-offset-ptheta)[j]*3600.))

    if 1 : # vs HA,DEC
        plt.figure("theta-ha-dec")
        x=HA
        y=DEC
        plt.scatter(table[x][ii],table[y][ii],c=z)
        print("rms={}".format(np.std(z)))
        plt.colorbar()
        plt.xlabel(x)
        plt.ylabel(y+" (<{:d})".format(decmax))
        
    if 1 : # vs HA,DEC
        plt.figure("theta-ha-dec-residuals")
        x=HA
        y=DEC
        print(table[x][ii].shape)
        print(table[y][ii].shape)
        print(z.shape)
        offset = np.median(z-ptheta)
        plt.scatter(table[x][ii],table[y][ii],c=z-ptheta-offset)
        print("rms={}".format(np.std(z-ptheta-offset)))
        plt.colorbar()
        plt.xlabel(x)
        plt.ylabel(y+" (<{:d})".format(decmax))
    if 1 : # vs RA,DEC
        plt.figure("theta-ra-dec-residuals")
        x=RA
        y=DEC
        print(table[x][ii].shape)
        print(table[y][ii].shape)
        print(z.shape)
        offset = np.median(z-ptheta)
        plt.scatter(table[x][ii],table[y][ii],c=z-ptheta-offset)
        print("rms={}".format(np.std(z-ptheta-offset)))
        plt.colorbar()
        plt.xlabel(x)
        plt.ylabel(y+" (<{:d})".format(decmax))
        

    if False :   
        selection=(np.abs(table[DEC][ii]-60.5)<0.5)&(table[HA][ii]<0)
        offset=np.mean(z[selection]-ptheta[selection])
        for j in np.where(selection)[0] :
            model=ptheta[j]+offset
            print("{:d}/{:08d} RA={:5.4f} Dec={:5.4f} HEXPOS=({},{},{},{},{},{}) Theta={:5.4f} model={:5.4f} Theta-model={:5.4f} deg".format(int(table["NIGHT"][ii][j]),int(table["EXPID"][ii][j]),table[RA][ii][j],table[DEC][ii][j],table["HEXPOS1"][ii][j],table["HEXPOS2"][ii][j],table["HEXPOS3"][ii][j],table["HEXPOS4"][ii][j],table["HEXPOS5"][ii][j],table["HEXPOS6"][ii][j],table["THETA"][ii][j],model,table["THETA"][ii][j]-model))

        for x in [RA,DEC,HA,"NIGHT","TTRUETBT","PARALLAC"] :
            title="selected-field-rotation-vs-"+x
            fig=plt.figure(title)
            first=True

            nightmin=int(np.min(table["NIGHT"][ii][selection]))
            nightmax=int(np.max(table["NIGHT"][ii][selection]))
            emin=int(np.min(table["EXPID"][ii][selection]))
            emax=int(np.max(table["EXPID"][ii][selection]))    
            label  = "measurements: nights {}-{}".format(nightmin,nightmax)
            if args.fit_offset :
                offset = np.mean(z-ptheta)
                print("file #{} offset={}".format(fid,np.mean(offset)))
                label  += " + offset"

            plt.scatter(table[x][ii][selection],(z-offset)[selection]*conversion,label=label,alpha=0.8)
            mlabel="model: polar misalignment + refraction"
            plt.plot(table[x][ii][selection],ptheta[selection]*conversion,mstyle,c="k",label=mlabel)
            legend_title='with ME={:d}" MA={:d}"'.format(int(me),int(ma))
            plt.legend(title=legend_title)
            plt.xlabel(x)
            if args.fit_offset :
                plt.ylabel("measured rotation angle on sky + offset ({})".format(unit))
            else :
                plt.ylabel("measured rotation angle on sky ({})".format(unit))
            plt.grid()
            #fig.savefig(title+".png")
            #print("wrote",title+".png")
    if 1 :
        plt.figure("theta-residuals-hist")
        bins=np.linspace(-0.1,0.1,100)
        offset = np.mean(z-ptheta)
        h0,bins=np.histogram(z-ptheta-offset,bins=bins)
        theta=bins[:-1]+(bins[1]-bins[0])/2.
        plt.errorbar(theta,h0,np.sqrt(h0),fmt='o-',label="DEC<{:d}".format(int(decmax)))
        plt.grid()
        plt.xlabel("Theta residual (deg)")
plt.show()
