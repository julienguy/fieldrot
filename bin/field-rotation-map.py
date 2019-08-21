#!/usr/bin/env python

import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import argparse

from fieldrot.rotlib import *



ha1d  = np.linspace(-88,88,60)
dec1d = np.linspace(-40,80,(8+4)*5)

ha   = np.tile(ha1d,(dec1d.size,1))
dec  = np.tile(dec1d,(ha1d.size,1)).T
theta = np.zeros(ha.shape)
dthetadt = np.zeros(ha.shape)

# latitude
# 31:57:50.5 DESI-3182-v4
lat = 31.96403 #deg
telescope_polar_axis_dha = 1.
telescope_polar_axis_ddec = 0.0335 
boresight_offset_dx = -0.017

# see rotlib tests for verification
# compute rotation matrix for polar misalignment
t,p = hadec2thetaphi(telescope_polar_axis_dha,90-telescope_polar_axis_ddec)
cross = np.cross(unit_vector(t,p),unit_vector(0,0))
norme = np.sqrt(np.sum(cross**2))
if norme > 0 :
    drot = rotation_matrix(cross/norme,np.arcsin(norme))
else :
    drot = np.eye(3)

dt  = 0.1  # time in minutes
dha = dt/(60.*24)*360. # now degrees
 
# set of points to measure numerically field rotation
tmp    = [0,np.pi/2.,np.pi,3*np.pi/2.]
x2d  = np.cos(tmp)
y2d  = np.sin(tmp)
for j in range(ha.shape[0]) :
    for i in range(ha.shape[1]) :

        for step in range(2) :
            # bore sight 
            ha_tmp,dec_tmp = xy2hadec(boresight_offset_dx,0,ha[j,i]+(dha*(step==1)),dec[j,i])
            ha1,dec1 = xy2hadec(x2d,y2d,ha_tmp,dec_tmp)
            # polar misalignment
            ha2,dec2 = rotation_hadec(ha1,dec1,drot)
            # refraction
            alt,az     = hadec2altaz(ha2,dec2,lat)
            if np.mean(alt)>30 :
                R = 79./3600.*np.tan(30.*D2R)/np.tan(alt*D2R)  # deg , refraction per point in field
                alt += R 
                ha2,dec2  = altaz2hadec(alt,az,lat)
                if step==0 :
                    theta[j,i] =  np.mean(measure_field_rotation(ha1,dec1,ha2,dec2,0,0)) # deg
                elif step==1 and theta[j,i]!=0 :
                    dthetadt[j,i] =  (np.mean(measure_field_rotation(ha1,dec1,ha2,dec2,0,0))-theta[j,i])/dt # deg per minute
                
theta[theta==0]=np.nan
dthetadt[dthetadt==0]=np.nan

title="theta"
fig=plt.figure(title)
p=plt.subplot(1,1,1)
p.set_title("field rotation theta (deg)")
plt.imshow(theta,origin=0,aspect="auto",extent=(ha1d[0],ha1d[-1],dec1d[0],dec1d[-1]))
#plt.imshow(ha,origin=0,aspect="auto",extent=(ha1d[0],ha1d[-1],dec1d[0],dec1d[-1]))
plt.xlabel("HA (deg)")
plt.ylabel("Dec (deg)")
plt.colorbar()
fig.savefig(title+".png")

title="dthetadt"
fig=plt.figure(title)
p=plt.subplot(1,1,1)
p.set_title("time derivative of field rotation theta (arcsec/minute)")
plt.imshow(dthetadt*3600,origin=0,aspect="auto",extent=(ha1d[0],ha1d[-1],dec1d[0],dec1d[-1]))
#plt.imshow(ha,origin=0,aspect="auto",extent=(ha1d[0],ha1d[-1],dec1d[0],dec1d[-1]))
plt.xlabel("HA (deg)")
plt.ylabel("Dec (deg)")
plt.colorbar()
fig.savefig(title+".png")
plt.show()


