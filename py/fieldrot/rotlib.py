#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys


D2R = np.pi/180.

def hadec2thetaphi(ha,dec) : # with ha dec in degrees
    theta = (90-dec)*D2R
    phi   = ha*D2R
    return theta,phi

def thetaphi2hadec(theta,phi) : # with ha dec in degrees
    dec = 90-theta/D2R
    ha  = (phi/D2R)%360.
    return ha,dec

def unit_vector(theta,phi) :
    ct=np.cos(theta)
    st=np.sin(theta)
    cp=np.cos(phi)
    sp=np.sin(phi)
    return np.array([st*cp,st*sp,ct])


def hadec2xy(ha,dec,cha,cdec) :
    """
    dec is up (y>0)
    ra is left (x<0)
    ha is right (x>0) (ha=lst-ra)
    CIN is top (y>0)
    CIW is right (x>0)
    CIE is left (x<0)
    CIS is bottom (y<0)
    """
    t,p = hadec2thetaphi(ha,dec)
    vec = unit_vector(t,p)
    t,p = hadec2thetaphi(cha,cdec)
    cp=np.cos(p)
    sp=np.sin(p)
    ct=np.cos(t)
    st=np.sin(t)
    rp=np.array([[cp,sp,0],[-sp,cp,0],[0,0,1]])
    rt=np.array([[ct,0,-st],[0,1,0],[+st,0,ct]])
    tmp = rt.dot(rp).dot(vec)
    x=tmp[1]/D2R
    y=-tmp[0]/D2R
    return x,y # deg
    

def xy2hadec(x,y,cha,cdec) :
    t,p = hadec2thetaphi(cha,cdec)
    cp=np.cos(p)
    sp=np.sin(p)
    ct=np.cos(t)
    st=np.sin(t)
    rp=np.array([[cp,-sp,0],[sp,cp,0],[0,0,1]])
    rt=np.array([[ct,0,st],[0,1,0],[-st,0,ct]])

    xr = x*D2R
    yr = y*D2R
    zr = np.sqrt(1-xr**2-yr**2)
    vec = np.array([-yr,xr,zr])
    vec = rp.dot(rt).dot(vec)
    t,p = vector2thetaphi(vec)
    ha,dec = thetaphi2hadec(t,p)
    return ha,dec

def vector2thetaphi(vec) :
    phi=np.arctan2(vec[1],vec[0])
    theta=np.arccos(vec[2])
    return theta,phi

def rotation_matrix(axis,theta) :
    # rotate around axis by an angle
    ct=np.cos(theta)
    st=np.sin(theta)
    u=axis
    uut=np.outer(u,u)
    uxx=np.array([[0,-u[2],u[1]],[u[2],0,-u[0]],[-u[1],u[0],0]])
    return ct*np.eye(axis.size)+st*uxx+(1-ct)*uut

def rotation_hadec(ha,dec,matrix) :
    t,p=hadec2thetaphi(ha,dec)
    v=unit_vector(t,p)
    v2=matrix.dot(v)
    t,p=vector2thetaphi(v2)
    return thetaphi2hadec(t,p)

def field_delta_in_deg(hacen,deccen,ha1,dec1,ha2,dec2) :
    cdec=np.cos(deccen*D2R)
    x1 =(ha1-hacen)*cdec
    y1 =(dec1-deccen)
    x2 =(ha2-hacen)*cdec
    y2 =(dec2-deccen)
    return np.sqrt((x1-x2)**2+(y1-y2)**2)


def field_radius_in_deg(hacen,deccen,ha1,dec1) :
    cdec=np.cos(deccen*D2R)
    x1 =(ha1-hacen)*cdec
    y1 =(dec1-deccen)
    return np.sqrt((x1**2+y1**2))
    
# measure the rotation angle when going from (ha1,dec1 to ha2,dec2)
# angle is positive for counterclockwise rotation
def field_rotation_angle_in_rad(ha1,dec1,ha2,dec2,dx=0,dy=0) :
    
    x1,y1=hadec2xy(ha1,dec1,ha1[0],dec1[0])
    x2,y2=hadec2xy(ha2,dec2,ha2[0],dec2[0])
    x1 -= np.mean(x1)
    y1 -= np.mean(y1)
    x2 -= np.mean(x2)
    y2 -= np.mean(y2)
    y2 += dy
    x2 += dx
    return np.arcsin((x1*y2-x2*y1)/np.sqrt((x1**2+y1**2)*(x2**2+y2**2)))

# measure the rotation angle when going from (ha1,dec1 to ha2,dec2)
# angle is positive for counterclockwise rotation
def field_rotation_angle_in_rad(ha1,dec1,ha2,dec2,dx=0,dy=0) :
    
    
    x1,y1=hadec2xy(ha1,dec1,ha1[0],dec1[0])
    if dx==0 and dy==0 :
        x2,y2=hadec2xy(ha2,dec2,ha2[0],dec2[0])
    else :
        offset_ha20,offset_dec20=xy2hadec(dx,dy,ha2[0],dec2[0])
        x2,y2=hadec2xy(ha2,dec2,offset_ha20,offset_dec20)
    
    x1 -= np.mean(x1)
    y1 -= np.mean(y1)
    x2 -= np.mean(x2)
    y2 -= np.mean(y2)
    return np.arcsin((x1*y2-x2*y1)/np.sqrt((x1**2+y1**2)*(x2**2+y2**2)))

def hadec2altaz(ha,dec,lat) : 
    cha  = np.cos(ha*D2R)
    sha  = np.sin(ha*D2R)
    cdec = np.cos(dec*D2R)
    sdec = np.sin(dec*D2R)
    clat = np.cos(lat*D2R)
    slat = np.sin(lat*D2R)
    
    x = - cha * cdec * slat + sdec * clat
    y = - sha * cdec
    z = cha * cdec * clat + sdec * slat
    r = np.sqrt(x**2 + y**2)
    az = np.arctan2(y,x)/D2R
    alt = np.arctan2(z,r)/D2R
    az = az%360
    if isinstance(az,float) :
        if az<0 :
            az += 360.
    else :
        az[az<0] += 360.
    return alt,az

def altaz2hadec(alt,az,lat) : 
    calt = np.cos(alt*D2R)
    salt = np.sin(alt*D2R)
    caz  = np.cos(az*D2R)
    saz  = np.sin(az*D2R)
    clat = np.cos(lat*D2R)
    slat = np.sin(lat*D2R)
    ha = np.arctan2( -saz*calt, -caz*slat*calt+salt*clat) / D2R
    ha = ha%360
    if isinstance(ha,float) :
        if ha<0 :
            ha += 360.
    else :
        ha[ha<0] += 360.
    dec = np.arcsin(slat*salt+clat*calt*caz)/ D2R
    return ha,dec

def field_rotation(ha,dec,ma=25/3600.,me=-103/3600.,lat=31.96403) :
    """
    Field rotation prediction given field center and polar axis misalignment

    Args :
      ha: hour angle in degrees of field center
      dec:  declination in degrees of field center
      me: degrees, in the northern hemisphere, positive ME means that the pole of the mounting is below the true (unrefracted) north (so elevation=latitude-me)
      ma: degrees, in the northern hemisphere, and positive MA means that the pole of the mounting is to the right of due north
      lat: degrees, latitude of telescope

    Returns:
      angle in degrees of stars rotation with respect to instrument.
      angle is increasing from X_tan to Y_tan,
                   where X_tan is along HA and increases with HA (LST-RA),
                   and Y_tan is along Dec and increases with Dec
     (angle goes counterclock-wise if X_tan points to the right and Y_tan points to the top)
    """

    # convert polar axis to HA and Dec
    axis_ha,axis_dec=altaz2hadec(alt=lat-me,az=ma,lat=lat)

    # now to theta phi
    axis_theta,axis_phi = hadec2thetaphi(axis_ha,axis_dec)

    # define a rotation matrix to move the polar axis to the north
    # vector product
    cross = np.cross(unit_vector(axis_theta,axis_phi),unit_vector(0,0))
    norme = np.sqrt(np.sum(cross**2))
    if norme > 0 :
        drot = rotation_matrix(cross/norme,np.arcsin(norme))
    else :
        drot = np.eye(3)

    # take a fiducial set of points in field of view
    phi = np.linspace(0,2*np.pi,10)
    x1  = np.cos(phi)
    y1  = np.sin(phi)

    aha   = np.atleast_1d(ha)
    adec  = np.atleast_1d(dec)
    angle = np.zeros(aha.shape)

    for i in range(aha.size) :
        # convert to ha and dec given the field center
        ha1,dec1 = xy2hadec(x1,y1,aha[i],adec[i])
        # rotate
        ha2,dec2   = rotation_hadec(ha1,dec1,drot)
        # convert to alt az
        alt,az     = hadec2altaz(ha2,dec2,lat)
        # apply refraction
        alt += 79./3600.*np.tan(30.*D2R)/np.tan(alt*D2R)  # deg , refraction per point in field
        # convert back to ha dec
        ha2,dec2   = altaz2hadec(alt,az,lat)
        # measure rotation
        angle[i] = np.mean(field_rotation_angle_in_rad(ha1,dec1,ha2,dec2,0,0))/D2R

    if isinstance(ha,np.ndarray) :
        return angle
    else :
        return angle[0]

#############################################################
# TESTS
#############################################################


def main() :
    ###################################
    # testing HA,Dec to tangent plane
    ###################################
    cha=12
    cdec=24
    x1 = np.random.uniform(size=3)-0.5
    y1 = np.random.uniform(size=3)-0.5
    ha1,dec1 = xy2hadec(x1,y1,cha,cdec)
    x2,y2 = hadec2xy(ha1,dec1,cha,cdec)
    print("x1=",x1)
    print("y1=",y1)
    print("x2=",x2)
    print("y2=",y2)

    assert(np.all(np.abs(x1-x2)<1e-6))
    assert(np.all(np.abs(y1-y2)<1e-6))



    cdec = 0.
    cha  = 30.
    plt.figure("HA Dec to tangent plane")
    x,y = hadec2xy(cha+0.,cdec+2.,cha,cdec)
    plt.plot(x,y,"o",label='dHA=0 dDec=2')
    print("x=",x)
    print("y=",y)
    plt.xlabel("x_t")
    plt.ylabel("y_t")
    
    x,y = hadec2xy(cha+1.,cdec+0.,cha,cdec)
    plt.plot(x,y,"x",label='dHA=1 (dRA=-1) dDec=0')
    print("x=",x)
    print("y=",y)
    r=3
    plt.xlim([-r,r])
    plt.ylim([-r,r])
    plt.grid()
    plt.legend()
    #plt.show(); sys.exit(12)

    ###################################
    # testing HA Dec to Alt Az
    ###################################
    lat = 31
    dec = 90.-25*np.ones(40)
    ha  = np.linspace(-100,100,40)%360.
    alt,az=hadec2altaz(ha,dec,lat)
    plt.figure("HA Dec to Alt Az")
    az[az>180] -= 360.
    plt.scatter(az,alt,c=(ha-360*(ha>180)))
    plt.plot(0,lat,"+",c="k")
    plt.colorbar(label="HA")
    plt.xlabel("WEST          AZ             EAST")              
    ha2,dec2=altaz2hadec(alt,az,lat)
    print(ha,dec)
    print(ha2,dec2)
    assert(np.all(np.abs(ha-ha2)<1e-6))
    assert(np.all(np.abs(dec-dec2)<1e-6))



    #########################################
    # testing measurement of rotation angle
    #########################################

    tmp    = np.linspace(0,2*np.pi,10)
    x1  = np.cos(tmp)
    y1  = np.sin(tmp)

    # counter clockwise rotation
    a   = 5.*D2R
    x2  = x1*np.cos(a)+y1*np.sin(-a)
    y2  = x1*np.sin(a)+y1*np.cos(a)

    plt.figure("measurement of rotation angle")
    plt.plot(x1,y1,"o")
    plt.plot(x2,y2,"x")


    for cha in np.linspace(0,360,11) :
        for cdec in np.linspace(0,80,10) :
            ha1,dec1 = xy2hadec(x1,y1,cha,cdec)
            ha2,dec2 = xy2hadec(x2,y2,cha,cdec)
            angle = field_rotation_angle_in_rad(ha1,dec1,ha2,dec2)
            print("dec=",cdec,"rotation angle result=",np.mean(angle/D2R),"expectation=",a/D2R)

    #######################################################
    # testing rotation angle due to polar misalignment
    #######################################################


    cha=0
    cdec=40.

    pole_ha  = -10.
    pole_dec = 90-1.
    px1,py1  = hadec2xy(pole_ha,pole_dec,cha,cdec)
    x1  = np.random.uniform(size=4)-0.5
    y1  = np.random.uniform(size=4)-0.5

    plt.figure("rotation angle due to polar misalignment")
    plt.subplot(1,3,1)
    plt.plot(x1,y1,"o",c="gray")
    plt.plot(px1,py1,"o",c="red")
    r=4*(90-cdec)
    plt.xlim([-r,r])
    plt.ylim([-r,r])
    plt.grid()

    pt,pp = hadec2thetaphi(pole_ha,pole_dec)
    pv    =  unit_vector(pt,pp)

    tdec = 90-pole_dec
    tha  = pole_ha
    cross = np.cross(unit_vector(0,0),unit_vector(-tdec*D2R,tha*D2R))
    norme = np.sqrt(np.sum(cross**2))
    drot=rotation_matrix(cross/norme,np.arcsin(norme))

    ha1,dec1 = xy2hadec(x1,y1,cha,cdec)
    ha2,dec2 = rotation_hadec(ha1,dec1,drot)
    angle    = field_rotation_angle_in_rad(ha1,dec1,ha2,dec2)
    print("mean field rotation angle=",np.mean(angle)/D2R,"deg")
    print("mean delta dec=",np.mean(dec2)-np.mean(dec1))
    print("mean delta ha=",np.mean(ha2)-np.mean(ha1))
    
    x2,y2    = hadec2xy(ha2,dec2,cha,cdec)

    ph2,pd2 = rotation_hadec(pole_ha,pole_dec,drot)
    print("HA,Dec of new pole=",ph2,pd2)
    px2,py2    = hadec2xy(ph2,pd2,cha,cdec)


    plt.subplot(1,3,2)
    plt.plot(x2,y2,"o",c="gray")
    plt.plot(px2,py2,"o",c="red")
    plt.plot(x1-np.mean(x1)+np.mean(x2),y1-np.mean(y1)+np.mean(y2),"x",c="k")

    plt.xlim([-r,r])
    plt.ylim([-r,r])
    plt.grid()


    ##############################
    # other calculation
    ##############################

    lat = 31
    
    print("pole ha dec=",pole_ha,pole_dec)
    
    palt,paz = hadec2altaz(pole_ha,pole_dec,lat)
    print("pole alt az=",palt,paz)
    
    # "in the northern hemisphere, positive MA means that the pole of the mounting is to the right of due north"
    # "in the northern hemisphere, positive ME means that the pole of the mounting is below the true (unrefracted) north" 
    
    ma = paz
    if ma > 180 : ma -= 360.
    me = -(palt-lat)
    print("ma,me=",ma,me)

    signe = 1.
    ha1  = pole_ha
    dec1 = pole_dec
    dha  = -me*np.sin(ha1*D2R)*np.tan(dec1*D2R) + ma*np.cos(ha1*D2R)*np.tan(dec1*D2R)
    ddec = -me*np.cos(ha1*D2R)                  - ma*np.sin(ha1*D2R)
    ha3        = ha1  + signe*dha
    dec3       = dec1 + signe*ddec
    print("HA,Dec of new pole= ",ha3,dec3)
    px3,py3 = hadec2xy(ha3,dec3,cha,cdec)
    
    ha1,dec1 = xy2hadec(x1,y1,cha,cdec)
    dha  = -me*np.sin(ha1*D2R)*np.tan(dec1*D2R) + ma*np.cos(ha1*D2R)*np.tan(dec1*D2R)
    ddec = -me*np.cos(ha1*D2R)                  - ma*np.sin(ha1*D2R)
    ha3  = ha1  + signe*dha
    dec3 = dec1 + signe*ddec
    x3,y3 = hadec2xy(ha3,dec3,cha,cdec)
    
    angle    = field_rotation_angle_in_rad(ha1,dec1,ha3,dec3)
    print("mean field rotation angle (other method, approximate)=",np.mean(angle)/D2R,"deg")
    print("mean delta dec (other method, approximate)=",np.mean(dec3)-np.mean(dec1))
    print("mean delta ha (other method, approximate)=",np.mean(ha3)-np.mean(ha1))
    
    plt.subplot(1,3,3)
    plt.plot(x3,y3,"o",c="gray")
    plt.plot(px3,py3,"o",c="red")
    
    plt.xlim([-r,r])
    plt.ylim([-r,r])
    plt.grid()
    
    plt.show()

    sys.exit(12)


    #plt.plot(x1,y1,"o")
    #plt.plot(x2,y2,"o")
    #plt.show()

    
if __name__ == '__main__' :

    main()

