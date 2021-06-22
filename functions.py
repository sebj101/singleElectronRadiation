import numpy as np
import constant

# Function definitions

def CalcAngCyclotronFreq(BField, KE):
    # Angular cyclotron frequency from magnetic field and electron KE    
    freq = (constant.COULOMBCHARGE * BField)/(constant.ERESTMASS + KE/(constant.CLIGHT*constant.CLIGHT))
    return freq

def CalcCyclotronFreq(BField, KE):
    # Angular cyclotron frequency from magnetic field and electron KE (in Joules)
    freq = (constant.COULOMBCHARGE * BField)/((constant.ERESTMASS + KE/(constant.CLIGHT*constant.CLIGHT)) * 2 * constant.PI)
    return freq

### In a "harmonic trap" ########
### Here we assume there is no energy loss as radiation
def BFieldHarmonicFromPos(zPos, l0, b0):
    field = b0 * (1 + zPos**2/l0**2)
    return field

def BFieldHarmonicFromTime(time, v0, l0, b0, bBkg, trapDepth):
    thetaBot = np.arcsin(np.sqrt(1 - trapDepth/bBkg))
    axFreq = v0 * np.sin(thetaBot) / l0
    zMax = l0 * 1/np.tan(thetaBot)    
    field = b0 * ( 1 + zMax*zMax/(2*l0*l0) - zMax*zMax/(2*l0*l0)*np.cos(2*axFreq*time) )
    return field

def XZPositionHarmonic(time, ke, v0, l0, b0, bBkg, trapDepth):
    angFreq = CalcAngCyclotronFreq(bBkg, ke)
    thetaBot = np.arcsin(np.sqrt(1 - trapDepth/bBkg))
    axFreq = v0 * np.sin(thetaBot) / l0
    zMax = l0 * 1/np.tan(thetaBot)
    freq = angFreq * (1 + zMax*zMax/(2*l0*l0) - zMax*zMax/(2*l0*l0)*np.cos(2*axFreq*time))
    gyroradius = (1. / np.sqrt(1 - (v0/constant.CLIGHT)**2 )) * constant.ERESTMASS * v0 / ( constant.COULOMBCHARGE * bBkg )
    arr = np.array([0., 0.])
    arr[0] = gyroradius * np.sin(angFreq * time)
    arr[1] = zMax * np.sin(axFreq * time) / 100.
    return arr

def XZVelocityHarmonic(time, ke, v0, l0, b0, bBkg, trapDepth):
    angFreq = CalcAngCyclotronFreq(bBkg, ke)
    thetaBot = np.arcsin(np.sqrt(1 - trapDepth/bBkg))
    axFreq = v0 * np.sin(thetaBot) / l0
    zMax = l0 * 1/np.tan(thetaBot)
    freq = angFreq * (1 + zMax*zMax/(2*l0*l0) - zMax*zMax/(2*l0*l0)*np.cos(2*axFreq*time))
    gyroradius = (1. / np.sqrt(1 - (v0/constant.CLIGHT)**2 )) * constant.ERESTMASS * v0 / ( constant.COULOMBCHARGE * bBkg )
    arr = np.array([0., 0.])
    arr[0] = gyroradius * angFreq * np.sin(angFreq * time)
    arr[1] = zMax * axFreq* np.cos(axFreq * time) / 100.
    return arr

def AngCyclFreqHarmonicFromTime(time, ke, v0, l0, b0, bBkg, trapDepth):
    premult = CalcAngCyclotronFreq(bBkg, ke)
    thetaBot = np.arcsin(np.sqrt(1 - trapDepth/bBkg))
    axFreq = v0 * np.sin(thetaBot) / l0
    zMax = l0 * 1/np.tan(thetaBot)
    freq = premult * (1 + zMax*zMax/(2*l0*l0) - zMax*zMax/(2*l0*l0)*np.cos(2*axFreq*time))
    return freq

def ReceiverFreqDoppler(vel, f0):
    # vel: Longitudinal velocity towards receiver
    # f0: Frequency at source
    f_r = np.sqrt( (1 - vel / constant.CLIGHT)/(1 + vel / constant.CLIGHT) ) * f0
    return f_r

#### Bathtub trap ####
def BFieldBathtubPlusBkgFromPos(zPos, l0, l1, b0, bBkg):
    # Bathtub approximation for trapping B field generated using two coils
    # zPos is position in field
    # l0 is measure of field gradient in curved region
    # l1 is width of the flat region
    # b0 is field scale
    # bBkg is background field
    
    field = 0.
    
    if zPos < -l1/2. :
        field = bBkg + b0 * (1 + ((zPos + l1/2)**2)/l0**2 )
    elif zPos >= -l1/2. and zPos <= l1/2. :
        field = bBkg + b0
    else :
        field = bBkg + b0 * (1 + ((zPos - l1/2)**2)/l0**2 )
        
    return field

def BFieldBathtubFromTime(time, v0, l0, l1, b0, bBkg, trapDepth):
    l0 = l0 / 100.
    l1 = l1 / 100.
    thetaBot = np.arcsin(np.sqrt(1 - trapDepth/bBkg))
    zMax = l0 / np.tan(thetaBot)
    omegaA = v0 * np.sin(thetaBot) / l0
    
    t1 = l1 / (v0 * np.cos(thetaBot))
    t2 = t1 + np.pi / omegaA
    t3 = t1 + t2
    T  = 2 * t2

    tRem = time % T

    field = 0.
    
    if tRem > 0. and tRem < t1 :
        field = b0
    elif tRem > t1 and tRem < t2 :
        field = b0 * ( 1 + zMax*zMax/(2*l0*l0) - zMax*zMax/(2*l0*l0) * np.cos(2*omegaA*(tRem-t1)) )
    elif tRem > t2 and tRem < t3 :
        field = b0 
    elif tRem > t3 and tRem < T :
        field = b0 * ( 1 + zMax*zMax/(2*l0*l0) - zMax*zMax/(2*l0*l0) * np.cos(2*omegaA*(tRem-t3)) )     

    return field

def AngCyclFreqBathtubFromTime(time, v0, l0, l1, b0, bBkg, trapDepth):
    l0 = l0 / 100.
    l1 = l1 / 100.
    thetaBot = np.arcsin(np.sqrt(1 - trapDepth/bBkg))
    zMax = l0 / np.tan(thetaBot)
    omegaA = v0 * np.sin(thetaBot) / l0
    gamma = 1 / np.sqrt( 1 - v0*v0/(constant.CLIGHT*constant.CLIGHT) )
    
    premult = constant.COULOMBCHARGE * b0 / (gamma * constant.ERESTMASS)
    
    t1 = l1 / (v0 * np.cos(thetaBot))
    t2 = t1 + np.pi / omegaA
    t3 = t1 + t2
    T  = 2 * t2

    tRem = time % T

    freq = 0.
    
    if tRem > 0. and tRem < t1 :
        freq = premult
    elif tRem > t1 and tRem < t2 :
        freq = premult * ( 1 + zMax*zMax/(2*l0*l0) - zMax*zMax/(2*l0*l0) * np.cos(2*omegaA*(tRem-t1)) )
    elif tRem > t2 and tRem < t3 :
        freq = premult 
    elif tRem > t3 and tRem < T :
        freq = premult * ( 1 + zMax*zMax/(2*l0*l0) - zMax*zMax/(2*l0*l0) * np.cos(2*omegaA*(tRem-t3)) ) 

    return freq

def ZPosVelBathtubFromTime(time, v0, l0, l1, b0, bBkg, trapDepth):
    l0 = l0 / 100.
    l1 = l1 / 100.
    thetaBot = np.arcsin(np.sqrt(1 - trapDepth/bBkg))
    vz0 = v0 * np.cos(thetaBot)
    omegaA = v0 * np.sin(thetaBot) / l0
    zMax = l0 / np.tan(thetaBot)
        
    t1 = l1 / (v0 * np.cos(thetaBot))
    t2 = t1 + np.pi / omegaA
    t3 = t1 + t2
    T  = 2 * t2

    tRem = time % T
    
    zPos = 0.
    zVel = 0.
    if tRem > 0. and tRem < t1 :
        zPos = vz0 * tRem - l1/2.
        zVel = vz0
    elif tRem > t1 and tRem < t2 :
        zPos = zMax * np.sin( omegaA*(tRem-t1) ) + l1/2.
        zVel = zMax * omegaA * np.cos( omegaA*(tRem - t1) )
    elif tRem > t2 and tRem < t3 :
        zPos = -vz0 * (tRem - t2) + l1/2.
        zVel = -vz0
    elif tRem > t3 and tRem < T :
        zPos = -zMax * np.sin( omegaA*(tRem-t3) ) - l1/2.
        zVel = -zMax * omegaA * np.cos( omegaA*(tRem - t3) )
    return zPos, zVel
