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

def BFieldBathtubPlusBkg(zPos, l0, l1, b0, bBkg):
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

def ZPositionHarmonic(time, v0, l0, b0, bBkg, trapDepth):
    thetaBot = np.arcsin(np.sqrt(1 - trapDepth/bBkg))
    axFreq = v0 * np.sin(thetaBot) / l0
    zMax = l0 * 1/np.tan(thetaBot)
    z = zMax * np.sin(axFreq * time)
    return z

def AngCyclFreqHarmonicFromTime(time, ke, v0, l0, b0, bBkg, trapDepth):
    premult = CalcAngCyclotronFreq(bBkg, ke)
    thetaBot = np.arcsin(np.sqrt(1 - trapDepth/bBkg))
    axFreq = v0 * np.sin(thetaBot) / l0
    zMax = l0 * 1/np.tan(thetaBot)
    freq = premult * (1 + zMax*zMax/(2*l0*l0) - zMax*zMax/(2*l0*l0)*np.cos(2*axFreq*time))
    return freq

def ZVelocityHarmonic(time, v0, l0, b0, bBkg, trapDepth):
    thetaBot = np.arcsin(np.sqrt(1 - trapDepth/bBkg))
    axFreq = v0 * np.sin(thetaBot) / l0
    zMax = l0 * 1/np.tan(thetaBot)
    v = zMax * axFreq* np.cos(axFreq * time) / 100.
    return v

def ReceiverFreqDoppler(vel, f0):
    # vel: Longitudinal velocity towards receiver
    # f0: Frequency at source
    f_r = np.sqrt( (1 - vel / constant.CLIGHT)/(1 + vel / constant.CLIGHT) ) * f0
    return f_r
