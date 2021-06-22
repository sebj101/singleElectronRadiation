#!/usr/bin/python3
# MagneticFieldPlots.py
# Make plots of the magnetic field using the bathtub approximation

import numpy as np
import matplotlib.pyplot as plt
import ROOT

def BFieldBathtub(x, p):
    # Bathtub approximation for B field generated using two coils
    # zPos is position in field
    # l0 is measure of field gradient in curved region
    # l1 is width of the flat region
    # b0 is field scale
    field = 0.
    zPos = x[0]
    l0 = p[0]
    l1 = p[1]
    b0 = p[2]
    
    if zPos < -l1/2. :
        field = b0 * (1 + ((zPos + l1/2)**2)/l0**2 )
    elif zPos >= -l1/2. and zPos <= l1/2. :
        field = b0
    else :
        field = b0 * (1 + ((zPos - l1/2)**2)/l0**2 )
        
    return field

def BFieldBathtubPlusBkg(x, p):
    # Bathtub approximation for B field generated using two coils
    # zPos is position in field
    # l0 is measure of field gradient in curved region
    # l1 is width of the flat region
    # b0 is field scale
    #
    field = 0.
    zPos = x[0]
    l0 = p[0]
    l1 = p[1]
    b0 = p[2]
    bBkg = p[3]
    
    if zPos < -l1/2. :
        field = bBkg + b0 * (1 + ((zPos + l1/2)**2)/l0**2 )
    elif zPos >= -l1/2. and zPos <= l1/2. :
        field = bBkg + b0
    else :
        field = bBkg + b0 * (1 + ((zPos - l1/2)**2)/l0**2 )
        
    return field

def BFieldHarmonic(x, p):
    # Harmonic trap using a single coil
    # zPos is position in field
    # l0 is measure of scale of curvature
    # b0 is field scale
    zPos = x[0]
    l0 = p[0]
    b0 = p[1]

    field = b0 * (1 + zPos*zPos/(l0*l0))

    return field

def BFieldHarmonicPlusBkg(x, p):
    # Harmonic trap using a single coil with a background field
    # zPos is position in field
    # l0 is measure of scale of curvature
    # b0 is field scale
    zPos = x[0]
    l0 = p[0]
    b0 = p[1]
    bBkg = p[2]

    field = b0 * (1 + zPos*zPos/(l0*l0)) + bBkg

    return field

def RunMagneticFieldPlots():
    BKGFIELD = 1. # Tesla
    
    f1 = ROOT.TF1("fBathtub; z [cm]; Total field - main field [mT]", BFieldBathtub, -15., 15., 3)
    f1.SetParameter(0, 10.)
    f1.SetParameter(1, 15.)
    f1.SetParameter(2, 1.)
    f1.SetNpx(1000)
    
    #c1 = ROOT.TCanvas("c1", "c1")
    #c1.cd()
    f1.Draw("lr")
    
    f2 = ROOT.TF1("fBathtubBkg; z [cm]; Total field [T]", BFieldBathtubPlusBkg, -40., 40., 4)
    f2.SetParameter(0, 4.5)
    f2.SetParameter(1, 20)
    f2.SetParameter(2, 0.0001)
    f2.SetParameter(3, 1.)
    f2.GetYaxis().SetRangeUser(0.9, 1)
    f2.SetNpx(1000)
    #f2.Draw("lr")

    fHarm = ROOT.TF1("fHarm; z [cm]; Total field [mT]", BFieldHarmonic, -2., 2., 2)
    fHarm.SetParameter(0, np.sqrt(1000.))
    fHarm.SetParameter(1, 1)
    fHarm.SetNpx(1000)
    #fHarm.Draw("lr")
    
    text = input()
    
# Run the code
if __name__ == "__main__":
    RunMagneticFieldPlots()
