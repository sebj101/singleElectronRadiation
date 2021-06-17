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

def RunMagneticFieldPlots():
    BKGFIELD = 1. # Tesla
    
    f = ROOT.TF1("fBathtub; z [cm]; Total field - main field [mT]", BFieldBathtub, -40., 40., 3)
    f.SetParameter(0, 4.5)
    f.SetParameter(1, 20)
    f.SetParameter(2, 0.1)
    f.SetNpx(1000)
    
    #c1 = ROOT.TCanvas("c1", "c1")
    #c1.cd()
    #f.Draw("lr")
    

    f2 = ROOT.TF1("fBathtubBkg; z [cm]; Total field [T]", BFieldBathtubPlusBkg, -40., 40., 4)
    f2.SetParameter(0, 100)
    f2.SetParameter(1, 20)
    f2.SetParameter(2, 0.01)
    f2.SetParameter(3, 1.)
    f2.GetYaxis().SetRangeUser(0.9, 1)
    f2.SetNpx(1000)
    f2.Draw("lr")
    
    text = input()
    
# Run the code
if __name__ == "__main__":
    RunMagneticFieldPlots()
