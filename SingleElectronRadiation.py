#!/usr/bin/python3
# SingleElectronRadiation.py
# Single electron in a magnetic field

import constant
import functions as func
import numpy as np
import matplotlib.pyplot as pyplot
import cv2

def RunSingleElectronSim():
    # Field constants for harmonic trap field
    BBKG = 1. # Tesla
    L0   = 20. # cm
    B0   = 1. # Tesla
    TRAPDEPTH = 0.004 # Tesla

    EKE = (18.6 * 10**3) * constant.COULOMBCHARGE # Electron KE in Joules
    
    nSteps = 500
    maxTime =  50000 / func.CalcCyclotronFreq(BBKG, EKE) 
    timeStepSize = maxTime / nSteps

    print(func.CalcCyclotronFreq(BBKG, EKE) )
    
    V0 = 0.263 * constant.CLIGHT

    TimeArray    = np.zeros(nSteps)
    BFieldArray  = np.zeros(nSteps)
    ZPosArray    = np.zeros(nSteps)
    AngFreqArray = np.zeros(nSteps)
    
    time = 0.
    
    for i in range(nSteps):
        #print(i)

        TimeArray[i]    = time
        ZPosArray[i]    = func.ZPositionHarmonic(time, V0, L0, B0, BBKG, TRAPDEPTH)
        BFieldArray[i]  = func.BFieldHarmonicFromTime(time, V0, L0, B0, BBKG, TRAPDEPTH)
        AngFreqArray[i] = func.AngCyclFreqHarmonicFromTime(time, EKE, V0, L0, B0, BBKG, TRAPDEPTH)
        
        time = time + timeStepSize

    # Now plot things
    fig1, ax1 = pyplot.subplots(nrows=3, ncols=1, figsize=[9, 15])
    f1graph0 = ax1[0].plot(TimeArray, ZPosArray)
    f1graph1 = ax1[1].plot(TimeArray, BFieldArray)
    f1graph2 = ax1[2].plot(TimeArray, AngFreqArray)
    ax1[0].set_title("Z position")
    ax1[0].set_xlabel('t [s]')
    ax1[0].set_ylabel('z [cm]')
    ax1[1].set_title("B field experienced by electron")
    ax1[1].set_xlabel('t [s]')
    ax1[1].set_ylabel('B [T]')    
    ax1[2].set_title("Electron angular frequency")
    ax1[2].set_xlabel('t [s]')
    ax1[2].set_ylabel(r'$\Omega_{c}\, [s^{-1}]$')
    
    pyplot.show()

# Run the code
if __name__ == "__main__":
    RunSingleElectronSim()
