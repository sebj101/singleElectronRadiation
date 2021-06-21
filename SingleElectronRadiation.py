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
    
    nSteps = 8961
    maxTime =  8 * 10**-6 #50000 / func.CalcCyclotronFreq(BBKG, EKE) 
    timeStepSize = maxTime / nSteps

    print(func.CalcCyclotronFreq(BBKG, EKE) )
    
    V0 = 0.263 * constant.CLIGHT

    TimeArray    = np.zeros(nSteps)
    BFieldArray  = np.zeros(nSteps)
    ZPosArray    = np.zeros(nSteps)
    ZVelArray    = np.zeros(nSteps)
    AngFreqArray = np.zeros(nSteps)
    CyclFreqArray = np.zeros(nSteps)
    ReceiverFreqArrayNoDoppler = np.zeros(nSteps)
    ReceiverFreqArrayDoppler = np.zeros(nSteps)
    
    time = 1 * 10**-9
    
    for i in range(nSteps):
        TimeArray[i]     = time
        ZPosArray[i]     = func.ZPositionHarmonic(time, V0, L0, B0, BBKG, TRAPDEPTH)
        ZVelArray[i]     = func.ZVelocityHarmonic(time, V0, L0, B0, BBKG, TRAPDEPTH)
        BFieldArray[i]   = func.BFieldHarmonicFromTime(time, V0, L0, B0, BBKG, TRAPDEPTH)
        AngFreqArray[i]  = func.AngCyclFreqHarmonicFromTime(time, EKE, V0, L0, B0, BBKG, TRAPDEPTH)
        CyclFreqArray[i] = func.AngCyclFreqHarmonicFromTime(time, EKE, V0, L0, B0, BBKG, TRAPDEPTH)/(2. * np.pi)
        
        # Assume there is detector above midpoint of CRES region
        # say d cm away from field axis
        # what frequencies does it see over time?
        
        # Calculate angle between electron and detector
        D = 5.
        eDetAngle = np.arctan(D / ZPosArray[i])
        ReceiverFreqArrayNoDoppler[i] = func.ReceiverFreqDoppler(0, func.AngCyclFreqHarmonicFromTime(time, EKE, V0, L0, B0, BBKG, TRAPDEPTH) / (2 * np.pi))
        ReceiverFreqArrayDoppler[i] = func.ReceiverFreqDoppler(ZVelArray[i] * np.cos(eDetAngle), func.AngCyclFreqHarmonicFromTime(time, EKE, V0, L0, B0, BBKG, TRAPDEPTH) / (2 * np.pi))
        
        time = time + timeStepSize

    # Now plot things
    fig1, ax1 = pyplot.subplots(nrows=3, ncols=2, figsize=[16, 15])
    f1graph0 = ax1[0][0].plot(TimeArray, ZPosArray)
    f1graph1 = ax1[1][0].plot(TimeArray, BFieldArray)
    f1graph2 = ax1[2][0].plot(TimeArray, AngFreqArray)
    f1graph3 = ax1[0][1].plot(TimeArray, ZVelArray)
    f1graph5 = ax1[2][1].plot(TimeArray, CyclFreqArray)
    ax1[0][0].set_title("Z position")
    ax1[0][0].set_xlabel('t [s]')
    ax1[0][0].set_ylabel('z [cm]')
    ax1[1][0].set_title("B field experienced by electron")
    ax1[1][0].set_xlabel('t [s]')
    ax1[1][0].set_ylabel('B [T]')    
    ax1[2][0].set_title("Electron angular frequency")
    ax1[2][0].set_xlabel('t [s]')
    ax1[2][0].set_ylabel(r'$\Omega_{c}\, [s^{-1}]$')
    ax1[0][1].set_title("Electron axial velocity")
    ax1[0][1].set_xlabel('t [s]')
    ax1[0][1].set_ylabel(r'$v_{z}\, [ms^{-1}]$')

    ax1[2][1].set_title("Electron frequency")
    ax1[2][1].set_xlabel('t [s]')
    ax1[2][1].set_ylabel(r'$\Omega_{c}\, [s^{-1}]$')
    
    fig2, ax2 = pyplot.subplots(nrows=1, ncols=1)
    f2graph0 = ax2.hist(ReceiverFreqArrayNoDoppler, bins=100)
    ax2.set_xlabel("Frequency [Hz]")
    ax2.set_ylabel("A. U.")

    fig3, ax3 = pyplot.subplots(nrows=1, ncols=1)
    f3graph0 = ax3.hist(ReceiverFreqArrayDoppler, bins=100)
    ax3.set_xlabel("Frequency [Hz]")
    ax3.set_ylabel("A. U.")
    
    pyplot.show()

# Run the code
if __name__ == "__main__":
    RunSingleElectronSim()
