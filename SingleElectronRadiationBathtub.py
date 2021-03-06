#!/usr/bin/python3
# SingleElectronRadiationBathtub.py
# Motion of an electron in bathtub magnetic trap

import sys
import constant
import functions as func
import numpy as np
import matplotlib.pyplot as pyplot

# Main functions
def RunSingleElectronBathtubSim():
    
    # Field constants for harmonic trap field
    BBKG = 1. # Tesla
    L0   = 35. # cm
    L1   = 10. # cm
    B0   = 1. # Tesla
    TRAPDEPTH = 0.004 # Tesla

    PITCHANGLE = 89. * np.pi / 180.
    
    DETPOS = np.array([0.3, 0.0]) # Position of receiver
    
    EKE = (18.6 * 10**3) * constant.COULOMBCHARGE # Electron KE in Joules
    GAMMA = (constant.ERESTMASS*constant.CLIGHT*constant.CLIGHT + EKE)/(constant.ERESTMASS*constant.CLIGHT*constant.CLIGHT)
    V0 = constant.CLIGHT * ( np.sqrt((GAMMA*GAMMA - 1)/GAMMA) )

    #1. / np.sqrt(1 - (V0/constant.CLIGHT)**2 )

    print("GAMMA =", GAMMA)
    print("BETA =", V0/constant.CLIGHT)

    l0 = L0 / 100.
    l1 = L1 / 100.
    
    THETABOT = func.CalcThetaBotMin(TRAPDEPTH, B0)
    
    ZMAX = l0 / np.tan(PITCHANGLE)
    OMEGAA = V0 * np.sin(PITCHANGLE) / l0

    AXFREQ = OMEGAA * ( 1 + l1 * np.tan(PITCHANGLE)/(np.pi*l0) )**-1
    
    t1 = l1 / (V0 * np.cos(PITCHANGLE))                                                               
    t2 = t1 + np.pi / OMEGAA
    t3 = t1 + t2
    T  = 2 * t2

    print('Cyclotron frequency =', func.CalcAngCyclotronFreq(BBKG, EKE) )
    print('omegaA =', OMEGAA)
    print('Axial frequency =', AXFREQ)
    print('Minimum pitch angle at bottom of trap=', THETABOT * 180. / np.pi)

    nSteps = 20000
    maxTime = 5*T
    timeStepSize = maxTime / nSteps
    
    TimeArray   = np.zeros(nSteps)
    BFieldArray = np.zeros(nSteps)
    AngCyclFreqArray = np.zeros(nSteps)
    AngCyclFreqArrayDoppler = np.zeros(nSteps)
    
    ZVelArray = np.zeros(nSteps)
    ZPosArray = np.zeros(nSteps)
    
    time = 1 * 10**-10

    maxAngFreq = -1 * sys.float_info.max
    minAngFreq = sys.float_info.max
    
    for i in range(nSteps):
        TimeArray[i] = time
        BFieldArray[i] = func.BFieldBathtubFromTime(time, V0, L0, L1, B0, BBKG, TRAPDEPTH)
        AngCyclFreqArray[i] = func.AngCyclFreqBathtubFromTime(time, V0, L0, L1, B0, BBKG, TRAPDEPTH)

        ZPosArray[i], ZVelArray[i] = func.ZPosVelBathtubFromTime(time, V0, L0, L1, B0, BBKG, TRAPDEPTH)

        pos = np.array([0., ZPosArray[i]])
        vel = np.array([0., ZVelArray[i]])
        electronDetVec = np.subtract(DETPOS, pos)
        electronDetVec_hat = electronDetVec / np.linalg.norm(electronDetVec)
        # Velocity in direction of receiver
        velReceiver = np.dot(vel, electronDetVec_hat)
        AngCyclFreqArrayDoppler[i] = func.ReceiverFreqDoppler(velReceiver, AngCyclFreqArray[i], 1.5)

        if AngCyclFreqArrayDoppler[i] > maxAngFreq :
            maxAngFreq = AngCyclFreqArrayDoppler[i]

        if AngCyclFreqArrayDoppler[i] < minAngFreq :
            minAngFreq = AngCyclFreqArrayDoppler[i]
        
        time = time + timeStepSize

    deltaOmega = maxAngFreq - minAngFreq
    print('Modulation index, h = ', deltaOmega/AXFREQ)
    
        
    # Now make the plots
    fig0, ax0 = pyplot.subplots(nrows=2, ncols=1, figsize=[6, 11])
    f0graph0 = ax0[0].plot(TimeArray, BFieldArray)
    f0graph1 = ax0[1].plot(TimeArray, AngCyclFreqArray)
    ax0[0].set_title('Magnetic field experienced')
    ax0[0].set_xlabel('t [s]')
    ax0[0].set_ylabel('B [T]')
    ax0[1].set_title('Magnetic field experienced')
    ax0[1].set_xlabel('t [s]')
    ax0[1].set_ylabel(r'$\Omega_{c}$ [radians $s^{-1}$]')

    fig1, ax1 = pyplot.subplots(nrows=2, ncols=1, figsize=[6, 11])
    f1graph0 = ax1[0].plot(TimeArray, ZPosArray)
    f1graph1 = ax1[1].plot(TimeArray, ZVelArray)
    ax1[0].set_title('Z position')
    ax1[0].set_xlabel('t [s]')
    ax1[0].set_ylabel('z [m]')
    ax1[1].set_title('Z velocity')
    ax1[1].set_xlabel('t [s]')
    ax1[1].set_ylabel('$v_{z}$ [m$s^{-1}$]')

    fig2, ax2 = pyplot.subplots(nrows=2, ncols=1, figsize=[6, 11])
    f2graph0 = ax2[0].hist(AngCyclFreqArray, bins=200)
    ax2[0].set_title('Angular frequency - no doppler effect')
    ax2[0].set_xlabel(r'$\Omega_{c}$ [radians $s^{-1}$]')
    ax2[0].set_ylabel('A. U.')
    #ax2[0].set_xlim([1.694e11, 1.696e11])
    f2graph0 = ax2[1].hist(AngCyclFreqArrayDoppler, bins=200)
    ax2[1].set_title('Angular frequency - doppler effect')
    ax2[1].set_xlabel(r'$\Omega_{c}$ [radians $s^{-1}$]') 
    ax2[1].set_ylabel('A. U.')
    #ax2[1].set_yscale('log')
    #ax2[1].set_xlim([1.694e11, 1.696e11])
    
    pyplot.show()
    
# Run the code
if __name__ == "__main__":
    RunSingleElectronBathtubSim()
