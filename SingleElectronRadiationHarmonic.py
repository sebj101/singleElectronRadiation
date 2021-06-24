#!/usr/bin/python3
# SingleElectronRadiation.py
# Single electron in a magnetic field

import constant
import functions as func
import numpy as np
import matplotlib.pyplot as pyplot

def RunSingleElectronSim():
    # Field constants for harmonic trap field
    BBKG = 1.004 # Tesla
    L0   = 0.2 # m
    B0   = 1. # Tesla
    TRAPDEPTH = 0.004 # Tesla
    
    DETPOS = np.array([.1, 0.]) # Position of receiver

    # Basic electron kinematics
    EKE = (30. * 10**3) * constant.COULOMBCHARGE # Electron KE in Joules
    GAMMA = (constant.ERESTMASS*constant.CLIGHT*constant.CLIGHT + EKE)/(constant.ERESTMASS*constant.CLIGHT*constant.CLIGHT)
    V0 = constant.CLIGHT * ( np.sqrt((GAMMA*GAMMA - 1)/GAMMA) )

    THETABOTMIN = func.CalcThetaBotMin(TRAPDEPTH, B0)
    PITCHANGLE = 89.43 * np.pi / 180.
    AXIALFREQUENCY = func.AxialFrequencyHarmonic(V0, L0, PITCHANGLE)
    ZMAX = L0/(np.tan(PITCHANGLE))
    
    print("GAMMA =", GAMMA)
    print("Cyclotron frequency =", func.CalcAngCyclotronFreq(B0, EKE) )
    print('Minimum pitch angle at trap bottom =', THETABOTMIN * 180. / np.pi)
    print('Axial frequency =', AXIALFREQUENCY)
    print('z_max =', ZMAX, 'm')
    
    
    # How many time steps do we want?
    nSteps = 20000
    maxTime =  1 * 10**-7 
    timeStepSize = maxTime / nSteps
    time = 1 * 10**-10
    
    TimeArray    = np.zeros(nSteps)
    BFieldArray  = np.zeros(nSteps)
    ZPosArray    = np.zeros(nSteps)
    ZVelArray    = np.zeros(nSteps)
    XPosArray    = np.zeros(nSteps)
    XVelArray    = np.zeros(nSteps)
    AngFreqArray = np.zeros(nSteps)
    CyclFreqArray = np.zeros(nSteps)
    ReceiverFreqArrayNoDoppler = np.zeros(nSteps)
    ReceiverFreqArrayDoppler   = np.zeros(nSteps)
    
    for i in range(nSteps):
        TimeArray[i] = time

        pos = func.XZPositionHarmonic(time, EKE, V0, L0, BBKG, PITCHANGLE)
        vel = func.XZVelocityHarmonic(time, EKE, V0, L0, BBKG, PITCHANGLE)
        pos[0] = 0.
        vel[0] = 0.
        XPosArray[i] = pos[0]
        XVelArray[i] = vel[0]
        ZPosArray[i] = pos[1]
        ZVelArray[i] = vel[1]
        
        BFieldArray[i]   = func.BFieldHarmonicFromTime(time, V0, L0, B0, PITCHANGLE)
        AngFreqArray[i]  = func.AngCyclFreqHarmonicFromTime(time, EKE, V0, L0, B0, BBKG, PITCHANGLE)
        CyclFreqArray[i] = func.AngCyclFreqHarmonicFromTime(time, EKE, V0, L0, B0, BBKG, PITCHANGLE)/(2. * np.pi)
        
        # Assume there is detector above midpoint of CRES region
        # say d cm away from field axis
        # what frequencies does it see over time?
        
        # Calculate electron velocity in direction of detector
        electronDetVec = np.subtract(DETPOS, pos)
        electronDetVec_hat = electronDetVec / np.linalg.norm(electronDetVec)

        # Velocity in direction of receiver
        velReceiver = np.dot(vel, electronDetVec_hat)
        vel_hat = vel / np.linalg.norm(vel)
        # Angle between source travel direction and receiver
        sourceReceiverAng = np.arccos( np.dot(vel_hat, electronDetVec_hat) )
        
        ReceiverFreqArrayNoDoppler[i] = func.ReceiverFreqDoppler(0, AngFreqArray[i], sourceReceiverAng)
        ReceiverFreqArrayDoppler[i] = func.ReceiverFreqDoppler(ZVelArray[i], AngFreqArray[i], sourceReceiverAng)
        #ReceiverFreqArrayDoppler[i] = func.ReceiverFreqDoppler(velReceiver, angFreqArray[i])
        
        time = time + timeStepSize

    # Now plot things
    fig1, ax1 = pyplot.subplots(nrows=2, ncols=2, figsize=[16, 11])
    f1graph0 = ax1[0][0].plot(TimeArray, ZPosArray)
    f1graph1 = ax1[0][1].plot(TimeArray, ZVelArray)
    f1graph2 = ax1[1][0].plot(TimeArray, XPosArray)
    f1graph3 = ax1[1][1].plot(TimeArray, XVelArray)
    ax1[0][0].set_title("Z position")
    ax1[0][0].set_xlabel('t [s]')
    ax1[0][0].set_ylabel('z [m]')
    ax1[0][1].set_title("Electron axial velocity")
    ax1[0][1].set_xlabel('t [s]')
    ax1[0][1].set_ylabel(r'$v_{z}$ [m$s^{-1}$]')
    ax1[1][0].set_title("X position")
    ax1[1][0].set_xlabel('t [s]')
    ax1[1][0].set_ylabel('x [m]')
    ax1[1][1].set_title("Electron x velocity")
    ax1[1][1].set_xlabel('t [s]')
    ax1[1][1].set_ylabel(r'$v_{x}$ [m$s^{-1}$]')

    fig2, ax2 = pyplot.subplots(nrows=3, ncols=1, figsize=[8, 15])    
    f2graph0 = ax2[0].plot(TimeArray, BFieldArray)
    f2graph1 = ax2[1].plot(TimeArray, AngFreqArray)
    f2graph2 = ax2[2].plot(TimeArray, CyclFreqArray)
    ax2[0].set_title("B field experienced by electron")
    ax2[0].set_xlabel('t [s]')
    ax2[0].set_ylabel('B [T]')    
    ax2[1].set_title("Electron angular frequency")
    ax2[1].set_xlabel('t [s]')
    ax2[1].set_ylabel(r'$\Omega_{c}$ [$s^{-1}$]')
    ax2[2].set_title("Electron frequency")
    ax2[2].set_xlabel('t [s]')
    ax2[2].set_ylabel(r'$\Omega_{c}$ [$s^{-1}$]')
    
    fig3, ax3 = pyplot.subplots(nrows=1, ncols=1)
    f3graph0 = ax3.hist(ReceiverFreqArrayNoDoppler, bins=100)
    ax3.set_title('Angular frequency - no doppler effect')
    ax3.set_xlabel('Frequency [Hz]')
    ax3.set_ylabel('A. U.')

    fig4, ax4 = pyplot.subplots(nrows=1, ncols=1)
    f4graph0 = ax4.hist(ReceiverFreqArrayDoppler, bins=100)
    ax4.set_title('Angular frequency - with doppler effect')
    ax4.set_xlabel('Frequency [Hz]')
    ax4.set_ylabel('A. U.')
    
    pyplot.show()

# Run the code
if __name__ == "__main__":
    RunSingleElectronSim()
