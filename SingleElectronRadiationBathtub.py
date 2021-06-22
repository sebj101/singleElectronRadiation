#!/usr/bin/python3
# SingleElectronRadiationBathtub.py
# Motion of an electron in bathtub magnetic trap

import constant
import functions as func
import numpy as np
import matplotlib.pyplot as pyplot

# Main functions
def RunSingleElectronBathtubSim():
    
    # Field constants for harmonic trap field
    BBKG = 1. # Tesla
    L0   = 10. # cm
    L1   = 15. # cm
    B0   = 1. # Tesla
    TRAPDEPTH = 0.004 # Tesla

    DETPOS = np.array([5., 5.]) # Position of receiver
    
    EKE = (18.6 * 10**3) * constant.COULOMBCHARGE # Electron KE in Joules
    V0 = 0.263 * constant.CLIGHT
    GAMMA = 1. / np.sqrt(1 - (V0/constant.CLIGHT)**2 )

    print("GAMMA = ", GAMMA)


    l0 = L0 / 100.
    l1 = L1 / 100.
    thetaBot = np.arcsin(np.sqrt(1 - TRAPDEPTH/BBKG))
    zMax = l0 / np.tan(thetaBot)
    omegaA = V0 * np.sin(thetaBot) / l0

    t1 = l1 / (V0 * np.cos(thetaBot))                                                               
    t2 = t1 + np.pi / omegaA
    t3 = t1 + t2
    T  = 2 * t2
    
    nSteps = 10000
    maxTime = T*5
    timeStepSize = maxTime / nSteps

    print("Cyclotron frequency = ", func.CalcCyclotronFreq(BBKG, EKE) )
    
    TimeArray    = np.zeros(nSteps)
    BFieldArray  = np.zeros(nSteps)

    time = 1 * 10**-9

    
    for i in range(nSteps):
        TimeArray[i] = time
        BFieldArray[i] = func.BFieldBathtubFromTime(time, V0, L0, L1, B0, BBKG, TRAPDEPTH)

        time = time + timeStepSize
        
    # Now make the plots
    fig0, ax0 = pyplot.subplots(nrows=2, ncols=1, figsize=[6, 11])
    f0graph0 = ax0[0].plot(TimeArray, BFieldArray)

    ax0[0].set_title('Magnetic field experienced')
    ax0[0].set_xlabel('t [s]')
    ax0[0].set_ylabel('B [T]')
    
    pyplot.show()
    
# Run the code
if __name__ == "__main__":
    RunSingleElectronBathtubSim()
