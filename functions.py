import numpy
import constant

# Function definitions

def CalcAngCyclotronFreq(BField, KE):
    # Angular cyclotron frequency from magnetic field and electron KE    
    freq = (constant.COULOMBCHARGE * BField)/(constant.ERESTMASS + KE/constant.CLIGHT)
    return freq

def CalcCylcotronFreq(BField, KE):
    # Angular cyclotron frequency from magnetic field and electron KE
    freq = (constant.COULOMBCHARGE * BField)/((constant.ERESTMASS + KE/(constant.CLIGHT*constant.CLIGHT)) * 2 * constant.PI)
    return freq


