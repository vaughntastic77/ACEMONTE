# WINDMONTE is licensed under GNU GPL v3, see COPYING.txt

import re
import numpy as np
import copy

""" 
This is where the test specific data reduction equations (DREs) should be included.  Replace everything in the eval() 
function with the appropriate DREs to map test inputs/measurands (data, testinfo) to test results/VOIs.  eval() should accept 
data and testinfo as inputs and return the test results in the same format as the 'data' variable.  Additional functions used 
inside eval() may be defined in this file as appropriate.
"""

# Air constants
GAM = 1.4
R = 287.052874  # J/kg K
MU0 = 1.716*10**(-5)  # Sutherland's constants
TMU = 273.15  # Sutherland's constants
SMU = 110.4  # Sutherland's constants

# Conversions
PSI2PA = 6894.7572931783
PSI2TORR = 51.7149325716
TORR2PA = 133.3223684211
IN2M = 0.0254

def eval(data,G):

    D = copy.deepcopy(data)

    for i in range(len(D)):
        # Wall static pressure
        D[i]['M(P)'] = np.sqrt((2/(GAM-1))*(((D[i]['P01']*PSI2PA)/D[i]['P'])**((GAM-1)/GAM) - 1))
        
        D[i]['T(P)'] = D[i]['T0']/(1 + ((GAM-1)/2)*(D[i]['M(P)']**2))

        D[i]['U(P)'] = D[i]['M(P)']*np.sqrt(GAM*R*D[i]['T(P)'])

        D[i]['Rho(P)'] = D[i]['P']/(R*D[i]['T(P)'])

        D[i]['Mu(P)'] = MU0*((TMU+SMU)/(D[i]['T(P)']+SMU))*(D[i]['T(P)']/TMU)**(3.0/2.0)

        D[i]['Re(P)'] = D[i]['Rho(P)']*D[i]['U(P)']/D[i]['Mu(P)']

        # Pitot stagnation pressure 
        Pratio = D[i]['P02']/(D[i]['P01']*PSI2TORR)
        Mach = M_P02P01(Pratio)
        D[i]['M(P02)'] = Mach

        D[i]['T(P02)'] = D[i]['T0']/(1 + ((GAM-1)/2)*(D[i]['M(P02)']**2))

        D[i]['U(P02)'] = D[i]['M(P02)']*np.sqrt(GAM*R*D[i]['T(P02)'])

        P02P1ratio = ((1-GAM+2*GAM*(Mach**2))/(GAM+1))*((((GAM+1)**2)*(Mach**2))/(4*GAM*(Mach**2)-2*(GAM-1)))**(GAM/(GAM-1))
        D[i]['P(P02)'] = D[i]['P02']*TORR2PA/P02P1ratio

        D[i]['Rho(P02)'] = D[i]['P(P02)']/(R*D[i]['T(P02)'])

        D[i]['Mu(P02)'] = MU0*((TMU+SMU)/(D[i]['T(P02)']+SMU))*(D[i]['T(P02)']/TMU)**(3.0/2.0)

        D[i]['Re(P02)'] = D[i]['Rho(P02)']*D[i]['U(P02)']/D[i]['Mu(P02)']

    return D

def Pr(M):
    PRout = ((GAM+1)*M**2/((GAM-1)*M**2+2))**(GAM/(GAM-1))*((GAM+1)/(2*GAM*M**2 - (GAM-1)))**(1/(GAM-1))
    return PRout

def M_P02P01(Pratio):
    Mach = 0
    if Pratio <= 1:
        err = 1
        Ma = 1  # initial left-hand value
        Mb = 10  # initial right-hand value
        tol =  1e-5

        Pra = Pr(Ma) - Pratio
        Prb = Pr(Mb) - Pratio

        while abs(err) > tol:
            Mach = (Ma + Mb)/2
            Prc = Pr(Mach) - Pratio
            if np.sign(Prc) == np.sign(Pra):
                err = Mach - Ma
                Ma = Mach
            elif np.sign(Prc) == np.sign(Prb):
                err = Mach - Mb
                Mb = Mach
    return Mach
