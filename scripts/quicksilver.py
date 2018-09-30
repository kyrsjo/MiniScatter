#!/usr/bin/env python3

# Quicksilver is a little package for beam dynamics computations;
# Twiss and track.
# Its name is a corruption of "MAD"->"MadHatter"->Quicksilver
# The goal is to have an easy to use collection of typical mathematical
# transformations done when computing charged particle dynamics.

import numpy as np

#Convert between parameters

#PLASMA LENSES
def computeGradient(current,radius):
    "Utility function to convert current [A] and radius [mm] to gradient [T/m]"
    mu0 = 4*np.pi*1e-7 # [Tm/A]
    gradient = mu0*current / (2*np.pi*(radius*1e-3)**2) # [T/m]
    return gradient

def computeCurrent(gradient,radius):
    "Utility function to convert gradient [T/m] and radius [mm] to current [A]"
    mu0 = 4*np.pi*1e-7 # [Tm/A]
    current = gradient * (2*np.pi*(radius*1e-3)**2) / mu0
    return current

def computePlasmaK(current, radius, P0_MeV):
    SI_c = 299792458 #[m/s]
    gradient = computeGradient(current,radius)
    return gradient*SI_c/(P0_MeV*1e6)
    #return 299.0*gradient/energy_MeV

#QUADRUPOLE MAGNETS
def computeClearK(current, P0_MeV):
    "Focusing strength of CLEAR quadrupoles in 1/m**2 given current [A]"
    EC_QD = 0.056;
    SI_c = 299792458 #[m/s]
    F_QD = SI_c * 1e-6 * EC_QD;
    return current * F_QD / P0_MeV;


#Functions defining the effect of various elements

def driftMatrix(L):
    "Produces a drift matrix given the element length L [m]"
    return np.asarray([[1.0, L],\
                       [0.0, 1.0]])
    return matrix

def focusThickMatrix(k,L):
    "Produces a thick focusing matrix given the magnet strength k [1/m**2] and length L [m]"
    if abs(k) < 1e-10:
        return driftMatrix(L)
    return np.asarray([[            np.cos(L*np.sqrt(k)), np.sin(L*np.sqrt(k))/np.sqrt(k)],\
                       [-np.sqrt(k)*np.sin(L*np.sqrt(k)), np.cos(L*np.sqrt(k))        ]])

def twissTransformMatrix(tm):
    """
    Convert a transfer matrix to a twiss transforming matrix.
    The twiss vector should be on the form (beta [m], alpha, gamma)
    """
    return np.asarray([[tm[0,0]**2      , -2.0*tm[0,0]*tm[0,1]           , tm[0,1]**2      ],\
                       [-tm[0,0]*tm[1,0], tm[0,0]*tm[1,1]+tm[0,1]*tm[1,0], -tm[0,1]*tm[1,1]],\
                       [tm[1,0]         , -2.0*tm[1,0]*tm[1,1]           , tm[1,1]**2]])

def getGamma(alpha,beta):
    return(1+alpha**2)/beta

def get_eps_g(eps_n,E0):
    """
    Compute the geometrical emittance [um given the normalized emittance [um]
    assuming electron beam and given the reference energy in MeV.
    """
    m0 = 511e-3 # [MeV/c**2]
    gamma_rel = E0/m0
    P0 = np.sqrt((E0-m0)*(E0+m0))
    beta_rel = np.sqrt((1.0-1.0/np.sqrt(gamma_rel))*(1.0+1.0/np.sqrt(gamma_rel)))

    return eps_n/(beta_rel*gamma_rel) #[um]

#Function for generating a beam distribution

def makeDist(COVAR, N, E0):
    """
    Generate a beam distribution given the optics parameters
    in the form expected by miniscatter; COVAR = (eps_n [um], beta[m], alpha) and E0 [MeV]
    """
    eps_g = get_eps_g(COVAR[0],E0)
    covar = np.asarray([[COVAR[1], -COVAR[2]                ],\
                        [-COVAR[2], (1+COVAR[2]**2)/COVAR[1]]])
    covar *= eps_g*1e-6

    return np.random.multivariate_normal([0.0,0.0],covar,N).T
