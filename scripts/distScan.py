#!/usr/bin/env python3

## Script to scan detector distance and plot detector-target distance vs. RMS for one given target

#THICK  = 5.00
THICK  = 0.010
PHYS   = "QGSP_FTFP_BERT__SS"
PHYS   = "QGSP_FTFP_BERT"
N      = 50000
ENERGY = 200

import subprocess

def runScatter(dist,seed):
    print ("running dist =",dist,"[mm], n =",N)
    cmd = ["./MiniScatter", "-t", str(THICK), "-p", PHYS, "-n", str(N), "-e", str(ENERGY), "-d", str(dist), "-s", str(seed)]
    cmdline = ""
    for c in cmd:
        cmdline += c + " "
    print ("Command line:", cmdline)
    runResults = subprocess.run(cmd, stdout=subprocess.PIPE )

    #print(runResults.stdout)
    cutoff=-1
    xave=xrms=yave=yrms=None
    for line in runResults.stdout.split(b'\n'):
        ls = line.decode()
        if ls.startswith("Above cutoff"):
            cutoff=0
            continue
        if cutoff>=0:
            cutoff+=1
            print(ls)
            if cutoff == 1:
                lss = ls.split()
                xave = float(lss[3])
                xrms = float(lss[7])
            if cutoff == 2:
                lss = ls.split()
                yave = float(lss[3])
                yrms = float(lss[7])
                break;
            if cutoff == 3:
                break
    print ()
    return (xave,xrms, yave,yrms)

import numpy as np
import matplotlib.pyplot as plt


D = np.linspace(0.01,20,10)
XRMS = np.zeros_like(D)
YRMS = np.zeros_like(D)

for i in range(len(D)):
    d = D[i]

    (xave,xrms, yave,yrms) = runScatter(d,i+1)
    XRMS[i] = xrms
    YRMS[i] = yrms

plt.plot(D,XRMS*1e3,'-+',label="Y")
plt.plot(D,YRMS*1e3,'-+',label="X")
plt.legend()
plt.xlabel("Distance [mm]")
plt.ylabel("RMS [um]")
plt.title("Scattering from t = "+str(THICK*1e3)+ " um aluminium,\n"+str(ENERGY) + " MeV electrons, "+PHYS)
plt.show()
