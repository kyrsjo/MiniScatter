#!/usr/bin/env python3

## Script to scan beam parameters at focal point, computing the twiss parameters modification (incl. emittance growth)
## by a scattering foil.

import numpy as np
import matplotlib.pyplot as plt

import miniScatterDriver

def TwissDrift(beta_in:float, gamma_in:float, alpha_in:float, L:float) -> tuple:
    beta_out  = beta_in + L**2 * gamma_in - 2*L*alpha_in
    gamma_out =                  gamma_in
    alpha_out =         -    L * gamma_in     + alpha_in
    return (beta_out, gamma_out, alpha_out)

def getGamma(beta_in:float, alpha_in:float):
    return (1.0 + alpha_in**2)/beta_in

#Parameters at the capillary
import sys
if len(sys.argv) == 3 or len(sys.argv) == 5:
    if ":" in sys.argv[1]:
        ls = sys.argv[1].split(":")
        eps_0x = float(ls[0])
        eps_0y = float(ls[1])
    else:
        eps_0x = float(sys.argv[1])
        eps_0y = float(sys.argv[1])
    alpha_0 = float(sys.argv[2])
else:
    eps_0x    = 1.0                   # [um]
    eps_0y    = 1.0                   # [um]
    alpha_0  = -5.0                   # [1]

beta_0   = np.logspace(np.log10(0.0001),np.log10(1.50),10) #[m]


#System parameters
L = 0.202 #[m] Capillary-to-foil distance
Lcap = 0.015 #[m] Capillary length
L = L + Lcap/2.0 #We actually care about the middle of the capillary

if len(sys.argv) == 5:
    THICK = float(sys.argv[3])
    MAT   = sys.argv[4]
    matstr = MAT
    if MAT == "G4_MYLAR":
        matstr = "Mylar"
    elif MAT == "G4_KAPTON":
        matstr = "Kapton"
    plottitle = "{0:.1f} $\\mu$m {1}, $\\alpha_0={2:.1f}$, $\\varepsilon_0 = {3:.1f}:{4:.1f}$ [$\\mu$m]".\
                                                  format(THICK*1000, matstr, alpha_0, eps_0x,eps_0y)
else:
    plottitle = "3.0 $\\mu$m Mylar, $\\alpha_0={0:.1f}$, $\\varepsilon_0 = {1:.1f}$ [$\\mu$m]".format(alpha_0, eps_0)
    THICK = 3e-3 #[mm]
    MAT   = "G4_MYLAR"
    #MAT   = "G4_Galactic"
    #MAT   ="G4_KAPTON"

#Simulation parameters
N         = 1000000
#DIST      = THICK
DIST      = L*1000 #[mm]
ENERGY    = 210.0
ZOFFSET   = -THICK
ZOFFSET_BACKTRACK = True
PHYS      = "QGSP_BERT__SS"
QUICKMODE = True

#Assume electron beam
gamma_rel = ENERGY/0.511
beta_rel  = np.sqrt(gamma_rel**2 - 1.0) / gamma_rel;

print ("plottile  = '"+plottitle+"'")
print ("gamma_rel =", gamma_rel)
print ("beta_rel  =", beta_rel)
#Make output folder
import os
dirpath="plots/beamParamScan_alpha0={0:.1f}_eps0={1:.1f}:{2:.1f}_THICK={3}_MAT={4}_Npoints={5}".\
                                     format(alpha_0,eps_0x,eps_0y,THICK,MAT,len(beta_0))

try:
    print()
    print("Making folder '"+dirpath+"'")
    os.mkdir(dirpath)
except FileExistsError:
    print("Aready there! nuking.")
    import shutil
    shutil.rmtree(dirpath)
    os.mkdir(dirpath)
print()

#Computed parameters
eps_x   = np.empty_like(beta_0)
beta_x  = np.empty_like(beta_0)
alpha_x = np.empty_like(beta_0)
sigma_x = np.empty_like(beta_0)

eps_y   = np.empty_like(beta_0)
beta_y  = np.empty_like(beta_0)
alpha_y = np.empty_like(beta_0)
sigma_y = np.empty_like(beta_0)

sigma_0x   = np.empty_like(beta_0)
sigma_0y   = np.empty_like(beta_0)
beta_foil  = np.empty_like(beta_0)
alpha_foil = np.empty_like(beta_0)

#Plot the beam development as function of s
default_size = plt.rcParams.get('figure.figsize')
double_width = np.asarray(default_size)
double_width[0] *= 2
plt.figure(1,figsize=double_width)
plt.subplots_adjust(left=0.06,right=0.98)
plt.figure(2,figsize=double_width)
plt.subplots_adjust(left=0.06,right=0.98)

sPos = np.linspace(-1.25*L,0.25*L,200)


for i in range(len(beta_0)):
    gamma_0 = getGamma(beta_0[i],alpha_0)
    (beta_foil[i],gamma_foil,alpha_foil[i]) = TwissDrift(beta_0[i],gamma_0,alpha_0,-L)

    COVAR = (eps_0x,beta_foil[i],alpha_foil[i],eps_0y,beta_foil[i],alpha_foil[i])

    sigma_0x[i] = np.sqrt(eps_0x*beta_0[i]*1e6/(gamma_rel*beta_rel))
    sigma_0y[i] = np.sqrt(eps_0y*beta_0[i]*1e6/(gamma_rel*beta_rel))
    print("Simulation", i+1, "of",len(beta_0))
    print("Simulating for beta_0 =", beta_0[i], "[m]")
    print("Sigma_0               =", sigma_0x[i], "and", sigma_0y[i], "[um]")

    miniScatterDriver.runScatter(THICK=THICK,MAT=MAT,DIST=DIST,PHYS=PHYS,N=N,\
                                 ENERGY=ENERGY,ZOFFSET=ZOFFSET,ZOFFSET_BACKTRACK=ZOFFSET_BACKTRACK,\
                                 COVAR=COVAR, QUICKMODE=QUICKMODE, SEED=i+1)
    emittances = miniScatterDriver.getData()

    eps_x[i]   = emittances[0]
    beta_x[i]  = emittances[1]
    alpha_x[i] = emittances[2]
    sigma_x[i] = np.sqrt(eps_x[i]*beta_x[i]*1e6/(gamma_rel*beta_rel))

    eps_y[i]   = emittances[3]
    beta_y[i]  = emittances[4]
    alpha_y[i] = emittances[5]
    sigma_y[i] = np.sqrt(eps_y[i]*beta_y[i]*1e6/(gamma_rel*beta_rel))

    plt.figure(1)
    beta0_s = TwissDrift(beta_0[i],gamma_0,alpha_0,sPos)[0]
    gamma_x = getGamma(beta_x[i],alpha_x[i])
    beta_xs = TwissDrift(beta_x[i],gamma_x,alpha_x[i],sPos)[0]
    l1 = plt.plot(sPos,np.sqrt(eps_0x*beta0_s*1e6/(gamma_rel*beta_rel)), ls="--" )[0]
    plt.plot(sPos,np.sqrt(eps_x[i]*beta_xs*1e6/(gamma_rel*beta_rel)),\
             label='{0:.3g}'.format(beta_0[i]*100), ls="-",color=l1.get_color())

    plt.figure(2)
    beta0_s = TwissDrift(beta_0[i],gamma_0,alpha_0,sPos)[0]
    gamma_y = getGamma(beta_y[i],alpha_y[i])
    beta_ys = TwissDrift(beta_y[i],gamma_y,alpha_y[i],sPos)[0]
    l1 = plt.plot(sPos,np.sqrt(eps_0y*beta0_s*1e6/(gamma_rel*beta_rel)), ls="--" )[0]
    plt.plot(sPos,np.sqrt(eps_y[i]*beta_ys*1e6/(gamma_rel*beta_rel)),\
             label='{0:.3g}'.format(beta_0[i]*100), ls="-",color=l1.get_color())

    print("Twiss (X):", emittances[0], emittances[1], emittances[2])
    print("Twiss (Y):", emittances[3], emittances[4], emittances[5])
    print("Sigma (X):",sigma_x[i], "[um]")
    print("Sigma (Y):",sigma_y[i], "[um]")
    print()

### Write h5 files ###
import h5py
outputfile_name = "beamParamScan.h5"
print ("Dumping data to file '"+ outputfile_name+ "'")
outputfile = h5py.File(dirpath + "/" + outputfile_name,'w')
outputfile.create_dataset("L",data=L)
outputfile.create_dataset("Lcap",data=Lcap)
outputfile.create_dataset("THICK",data=THICK)
outputfile.create_dataset("MAT",data=MAT)
outputfile.create_dataset("N",data=N)
outputfile.create_dataset("ENERGY",data=N)
outputfile.create_dataset("PHYS",data=PHYS)

outputfile.create_dataset("eps_0x",data=eps_0x)
outputfile.create_dataset("eps_0y",data=eps_0y)
outputfile.create_dataset("alpha_0",data=alpha_0)
outputfile.create_dataset("beta_0",data=beta_0)
outputfile.create_dataset("sigma_0x",data=sigma_0x)
outputfile.create_dataset("sigma_0y",data=sigma_0y)

outputfile.create_dataset("alpha_foil",data=alpha_foil)
outputfile.create_dataset("beta_foil",data=beta_foil)

outputfile.create_dataset("eps_x",data=eps_x)
outputfile.create_dataset("beta_x",data=beta_x)
outputfile.create_dataset("alpha_x",data=alpha_x)
outputfile.create_dataset("sigma_x",data=sigma_x)

outputfile.create_dataset("eps_y",data=eps_y)
outputfile.create_dataset("beta_y",data=beta_y)
outputfile.create_dataset("alpha_y",data=alpha_y)
outputfile.create_dataset("sigma_y",data=sigma_y)

outputfile.close()
print("Done.")

### Finish plots ###
print("Plotting..")

import matplotlib.transforms as mtransforms

plt.figure(1)
plt.title(plottitle)
plt.axvline(-L, ls='-',color='blue')
trans = mtransforms.blended_transform_factory(plt.gca().transData, plt.gca().transAxes)
plt.fill_between([-Lcap/2.0, Lcap/2.0], 0, 1, facecolor='gray', alpha=0.25, transform=trans)
plt.xlabel('s [m]')
plt.ylabel('$\\sigma_x$ [$\\mu$m]')
plt.legend(loc=0,ncol=4,title="$\\beta_0$ [cm]")
plt.savefig(dirpath + "/" + "sigma-s.png")
plt.ylim(0,200)
plt.savefig(dirpath + "/" + "sigmax-s-ZOOM.png")
plt.ylim(0,100)
plt.savefig(dirpath + "/" + "sigmax-s-ZOOM2.png")

plt.figure(2)
plt.title(plottitle)
plt.axvline(-L, ls='-',color='blue')
trans = mtransforms.blended_transform_factory(plt.gca().transData, plt.gca().transAxes)
plt.fill_between([-Lcap/2.0, Lcap/2.0], 0, 1, facecolor='gray', alpha=0.25, transform=trans)
plt.xlabel('s [m]')
plt.ylabel('$\\sigma_y$ [$\\mu$m]')
plt.legend(loc=0,ncol=4,title="$\\beta_0$ [cm]")
plt.savefig(dirpath + "/" + "sigma-s.png")
plt.ylim(0,200)
plt.savefig(dirpath + "/" + "sigmay-s-ZOOM.png")
plt.ylim(0,100)
plt.savefig(dirpath + "/" + "sigmay-s-ZOOM2.png")

#Used for a lot of vertical lines
minEpsIdx_x = np.argmin(eps_x)
minEpsIdx_y = np.argmin(eps_y)

plt.figure()
plt.title(plottitle)
l_x = plt.plot(beta_0*100,eps_x, label="$\\varepsilon_x$")[0]
l_y = plt.plot(beta_0*100,eps_y, label="$\\varepsilon_y$")[0]
plt.xlabel("$\\beta_{0,capillary}$ [cm]")
plt.ylabel("$\\varepsilon_N$ [$\\mu$m]")
plt.axhline(eps_0x,ls='--', label="$\\varepsilon_{x,0}$", color=l_x.get_color())
plt.axhline(eps_0y,ls='--', label="$\\varepsilon_{y,0}$", color=l_y.get_color())
plt.legend(loc=0)
plt.axvline(beta_0[minEpsIdx_x]*100,color=l_x.get_color(),ls='--')
plt.axvline(beta_0[minEpsIdx_y]*100,color=l_y.get_color(),ls='--')
plt.savefig(dirpath + "/" + "epsN-beta0.png")
plt.ylim(0,max(eps_x[-1],eps_y[-1])*1.2)
plt.savefig(dirpath + "/" + "epsN-beta0-ZOOM.png")

plt.figure()
ax1 = plt.gca()
plt.title(plottitle)
l_bf = plt.plot(beta_0*100,beta_foil*100)[0]
plt.ylabel("$\\beta_{foil} [cm]$")
plt.xlabel("$\\beta_{0,capillary}$ [cm]")
plt.axvline(beta_0[minEpsIdx_x]*100,color=l_x.get_color(),ls='--')
plt.axvline(beta_0[minEpsIdx_y]*100,color=l_y.get_color(),ls='--')
ax2=plt.twinx()
l_af = ax2.plot(beta_0*100,alpha_foil, ls="--")[0]
ax2.set_ylabel("$\\alpha_{foil}$")
plt.subplots_adjust(right=0.82)
plt.figlegend(labels=("$\\beta$","$\\alpha$"), handles=(l_bf, l_af), loc='upper center',bbox_to_anchor=(0.5, 0.88))
plt.savefig(dirpath + "/" + "twissFoil-beta0.png")
ax1.set_ylim(0,beta_foil[-1]*100*1.5)
ax2.set_ylim(0,alpha_foil[np.argmin(beta_foil)]*1.5)
plt.savefig(dirpath + "/" + "twissFoil-beta0-ZOOM.png")


plt.figure()
plt.title(plottitle)
l_bx=plt.plot(beta_0*100,beta_x*100)[0]
l_by=plt.plot(beta_0*100,beta_y*100)[0]
plt.ylabel("$\\beta_{capillary} [cm]$")
plt.xlabel("$\\beta_{0,capillary}$ [cm]")
plt.axvline(beta_0[minEpsIdx_x]*100,color=l_bx.get_color(),ls='--')
plt.axvline(beta_0[minEpsIdx_y]*100,color=l_by.get_color(),ls='--')
ax2=plt.twinx()
l_ax=ax2.plot(beta_0*100,alpha_x, ls="--")[0]
l_ay=ax2.plot(beta_0*100,alpha_y, ls="--")[0]
ax2.set_ylabel("$\\alpha_{capillary}$")
plt.subplots_adjust(right=0.85)
plt.figlegend(labels=("$\\beta_x$","$\\beta_y$","$\\alpha_x$","$\\alpha_y$"),\
              handles=(l_bx,l_by, l_ax,l_ay),\
              loc='upper center',bbox_to_anchor=(0.5, 0.88),\
              ncol=2)
plt.savefig(dirpath + "/" + "twissCap-beta0.png")


plt.figure()
plt.title(plottitle)
l_x = plt.plot(beta_0*100,sigma_x, label="$\\sigma_x$")[0]
l_y = plt.plot(beta_0*100,sigma_y, label="$\\sigma_y$")[0]
plt.axvline(beta_0[minEpsIdx_x]*100,color=l_x.get_color(),ls='--')
plt.axvline(beta_0[minEpsIdx_y]*100,color=l_y.get_color(),ls='--')
plt.ylabel("$\\sigma_{capillary}$ [$\\mu$m]")
plt.xlabel("$\\beta_{0,capillary}$ [cm]")
ylim_old = plt.ylim()
plt.ylim(0.0,ylim_old[1])
plt.legend()
plt.savefig(dirpath + "/" + "sigmaCap-beta0.png")


plt.figure()
plt.title(plottitle)
l_x = plt.plot(sigma_0x, sigma_x, label="$\\sigma_x$")[0]
l_y = plt.plot(sigma_0y, sigma_y, label="$\\sigma_y$")[0]
plt.axvline(sigma_0x[minEpsIdx_x],color=l_x.get_color(),ls='--')
plt.axvline(sigma_0y[minEpsIdx_y],color=l_y.get_color(),ls='--')
plt.ylabel("$\\sigma_{capillary}$ [$\\mu$m]")
plt.xlabel("$\\sigma_{0,capillary}$ [$\\mu$m]")
ylim_old = plt.ylim()
plt.ylim(0.0,ylim_old[1])
xlim_old = plt.xlim()
plt.xlim(0.0,xlim_old[1])
plt.plot([0.0,xlim_old[1]],[0.0,ylim_old[1]],color='green',ls='--', label="$\\sigma_{capillary} = \\sigma_{0,capillary}$")
plt.legend()
plt.savefig(dirpath + "/" + "sigmaCap-sigmaCap0.png")

#plt.show()
