# Source this script to setup an LCG environment on lxplus7.cern.ch
# capable of building and running MiniScatter
# K. Sjobak, 11/02-2019

export LCGENV_PATH=/cvmfs/sft.cern.ch/lcg/releases

#Put all the commands in a script `setupLCG.sh` which is then ran,
# since running the commands breaks the Python2 environment which LCGenv is running on
/cvmfs/sft.cern.ch/lcg/releases/LCG_94python3/lcgenv/1.3.6/x86_64-centos7-gcc7-opt/lcgenv -p LCG_94python3 x86_64-centos7-gcc7-opt gcc     >  setupLCG_runme.sh

/cvmfs/sft.cern.ch/lcg/releases/LCG_94python3/lcgenv/1.3.6/x86_64-centos7-gcc7-opt/lcgenv -p LCG_94python3 x86_64-centos7-gcc7-opt CMake   >> setupLCG_runme.sh

/cvmfs/sft.cern.ch/lcg/releases/LCG_94python3/lcgenv/1.3.6/x86_64-centos7-gcc7-opt/lcgenv -p LCG_94python3 x86_64-centos7-gcc7-opt hdf5    >> setupLCG_runme.sh

/cvmfs/sft.cern.ch/lcg/releases/LCG_94python3/lcgenv/1.3.6/x86_64-centos7-gcc7-opt/lcgenv -p LCG_94python3 x86_64-centos7-gcc7-opt Geant4  >> setupLCG_runme.sh
echo "source geant4.sh"                                                                                                                    >> setupLCG_runme.sh 

/cvmfs/sft.cern.ch/lcg/releases/LCG_94python3/lcgenv/1.3.6/x86_64-centos7-gcc7-opt/lcgenv -p LCG_94python3 x86_64-centos7-gcc7-opt ROOT    >> setupLCG_runme.sh

source setupLCG_runme.sh

#Cleanup; comment out for debugging
rm setupLCG_runme.sh
