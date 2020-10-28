# Source this script to setup an LCG environment on lxplus7.cern.ch
# capable of building and running MiniScatter
# K. Sjobak, 11/02-2019

export LCGENV_PATH=/cvmfs/sft.cern.ch/lcg/releases

LCGENV=/cvmfs/sft.cern.ch/lcg/releases/LCG_98python3/lcgenv/1.3.14/x86_64-centos7-gcc10-opt/lcgenv
LCGREL=LCG_97python3
#LCGPLAT=x86_64-centos7-gcc9-opt
LCGPLAT=x86_64-centos7-clang8-opt

#Put all the commands in a script `setupLCG.sh` which is then ran,
# since running the commands breaks the Python2 environment which LCGenv is running on (AT LEAST IN LCG94)
echo "" > setupLCG_runme.sh

#$LCGENV -p $LCGREL $LCGPLAT gcc    >> setupLCG_runme.sh

echo "# CMAKE:" >> setupLCG_runme.sh
$LCGENV -p $LCGREL $LCGPLAT CMake   >> setupLCG_runme.sh
echo "" >> setupLCG_runme.sh
echo "" >> setupLCG_runme.sh

echo "# CMAKE:" >> setupLCG_runme.sh
$LCGENV -p $LCGREL $LCGPLAT hdf5    >> setupLCG_runme.sh
echo "" >> setupLCG_runme.sh
echo "" >> setupLCG_runme.sh

echo "# GEANT4:" >> setupLCG_runme.sh
$LCGENV -p $LCGREL $LCGPLAT XercesC >> setupLCG_runme.sh
$LCGENV -p $LCGREL $LCGPLAT Geant4  >> setupLCG_runme.sh
echo "source geant4.sh"             >> setupLCG_runme.sh # From PATH
echo "" >> setupLCG_runme.sh
echo "" >> setupLCG_runme.sh

echo "# ROOT:" >> setupLCG_runme.sh
$LCGENV -p $LCGREL $LCGPLAT ROOT    >> setupLCG_runme.sh
echo "" >> setupLCG_runme.sh
echo "" >> setupLCG_runme.sh

#RUN IT
source setupLCG_runme.sh

#Cleanup; comment out for debugging
rm setupLCG_runme.sh
