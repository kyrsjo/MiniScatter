# MiniScatter

Simple Geant4+Root program for simulating the distribution of the scattered particles after passing through a block of copper.

To build: In a separate directory, run:
> cmake -DGeant4_dir=/PATH/TO/geant4.10.01.p02-install/lib64/Geant4-10.1.2/ PATH/TO/THIS/FOLDER

> make

Finally, create a sub folder "plots", where the root file will be stored. Then run some events with
> /run/beamOn NUMBER
