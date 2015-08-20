# MiniScatter

Simple Geant4+Root program for simulating the distribution of the scattered particles after passing through a block of copper.

To build: In a separate directory, run:
> cmake -DGeant4_dir=/PATH/TO/geant4.10.01.p02-install/lib64/Geant4-10.1.2/ PATH/TO/THIS/FOLDER

If wanted, the flag 
> -DCMAKE_BUILD_TYPE=Debug
can be added to cmake in order to make debugging possible.

Then compile using:
> make

Finally, create a sub folder "plots", where the root file will be stored.

The program execution is normally controlled via command line arguments.
These set the target thickness, the physics list, the number of events to generate, and wether to display a GUI after running the specified number of events.
To see the available options, run:
> ./MiniScatter -h
To generate 100000 events, using physics list QGSP_BERT_HP, on a target 1.0 mm thick, run:
> ./MiniScatter -t 1.0 -p QGSP_BERT_HP -n 100000
