# MiniScatter's command line interface

The basic user interface of MiniScatter is through the command line.
This allows the user to specify the simulation parameters such as target thickness, target material, and incident beam properties as command line flags.
An altertative is to use the [Python interface](PyInterface.m), this interface is a wrapper which behind the scences is calling MiniScatter through the command line interface.
It is therefore useful to have some understanding of the command line interface even if one intends to use the Python interface.

## Quick example
As an example, in order to run a simulation with 1000 incident electrons at 200 MeV in a pencil-beam configuration (i.e. all particles travelling in parallel and starting at the same point) with a 20 mm thick water target, one could use the command:
```
./MiniScatter -n 1000 -e 200 -b e- -t 20 -m G4_WATER
```
This initializes Geant4, runs the simulation, and provides various output on the terminal.
Furthermore, additional output is given through the file `plots/output.root`, which can be inspected with e.g. `TBrowser`:
```
$ cd plots
$ root output.root
<ROOT loading>
[1] new TBrowser();
<TBrowser window pops up, allowing you to explore the file>
[2].q
<ROOT exits>
```

## Built-in documentation
As the list of available options is changing often, a built-in documentation is maintained instead of a separate manual.
This is accessed by running MiniScatter with the `-h` flag, as shown below.

```
$ ./MiniScatter -h
Welcome to MiniScatter!

Usage/options:
-t <double> : Target thickness [mm],  default/current value = 1
-m <string> : Target material name,   default/current       = 'G4_Al'
 Valid choices: 'G4_Al', 'G4_C', 'G4_Cu', 'G4_Pb', 'G4_Ti', 'G4_Si', 'G4_W', 'G4_U', 'G4_MYLAR', 'G4_KAPTON', 'G4_STAINLESS-STEEL', 'G4_WATER', 'G4_Galactic', 'Sapphire'
 Also possible: 'gas::pressure'  where 'gas' is 'H_2', 'He', 'N_2', 'Ne', or 'Ar', and pressure is given in mbar (T=300K is assumed).
-d <double> : Detector distance [mm], default/current value = 50
-a <double> : Detector angle [deg],   default/current value = 0
-w <double> : World size X/Y [mm],    default/current value = 0
-p <string> : Physics list name,      default/current       = 'QGSP_FTFP_BERT
-n <int>    : Run a given number of events automatically
-e <double> : Beam energy [MeV],      default/current value = 200
-b <string> : Particle type,          default/current value = e-
 This accepts standard Geant4 particle types (see /gun/List for all of them),
 typcial ones are 'e-' and 'proton'.
 Ions can also be specified as 'ion::Z,A' where Z and A are the nucleus charge and mass number.
-x <double> : Beam offset (x) [mm],   default/current value = 0
-z (*)<double> : Beam offset (z) [mm],   default/current value = 0, doBacktrack = false
 If set to 0.0, start at half the buffer distance. Note that target always at z=0.
 If a '*' is prepended, the distribution is to be generated at z=0,
 then backtracked to the given z value (which may be 0.0)
-c epsN[um]:beta[m]:alpha(::epsN_Y[um]:betaY[m]:alphaY) : 
 Set realistic beam distribution (on target surface); 
 if optional part given then x,y are treated separately
--beamRcut <double> : Radial cutoff for the beam distribution.
 If given alone, generate a circular uniform distribution.
 If given together with -c, generate a multivariate gaussian with all particles starting within the given radius.
 Default/current value = 0
-s <int>    : Set the initial seed,   default/current value = 123
-g : Use a GUI
-q : Quickmode, skip most post-processing and plots, default/current value = false
-r : miniROOTfile, write small root file with only anlysis output, no TTrees, default/current value = false
-f <string> : Output filename,        default/current value = output
-o <string : Output folder,           default/current value = plots
--cutoffEnergyfraction : Minimum of beam energy to require for 'cutoff' plots, default/current value = 0.95
--cutoffRadius         : Maximum radius on target to require for 'cutoff' plots, default/current value = 1 [mm]
--edepDZ               : Z bin width for energy deposit histograms default/current value = 0 [mm]
--magnet (*)pos:type:length:gradient(:type=val1:specific=val2:arguments=val3) :  Create a magnet of the given type at the given position. 
 If a '*' is prepended the position (<double> [mm]), the position is the start of the active element relative to the end of the target; otherwize it is the z-position of the middle of the element.
 The gradient (<double> [T/m]) is the focusing gradient of the device.
 The length <double> [mm] is the total length of the volumes used by the device.
 The type-specific arguments are given as key=value pairs.
 Accepted types and their arguments:
  'PLASMA1':
     radius:    Capillary radius (<double> [mm])
     totalAmps: Flag (<True/False>) to interpret the gradient parameter as the total current [A] instead of in [T/m].
     width:     Capillay crystal width (<double> [mm])
     height:    Capillay crystal height (<double> [mm])
  'COLLIMATOR1':
     radius:    Channel radius (<double> [mm])
     width:     Absorber width (<double> [mm])
     height:    Absorber height (<double> [mm])
Currently have the following magnet setups:


Note that if both -g and -n is used, the events are ran before the GUI is opened.
One may also use one or more arguments which does not include a '-n' -- these are forwarded untouched to Geant4
The first argument not in the form '-char' is interpreted as a macro to run. Don't use vis.mac, it will crash.
Extra arguments are not compatible with -g
```

A further advantage of providing the description of the flags in this form is that it allows the user see the default or current value of each flags.
Please note that since the command line arguments are parsed sequentially, the values displayed by `-h` are modified by the flags that precede the `-h`, i.e. if you run `./MiniScatter -b proton -h`, the output will now state: `-b <string> : Particle type,          default/current value = proton`


## Geant4 macros
It is sometimes useful to run Geant4 macros, calling the built-in command interface.
To do this, add the name of the macro file to the end of argument list, i.e. `./MiniScatter -n 10 verbose.mac`.

## GUI
Geant4 comes with a built-in GUI, which can be accessed by running MiniScatter with the -g flag.
This allows inspecting the 3D geometry as well as the particle tracks, and run Geant4 commands such as `/run/beamOn 10`.

![Screenshot of GUI for `./MiniScatter -g -t 10`](MiniScatterGUI.png)

Note that Geant4 macros cannot be specified on the MiniScatter command line when using the GUI;
to run macros such as `gammaFilter.mac` do it from inside the started GUI before running the `beamOn` command by directly calling `/control/execute`, e.g. `/control/execute gammaFilter.mac`.
