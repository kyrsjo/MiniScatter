# MiniScatter

MiniScatter is a straight-forward program for simulating what happens a particle beam after passing through a target of some material and thickness, taking the effect of magnetic fields into account.
It is based on Geant4 [1], and simulation output is handled via ROOT [2].
The geometry is created based on user input, built up from a range of built-in and configurable elements.
No config files are used; the program is configured by the command line used to run it.

In addition to MiniScatter itself, a set of Python libraries for interacting with MiniScatter are also provided.
This makes it very convenient to use MiniScatter for running scans and extracting the output.
Furthermore, the scan library provides behind-the-scenes saving of results, meaning that if you have previously ran a simulation and then restart the Jupyter notebook, running the same python call again will instantly return the same results as when running the simulation in the first round.
This is implemented using ROOT and HDF5.

## Installing MiniScatter

*Note: If you have access to CERN computing resources, you can run MiniScatter through SWAN or LxPlus instead of installing MiniScatter on your local machine.
This is useful for testing, teaching, etc. however the performance is most likely slower.
For more information, see [SWAN](SWAN.md).*

To install MiniScatter, you must first install Geant4 and ROOT.

For Geant4, build and install as described on the framework's webpage [3].
Assuming that it is installed in e.g. `~/code/geant4/geant4.10.04.p02-install`, then load Geant4 in the current shell by running `source ~/code/geant4/geant4.10.04.p02-install/bin/geant4.sh`.
If you intend to use MiniScatter with Jupyter [4], do this in the shell you plan to start Jupyter before launching `jupyter-notebook`.
You must also remember to load Geant4 before compiling and running MiniScatter.

For ROOT, simply install it using your distribution's package manager, e.g. `dnf` on Fedora.

Then, clone MiniScatter into a new folder, e.g. from ~/code run `git clone https://github.com/kyrsjo/MiniScatter.git`.
Enter `~/code/MiniScatter`, and create a build folder `~/code/MiniScatter/build`.
In this folder, first configure the build system using CMAKE [5] by running `cmake ../.`, and then compile using `make -j N` where N is the number of CPUs to use; if for example you are running on a 16-core machine use `make -j 16`.
If all goes well, you have now built the executable!

Note that you may also change various build options, for example to enable debugging symbols to be written to the executable.
To do this, run `ccmake .` in your build folder, change what you need to change in the menu, then reconfigure and regenerate the makefiles.
Finally, run `make` again.

## Running MiniScatter from the command line

[More information](CommandLineUse.md)

The simplest way to use MiniScatter is to launch simulations from the command line.
Furthermore, the `-h` option will show you a list of available options and their default values, and `-g` will open the standard Geant4 GUI which is useful to check what the geometry looks like and how typical events look like.
In addition to the "minus-options", Geant4 macros can also be used; as an example the included macro `verbose.mac` makes Geant4 print extra information about the particle tracks and the interactions as the simulation progresses.
These macros can be ran by specifying them at the end of the command line, for example as `./MiniScatter -n 10 -- verbose.mac`.

Note that some options, like `--magnet` which creates a magnet or collimator with a given set of parameters, can be specified more than once.

## Running MiniScatter via Python

[More information](PyInterface.md)

Two Python libraries for running MiniScatter are provided:
 * `miniScatterDriver.py`  : Runs a single miniScatter simulation, or extracts and returns the results from a .root file written by MiniScatter.
 To see the available options, please read the sources (i.e. the definitions of the functions).
 In general, it is a quite thin "shim" over the command line interface, so see `./MiniScatter -h` for a description of the options.
 * `miniScatterScanner.py` : Runs a scan over several miniScatter simulations (using `miniScatterDriver`).
 Returns a standard set of observables as numpy arrays as a function of the scanned variable, as well as a requested set of ROOT histograms from each scan point.
 This code has the capability to run the simulations in parallel, in which case it builds a job queue with a given number of CPUs which then "eats" the queue untill it is empty.
 For a full list of available options, please see the source.

For examples of using the Python interface as well as running the code via Jupyter examples, please see in the folder `examples`.

## References

[1] : https://geant4.web.cern.ch/

[2] : https://root.cern.ch/

[3] : http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/InstallationGuide/html/

[4] : http://jupyter.org/

[5] : https://cmake.org/
