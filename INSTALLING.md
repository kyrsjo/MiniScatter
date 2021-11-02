## Installing MiniScatter

To install MiniScatter, you must first install Geant4 and ROOT.
Both should be compiled with support for C++17.
Then MiniScatter itself can be compiled from source.

### Geant4
For Geant4, build and install as described on the framework's webpage [4].
A tested set of CMAKE flags for `Geant4.10.07.p02` is:
* `GEANT4_BUILD_CXXSTD=17`
* `GEANT4_BUILD_MULTITHREADED=ON`
* `GEANT4_INSTALL_DATA=ON`
* `GEANT4_USE_GDML=ON`
* `GEANT4_USE_OPENGL_X11=ON`
* `GEANT4_USE_QT=ON`
* `GEANT4_USE_RAYTRACER_X11=ON`
* `GEANT4_USE_XM=ON`

This enables C++17 and the optional GUI; the GUI needs the `GEANT4_USE_OPENGL_X11=ON` and `GEANT4_USE_QT=ON`. You probably also want to set `CMAKE_INSTALL_PREFIX` and `GEANT4_INSTALL_DATADIR`.

Assuming that it is installed in e.g. `~/code/geant4/geant4.10.04.p02-install`, then load Geant4 in the current shell by running `source ~/code/geant4/geant4.10.04.p02-install/bin/geant4.sh`.
If you intend to use MiniScatter with Jupyter [5], do this in the shell you plan to start Jupyter before launching `jupyter-notebook`.
You must also remember to load Geant4 before compiling and running MiniScatter.

On Fedora we have also seen that you need the following packages intalled in order to make the GUI work:
 * glui-devel
 * libglu
 * libglu-devel
 * libXmu-devel
 * qt-devel
 * glx-utils

Note that if your `.bashrc` (etc.) is loading a modified environment (such as Anaconda's default installation), this may lead to nonstandard version of e.g. QT libraries being used instead of the ones installed with your distribution's package manager.
To be safe, make sure anaconda (etc.) is sourced when installing Geant4.
In order to get rid of this problem, one must first remove the anaconda `source` command from `.bashrc`, then completely log the user out and in again (or reboot), before running cmake/make/make install for Geant4.

### ROOT
For ROOT, simply install it using your distribution's package manager, e.g. `dnf` on Fedora.
If compiling from source, these are appropriate flags in addition to `CMAKE_INSTALL_PREFIX`:
* `CMAKE_CXX_STANDARD=17`

Note that if using Jupyter etc. from a virtualenv or Anaconda, this should probably be loaded before compiling ROOT, so that ROOT gets linked to the correct Python libraries.

### MiniScatter
Then, clone MiniScatter into a new folder, e.g. from ~/code run `git clone https://github.com/kyrsjo/MiniScatter.git`.
Enter `~/code/MiniScatter`, and create a build folder `~/code/MiniScatter/build`.
In this folder, first configure the build system using CMAKE [6] by running `cmake ../.`, and then compile using `make -j N` where N is the number of CPUs to use; if for example you are running on a 16-core machine use `make -j 16`.
If all goes well, you have now built the executable!

Note that you may also change various build options, for example to enable debugging symbols to be written to the executable.
To do this, run `ccmake .` in your build folder, change what you need to change in the menu, then reconfigure and regenerate the makefiles.
Finally, run `make` again.
