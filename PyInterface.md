# MiniScatter's Python Interface

As a layer on top of the (command line interface)[CommandLineUse.md], a Python interface has been provided.
This allows to run MiniScatter (once or in scans), and to extract data from these runs.

The Python libraries are contained in the the `scripts` folder, in the files [miniScatterDriver.py](scripts/miniScatterDriver.py) and [miniScatterScanner.py](scripts/miniScatterScanner.py).
When building MiniScatter, these libraries gets symlinked into the build folder.
This is done by symlink instead of file copy, so that one does not loose changes the script if modifying the "wrong" copy of the same source code.

Note that the libraries are written for Python3, and will not work with Python2.

For examples on how to use the modules, see the Jupyter notebooks in the examples folder.

## Importing the modules
One should always import the modules from the build folder, not the scripts folder in the sources, as they search for the MiniScatter binary in their installation path.
Thus, one can do:
```
#Setup MiniScatter
import sys
MiniScatter_path="../build/."
sys.path.append(MiniScatter_path)

import miniScatterDriver
import miniScatterScanner
```
After this, MiniScatter can be ran.

## The driver module
This module, found in [miniScatterDriver.py](scripts/miniScatterDriver.py), is used for running the simulation and for loading the `.root` file data.
Three functions are provided:

### `runScatter(simSetup,quiet=True)`
This function is used to run MiniScatter.
The first argument is a Python map describing the values of various command line arguments, and can as an example be used as:
```
mySimSetup = {}
mySimSetup["N"]        = 1000
mySimSetup["ENERGY"]   = 200
mySimSetup["BEAM"]     = "e-"
mySimSetup["THICK"]    = 20
mySimSetup["MATERIAL"] = "G4_WATER"

miniScatterDriver.runScatter(mySimSetup, True)
```

To see the full list of the various possible options, please see the source code of the file.
Note that keys/values are for the most part mapping directly to the various command line options found by running `./MiniScatter -h`.

### `getData(filename="plots/output.root", quiet=False, getRaw=False, getObjects=None)`
This function is used to load the output ROOT file produced by running MiniScatter.

By default, it loads the default output file (as if standing in the build directory),
and returns a list of twiss paramters for different detector planes,
and the number of each type of particle at each detector plane.
Optionally, a list of names in the ROOT file can be provided (i.e. histograms), and these objects are cloned and returned. 
Optionally, if getRaw is true, a pointer to the ROOT TFile object is also returned; it is then up to the callee to close the root file.

To see a list of the available names inside the ROOT file, use `TBrowser` or the command `rootls FILENAME.root`.

### `getData_tryLoad(simSetup, quiet=False, getRaw=False, getObjects=None, tryload=True)`
This function is a combination of `runScatter` and `getData`; it is used to only run the (potentially time-consuming) simulation if necessary.
This is very useful e.g. in Jupyter notebooks, as a cell containing a simulation will take a significant ammount of time to run the first time, and then the next time (e.g. after reloading the notebook) it will load very quickly, without any modification such as commenting out of the call to runScatter.

It works by checking if the output file name described by the simSetup is there, and if it is, getLoad is called to load the data.
If it is not there, it calls `runScatter` and then `getData`.
The return data is the same as what is passed from `getData`.

## The scanner module
This module, found in [miniScatterScanner.py](scripts/miniScatterScanner.py), is used for repeatedly running the simulation with different parameters.
If available, it can take advantage of multi-core parallelism, scheduling a number of jobs to run in parallel until the list of parameters is exhausted.
Furthermore, similarly to `getData_tryLoad` a caching mechanism is included, which means that on re-run the pre-computed simulation data can be loaded from file.

A single function, `ScanMiniScatter`, is provided for this.
It takes 3 mandatory arguments: `scanVar`, `scanVarRange` and `baseSimSetup`; the last one is a simSetup map like the one described above except that it must be missing the key for the variable to be scanned.
The two other ones are a string describing which variable to scan, and a list of the values that should be used with this variable.

The rest of the options are optional and describe things like how many CPUs to use, a comment to insert into the filename, in which folder to run MiniScatter, etc.
One important option is `getObjects`; this is the list of ROOT objects to get from each simulation file, and objects that are not included in this file are deleted, only possible to recover by re-running the scan.

