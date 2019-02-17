# Running MiniScatter on CERN SWAN

MiniScatter can be ran on CERN's SWAN system [1], which provides a ready-to-go Linux environment in a virtual machine with Geant4, ROOT, Jupyter, and Python3.

## Prerequisites
In order to use this, you need a CERN account.
If you do not have that, you can still install and run MiniScatter on your local machine.

Furthermore, you need:
 * Access to CERNbox / EOS
 * Access to SWAN

## First-time setup
The easiest way to setup MiniScatter is to login to lxplus7 with SSH and compile it from there.
This can be done by simply opening a terminal and running `ssh USERNAME@lxplus7.cern.ch` if you are on a Linux or Mac machine;
from Windows a 3rd party program such as Putty [2] is required.

Once connected to lxplus7, download and compile MiniScatter as follows:
 * Change folder to your CERNBOX:
   `cd /eos/user/first_letter_in_username/username/`
   (remember to use your actual username and the first letter in the username in this command)
 * Clone the GIT repository:
   `git clone https://github.com/kyrsjo/MiniScatter.git`<sup>1</sup>
 * Enter the MiniScatter folder:
   `cd MiniScatter`
 * Load the LCG computing environment:
   `source setupLCG.sh`
 * Compile it:
   ```
   mkdir build_SWAN
   cmake ../.
   make -j
   ```

Note that if a SSH client cannot be installed on your local machine, it is possible to access a terminal once logged into SWAN.
However this browser-based terminal will likely be significantly less confortable than a "real" SSH client.
To do so, click on the  `>_` symbol in upper right part of screen, then follow the instructions above starting with cloning the GIT repository. 

Logging in to LxPlus7 can also be a good way for debugging issues with MiniScatter, as it provides an already setup environment including the supporting libraries<sup>2</sup>.
As an example, if you are standing in the `build_SWAN` folder you can run MiniScatter as follows:
`./MiniScatter -n 100`
Much output but no errors should be produced.
For more information, see [Command Line Use](CommandLineUse.md).

## Using MiniScatter with SWAN
Assuming that you have compiled MiniScatter in CERNBOX as described above and the compilation was sucessfull, you are now ready to use MiniScatter through it's Python interface.

 * Spin up a SWAN virtual machine.
   In the configuration dialog, make sure to select a software stack and platform which matches the one used in `setupLCG.sh` when compiling; currently this is:
   - Software stack = `94 Python3`
   - Platform = `x86_64_centos7-gcc7-opt`
   - Environment script = `geant4.sh`
 * Switch to tab `CERNBox`, then open the MiniScatter folder
 * Open one notebooks in the `examples` folder, and try to execute it step-by-step

It should run and produce data and plots. You may similarly create new notebooks in any folder; the only prerequisite is that the MiniScatter python libraries are loaded.
To do this, have a cell containing the following near the top of your notebook
```
#Setup MiniScatter
import sys
MiniScatter_path="../build/."
sys.path.append(MiniScatter_path)

import miniScatterDriver
import miniScatterScanner
````
where MiniScatter_path should point to the folder where you built the MiniScatter binary.
Please note that you may need to change some of the examples to point to the right folder, i.e. if you built MiniScatter in `build_SWAN`, the path setting should be `MiniScatter_path="../build_SWAN/."`.

For further information, see the documentation on the [Python Interface](PyInterface.md).

## References

[1] https://swan.cern.ch

[2] https://www.putty.org/

## Footnotes

1: You can of course also clone MiniScatter from a fork or using the SSH url, if you are logged into GitHub on LxPlus7.

2: But remember to `source setupLCG.sh`!
