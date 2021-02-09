#!/usr/bin/env python3

"""
This file is part of MiniScatter.

MiniScatter is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MiniScatter is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MiniScatter.  If not, see <https://www.gnu.org/licenses/>.
"""

import subprocess
import os
import ROOT
import ROOT.TFile, ROOT.TVector
import datetime

def runScatter(simSetup, quiet=False,allOutput=False, logName=None, onlyCommand=False):
    "Run a MiniScatter simulation, given the parameters that are described by running './MiniScatter -h'. as the map simSetup."

    if quiet and allOutput:
        raise AssertionError("Setting both 'quiet' and 'alloutput' makes no sense")

    for key in simSetup.keys():
        if not key in ("THICK", "MAT", "PRESS", "DIST", "ANG", "TARG_ANG", "WORLDSIZE", "PHYS", "PHYS_CUTDIST",\
                       "N", "ENERGY", "ENERGY_FLAT",\
                       "BEAM", "XOFFSET", "ZOFFSET", "ZOFFSET_BACKTRACK",\
                       "COVAR", "BEAM_RCUT", "SEED", \
                       "OUTNAME", "OUTFOLDER", "QUICKMODE", "MINIROOT",\
                       "CUTOFF_ENERGYFRACTION", "CUTOFF_RADIUS", "EDEP_DZ", "ENG_NBINS"):
            if key.startswith("MAGNET"):
                continue
            raise KeyError("Did not expect key {} in the simSetup".format(key))

    cmd = ["./MiniScatter"]

    if "THICK" in simSetup:
        cmd += ["-t", str(simSetup["THICK"])]
    if "MAT" in simSetup:
        if "PRESS" in simSetup:
            cmd += ["-m", simSetup["MAT"]+'::'+str(simSetup["PRESS"])]
        else:
            cmd += ["-m", simSetup["MAT"]]
    else:
        if "PRESS" in simSetup:
            raise ValueError("Found PRESS="+str(simSetup["PRESS"]) + " but no MAT. This makes no sense.")

    if "DIST" in simSetup:
        if simSetup["DIST"] == "NONE":
            cmd += ["-d", "NONE"]
        elif type(simSetup["DIST"]) == float:
            cmd += ["-d", str(simSetup["DIST"])]
        else: # It's a list of distances
            distStr = ""
            for d in simSetup["DIST"]:
                distStr = distStr + str(d) + ":"
            distStr = distStr[0:-1]
            print("distStr=",distStr)
            cmd += ["-d", distStr]

    if "ANG" in simSetup:
        cmd += ["-a", str(simSetup["ANG"])]

    if "TARG_ANG" in simSetup:
        cmd += ["-A", str(simSetup["TARG_ANG"])]

    if "WORLDSIZE" in simSetup:
        cmd += ['-w', str(simSetup["WORLDSIZE"])]

    if "PHYS" in simSetup:
        cmd += ["-p", simSetup["PHYS"]]

    if "PHYS_CUTDIST" in simSetup:
        cmd += ["--physCutoffDist", str(simSetup["PHYS_CUTDIST"])]

    if "N" in simSetup:
        cmd += ["-n", str(simSetup["N"])]

    if "ENERGY" in simSetup:
        cmd += ["-e", str(simSetup["ENERGY"])]

    if "ENERGY_FLAT" in simSetup:
        cmd += ["--energyDistFlat", str(simSetup["ENERGY_FLAT"][0])+':'+str(simSetup["ENERGY_FLAT"][1]) ]

    if "BEAM" in simSetup:
        cmd += ["-b", str(simSetup["BEAM"])]

    if "XOFFSET" in simSetup:
        cmd += ["-x", str(simSetup["XOFFSET"])]

    if "ZOFFSET" in simSetup:
        if (not "ZOFFSET_BACKTRACK" in simSetup) or (simSetup["ZOFFSET_BACKTRACK"] == False):
            cmd += ["-z", str(simSetup["ZOFFSET"])]
        elif simSetup["ZOFFSET_BACKTRACK"] == True:
            cmd += ["-z", "*"+str(simSetup["ZOFFSET"])]
        else:
            raise ValueError("ZOFFSET_BACKTRACK=" +\
                             str(simSetup["ZOFFSET_BACKTRACK"]) +\
                             " is inconsistent with ZOFFSET="+str(simSetup["ZOFFSET"]))
    else:
        if "ZOFFSET_BACKTRACK" in simSetup:
            raise ValueError("ZOFFSET_BACKTRACK present but not ZOFFSET?")

    if "COVAR" in simSetup:
        if len(simSetup["COVAR"]) == 3:
            cmd += ["-c", str(simSetup["COVAR"][0]) + ":" + str(simSetup["COVAR"][1]) + ":" + str(simSetup["COVAR"][2])]
        elif len(simSetup["COVAR"]) == 6:
            cmd += ["-c", str(simSetup["COVAR"][0]) + ":" + str(simSetup["COVAR"][1]) + ":" + str(simSetup["COVAR"][2])+"::" \
                    + str(simSetup["COVAR"][3]) + ":" + str(simSetup["COVAR"][4]) + ":" + str(simSetup["COVAR"][5])]
        else:
            raise ValueError("Expected len(COVAR) == 3 or 6")

    if "BEAM_RCUT" in simSetup:
        cmd += ["--beamRcut", str(simSetup["BEAM_RCUT"])]

    if "SEED" in simSetup:
        cmd += ["-s", str(simSetup["SEED"])]

    if "OUTNAME" in simSetup:
        cmd += ["-f", simSetup["OUTNAME"]]

    if "OUTFOLDER" in simSetup:
        cmd += ["-o", simSetup["OUTFOLDER"]]

    if "QUICKMODE" in simSetup:
        if simSetup["QUICKMODE"] == True:
            cmd += ["-q"]
        else:
            assert simSetup["QUICKMODE"] == False

    if "MINIROOT" in simSetup:
        if simSetup["MINIROOT"] == True:
            cmd += ["-r"]
        else:
            assert simSetup["MINIROOT"] == False

    if "CUTOFF_ENERGYFRACTION" in simSetup:
        cmd += ["--cutoffEnergyFraction", str(simSetup["CUTOFF_ENERGYFRACTION"])]

    if "CUTOFF_RADIUS" in simSetup:
        cmd += ["--cutoffRadius", str(simSetup["CUTOFF_RADIUS"])]

    if "EDEP_DZ" in simSetup:
        cmd += ["--edepDZ", str(simSetup["EDEP_DZ"])]

    if "ENG_NBINS" in simSetup:
        cmd += ["--engNbins", str(simSetup["ENG_NBINS"])]

    if "MAGNET" in simSetup:
        for mag in simSetup["MAGNET"]:
            mag_cmd = ""
            if "mag_pos_relative" in mag and mag["mag_pos_relative"] == True:
                mag_cmd += "*"
            mag_cmd += str(float(mag["pos"]))       + ":"
            mag_cmd += str(mag["type"])             + ":"
            mag_cmd += str(float(mag["length"]))    + ":"
            mag_cmd += str(float(mag["gradient"]))
            for k,v in mag["keyval"].items():
                mag_cmd += ":" + str(k)+"="+str(v)
            cmd += ["--magnet", mag_cmd]

    cmdline = ""
    for c in cmd:
        cmdline += c + " "

    runFolder = os.path.dirname(os.path.abspath(__file__))

    if logName is None:
        logName = os.path.join(runFolder,"MiniScatterLog_" + datetime.datetime.now().isoformat()+".txt")
    logFile = open(logName, 'w')

    if not quiet:
        print ("Running command line: '" + cmdline[:-1] + "'")
        if not onlyCommand:
            print ("RunFolder = '" + runFolder + "'")
            print ("logName   = '" + logName   + "'")

    # runResults = subprocess.run(cmd, close_fds=True, stdout=subprocess.PIPE, cwd=runFolder)
    # #print (runResults)

    #Inspired by https://stackoverflow.com/questions/52545512/realtime-output-from-a-shell-command-in-jupyter-notebook
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, close_fds=True, cwd=runFolder, universal_newlines=True)
    linebuff = ""
    spinnerState = None # 0=/, 1=-, 2=\, 3=| (cycle); None: Last printout was an event, so issue a newline not a carriage return

    for line in iter(process.stdout.readline, ''):
        if allOutput:
            print(line.rstrip())
            logFile.write(line.rstrip()+"\n")

        elif not quiet and not onlyCommand:
            # Collect data until we have a whole line,
            # then print it if AND ONLY IF it is a progress report.
            ls = line.split('\n')
            lsl = len(ls)
            assert lsl <= 2
            linebuff += ls[0]
            if len(ls) > 1:
                logFile.write(ls[0]+'\n')
                logFile.flush()

                if linebuff.startswith("Event#"):
                    if not spinnerState is None:
                        print('',end='\n')
                    print(linebuff)
                    spinnerState = None

                else:
                    if spinnerState is None:
                        spinnerState = 0
                    print('',end='\r')
                    spinnerState = (spinnerState+1) % 4
                    if spinnerState == 0:
                        print ('/', end='')
                    elif spinnerState == 1:
                        print ('-', end='')
                    elif spinnerState == 2:
                        print ("\\", end='')
                    elif spinnerState == 3:
                        print ('|', end='')
                    else:
                        print("WTF?")
                linebuff = ls[1]
    process.stdout.close()
    returncode = process.wait()
    logFile.close()

    if not quiet:
        print('',end='\n')
        print ("Done!")
    
    if returncode != 0:
        print("Errors encountered during simulation")
        print("Please check logfile '"+logName+"'")
        print("Command line: '" +cmdline[:-1] +"'")

        raise SimulationError(returncode,logName,cmdline)

def getData(filename="plots/output.root", quiet=False, getRaw=False, getObjects=None):
    """
    Collects data from the ROOT file, and optionally returns the file for looping over the ttrees.
    If the file is returned, the caller is responsible for closing it.
    """
    dataFile = ROOT.TFile(filename,'READ')

    #Load TWISS data
    twiss = {}
    for key in dataFile.GetListOfKeys():
        keyName = key.GetName()
        if not keyName.endswith("_TWISS"):
            continue #Skip this one

        xy = None
        if '_x_' in keyName:
            xy = 'x'
        elif '_y_' in keyName:
            xy = 'y'
        else:
            raise AssertionError("Could not determine 'xy' for TwissDet '"+key+"'")

        key_short = keyName.replace('_'+xy+'_','')
        if '_TWISS' in key_short:
            key_short = key_short.replace('_TWISS','')
        else:
            #Some are named like ..._x_TWISS e.g. init_x_TWISS, and one arm of the "_x_" is already deleted.
            key_short = key_short.replace('TWISS','')

        if not key_short in twiss:
            twiss[key_short] = {}

        twissData = dataFile.Get(keyName)
        twiss[key_short][xy] = {'eps':twissData[0], 'beta':twissData[1], 'alpha':twissData[2]}
        if len(twissData) > 3:
            twiss[key_short][xy]['posAve'] = twissData[3]
            twiss[key_short][xy]['angAve'] = twissData[4]
        if len(twissData) > 4:
            twiss[key_short][xy]['posVar'] = twissData[5]
            twiss[key_short][xy]['angVar'] = twissData[6]
            twiss[key_short][xy]['coVar']  = twissData[7]

    #Load the NumPart data
    numPart = {}
    for key in dataFile.GetListOfKeys():
        keyName = key.GetName()
        if not keyName.endswith("_ParticleTypes_PDG"):
            continue #Skip
        keyBase = keyName[0:-18]
        if not (dataFile.GetListOfKeys().Contains(keyBase+"_ParticleTypes_numpart")):
            raise AssertionError("Found PDG IDs but not numpart for '"+keyBase+"'")
        numPart[keyBase] = {}

        numPart_PDG = dataFile.Get(keyBase+"_ParticleTypes_PDG")
        numPart_num = dataFile.Get(keyBase+"_ParticleTypes_numpart")

        for i in range(len(numPart_PDG)):
            numPart[keyBase][int(numPart_PDG[i])] = int(numPart_num[i])

    #Load the requested objects
    objects = None
    if getObjects:
        objects = {}
        for objName in getObjects:
            if not dataFile.GetListOfKeys().Contains(objName):
                dataFile.ls()
                dataFile.Close()
                raise KeyError("Object {} not found in file {}".format(objName,filename))
            #Clone the object -- the name can be changed later, but make it unique
            objects[objName] = dataFile.Get(objName).Clone(objName+"-localClone")
            # make it independent of the datafile TFile (histograms only)
            try:
                objects[objName].SetDirectory(0)
            except AttributeError:
                pass

    if getRaw:
        return (twiss, numPart, objects, dataFile)
    else:
        dataFile.Close()
        return(twiss, numPart, objects)

def getData_tryLoad(simSetup, quiet=False, getRaw=False, getObjects=None, tryload=True, allOutput=False):
    """
    Checks if the ROOT file given by the parameters in simsetup exists;
    if it does then load.

    If it does not exist, run the simulation then load.

    This is quite practical for e.g. Jupyter, since when using it re-running a notebook
    will quickly load already computed data without having to comment out
    the call to runScatter().
    """

    ROOTfilename = 'output.root'
    if "OUTNAME" in simSetup:
        ROOTfilename = simSetup["OUTNAME"]+".root"
    if "OUTFOLDER" in simSetup:
        ROOTfilename = os.path.join(simSetup["OUTFOLDER"],ROOTfilename)
        if not os.path.isdir(simSetup["OUTFOLDER"]):
            os.mkdir(simSetup["OUTFOLDER"])
    else:
        runFolder = os.path.dirname(os.path.abspath(__file__))
        ROOTfilename = os.path.join(runFolder,"plots",ROOTfilename)

    logName = ROOTfilename[:-5]+".txt"

    if not os.path.exists(ROOTfilename):
        print("Did not find any pre-computed data at '"+ROOTfilename+"', computing now.")
        runScatter(simSetup, quiet=quiet, allOutput=allOutput, logName=logName)
    elif tryload==False:
        if not quiet:
            print("TryLoad is False, computing now.")
        runScatter(simSetup, quiet=quiet, allOutput=allOutput, logName=logName)
    else:
        if not quiet:
            print("Found a file at '"+ROOTfilename+"', loading!")

    return getData(ROOTfilename, quiet=quiet, getRaw=getRaw, getObjects=getObjects)

class SimulationError(Exception):
    "Base class for simulation crashes"
    def __init__(self,returncode,logName,cmdline):
        self.returncode = returncode
        self.logName = logName
        self.cmdline = cmdline
        super().__init__("Simulation crashed with returncode = " + str(self.returncode) + ", see log file '" +self.logName + "' for more info.")

    returncode  = None
    logName = None
    returncode  = None
