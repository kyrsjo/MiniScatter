#!/usr/bin/env python3

import subprocess
import os
import ROOT
import ROOT.TFile, ROOT.TVector

def runScatter(simSetup, quiet=False):
    "Run a MiniScatter simulation, given the parameters that are described by running './MiniScatter -h'. as the map simSetup."

    for key in simSetup.keys():
        if not key in ("THICK", "MAT", "PRESS", "DIST", "ANG", "WORLDSIZE", "PHYS",\
                       "N", "ENERGY", "BEAM", "XOFFSET", "ZOFFSET", "ZOFFSET_BACKTRACK",\
                       "COVAR", "BEAM_RCUT", "SEED", \
                       "OUTNAME", "OUTFOLDER", "QUICKMODE", "MINIROOT",\
                       "CUTOFF_ENERGYFRACTION", "CUTOFF_RADIUS", "EDEP_DZ"):
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
        cmd += ["-d", str(simSetup["DIST"])]

    if "ANG" in simSetup:
        cmd += ["-a", str(simSetup["ANG"])]

    if "WORLDSIZE" in simSetup:
        cmd += ['-w', str(simSetup["WORLDSIZE"])]

    if "PHYS" in simSetup:
        cmd += ["-p", simSetup["PHYS"]]

    if "N" in simSetup:
        cmd += ["-n", str(simSetup["N"])]

    if "ENERGY" in simSetup:
        cmd += ["-e", str(simSetup["ENERGY"])]

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
        if "ZOFFSET_BACKTRACK" in simSetup or simSetup["ZOFFSET_BACKTRACK"] == False:
            raise ValueError("ZOFFSET_BACKTRACK=" +\
                             str(simSetup["ZOFFSET_BACKTRACK"]) +\
                             " is inconsistent with ZOFFSET="+str(simSetup["ZOFFSET"]))

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

    if "MAGNET" in simSetup:
        for mag in simSetup["MAGNET"]:
            mag_cmd = ""
            if mag["mag_pos_relative"] == True:
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
    if not quiet:
        print ("Running command line: '" + cmdline[:-1] + "'")
        print ("RunFolder = '"+ runFolder +"'")

    runResults = subprocess.run(cmd, close_fds=True, stdout=subprocess.PIPE, cwd=runFolder)
    #print (runResults)
    if not quiet:
        print ("Done!")

#Names of the planes in which the twiss parameters / number of particles of each type
# have been extracted
twissDets    = ("init","target_exit","target_exit_cutoff","tracker","tracker_cutoff")
numPartDets = ("tracker", "tracker_cutoff", "target", "target_cutoff")

def getData(filename="plots/output.root", quiet=False, getRaw=False, getObjects=None):
    """
    Collects data from the ROOT file, and optionally returns the file for looping over the ttrees.
    If the file is returned, the caller is responsible for closing it.
    """
    dataFile = ROOT.TFile(filename,'READ')

    twiss = {}
    for det in twissDets:
        twiss[det] = {}
        for pla in ("x","y"):
            dataName = det + "_" + pla + "_TWISS"
            if not dataFile.GetListOfKeys().Contains(dataName):
                raise KeyError("Object {} not found in file".format(dataName,filename))
            twissData = dataFile.Get(dataName)
            twiss[det][pla] = {'eps':twissData[0], 'beta':twissData[1], 'alpha':twissData[2]}
            if len(twissData) > 3:
                twiss[det][pla]['posAve'] = twissData[3]
                twiss[det][pla]['angAve'] = twissData[4]

    numPart = {}
    for det in numPartDets:
        numPart[det] = {}
        if not dataFile.GetListOfKeys().Contains(det+"_ParticleTypes_PDG") or \
           not dataFile.GetListOfKeys().Contains(det+"_ParticleTypes_numpart"):
            if not quiet:
                print("No particles found for det={}".format(det))
            continue
        numPart_PDG = dataFile.Get(det+"_ParticleTypes_PDG")
        numPart_num = dataFile.Get(det+"_ParticleTypes_numpart")
        for i in range(len(numPart_PDG)):
            numPart[det][int(numPart_PDG[i])] = numPart_num[i]

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
