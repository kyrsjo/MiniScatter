#!/usr/bin/env python3

import subprocess
import ROOT
import ROOT.TFile, ROOT.TVector

def runScatter(simSetup, quiet=False):
    "Run a MiniScatter simulation, given the parameters that are described by running './MiniScatter -h'. as the map simSetup."

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
            print ("Found PRESS="+str(simSetup["PRESS"]) + " but no MAT. This makes no sense.")
            exit(1)

    if "DIST" in simSetup:
        cmd += ["-d", str(simSetup["DIST"])]

    if "ANG" in simSetup:
        cmd += ["-a", str(simSetup["ANG"])]

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
            print ("ZOFFSET_BACKTRACK=",simSetup["ZOFFSET_BACKTRACK"], "is inconsistent with ZOFFSET=",simSetup["ZOFFSET"])
            exit(1)
    else:
        if "ZOFFSET_BACKTRACK" in simSetup or simSetup["ZOFFSET_BACKTRACK"] == False:
            print ("ZOFFSET_BACKTRACK=",simSetup["ZOFFSET_BACKTRACK"], "is inconsistent with ZOFFSET=",simSetup["ZOFFSET"])
            exit(1)

    if "COVAR" in simSetup:
        if len(simSetup["COVAR"]) == 3:
            cmd += ["-c", str(simSetup["COVAR"][0]) + ":" + str(simSetup["COVAR"][1]) + ":" + str(simSetup["COVAR"][2])]
        elif len(simSetup["COVAR"]) == 6:
            cmd += ["-c", str(simSetup["COVAR"][0]) + ":" + str(simSetup["COVAR"][1]) + ":" + str(simSetup["COVAR"][2])+"::" \
                    + str(simSetup["COVAR"][3]) + ":" + str(simSetup["COVAR"][4]) + ":" + str(simSetup["COVAR"][5])]
        else:
            print("Expected len(COVAR) == 3 or 6")
            exit(1)
    if "SEED" in simSetup:
        cmd += ["-s", str(simSetup["SEED"])]

    if "OUTNAME" in simSetup:
        cmd += ["-f", simSetup["OUTNAME"]]

    if "QUICKMODE" in simSetup:
        cmd += ["-q"]

    if "MINIROOT" in simSetup:
        cmd += ["-r"]

    if "CUTOFF_ENERGYFRACTION" in simSetup:
        cmd += ["--cutoffEnergyFraction", str(simSetup["CUTOFF_ENERGYFRACTION"])]

    if "CUTOFF_RADIUS" in simSetup:
        cmd += ["--cutoffRadius", str(simSetup["CUTOFF_RADIUS"])]

    cmdline = ""
    for c in cmd:
        cmdline += c + " "
    if not quiet:
        print ("Running command line: '" + cmdline[:-1] + "'")
    runResults = subprocess.run(cmd, close_fds=True, stdout=subprocess.PIPE)
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
        for p in ("x","y"):
            dataName = det + "_" + p + "_TWISS"
            if not dataFile.GetListOfKeys().Contains(dataName):
                raise KeyError("Object {} not found in file".format(dataName,filename))
            twissData = dataFile.Get(dataName)
            twiss[det][p] = {'eps':twissData[0], 'beta':twissData[1], 'alpha':twissData[2]}

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
            objects[objName].SetDirectory(0) # make it independent of the datafile TFile

    if getRaw:
        return (twiss, numPart, objects, dataFile)
    else:
        dataFile.Close()
        return(twiss, numPart, objects)
