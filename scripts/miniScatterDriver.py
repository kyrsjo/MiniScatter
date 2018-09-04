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

    cmdline = ""
    for c in cmd:
        cmdline += c + " "
    if not quiet:
        print ("Running command line: '" + cmdline[:-1] + "'")
    runResults = subprocess.run(cmd, close_fds=True, stdout=subprocess.PIPE)
    #print (runResults)
    if not quiet:
        print ("Done!")

def getData(filename="plots/output.root", quiet=False, getRaw=False):
    """
    Collects twiss parameters from the ROOT file, and optionally returns the file for looping over the ttrees.
    If the file is returned, the caller is responsible for closing it.
    """
    data = ROOT.TFile(filename)

    x_init = data.Get("initPhasespaceX_TWISS")
    y_init = data.Get("initPhasespaceY_TWISS")

    x_final = data.Get("trackerPhasespaceX_cutoff_TWISS")
    y_final = data.Get("trackerPhasespaceY_cutoff_TWISS")

    if not quiet:
        print("Got initial parameters:")
        print("X :", x_init[0],x_init[1],x_init[2])
        print("Y :", y_init[0],y_init[1],y_init[2])

    numPart = {}
    for det in ("tracker", "tracker_cutoff", "target", "target_cutoff"):
        numPart[det] = {}
        if not data.GetListOfKeys().Contains(det+"_ParticleTypes_PDG") or \
           not data.GetListOfKeys().Contains(det+"_ParticleTypes_numpart"):
            print("No particles found for det={}".format(det))
            continue
        numPart_PDG = data.Get(det+"_ParticleTypes_PDG")
        numPart_num = data.Get(det+"_ParticleTypes_numpart")
        for i in range(len(numPart_PDG)):
            numPart[det][int(numPart_PDG[i])] = numPart_num[i]

    final_tup = (x_final[0],x_final[1],x_final[2],y_final[0],y_final[1],y_final[2])

    if getRaw:
        return (final_tup, numPart, data)
    else:
        data.Close()
        return(final_tup, numPart)
