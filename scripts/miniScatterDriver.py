#!/usr/bin/env python3

import subprocess
import ROOT
import ROOT.TFile, ROOT.TVector

def runScatter(THICK:float=None, MAT:str=None, DIST:float=None, ANG:float=None, PHYS:str=None, \
               N:int=None, ENERGY:float=None, BEAM:str=None, \
               XOFFSET:float=None, ZOFFSET:float=None, ZOFFSET_BACKTRACK:bool=None, COVAR:tuple=None, \
               SEED:int=None, OUTNAME:str=None, QUICKMODE:bool=None, quiet=False) -> None:
    #Parameters are descibed by running "./Miniscatter -h"

    cmd = ["./MiniScatter"]

    if THICK != None:
        cmd += ["-t", str(THICK)]

    if MAT != None:
        cmd += ["-m", MAT]

    if DIST != None:
        cmd += ["-d", str(DIST)]

    if ANG != None:
        cmd += ["-a", str(ANG)]

    if PHYS != None:
        cmd += ["-p", PHYS]

    if N != None:
        cmd += ["-n", str(N)]

    if ENERGY != None:
        cmd += ["-e", str(ENERGY)]

    if BEAM != None:
        cmd += ["-b", str(BEAM)]

    if XOFFSET != None:
        cmd += ["-x", str(XOFFSET)]

    if ZOFFSET != None:
        if ZOFFSET_BACKTRACK == None or ZOFFSET_BACKTRACK == False:
            cmd += ["-z", str(ZOFFSET)]
        elif ZOFFSET_BACKTRACK == True:
            cmd += ["-z", "*"+str(ZOFFSET)]
        else:
            print ("ZOFFSET_BACKTRACK=",ZOFFSET_BACKTRACK, "is inconsistent with ZOFFSET=",ZOFFSET)
            exit(1)
    else:
        if not (ZOFFSET_BACKTRACK == None or ZOFFSET_BACKTRACK == False):
            print ("ZOFFSET_BACKTRACK=",ZOFFSET_BACKTRACK, "is inconsistent with ZOFFSET=",ZOFFSET)
            exit(1)

    if COVAR !=None:
        if len(COVAR) == 3:
            cmd += ["-c", str(COVAR[0]) + ":" + str(COVAR[1]) + ":" + str(COVAR[2])]
        elif len(COVAR) == 6:
            cmd += ["-c", str(COVAR[0]) + ":" + str(COVAR[1]) + ":" + str(COVAR[2]) +"::" \
                    + str(COVAR[3]) + ":" + str(COVAR[4]) + ":" + str(COVAR[5])]
        else:
            print("Expected len(COVAR) == 3 or 6")
            exit(1)
    if SEED != None:
        cmd += ["-s", str(SEED)]

    if OUTNAME != None:
        cmd += ["-f", OUTNAME]

    if QUICKMODE != None:
        cmd += ["-q"]

    cmdline = ""
    for c in cmd:
        cmdline += c + " "
    if not quiet:
        print ("Running command line: '" + cmdline[:-1] + "'")
    runResults = subprocess.run(cmd, close_fds=True, stdout=subprocess.PIPE)
    #print (runResults)
    if not quiet:
        print ("Done!")

def getData(filename="plots/output.root", quiet=False):
    data = ROOT.TFile(filename)
    
    x_init = data.Get("initPhasespaceX_TWISS")
    y_init = data.Get("initPhasespaceY_TWISS")

    x_final = data.Get("trackerPhasespaceX_cutoff_TWISS")
    y_final = data.Get("trackerPhasespaceY_cutoff_TWISS")

    #print("Got initial parameters:")
    if not quiet:
        print("X :", x_init[0],x_init[1],x_init[2])
        print("Y :", y_init[0],y_init[1],y_init[2])

    final_tup = (x_final[0],x_final[1],x_final[2],y_final[0],y_final[1],y_final[2])

    data.Close()
    
    return(final_tup)
