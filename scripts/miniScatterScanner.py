import miniScatterDriver

import numpy as np

import h5py

from queue import Queue
import threading

import sys
import traceback

import os

SEED = 1

#Used for particle-type countings, what to return
PDG_keep    = (11,-11,22,2212,2112,'other');

m_twiss = 0.511 #[MeV/c^2], assumed mass of the particles used for TWISS computation, assume electrons.

def ScanMiniScatter(scanVar,scanVarRange,baseSimSetup, \
                    NUM_THREADS=4, tryLoad=False, COMMENT=None, QUIET=True, \
                    detailedAnalysisRoutine=None, detailedAnalysisRoutine_names=None, \
                    cleanROOT=True):
    """
    This routine is built to scan arbitrary parameters with MiniScatter.
    It can cache the results in HDF5-files with long and difficult names, as well as call detailed analysis routines.
    """

    global SEED # Updated every time one does a scan

    ### Create the output arrays ###

    twiss = {}
    for det in miniScatterDriver.twissDets:
        twiss[det] = {}
        for p in ('x','y'):
            twiss[det][p] = {}
            for t in ('eps','beta','alpha','sigma'):
                twiss[det][p][t] = np.zeros_like(scanVarRange)

    numPart = {}
    for det in miniScatterDriver.numPartDets:
        numPart[det] = {}
        for pdg in PDG_keep:
            numPart[det][pdg] = np.zeros_like(scanVarRange)

    analysis_output = None
    if detailedAnalysisRoutine:
        analysis_output = {}
        assert type(detailedAnalysisRoutine_names)==list
        for name in detailedAnalysisRoutine_names:
            analysis_output[name] = np.zeros_like(scanVarRange)
    else:
        assert detailedAnalysisRoutine_names == None

    #Later we have assumed electron beam for emittance calculation
    assert baseSimSetup["BEAM"] == "e+" or baseSimSetup["BEAM"] == "e-"
    assert "ENERGY" in baseSimSetup or scanVar=="ENERGY"

    #Sanity check
    if scanVar in baseSimSetup:
        print("Please do not put scanVar in the baseSimSetup!")
        raise ValueError("Found scanVar in the baseSimSetup")

    #Loading a pre-ran simulation?
    loadFileName = "SaveSim_{}_{}.h5".format(scanVar,COMMENT)
    print("LoadFile filename and status: '" + loadFileName + "'", tryLoad)
    if tryLoad:
        print("Loading...")
        try:
            loadFile = h5py.File(loadFileName,mode='r')

            #Check that the file is workable
            scanVar_loaded = loadFile.attrs["scanVarName"]
            if scanVar != scanVar_loaded:
                print("Scan variables did not match even tough the filename did")
                print("please run with tryLoad=False to recompute.")
                loadFile.close()
                raise ValueError("scanVar did not match loaded file")

            for key in baseSimSetup.keys():
                if not (key in loadFile.attrs):
                    print ("Key '{}' found in baseSimSetup but not in the file.".format(key))
                    print ("Please run with tryLoad=False to recompute.")
                    loadFile.close()
                    raise ValueError("Found key {} in baseSimSetup but not in the file".format(key))

                if not loadFile.attrs[key] == baseSimSetup[key]:
                    print ("Value of key '{}' found in baseSimSetup".format(key))
                    print (" did not match what was found in the file.")
                    print ("Values: baseSimSetup={}, file={}".format(baseSimSetup[key],loadFile.attrs[key]))
                    print ("Please run with tryLoad=False to recompute.")
                    loadFile.close()
                    raise ValueError("Value of key {} did not match with baseSimSetup and file".format(key))

            for key in loadFile.attrs:
                if key == scanVar or key=='scanVarName':
                    continue
                if not key in baseSimSetup.keys():
                    print("Key {} found in file but not in baseSimSetup.".format(key))
                    print ("Please run with tryLoad=False to recompute.")
                    loadFile.close()
                    raise ValueError("Key {} found in file but not in baseSimSetup".format(key))

            scanVarRange_loaded = np.asarray(loadFile.attrs[scanVar])
            if len(scanVarRange) != len(scanVarRange_loaded) or \
               not np.all(np.equal(scanVarRange, scanVarRange_loaded)):
                print ("Scan variable ranges did not match, run with tryLoad=False to recompute.")
                print ("Now :", scanVarRange)
                print ("File:", scanVarRange_loaded)
                # for i in range(min(len(scanVarRange),len(scanVarRange_loaded))):
                #     print(scanVarRange[i], scanVarRange_loaded[i],\
                #           scanVarRange[i]-scanVarRange_loaded[i], \
                #           (scanVarRange[i]-scanVarRange_loaded[i])==0.0)
                loadFile.close()
                raise ValueError("ScanVar range did not match with loaded file")

            print("Scan variable ranges match, let's load!")

            for det in miniScatterDriver.twissDets:
                for p in ('x','y'):
                    for t in ('eps','beta','alpha','sigma'):
                        twiss[det][p][t] = np.asarray(loadFile["twiss_"+det+"_"+p+"_"+t])

            for det in miniScatterDriver.numPartDets:
                for pdg in PDG_keep:
                    arrayName = "numPart_"+det+"_"+str(pdg)
                    if not arrayName in loadFile:
                        print("Could not find numPart for {}, PDG={} in the file. Please recompute.".format(det,pdg))
                        loadFile.close()
                        raise ValueError("NumPart for det={}, PDG={} was not found in the file".format(det,pdg))

                    numPart[det][pdg] = np.asarray(loadFile[arrayName])

            if detailedAnalysisRoutine:
                for name in detailedAnalysisRoutine_names:
                    nameMangle = "ANALYSIS_"+name
                    if not (nameMangle) in loadFile:
                        print("Could not find '"+name+"' in the file. Please recompute.")
                        loadFile.close()
                        return
                    analysis_output[name] = np.asarray(loadFile["ANALYSIS_"+name])
            loadFile.close()

            print("Loaded! That was fast.")
            return (twiss, numPart, analysis_output)

        except OSError:
            print("File not found. Computing...")

    ### Build the job queue ###

    def computeOnePoint(var,i,lock):
        with lock:
            #print ("pressure = {0} [mbar] ({1}/{2})".format(p,i+1,len(pressures)))
            print("{} = {} ({}/{})".format(scanVar, var, i+1, len(scanVarRange)))

        filenameROOT = 'output_'+scanVar+"="+str(var)
        simSetup = baseSimSetup.copy()
        simSetup[scanVar]   = var
        simSetup["SEED"]    = SEED + i
        simSetup["OUTNAME"] = filenameROOT
        miniScatterDriver.runScatter(simSetup,quiet=QUIET)

        filenameROOTfile = "plots/"+filenameROOT+".root"
        badSim=False
        if os.path.isfile(filenameROOTfile):
            (twiss_singlesim,numPart_singleSim,datafile) = \
                miniScatterDriver.getData(filename=filenameROOTfile,quiet=QUIET,getRaw=True)
        else:
            with lock:
                print("Did not find file '{}', simulation crashed?".format(filenameROOTfile))
            badSim=True

        #Fill the emittance arrays
        if not badSim:
            gamma_rel = simSetup["ENERGY"]/m_twiss
            beta_rel  = np.sqrt(gamma_rel**2 - 1.0) / gamma_rel;
            for det in miniScatterDriver.twissDets:
                for p in ('x','y'):
                    for t in ('eps','beta','alpha'):
                        twiss[det][p][t][i] = twiss_singlesim[det][p][t]
                    twiss[det][p]['sigma'][i] = \
                        np.sqrt( twiss[det][p]['eps'][i] * twiss[det][p]['beta'][i]*1e6/  \
                                 (gamma_rel*beta_rel)                                     )

            #Fill the NumPart array
            for detDictKey in numPart_singleSim.keys():
                for pdg in numPart_singleSim[detDictKey].keys():
                    if pdg in PDG_keep:
                        numPart[detDictKey][pdg][i] = numPart_singleSim[detDictKey][pdg]
                    else:
                        numPart[detDictKey]['other'][i] += numPart_singleSim[detDictKey][pdg]
                        with lock:
                            print("Found pdg={} for detector={}".format(pdg,detDictKey))

            #Do special analysis over the TTrees
            if detailedAnalysisRoutine:
                #Put the call to the external routine in a try/catch,
                # so that the thread will actually finish correctly in case of a user error.
                try:
                    detailedData = detailedAnalysisRoutine(datafile)
                    for k in detailedData.keys():
                        analysis_output[k][i]=detailedData[k]
                except Exception as err:
                    with lock:
                        traceback.print_tb(err.__traceback__)

            if cleanROOT:
                datafile.Close()
                os.remove(filenameROOTfile)
                if not QUIET:
                    with lock:
                        print("Deleting '{}'.".format(filenameROOTfile))

    def threadWorker(jobQueue_local,lock):
        while not jobQueue_local.empty():
            (var,i) = jobQueue_local.get()
            computeOnePoint(var,i,lock)
            jobQueue_local.task_done()

    jobQueue = Queue(0)
    i=0 # Array index where the jobs should write their data
    for var in scanVarRange:
        jobQueue.put((var,i))
        i = i+1
    lock = threading.Lock() # Shared lock to synchronize output
    for i in range(NUM_THREADS):
        worker = threading.Thread(target=threadWorker, args=(jobQueue,lock))
        worker.start()

    jobQueue.join()
    SEED = SEED+i #Prepare the SEED for the next run of simulations

    print("Simulation complete, saving data to h5 for later retrival.")
    #Write out the data
    saveFile = h5py.File(loadFileName,mode="w")

    saveFile.attrs["scanVarName"] = scanVar
    saveFile.attrs[scanVar]       = scanVarRange
    for key in baseSimSetup.keys():
        if key in saveFile.attrs:
            saveFile.close()
            raise RuntimeError("Attribute '{}' already written to file? Probably a bug!".format(key))
        value = baseSimSetup[key]
        saveFile.attrs[key] = value

    for det in miniScatterDriver.twissDets:
        for p in ('x','y'):
            for t in ('eps','beta','alpha','sigma'):
                saveFile["twiss_"+det+"_"+p+"_"+t] = twiss[det][p][t]

    for det in miniScatterDriver.numPartDets:
        for pdg in PDG_keep:
            arrayName = "numPart_"+det+"_"+str(pdg)
            saveFile[arrayName] = numPart[det][pdg]

    if detailedAnalysisRoutine:
        for name in detailedAnalysisRoutine_names:
            nameMangle = "ANALYSIS_"+name
            saveFile["ANALYSIS_"+name] = analysis_output[name]

    saveFile.close()

    return (twiss, numPart, analysis_output)
