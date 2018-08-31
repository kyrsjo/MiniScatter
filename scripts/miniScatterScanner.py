import miniScatterDriver

import numpy as np

import h5py

from queue import Queue
import threading

import sys
import traceback

import os

SEED = 1

PDG_keep = (11,-11,22,2212,2112,'other');

def ScanMiniScatter(scanVar,scanVarRange,baseSimSetup, \
                    NUM_THREADS=4, tryLoad=False, COMMENT=None, QUIET=True, \
                    detailedAnalysisRoutine=None, detailedAnalysisRoutine_names=None, \
                    cleanROOT=True):
    """
    This routine is built to scan arbitrary parameters with MiniScatter.
    It can cache the results in HDF5-files with long and difficult names, as well as call detailed analysis routines.
    """

    global SEED # Updated every time one does a scan

    #Output arrays
    eps_x  = np.zeros_like(scanVarRange)
    eps_y  = np.zeros_like(scanVarRange)

    beta_x  = np.zeros_like(scanVarRange)
    beta_y  = np.zeros_like(scanVarRange)

    alpha_x = np.zeros_like(scanVarRange)
    alpha_y = np.zeros_like(scanVarRange)

    sigma_x = np.zeros_like(scanVarRange)
    sigma_y = np.zeros_like(scanVarRange)

    numPart = {}
    for pdg in PDG_keep:
        numPart[pdg] = np.zeros_like(scanVarRange)

    analysis_output = {}
    if detailedAnalysisRoutine:
        assert type(detailedAnalysisRoutine_names)==list
        for name in detailedAnalysisRoutine_names:
            analysis_output[name] = np.zeros_like(scanVarRange)
    else:
        assert detailedAnalysisRoutine_names == None

    #Later we have assumed electron beam for emittance calculation
    assert baseSimSetup["BEAM"] == "e+" or baseSimSetup["BEAM"] == "e-"
    assert "ENERGY" in baseSimSetup or scanVar=="ENERGY"

    #Loading a pre-ran simulation?
    loadFileName = ""
    for key in sorted(baseSimSetup.keys()):
        value = baseSimSetup[key]
        loadFileName += str(key) + "=" + str(value) + "--"
    loadFileName += "SCANVAR=" + scanVar + '=['
    for var in scanVarRange:
        loadFileName += str(var) + ','
    loadFileName = loadFileName[:-1] + "]--"
    loadFileName += "COMMENT=" + COMMENT
    loadFileName = "SaveSim_" + loadFileName +".h5"

    print("LoadFile filename and status: '" + loadFileName + "'", tryLoad)
    if tryLoad:
        print("Loading...")
        try:
            loadFile = h5py.File(loadFileName,mode='r')

            scanVar_loaded = loadFile["scanVarName"]
            if scanVar != scanVar_loaded:
                print("Scan variables did not match even tough the filename did")
                print("please run with tryLoad=False to recompute.")
                loadFile.close()
                return

            scanVarRange_loaded = np.asarray(loadFile[scanVar])
            if len(scanVarRange)==len(scanVarRange_loaded) and np.all(np.equal(scanVarRange, scanVarRange_loaded)):
                print("Scan variable ranges match, let's load!")

                eps_x = np.asarray(loadFile["eps_x"])
                eps_y = np.asarray(loadFile["eps_y"])

                alpha_x = np.asarray(loadFile["alpha_x"])
                alpha_y = np.asarray(loadFile["alpha_y"])

                sigma_x = np.asarray(loadFile["sigma_x"])
                sigma_y = np.asarray(loadFile["sigma_y"])

                numPart = {}
                for pdg in PDG_keep:
                    if not "numPart_"+pdg in loadFile:
                        print("Could not find numPart for PDG={} in the file. Please recompute.".format(pdg))
                        loadFile.close()
                        return
                    numPart[pdg] = loadFile["numPart_"+pdg]

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
                return (eps_x,eps_y, beta_x,beta_y, alpha_x,alpha_y, sigma_x,sigma_y, numPart, analysis_output)

            else:
                print ("Scan variable ranges did not match, run with tryLoad=False to recompute.")
                loadFile.close()

                return
            return

        except OSError:
            print("File not found. Computing...")

    ## Build the job queue

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
            (emittances,numPart_singleSim,datafile) = miniScatterDriver.getData(filename=filenameROOTfile,quiet=QUIET,getRaw=True)
        else:
            with lock:
                print("Did not find file '{}', simulation crashed?".format(filenameROOTfile))
            badSim=True
            emittances = (float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'))
            numPart_singleSim = {}

        #Fill the raw emittance arrays
        gamma_rel = simSetup["ENERGY"]/0.511 #assume electron beam!
        beta_rel  = np.sqrt(gamma_rel**2 - 1.0) / gamma_rel;

        eps_x[i]   = emittances[0]
        beta_x[i]  = emittances[1]
        alpha_x[i] = emittances[2]
        sigma_x[i] = np.sqrt(eps_x[i]*beta_x[i]*1e6/(gamma_rel*beta_rel))
        eps_y[i]   = emittances[3]
        beta_y[i]  = emittances[4]
        alpha_y[i] = emittances[5]
        sigma_y[i] = np.sqrt(eps_y[i]*beta_y[i]*1e6/(gamma_rel*beta_rel))

        #Fill the NumPart array
        for pdg in numPart_singleSim.keys():
            if pdg in PDG_keep:
                numPart[pdg][i] = numPart_singleSim[pdg]
            else:
                print(pdg, numPart_singleSim)
                numPart['other'][i] += numPart_singleSim[pdg]
                with lock:
                    print("Found pdg={}".format(pdg))

        #Do special analysis over the TTrees
        if detailedAnalysisRoutine and not badSim:
            #Put the call to the external routine in a try/catch,
            # so that the thread will actually finish correctly in case of a user error.
            try:
                detailedData = detailedAnalysisRoutine(datafile)
                for k in detailedData.keys():
                    analysis_output[k][i]=detailedData[k]
            except Exception as err:
                with lock:
                    traceback.print_tb(err.__traceback__)

        if cleanROOT and not badSim:
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
    saveFile["scanVarName"] = scanVar
    saveFile[scanVar] = scanVarRange

    saveFile["eps_x"] = eps_x
    saveFile["eps_y"] = eps_y

    saveFile["alpha_x"] = alpha_x
    saveFile["alpha_y"] = alpha_y

    saveFile["sigma_x"] = sigma_x
    saveFile["sigma_y"] = sigma_y

    for pdg in PDG_keep:
        saveFile["numPart_"+str(pdg)] = numPart[pdg]

    if detailedAnalysisRoutine:
        for name in detailedAnalysisRoutine_names:
            nameMangle = "ANALYSIS_"+name
            saveFile["ANALYSIS_"+name] = analysis_output[name]

    saveFile.close()

    if detailedAnalysisRoutine:
        return (eps_x,eps_y, beta_x,beta_y, alpha_x,alpha_y, sigma_x,sigma_y, numPart, analysis_output)
    else:
        return (eps_x,eps_y, beta_x,beta_y, alpha_x,alpha_y, sigma_x,sigma_y, numPart)
