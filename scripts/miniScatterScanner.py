import miniScatterDriver

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

import numpy as np

import h5py
import ROOT

from queue import Queue
import threading

import sys
import traceback

import os

SEED = 1

#Used for particle-type countings and twiss, what to return
PDG_keep    = (11,-11,22,2212,2112,'other')
twissVar    = ('eps','beta','alpha', 'posAve','angAve', 'posVar','angVar','coVar')
def ScanMiniScatter(scanVar,scanVarRange,baseSimSetup, \
                    NUM_THREADS=4, tryLoad=False, COMMENT=None, QUIET=True, \
                    detailedAnalysisRoutine=None, detailedAnalysisRoutine_names=None, \
                    cleanROOT=True, getObjects=None, sameSeedValue=None, tmpFolder=None):
    """
    This routine is built to scan arbitrary parameters with MiniScatter.
    It can cache the results in HDF5-files with long and difficult names, as well as call detailed analysis routines.
    """

    global SEED # Updated every time one does a scan

    if COMMENT==None:
        raise ValueError("COMMENT argument is required.")

    # Initialize the scanVar in case of magnet
    scanVarMagnet = None
    if scanVar.startswith("MAGNET"):
        scanVarMagnet = scanVar.split(".")
        assert len(scanVarMagnet) > 1

    ### Create the output arrays ###
    
    twiss   = {} #Initialize on load
    numPart = {} #Initialize on load

    analysis_output = None
    if detailedAnalysisRoutine:
        analysis_output = {}
        assert type(detailedAnalysisRoutine_names)==list
        for name in detailedAnalysisRoutine_names:
            analysis_output[name] = np.zeros(len(scanVarRange),dtype=float)
    else:
        assert detailedAnalysisRoutine_names == None

    objects = None
    objectNames = None #Some objects can use SetName(), but we can still Write() them with a separate name...
    if getObjects:
        objects = {}
        objectNames = {}
        assert type(getObjects) == list
        for objName in getObjects:
            objects[objName]     = [None]*len(scanVarRange) #One entry per scanVarRange value
            objectNames[objName] = [None]*len(scanVarRange)

    #Later we have assumed electron beam for emittance calculation
    #assert baseSimSetup["BEAM"] == "e+" or baseSimSetup["BEAM"] == "e-" #(I think this was actually fixed?)
    assert "ENERGY" in baseSimSetup or scanVar=="ENERGY"

    #Sanity check
    if scanVar in baseSimSetup:
        print ("Please do not put scanVar in the baseSimSetup!")
        raise ValueError("Found scanVar in the baseSimSetup")

    if scanVarMagnet != None:
        assert "MAGNET" in baseSimSetup

        if not len(scanVarMagnet[0]) > 6:
            raise ValueError("When parsing scanvar MAGNETi.arg; expected"+\
                             " an integer i (from 0 to len(baseSimSetup['MAGNET']) indicating which magnet.")
        whichMagnet   = int(scanVarMagnet[0][6:])

        if len(scanVarMagnet) == 2 and scanVarMagnet[1] in baseSimSetup["MAGNET"][whichMagnet]:
            print ("Please do not put scanVar in the baseSimSetup['MAGNET'][{}]".\
                   format(whichMagnet))
            raise ValueError("Found scanVar in the baseSimSetup['MAGNET'][{}]".\
                             format(whichMagnet))

        elif len(scanVarMagnet) == 3 and scanVarMagnet[2] in baseSimSetup["MAGNET"][whichMagnet]["keyval"]:
            print ("Please do not put scanVar in the baseSimSetup['MAGNET'][{}]['keyval']".\
                   format(whichMagnet))
            raise ValueError("Found scanVar in the baseSimSetup['MAGNET'][{}]['keyval']".\
                             format(whichMagnet))

    #Loading a pre-ran simulation?
    loadFileName = "SaveSim_{}_{}.h5".format(scanVar,COMMENT)
    print ("LoadFile filename and status: '" + loadFileName + "'", tryLoad)
    if tryLoad:
        print ("Loading...")
        try:
            loadFile = h5py.File(loadFileName,mode='r')

            #Check that the file is workable
            scanVar_loaded = loadFile.attrs["scanVarName"]
            if scanVar != scanVar_loaded:
                print ("Scan variables did not match even tough the filename did")
                print ("please run with tryLoad=False to recompute.")
                loadFile.close()
                raise ValueError("scanVar did not match loaded file")

            # Check that all the keys in the basesimsetup are in the loadFile
            for key in baseSimSetup.keys():
                if key == "MAGNET":
                    #Handle this separately
                    continue
                if not (key in loadFile.attrs):
                    print ("Key '{}' found in baseSimSetup but not in the file.".format(key))
                    print ("Please run with tryLoad=False to recompute.")
                    loadFile.close()
                    raise ValueError("Found key {} in baseSimSetup but not in the file".format(key))

                if key == "COVAR":
                    if not (np.all(loadFile.attrs[key] == baseSimSetup[key])):
                        print ("Value of key '{}' found in baseSimSetup".format(key))
                        print (" did not match what was found in the file.")
                        print ("Values: baseSimSetup={}, file={}".format(baseSimSetup[key],loadFile.attrs[key]))
                        print ("Please run with tryLoad=False to recompute.")
                        loadFile.close()
                        raise ValueError("Value of key {} did not match with baseSimSetup and file".format(key))
                    continue

                if not loadFile.attrs[key] == baseSimSetup[key]:
                    print ("Value of key '{}' found in baseSimSetup".format(key))
                    print (" did not match what was found in the file.")
                    print ("Values: baseSimSetup={}, file={}".format(baseSimSetup[key],loadFile.attrs[key]))
                    print ("Please run with tryLoad=False to recompute.")
                    loadFile.close()
                    raise ValueError("Value of key {} did not match with baseSimSetup and file".format(key))

            if "MAGNET" in baseSimSetup.keys():
                magnetCounter = 0
                for mag in baseSimSetup["MAGNET"]:

                    for key,val in mag.items():

                        if key == "keyval":
                            #another level
                            for key2,val2 in val.items():
                                mkey = "MAGNET"+str(magnetCounter)+".keyval."+key2
                                if not (mkey in loadFile.attrs):
                                    print (("Key '{}' found in baseSimSetup['MAGNET'][{}]['keyval'] "+\
                                           "but not in the file.").format(mkey,magnetCounter))
                                    print ("Please run with tryLoad=False to recompute.")
                                    loadFile.close()
                                    raise ValueError(("Found key {} in "+\
                                                      "baseSimSetup['MAGNET'][{}]['keyval'] "+\
                                                      "but not in the file").\
                                                     format(mkey,magnetCounter))
                                if not (loadFile.attrs[mkey] == val2):
                                    print (("Value of key '{}' found in "+\
                                            "baseSimSetup['MAGNET'][{}]['keyval']").\
                                           format(mkey,magnetCounter))
                                    print (" did not match what was found in the file.")
                                    print ("Values: baseSimSetup={}, file={}".\
                                           format(val2,loadFile.attrs[mkey]))
                                    print ("Please run with tryLoad=False to recompute.")
                                    loadFile.close()
                                    raise ValueError(("Value of key {} did not match "+\
                                                      "with baseSimSetup['MAGNET'][{}]['keyval'] "+\
                                                      "and file").format(mkey,magnetCounter))
                            #go on to the next key!
                            continue

                        mkey = "MAGNET"+str(magnetCounter)+"."+key
                        if not (mkey in loadFile.attrs):
                            print (("Key '{}' found in baseSimSetup['MAGNET'][{}] "+\
                                    "but not in the file.").\
                                   format(mkey,magnetCounter))
                            print ("Please run with tryLoad=False to recompute.")
                            loadFile.close()
                            raise ValueError(("Found key {} in baseSimSetup['MAGNET'][{}] "+\
                                              "but not in the file").format(mkey,magnetCounter))

                        if not (loadFile.attrs[mkey] == val):
                            print ("Value of key '{}' found in baseSimSetup['MAGNET'][{}]".\
                                   format(mkey,magnetCounter))
                            print (" did not match what was found in the file.")
                            print ("Values: baseSimSetup={}, file={}".format(val,loadFile.attrs[mkey]))
                            print ("Please run with tryLoad=False to recompute.")
                            loadFile.close()
                            raise ValueError(("Value of key {} did not match with "+\
                                              "baseSimSetup['MAGNET'][{}] "+\
                                              "and file").format(mkey,magnetCounter))

                    #END loop over a magnet's items

                    magnetCounter += 1
                #END loop over magnets

            # Check that all the keys in the loadFile are in the basesimsetup
            for key in loadFile.attrs:
                if key == scanVar or key=='scanVarName' or key=='objectsFileName':
                    continue
                if key.startswith("MAGNET"):
                    mkeyS = key.split(".")
                    magnetCounter = int(mkeyS[0][6:])
                    if len(mkeyS) > 2:
                           assert mkeyS[1] == "keyval"
                           if not mkeyS[2] in baseSimSetup["MAGNET"][magnetCounter]["keyval"].keys():
                               print (("Key '{}' (sub-part '{}') found in file "+\
                                       "but not in baseSimSetup['MAGNET'][{}]['keyval]'").\
                                      format(key,mkeyS[2],magnetCounter))
                               print ("Please run with tryLoad=False to recompute.")
                               loadFile.close()
                               raise ValueError(("Key {} found in file but not in "+\
                                                 "baseSimSetup['MAGNET'][{}]['keyval']").\
                                                format(key,magnetCounter))
                    else:
                        if not mkeyS[1] in baseSimSetup["MAGNET"][magnetCounter].keys():
                           print (("Key '{}' (sub-part '{}') found in file "+\
                                   "but not in baseSimSetup['MAGNET'][{}]").\
                                  format(key,mkeyS[1],magnetCounter))
                           print ("Please run with tryLoad=False to recompute.")
                           loadFile.close()
                           raise ValueError(("Key {} found in file but not in "+\
                                             "baseSimSetup['MAGNET'][{}]").\
                                            format(key,magnetCounter))
                    continue
                if not key in baseSimSetup.keys():
                    print ("Key {} found in file but not in baseSimSetup.".format(key))
                    print ("Please run with tryLoad=False to recompute.")
                    loadFile.close()
                    raise ValueError("Key {} found in file but not in baseSimSetup".format(key))

            scanVarRange_loaded = np.asarray(loadFile.attrs[scanVar])
            scanVarRange_safe = scanVarRange
            if type(scanVarRange[0]) == str:
                #Python3, it's unicode. HDF5 dislikes unicode arrays
                scanVarRange_safe = []
                for scanVarRange_iter in scanVarRange:
                    scanVarRange_safe.append( scanVarRange_iter.encode('ascii') )

            scanVarRange_sameLen = len(scanVarRange) != len(scanVarRange_loaded)
            scanVarRange_sameContent = False
            if scanVarRange_sameLen:
                scanVarRange_sameContent = True
                for i in range(len(scanVarRange_safe)):
                    if scanVarRange_safe[i] != scanVarRange_loaded[i]:
                        scanVarRange_sameContent = False

            if scanVarRange_sameContent:
                print ("Scan variable ranges did not match, run with tryLoad=False to recompute.")
                print ("Now :", scanVarRange)
                print ("File:", scanVarRange_loaded)
                # for i in range(min(len(scanVarRange),len(scanVarRange_loaded))):
                #     print (  scanVarRange[i], scanVarRange_loaded[i],\
                #              scanVarRange[i]-scanVarRange_loaded[i], \
                #            ( scanVarRange[i]-scanVarRange_loaded[i])==0.0)
                loadFile.close()
                raise ValueError("ScanVarRange did not match with loaded file")

            print ("Scan variable ranges match, let's load!")

            #Load Twiss data
            for key in loadFile:
                #key format: twiss_DET_XY_VAR, where DET is e.g. 'tracker_cutoff', XY is 'x' or 'y', and VAR is one of twissVar
                if not key.startswith('twiss_'):
                    continue

                xy = None
                if '_x_' in key:
                    xy = 'x'
                elif '_y_' in key:
                    xy = 'y'
                else:
                    loadFile.close()
                    raise AssertionError("Found neither _x_ or _y_ in key='"+key+"' while loading twiss data from file")
                
                key_short = key[6:]
                key_split = key_short.split("_"+xy+"_")
                assert len(key_split) == 2, "Got: '" + key + ", " + key_short + ", " + str(key_split)
                det = key_split[0]
                var = key_split[1]

                if not det in twiss:
                    twiss[det] = {}
                    twiss[det]['x'] = {}
                    twiss[det]['y'] = {}

                if(var in twiss[det][xy]):
                    loadFile.close()
                    raise KeyError('Variable {} was already loaded for output plane {}, xy={}'.format(var,det,xy))

                twiss[det][xy][var] = np.asarray(loadFile[key])

            #TWISS: Check that we got everything (but not more) for every output plane and and x/y
            for det in twiss:
                for xy in twiss[det]:
                    for var in twissVar:
                        if not var in twiss[det][xy]:
                            loadFile.close()
                            raise KeyError("Did not successfully load variable '"+var+"' from det='"+det+"', xy='"+xy+"'")
                    for var in twiss[det][xy]:
                        if not var in twissVar:
                            loadFile.close()
                            raise KeyError("Loaded unknown variable '"+var+"' from det='"+det+"', xy='"+xy+"'")

            #Load numPart data
            for key in loadFile:
                #key format: numPart_DET_PDG, where DET is e.g. target_cutoff and VAR is one of the ID's in PDG_keep
                if not key.startswith('numPart_'):
                    continue

                key_short = key[8:]

                PDG_wasfound = False
                for pdg in PDG_keep:
                    if not key_short.endswith("_"+str(pdg)):
                        continue
                    PDG_wasfound = True
                    det = key_short[:-(len(str(pdg))+1)] #Remove the PDG id from the key_short to find the detector name
                    if not det in numPart:
                        #Initialize numpart[det]
                        numPart[det] = {}
                        for p in PDG_keep:
                            numPart[det][p] = None

                    numPart[det][pdg]=np.asarray(loadFile[key])
                    break
                if not PDG_wasfound:
                    loadFile.close()
                    raise KeyError("numPart PDG not recognized in key='{}' from the file".format(key))

            #NumPart: Check that we got everything for every output plane
            for det in numPart:
                for PDG in PDG_keep:
                    if not PDG in numPart[det]:
                        loadFile.close()
                        raise AssertionError("Internal error, PDG={} not initialized for numPart[{}]".format(PDG,det))
                    if numPart[det] is None:
                        loadFile.close()
                        raise KeyError("Did not successfully load PDG={} for det={}".format(PDG,det))

            if detailedAnalysisRoutine:
                for name in detailedAnalysisRoutine_names:
                    nameMangle = "ANALYSIS_"+name
                    if not (nameMangle) in loadFile:
                        print ("Could not find '"+name+"' in the file. Please recompute.")
                        loadFile.close()
                        return
                    analysis_output[name] = np.asarray(loadFile["ANALYSIS_"+name])

            if getObjects:
                if not "objectsFileName" in loadFile.attrs:
                    raise ValueError("No objectsFileName in loadFile={}".format(loadFileName))
                objectsFileName = str(loadFile.attrs["objectsFileName"])
                objectsFile = ROOT.TFile(objectsFileName, 'READ')

                #Load the objects from file
                for objName in getObjects:
                    for i in range(len(scanVarRange)):
                        thisObjName = objName + "_" + scanVar + "_" + str(i) + "_" + COMMENT + "-fileClone"
                        if not objectsFile.GetListOfKeys().Contains(thisObjName):
                            objectsFile.Close()
                            loadFile.close()
                            raise KeyError("Did not find key '{}' in ROOT file '{}'".format(thisObjName,objectsFileName))

                        thisObj = objectsFile.Get(thisObjName)
                        objects[objName][i] = thisObj.Clone(thisObj.GetName()+ "-clone")
                        try:
                            objects[objName][i].SetDirectory(0)
                        except AttributeError:
                            #Doesn't work on TVector etc.
                            print("WARNING while loading aux ROOT file '{}': Can't SetDirectory() on object '{}'.".format(objectsFile,objName))

                objectsFile.Close()

                #Fix the names
                for objName in getObjects:
                    assert len(objects[objName]) == len(scanVarRange)
                    for i in range(len(scanVarRange)):
                        thisObjName = objects[objName][i].GetName()
                        try:
                            # remove '-fileClone-clone'
                            objects[objName][i].SetName(thisObjName[:-16])
                        except AttributeError:
                            #Doesn't work on TVectorT etc.
                            print("WARNING while loading aux ROOT file '{}': Can't SetName() on object '{}'.".format(objectsFile,objName))
                print ("Auxillary ROOT file {} loaded.".format(objectsFileName))

            loadFile.close()
            print ("Loaded! That was fast.")
            return (twiss, numPart, objects, analysis_output)

        except OSError:
            print ("File not found. Computing...")

    ### Build the job queue ###
    def computeOnePoint(var,simIdx,lock):
        with lock:
            print ("{} = {} ({}/{})".format(scanVar, var, simIdx+1, len(scanVarRange)))

        #Run the simulation -- this MUST be done in parallel
        filenameROOT = 'output_'+scanVar+"="+str(var)
        simSetup = baseSimSetup.copy()
        if "MAGNET" in baseSimSetup.keys():
            # It needs to be a deep copy!
            simSetup["MAGNET"] = []
            for basemag in baseSimSetup["MAGNET"]:
                simSetup["MAGNET"].append(basemag.copy())
                if "keyval" in basemag:
                    simSetup["MAGNET"][-1]["keyval"] = basemag['keyval'].copy()

        if scanVar.startswith("MAGNET"):
            scanVarMag = scanVar.split(".")
            magnetCounter = int(scanVarMag[0][6:])
            if len(scanVarMag) == 2:
                simSetup["MAGNET"][magnetCounter][scanVarMag[1]] = var
            elif len(scanVarMag) == 3:
                assert scanVarMag[1] == 'keyval'
                simSetup["MAGNET"][magnetCounter]['keyval'][scanVarMag[2]] = var
            else:
                raise ValueError("Expected len(scanVarMag) == 2 or 3")
        else:
            simSetup[scanVar]   = var

        if sameSeedValue == None:
            simSetup["SEED"]    = SEED + simIdx
        else:
            assert type(sameSeedValue) == int
            simSetup["SEED"]    = sameSeedValue
            print("WARNING: Running all simulations with same seed = ", simSetup["SEED"])

        simSetup["OUTNAME"] = filenameROOT

        assert not "OUTFOLDER" in simSetup
        if tmpFolder is not None:
            simSetup["OUTFOLDER"] = tmpFolder

        ##RUN THE SIMULATION!
        miniScatterDriver.runScatter(simSetup,quiet=QUIET, onlyCommand=True)
        ##SIMULATION COMPLETE!

        filenameROOTfile = None
        if tmpFolder is None:
            runFolder = os.path.dirname(os.path.abspath(__file__))
            filenameROOTfile = os.path.join(runFolder,"plots",filenameROOT+".root")
        else:
            filenameROOTfile = os.path.join(tmpFolder,filenameROOT+".root")

        badSim=False

        ## Extract the data
        # grab the lock since root histograms tend to have identical standard names
        # when fresh off the file => occational crashes.
        # Also we are manpulating the twiss and numPart arrays.
        with lock:
            if os.path.isfile(filenameROOTfile):
                #Always getRaw, and handle the cleanup here in Scanner.
                (twiss_singleSim, numPart_singleSim, objects_singleSim, datafile) = \
                    miniScatterDriver.getData(filename=filenameROOTfile,quiet=QUIET,getRaw=True,getObjects=getObjects)
            else:
                #with lock:
                print ("Did not find file '{}', simulation crashed?".format(filenameROOTfile))
                badSim=True

            
            if not badSim:

                #Fill the twiss array
                for det in twiss_singleSim:
                    #Initialize the detector if needed
                    if not det in twiss:
                        twiss[det] = {}
                        twiss[det]['x'] = {}
                        twiss[det]['y'] = {}
                        for var in twissVar:
                            twiss[det]['x'][var] = np.zeros(len(scanVarRange),dtype=float)
                            twiss[det]['y'][var] = np.zeros(len(scanVarRange),dtype=float)
                    #Load data from this simulation
                    for xy in ('x','y'):
                        for var in twiss_singleSim[det][xy]:
                            twiss[det][xy][var][simIdx] = twiss_singleSim[det][xy][var]

                #Fill the NumPart array
                for det in numPart_singleSim:
                    #Initialize the detector if needed
                    if not det in numPart:
                        numPart[det] = {}
                        for pdg in PDG_keep:
                            numPart[det][pdg] = np.zeros(len(scanVarRange),dtype=int)
                    #Load data from this simulation
                    for pdg in numPart_singleSim[det].keys():
                        if pdg in PDG_keep:
                            numPart[det][pdg][simIdx] = numPart_singleSim[det][pdg]
                        else:
                            numPart[det]['other'][simIdx] += numPart_singleSim[det][pdg]
                            print ("Found pdg={} for detector={}".format(pdg,det))

                #File the objects into the appropriate positions
                if getObjects:
                    for objName in getObjects:
                        thisObj = objects_singleSim[objName]
                        #thisObjName = thisObj.GetName()
                        thisObjName = objName + "_" + scanVar + "_" + str(simIdx) + "_" + COMMENT
                        try:
                            thisObj.SetName(thisObjName)
                        except AttributeError:
                            #Doesn't work on TVectorT etc.
                            print("WARNING while retrieving object from single-sim ROOT file {}: Can't SetName() on object '{}'.".format(filenameROOTfile,objName))

                        objects[objName][simIdx]     = thisObj
                        objectNames[objName][simIdx] = thisObjName

                #If requested: do special analysis over the TTrees
                if detailedAnalysisRoutine:
                    #Put the call to the external routine in a try/catch,
                    # so that the thread will actually finish correctly in case of a user error.
                    try:
                        detailedData = detailedAnalysisRoutine(datafile)
                        for k in detailedData.keys():
                            analysis_output[k][simIdx]=detailedData[k]
                    except Exception as err:
                        traceback.print_tb(err.__traceback__)

                # Cleanup the ROOT file
                datafile.Close()
                if cleanROOT:
                    os.remove(filenameROOTfile)
                    if not QUIET:
                        #with lock:
                        print ("Deleting '{}'.".format(filenameROOTfile))

    def threadWorker(jobQueue_local,lock):
        while not jobQueue_local.empty():
            (var,simIdx) = jobQueue_local.get()
            computeOnePoint(var,simIdx,lock)
            jobQueue_local.task_done()

    jobQueue = Queue(0)
    simIdx=0 # Array index where the jobs should write their data
    for var in scanVarRange:
        jobQueue.put((var,simIdx))
        simIdx = simIdx+1
    lock = threading.Lock() # Shared lock to synchronize output
    for i in range(NUM_THREADS):
        worker = threading.Thread(target=threadWorker, args=(jobQueue,lock))
        worker.start()

    jobQueue.join()
    if sameSeedValue != None:
        SEED = SEED+i #Prepare the SEED for the next run of simulations

    print ("Simulation complete, saving data to h5 for later retrival.")
    #Write out the data
    saveFile = h5py.File(loadFileName,mode="w")

    saveFile.attrs["scanVarName"] = scanVar
    if type(scanVarRange[0]) == str:
        #Python3, it's unicode. HDF5 dislikes unicode arrays
        scanVarRange_ascii = []
        for scanVarRange_iter in scanVarRange:
            scanVarRange_ascii.append( scanVarRange_iter.encode('ascii') )
        saveFile.attrs[scanVar] = scanVarRange_ascii
    else:
        saveFile.attrs[scanVar] = scanVarRange
    for key in baseSimSetup.keys():
        if key in saveFile.attrs:
            saveFile.close()
            raise RuntimeError("Attribute '{}' already written to file? Probably a bug!".format(key))

        if key == "MAGNET":
            magnetCounter = 0
            for mag in baseSimSetup["MAGNET"]:
                for key2 in mag:
                    if key2 == "keyval":
                        for key3 in mag[key2].keys():
                            mkey = "MAGNET"+str(magnetCounter)+".keyval."+key3
                            if mkey in saveFile.attrs:
                                saveFile.close()
                                raise RuntimeError\
                                    ("Attribute '{}' already written to file? Probably a bug!".\
                                     format(mkey))
                            value = mag[key2][key3]
                            saveFile.attrs[mkey] = value
                    else:
                        mkey = "MAGNET"+str(magnetCounter)+"."+key2
                        if mkey in saveFile.attrs:
                            saveFile.close()
                            raise RuntimeError\
                                ("Attribute '{}' already written to file? Probably a bug!".\
                                 format(mkey))
                        value = mag[key2]
                        saveFile.attrs[mkey] = value
                magnetCounter += 1
        else:
            value = baseSimSetup[key]
            saveFile.attrs[key] = value

    for det in twiss:
        for xy in twiss[det]:
            for var in twiss[det][xy]:
                saveFile["twiss_"+det+"_"+xy+"_"+var] = twiss[det][xy][var]

    for det in numPart:
        for pdg in numPart[det]: #Save'em all, don't rely on PDGkeep
            saveFile["numPart_"+det+"_"+str(pdg)] = numPart[det][pdg]

    if detailedAnalysisRoutine:
        for name in detailedAnalysisRoutine_names:
            nameMangle = "ANALYSIS_"+name
            saveFile["ANALYSIS_"+name] = analysis_output[name]

    if getObjects:
        objectsFileName = loadFileName[:-3] + ".root" #remove the .h5
        saveFile.attrs["objectsFileName"] = objectsFileName
        objectsFile = ROOT.TFile(objectsFileName,'RECREATE')
        for objName in getObjects:
            for i in range(len(scanVarRange)):
                #thisObjName = objects[objName][i].GetName() + "-fileClone"
                thisObjName = objectNames[objName][i] + "-fileClone"
                objCopy = objects[objName][i].Clone(thisObjName)
                objCopy.Write(thisObjName)
        objectsFile.Close()

    saveFile.close()

    return (twiss, numPart, objects, analysis_output)
