/* MiniScatter.cc:
 * Helga Holmestad 2015-2019
 * Kyrre Sjobak    2015-
 * Eric Fackelman  2022
 *
 * This file is part of MiniScatter.
 *
 *  MiniScatter is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MiniScatter is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MiniScatter.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "DetectorConstruction.hh"
#include "MagnetSensorWorldConstruction.hh"
#include "VirtualTrackerWorldConstruction.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"

#include "G4PhysListFactory.hh"
#include "G4ParallelWorldPhysics.hh"

#include "G4Version.hh"
#include "G4Exception.hh"
#if G4VERSION_NUMBER >= 1060
//These macros were removed in version 10.6;
// from what I can see it looks like they are are no longer needed
#define G4VIS_USE
#define G4UI_USE
#endif

#include "MiniScatterVersion.hh"

#include "RootFileWriter.hh"

#include "G4SystemOfUnits.hh"
#include "G4String.hh"
#include <string> //C++11 std::stoi

#include "TROOT.h"


#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//#include <unistd.h> //getopt()
#include <getopt.h> // Long options to getopt (GNU extension)

void printHelp(G4double target_thick,
               G4String target_material,
               G4String background_material,
               std::vector<G4double>* detector_distances,
               G4double detector_angle,
               G4double target_angle,
               G4double world_size,
               G4String physListName,
               G4double physCutoffDist,
               G4int    numEvents,
               G4double beam_energy,
               G4double beam_eFlat_min,
               G4double beam_eFlat_max,
               G4String beam_type,
               G4double beam_offset,
               G4double beam_zpos,
               G4double beam_angle,
               G4String covarianceString,
               G4double beam_rCut,
               G4bool   doBacktrack,
               G4String beam_loadfile,
               G4int    rngSeed,
               G4String filename_out,
               G4String foldername_out,
               G4bool   quickmode,
               G4bool   anaScatterTest,
               G4bool   miniROOTfile,
               G4double cutoff_energyFraction,
               G4double cutoff_radius,
               G4double edep_dens_dz,
               G4int    engNbins,
               G4double histPosLim,
               G4double histAngLim,
               std::vector<G4String> &magnetDefinitions);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

    // Avoid problem (segfault)
    // "Error in <UnknownClass::InitInterpreter()>: LLVM SYMBOLS ARE EXPOSED TO CLING! This will cause problems; please hide them or dlopen() them after the call to TROOT::InitInterpreter()!"
    // when running with GUI.
    gROOT->Reset();

    G4cout << "**********************************************************" << G4endl
           << "**********************************************************" << G4endl
           << G4endl;
    G4cout << " Welcome to MiniScatter, release version " << miniscatter_version << "!" << G4endl
           << " Release dated: " << miniscatter_date << G4endl
           << G4endl;
    G4cout << " MiniScatter is controlled with command line arguments." << G4endl
           << " Run MiniScatter --help to see all possible arguments,"
           << "  a description of each of them, and their default/current values." << G4endl
           << G4endl;
    G4cout << " If publishing results obtained with this code, " << G4endl
           << "  please cite K.Sjobak and H.Holmestad, 'MINISCATTER, A SIMPLE GEANT4 WRAPPER'," << G4endl
           << "  https://dx.doi.org/10.18429/JACoW-IPAC2019-WEPTS025" << G4endl
           << G4endl;

    G4cout << "**********************************************************" << G4endl
           << "**********************************************************" << G4endl
           << G4endl;

    //Parse command line arguments
    int getopt_char;
    int getopt_idx;

    G4double target_thick        = 1.0;           // Target thickness [mm]; 0.0 for no target slab (only magnets)
    G4String target_material     = "G4_Al";       // Name of target material to use
    G4double target_angle        = 0.0;           // Target angle around y-axis [deg]

    G4String background_material = "G4_Galactic"; // Background material

    std::vector<G4double> detector_distances;     // Detector distances at x=y=0  [mm]
    detector_distances.push_back(50.0);
    G4double detector_angle      = 0.0;           // Detectector angle around y-axis [deg]

    G4double world_size          = 500.0;        // World full size X/Y [mm]
                                                 // (not half-size as taken as argument by G4Box etc.)

    G4double beam_energy = 200;                   // Beam kinetic energy [MeV]
    G4double beam_eFlat_min = -1.0;               // For flat-spectrum energy distribution,
    G4double beam_eFlat_max = -1.0;               //  set these to min/max, both > 0.0.

    G4String beam_type   = "e-";                  // Beam particle type
    G4double beam_offset = 0.0;                   // Beam offset (x) [mm]
    G4double beam_zpos   = 0.0;                   // Beam offset (z) [mm]
    G4bool   doBacktrack = false;                 // Backtrack to the z-position?
    G4double beam_angle  = 0.0;                   // Beam angle [deg]
    G4String covarianceString = "";               // Beam covariance matrix parameters
    G4double beam_rCut = 0.0;                     // Beam distribution radial cutoff

    G4String beam_loadFile = "";                  // Filename to load beam distribution from

    G4String physListName = "QGSP_FTFP_BERT";     // Name of physics list to use
    G4double physCutoffDist = 0.1;                // Default physics cutoff distance [mm]

    G4int    numEvents    = 0;                    // Number of events to generate

    G4bool   useGUI         = false;              // GUI on/off
    G4bool   quickmode      = false;              // Don't make slow plots
    G4bool   anaScatterTest = false;              // Compute analytical multiple scattering,
                                                  // for comparison with simulation output.
    G4String filename_out   = "output";           // Output filename
    G4String foldername_out = "plots";            // Output foldername
    G4bool   miniROOTfile   = false;              // Write small root file
                                                  // (only analysis output, no TTrees)

    G4int    rngSeed        = 0;                  // RNG seed

    G4double cutoff_energyFraction = 0.95;        // [fraction]
    G4double cutoff_radius         = world_size;  // [mm]

    G4double edep_dens_dz          = 0.0;         // Z bin width for energy deposit histograms [mm]
    G4int    engNbins              = 0;           // Number of bins for the 1D energy histograms

    G4double histPosLim            = 10.0;        // Position histogram Limit, default 10
    G4double histAngLim            = 5.0;         // Angle histogram Limit, default 5

    std::vector<G4String> magnetDefinitions;

    static struct option long_options[] = {
                                           {"thick",                 required_argument, NULL, 't'  },
                                           {"mat",                   required_argument, NULL, 'm'  },
                                           {"backgroundMaterial",    required_argument, NULL, 1600 },
                                           {"dist",                  required_argument, NULL, 'd'  },
                                           {"ang",                   required_argument, NULL, 'a'  },
                                           {"targetAngle",           required_argument, NULL, 'A'  },
                                           {"worldsize",             required_argument, NULL, 'w'  },
                                           {"ang",                   required_argument, NULL, 'a'  },
                                           {"phys",                  required_argument, NULL, 'p'  },
                                           {"physCutoffDist",        required_argument, NULL, 1400 },
                                           {"numEvents",             required_argument, NULL, 'n'  },
                                           {"energy",                required_argument, NULL, 'e'  },
                                           {"energyDistFlat",        required_argument, NULL, 1300 },
                                           {"beam",                  required_argument, NULL, 'b'  },
                                           {"xoffset",               required_argument, NULL, 'x'  },
                                           {"zoffset",               required_argument, NULL, 'z'  },
                                           {"beamAngle",             required_argument, NULL, 1500 },
                                           {"covar",                 required_argument, NULL, 'c'  },
                                           {"beamRcut",              required_argument, NULL, 1200 },
                                           {"beamFile",              required_argument, NULL, 1700 },
                                           {"outname",               required_argument, NULL, 'f'  },
                                           {"outfolder",             required_argument, NULL, 'o'  },
                                           {"seed",                  required_argument, NULL, 's'  },
                                           {"help",                  no_argument,       NULL, 'h'  },
                                           {"gui",                   no_argument,       NULL, 'g'  },
                                           {"quickmode",             no_argument,       NULL, 'q'  },
                                           {"anaScatterTest",        no_argument,       NULL, 1004 },
                                           {"miniroot",              no_argument,       NULL, 'r'  },
                                           {"cutoffEnergyFraction",  required_argument, NULL, 1000 },
                                           {"cutoffRadius",          required_argument, NULL, 1001 },
                                           {"edepDZ",                required_argument, NULL, 1002 },
                                           {"engNbins",              required_argument, NULL, 1003 },
                                           {"histPosLim",            required_argument, NULL, 1005 }, 
                                           {"histAngLim",            required_argument, NULL, 1006 },
                                           {"magnet",                required_argument, NULL, 1100 },
                                           {"object",                required_argument, NULL, 1100 }, //synonymous with --magnet
                                           {0,0,0,0}
    };

    while ( (getopt_char = getopt_long(argc,argv, "t:m:d:a:A:w:p:n:e:b:x:z:c:f:o:s:hgqr", long_options, &getopt_idx)) != -1) {
        switch(getopt_char) {
        case 'h': //--help ; Print help
            printHelp(target_thick,
                      target_material,
                      background_material,
                      &detector_distances,
                      detector_angle,
                      target_angle,
                      world_size,
                      physListName,
                      physCutoffDist,
                      numEvents,
                      beam_energy,
                      beam_eFlat_min,
                      beam_eFlat_max,
                      beam_type,
                      beam_offset,
                      beam_zpos,
                      beam_angle,
                      covarianceString,
                      beam_rCut,
                      doBacktrack,
                      beam_loadFile,
                      rngSeed,
                      filename_out,
                      foldername_out,
                      quickmode,
                      anaScatterTest,
                      miniROOTfile,
                      cutoff_energyFraction,
                      cutoff_radius,
                      edep_dens_dz,
                      engNbins,
                      histPosLim,
                      histAngLim,
                      magnetDefinitions);
            exit(1);
            break;

        case 'g': //--gui ; Use the GUI
            useGUI = true;
            break;

        case 't': //--thick ; Target thickness
            try {
                target_thick = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading target thickness" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
            break;

        case 'd': //--dist ; Detector distance(s)
            { //Scope to avoid spilling temp variables
                detector_distances.clear();
                G4String detectorString = std::string(optarg);
                if (detectorString == "NONE") break; //Just clear it

                //Split by :
                str_size startPos = 0;
                str_size endPos = 0;
                do {
                    endPos          = detectorString.index(":",startPos);
                    G4String detStr = detectorString(startPos,endPos-startPos);

                    G4double dist;
                    try {
                        dist = std::stod(detStr);
                    }
                    catch (const std::invalid_argument& ia) {
                        G4cout << "Invalid argument when reading detector distance" << G4endl
                            << "Got: '" << detStr << "'" << G4endl
                            << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                        exit(1);
                    }
                    detector_distances.push_back(dist);

                    startPos = endPos+1;

                } while (endPos != std::string::npos);

            }
            break;

        case 'a': //--angle ; Detector angle
            try {
                detector_angle = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading detector angle" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
            break;

        case 'A': //--targetAngle ; Target angle
            try {
                target_angle = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading target angle" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
            break;

        case 'w': //--worldsize ; World size
            try {
                world_size = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading world size" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
            break;

        case 'p': //--phys ; Named physics list
            physListName = G4String(optarg);
            break;

        case 1400: //--physCutoffDist ; Physics cutoff distance (--physCutoffDist)
            try {
                physCutoffDist = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading physCutoffDist" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
            break;

        case 'm': //--mat ; Target material
            target_material = G4String(optarg);
            break;

        case 1600: //--backgroundMaterial ; Background material
            background_material = G4String(optarg);
            break;

        case 'n': //--numEvents ; Number of events
            try {
                numEvents = std::stoi(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading number of events" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected an integer!" << G4endl;
                exit(1);
            }
            break;

        case 'e': //--energy ; The (reference) beam kinetic energy
            try {
                beam_energy = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading beam kinetic energy" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
            break;

        case 1300: {//--energyDistFlat ; beam energy (flat distribution)
            G4String edist_str = G4String(optarg);

            str_size startPos = 0;
            str_size endPos   = edist_str.index(":",startPos);
            if (endPos == std::string::npos) {
                G4cout << " Error while searching for ':' in edist_str = "
                       << edist_str << "', did not find?" << G4endl;
                exit(1);
            }

            try {
                beam_eFlat_min = std::stod(string(edist_str(0,endPos)));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading minimum energy '" << edist_str(0,endPos) << "' for eFlat." << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }

            startPos = endPos+1;
            endPos   = edist_str.index(":",startPos);
            if (endPos != std::string::npos) {
                G4cout << " Error while searching for ':' in edist_str = "
                       << edist_str << "', found a second one?" << G4endl;
                exit(1);
            }

            try {
                beam_eFlat_max = std::stod(string(edist_str(startPos,endPos-startPos)));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading maximum energy '" << edist_str(startPos,endPos-startPos) << "' for eFlat." << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }

            if (beam_eFlat_min < 0 or beam_eFlat_max <= 0 or beam_eFlat_min >= beam_eFlat_max) {
                G4cout << "Invalid eFlat_min/_max" << G4endl;
                exit(1);
            }

            }
            break;

        case 'b': //--beam ; Beam particle type
            beam_type = G4String(optarg);
            break;

        case 'x': //--xoffset / Beam offset (x)
            try {
                beam_offset = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading beam offset (x)" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
            break;

        case 'z': //--zoffset ; Beam offset (z)
            if (strlen(optarg) > 0 && optarg[0] == '*') {
                doBacktrack = true;
            }
            try {
                if (doBacktrack) {
                    beam_zpos = std::stod(string(optarg).substr(1));
                }
                else {
                    beam_zpos = std::stod(string(optarg));
                }
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading beam offset (z)" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation and prepended '*' is accepted)" << G4endl;
                exit(1);
            }
            break;

        case 1500: //--beamAngle ; Initial angle of beam particles
            try {
                beam_angle = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading beamAngle" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
            break;
        
        case 'c': //--covar ; Beam covariance matrix from Twiss parameters
            //-c epsN[um]:beta[m]:alpha(::epsN_Y[um]:betaY[m]:alphaY)
            covarianceString = G4String(optarg);
            break;

        case 1200: //--beamRcut ; Beam radial cutoff [mm]
            try {
                beam_rCut = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading beam_rCut" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
            break;

        case 1700: //--beamFile ; Input filename to load beam from
            beam_loadFile = G4String(optarg);
            break;

        case 'f': //--outname ; Output filename
            filename_out = G4String(optarg);
            break;

        case 'o': //--outname ; Output foldername
            foldername_out = G4String(optarg);
            break;

        case 'r': //--miniroot ; Mini ROOT file
            miniROOTfile = true;
            break;

        case 's': //--seed ; RNG seed
            try {
                rngSeed = std::stoi(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading rngSeed" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected an integer!" << G4endl;
                exit(1);
            }
            break;

        case 'q': //--quickmode ; Quick mode (skip most plots)
            quickmode = true;
            break;

        case 1004: //--anaScatterTest ; Make analytical plots for scattering
            anaScatterTest = true;
            break;

        case 1000: //--cutoffEnergyFraction ; Cutoff energy fraction
            try {
                cutoff_energyFraction = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading cutoff_energyFraction" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
            
            if (cutoff_energyFraction > 1.0 || cutoff_energyFraction < 0) {
                G4ExceptionDescription errormessage;
                errormessage << "Expected an energy cut off fraction between 0.0 and 1.0, but " << cutoff_energyFraction << " was given!";
                G4Exception("MiniScatter.cc Flag Check: --cutoff_energyFraction","MSFlagCheck1000",FatalException,errormessage);
            }
            break;

        case 1001: //--cutoffRadius ; Cutoff radius [mm]
            try {
                cutoff_radius = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading cutoff_radius" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
            break;

        case 1002: //--edepDZ ; Z bin width for energy deposit histograms [mm]
            try {
                edep_dens_dz = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading edep_dens_dz" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
            break;

        case 1003: //--engNbins ; Number of bins for the energy 1D histograms
            if (engNbins != 0) {
                G4cout << "Can only set engNbins once." << G4endl;
            }

            try {
                engNbins = std::stoi(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading engNbins" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected an integer!" << G4endl;
                exit(1);
            }

            if (engNbins < 0){
                G4cout << "engNbins must be >= 0" << G4endl;
            }
            break;
            
        case 1005: //--histPosLim ; Position Histogram Limit change
            try {
                histPosLim = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading histPosLim" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number!" << G4endl;
                exit(1);
            }
            break;

        case 1006: //--histAngLim ; Angle Histogram Limit change
            try {
                histAngLim = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading histAngLim" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number!" << G4endl;
                exit(1);
            }
            break;

        case 1100: //--object/--magnet ; Definition of an Object, also known as a Magnet
            magnetDefinitions.push_back(string(optarg));
            break;
        
        default: // WTF?
            G4cout << "Got an unknown getopt_char '" << char(getopt_char) << "' ("<< getopt_char<<")"
                   << " when parsing command line arguments." << G4endl;
            exit(1);
        }
    }

    //Copy remaining arguments to array that is passed to Geant4
    int argc_effective = argc-optind+1;
    char** argv_effective = new char*[argc_effective];
    argv_effective[0] = argv[0]; //First arg is always executable name
    for (int i = optind; i < argc;i++){
        argv_effective[i+1-optind] = argv[i];
    }

    //Print the gotten/default arguments
    printHelp(target_thick,
              target_material,
              background_material,
              &detector_distances,
              detector_angle,
              target_angle,
              world_size,
              physListName,
              physCutoffDist,
              numEvents,
              beam_energy,
              beam_eFlat_min,
              beam_eFlat_max,
              beam_type,
              beam_offset,
              beam_zpos,
              beam_angle,
              covarianceString,
              beam_rCut,
              doBacktrack,
              beam_loadFile,
              rngSeed,
              filename_out,
              foldername_out,
              quickmode,
              anaScatterTest,
              miniROOTfile,
              cutoff_energyFraction,
              cutoff_radius,
              edep_dens_dz,
              engNbins,
              histPosLim,
              histAngLim,
              magnetDefinitions);

    G4cout << "Status of other arguments:" << G4endl
           << "numEvents         =  " << numEvents << G4endl
           << "useGUI            =  " << (useGUI==true ? "yes" : "no") << G4endl;
    if (covarianceString != "") {
        G4cout << "Covariance string = '" << covarianceString << "'" << G4endl;
    }
    G4cout << "Arguments which are passed on to Geant4:" << G4endl;
    for (int i = 0; i < argc_effective; i++) {
        G4cout << i << " '" << argv_effective[i] << "'" << G4endl;
    }
    G4cout << G4endl;

    G4cout << "Starting Geant4..." << G4endl << G4endl;

    G4RunManager* runManager = new G4RunManager;

    //Set the initial seed
    if (rngSeed == 0) {
        G4Random::setTheSeed(time(NULL));
    }
    else {
        G4Random::setTheSeed(rngSeed);
    }

    // ** Set mandatory initialization classes **

    // Physics
    G4int verbose=0;
    G4PhysListFactory plFactory;
    G4VModularPhysicsList* physlist = plFactory.GetReferencePhysList(physListName);
    if (physlist==NULL) {
        G4cerr << "Bad physics list!" << G4endl;
        G4cerr << G4endl;

        G4cerr << "Possiblities:" << G4endl;
        const std::vector<G4String>& listnames_hadr =  plFactory.AvailablePhysLists();
        for (auto l : listnames_hadr) {
            G4cerr << "'" << l << "'" << G4endl;
        }
        G4cerr << G4endl;

        G4cerr << "EM options:" << G4endl;
        const std::vector<G4String>& listnames_em =  plFactory.AvailablePhysListsEM();
        for (auto l : listnames_em) {
            G4cerr << "'" << l << "'" << G4endl;
        }
        G4cerr << G4endl;

        exit(1);
    }
    physlist->SetVerboseLevel(verbose);
    runManager->SetUserInitialization(physlist);
    //physlist->SetDefaultCutValue( 0.00001*mm);
    physlist->SetDefaultCutValue( physCutoffDist*mm);

    // Geometry and detectors
    
    G4double world_min_length_detectors = 
        VirtualTrackerWorldConstruction::ComputeMaxAbsZ(&detector_distances, detector_angle, world_size) * 2.0;

    G4double world_min_length_beam = fabs(PrimaryGeneratorAction::GetDefaultZpos(target_thick)) * 2.0;
    if (beam_zpos != 0.0) {
        world_min_length_beam = fabs(beam_zpos) * 2.0;
    }

    G4double world_min_length = world_min_length_detectors;
    if (world_min_length_beam > world_min_length) {
        world_min_length = world_min_length_beam;
    }

    DetectorConstruction* physWorld = new DetectorConstruction(target_thick,
                                                               target_material,
                                                               target_angle,
                                                               background_material,
                                                               world_size,
                                                               world_min_length,
                                                               magnetDefinitions);

    MagnetSensorWorldConstruction* magnetSensorWorld =
        new MagnetSensorWorldConstruction("MagnetSensorWorld",physWorld);
    physWorld->RegisterParallelWorld(magnetSensorWorld);
    physlist->RegisterPhysics(new G4ParallelWorldPhysics("MagnetSensorWorld"));

    VirtualTrackerWorldConstruction* virtualTrackerWorld = 
        new VirtualTrackerWorldConstruction("VirtualTrackerWorld",physWorld, &detector_distances, detector_angle);
    physWorld->RegisterParallelWorld(virtualTrackerWorld);
    physlist->RegisterPhysics(new G4ParallelWorldPhysics("VirtualTrackerWorld"));

    runManager->SetUserInitialization(physWorld);

    // ** Set user action classes **

    PrimaryGeneratorAction* gen_action = new PrimaryGeneratorAction(physWorld,
                                                                    beam_energy,
                                                                    beam_type,
                                                                    beam_offset,
                                                                    beam_zpos,
                                                                    doBacktrack,
                                                                    beam_angle,
                                                                    covarianceString,
                                                                    beam_rCut,
                                                                    rngSeed,
                                                                    beam_eFlat_min,
                                                                    beam_eFlat_max,
                                                                    beam_loadFile,
                                                                    numEvents);
    runManager->SetUserAction(gen_action);
    //
    RunAction* run_action = new RunAction;
    runManager->SetUserAction(run_action);
    //
    EventAction* event_action = new EventAction(run_action);
    runManager->SetUserAction(event_action);

    // ** Final initializations **

    //Initialize G4 kernel
    runManager->Initialize();
    //Initialize magnetic fields
    physWorld->PostInitialize();

    //Configure ROOT output
    RootFileWriter::GetInstance()->setFilename(filename_out);
    RootFileWriter::GetInstance()->setFoldername(foldername_out);
    RootFileWriter::GetInstance()->setQuickmode(quickmode);
    RootFileWriter::GetInstance()->setanaScatterTest(anaScatterTest);
    RootFileWriter::GetInstance()->setMiniFile(miniROOTfile);
    RootFileWriter::GetInstance()->setBeamEnergyCutoff(cutoff_energyFraction);
    RootFileWriter::GetInstance()->setPositionCutoffR(cutoff_radius);
    RootFileWriter::GetInstance()->setEdepDensDZ(edep_dens_dz);
    RootFileWriter::GetInstance()->setEngNbins(engNbins); // 0 = auto
    RootFileWriter::GetInstance()->setNumEvents(numEvents); // May be 0
    RootFileWriter::GetInstance()->setRNGseed(rngSeed);
    RootFileWriter::GetInstance()->setHistPosLim(histPosLim);
    RootFileWriter::GetInstance()->setHistAngLim(histAngLim);

#ifdef G4VIS_USE
    // Initialize visualization
    G4VisManager* visManager = new G4VisExecutive;
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();
#endif

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    // ** Start the run **

    if (argc_effective != 1) { // batch mode
        if (useGUI) {
            G4cout << "UseGUI is not compatible with batch mode!" << G4endl;
            exit(1);
        }
        G4String command = "/control/execute ";
        G4String fileName = argv_effective[1];
        UImanager->ApplyCommand(command+fileName);
    }
    else if (useGUI == true) {  // interactive mode : define UI session
#ifdef G4UI_USE
        G4UIExecutive* ui = new G4UIExecutive(argc_effective,argv_effective);
#ifdef G4VIS_USE
        UImanager->ApplyCommand("/control/execute vis.mac");
#endif
        if (numEvents > 0) {
            G4cout << G4String("'/run/beamOn ") + std::to_string(numEvents) << "'" << G4endl;
            UImanager->ApplyCommand(G4String("/run/beamOn ") + std::to_string(numEvents));
        }

        if (ui->IsGUI())
            ui->SessionStart(); //Returns when GUI is closed.
        delete ui;
#else
        G4cout << "ERROR in initialization: GUI is turned on (-g), but G4UI_USE is not set!" << G4endl;
        exit(1);
#endif
    }

    //Run given number of events
    if (useGUI==false and numEvents > 0) {
        G4cout << G4String("'/run/beamOn ") + std::to_string(numEvents) << "'" << G4endl;
        UImanager->ApplyCommand(G4String("/run/beamOn ") + std::to_string(numEvents));
    }

    G4cout <<"Done." << G4endl;

    // ** Job termination and cleanup **

    // Free the store: user actions, physics_list and detector_description are
    //                 owned and deleted by the run manager, so they should not
    //                 be deleted in the main() program !
#ifdef G4VIS_USE
    delete visManager;
#endif

    // Delete runManager
    delete runManager;
    G4cout << "runManager deleted" << G4endl;

    //Cleanup memory
    delete[] argv_effective;
    argv_effective = NULL;

    delete magnetSensorWorld;
    magnetSensorWorld = NULL;

    G4cout << "MiniScatter is finished" << G4endl;

    return 0;
}

//--------------------------------------------------------------------------------

void printHelp(G4double target_thick,
               G4String target_material,
               G4String background_material,
               std::vector<G4double>* detector_distances,
               G4double detector_angle,
               G4double target_angle,
               G4double world_size,
               G4String physListName,
               G4double physCutoffDist,
               G4int    numEvents,
               G4double beam_energy,
               G4double beam_eFlat_min,
               G4double beam_eFlat_max,
               G4String beam_type,
               G4double beam_offset,
               G4double beam_zpos,
               G4double beam_angle,
               G4String covarianceString,
               G4double beam_rCut,
               G4bool   doBacktrack,
               G4String beam_loadFile,
               G4int    rngSeed,
               G4String filename_out,
               G4String foldername_out,
               G4bool   quickmode,
               G4bool   anaScatterTest,
               G4bool   miniROOTfile,
               G4double cutoff_energyFraction,
               G4double cutoff_radius,
               G4double edep_dens_dz,
               G4int    engNbins,
               G4double histPosLim,
               G4double histAngLim,
               std::vector<G4String> &magnetDefinitions) {
            G4cout << "Usage/options, long and short forms:" << G4endl<<G4endl;

            G4cout << " --thick/-t <double>" << G4endl
                   << "\t Target thickness [mm]"
                   << "\t Set thickness to 0.0 for no target (only objects)" << G4endl
                   << "\t Default/current value = " << target_thick << G4endl <<G4endl;

            G4cout << " --mat/-m <string>" << G4endl
                   << "\t Target material name" << G4endl
                   << "\t Valid choices: 'G4_Al', 'G4_Au','G4_C', 'G4_Cu', 'G4_Pb', 'G4_Ti', 'G4_Si', 'G4_W', 'G4_U', 'G4_Fe'," << G4endl
                   << "\t                'G4_MYLAR', 'G4_KAPTON', 'G4_STAINLESS-STEEL', 'G4_WATER', 'G4_SODIUM_IODIDE', 'G4_Galactic', 'G4_AIR'," << G4endl
                   << "\t                'Sapphire', 'ChromoxPure', 'ChromoxScreen'." << G4endl
                   << "\t Also possible: 'gas::pressure' " << G4endl
                   << "\t                where 'gas' is 'H_2', 'He', 'N_2', 'Ne', or 'Ar'," << G4endl
                   << "\t                and pressure is given in mbar (T=300K is assumed)."  << G4endl
                   << "\t Default/current value = '"
                   << target_material << "'" << G4endl << G4endl;

            G4cout << "--backgroundMaterial <string>" << G4endl
                   << "\t The name of the background material (world volume) to use." << G4endl
                   << "\t Default/current value = '"
                   << background_material << "'" << G4endl << G4endl;

            G4cout << " --dist/-d <double>(:<double>:<double>:...) or NONE" << G4endl
                   << "\t Detector distance(s) [mm]" << G4endl
                   << "\t Default/current value = ";
            if (detector_distances->size() == 0) {
                G4cout << "NONE" << G4endl<<G4endl;
            }
            else {
                for (int i = 0; i < ((int)detector_distances->size())-1; i++) {
                    G4cout << detector_distances->at(i) << ":";
                }
                G4cout << detector_distances->back() << G4endl<<G4endl;
            }

            G4cout << " --ang/-a <double>" << G4endl
                   << "\t Detector angle [deg]" << G4endl
                   << "\t Default/current value = "
                   << detector_angle << G4endl << G4endl;

            G4cout << " --targetAngle/-A <double>" << G4endl
                   << "\t Target angle [deg]" << G4endl
                   << "\t Default/current value = "
                   << target_angle << G4endl << G4endl;

            G4cout << " --worldsize/-w <double>" << G4endl
                   << "\t World size X/Y [mm]" << G4endl
                   << "\t Default/current value = "
                   << world_size << G4endl << G4endl;

            G4cout << " --phys/-p <string>" << G4endl
                   << "\t Physics list name" << G4endl
                   << "\t Default/current value = '"
                   << physListName << "'" << G4endl << G4endl;

            G4cout << " --physCutoffDist <double>" << G4endl
                   << "\t Standard physics cutoff distance [mm]" << G4endl
                   << "\t Default/current value = "
                   << physCutoffDist << G4endl << G4endl;

            G4cout << " --numEvents/-n <int>" << G4endl
                   << "\t Run a given number of events automatically" << G4endl
                   << "\t Default/current value = "
                   << numEvents << G4endl << G4endl;

            G4cout << " --energy/-e <double>" << G4endl
                   << "\t Beam kinetic energy [MeV]" << G4endl
                   << "\t Default/current value = "
                   << beam_energy << G4endl << G4endl;

            G4cout << " --energyDistFlat <double>[MeV]:<double>[MeV]" << G4endl
                   << "\t Use a flat kinetic energy distribution, between the two limits." << G4endl
                   << "\t If --energyDistFlat is used, --energy/-e will just be the reference energy." << G4endl
                   << "\t Default/current min/max : "
                   << beam_eFlat_min << ", " << beam_eFlat_max << G4endl << G4endl;

            G4cout << " --beam/-b <string>" << G4endl
                   << "\t Particle type" << G4endl
                   << "\t This accepts standard Geant4 particle types (see /gun/List for all of them)," << G4endl
                   << "\t typcial examples are 'e-', 'proton', 'gamma'." << G4endl
                   << "\t Ions can also be specified as 'ion::Z;A' where Z and A are the nucleus charge and mass number." << G4endl
                   << "\t Default/current value = '" << beam_type << "'" << G4endl << G4endl;

            G4cout << " --xoffset/-x <double>" << G4endl
                   << "\t Beam offset (x) [mm]" << G4endl
                   << "\t Default/current value = " << beam_offset << G4endl << G4endl;

            G4cout << " --zoffset/-z (*)<double>" << G4endl
                   << "\t Beam offset (z) [mm]" << G4endl
                   << "\t If set to 0.0, start at half the buffer distance. Note that target always at z=0." << G4endl
                   << "\t If a '*' is prepended, the distribution is to be generated at z=0," << G4endl
                   << "\t then backtracked to the given z value (which may be 0.0)" << G4endl
                   << "\t Default/current value = " << beam_zpos
                   << ", doBacktrack = " << (doBacktrack?"true":"false") << G4endl << G4endl;

            G4cout << " --beamAngle <double>" << G4endl
                   << "\t Initial beam angle [deg]" << G4endl
                   << "\t If set to nonzero value, rotate beam generation point and axis around x=y=z=0 by this amount." << G4endl
                   << "\t The beam will still be aimed at x=y=z=0." << G4endl
                   << "\t The beam will be generated in the x/z plane, y=y'=0. Positive angle => positive x." << G4endl
                   << "\t Angles with absolute value >= 90 degrees are not accepted." << G4endl
                   << "\t Currently not compatible with --xoffset, --covar, --beamRcut, or backtracking from * in --zoffset." << G4endl
                   << "\t Note that the absolute value of --zoffset (which is generally negative)" << G4endl
                   << "\t will mean the distance from x=y=z=0 to the starting point."
                   << "\t Default/current value = " << beam_angle << G4endl << G4endl;

            G4cout << " --covar/-c epsN[um]:beta[m]:alpha(::epsN_Y[um]:betaY[m]:alphaY)" << G4endl
                   << "\t Set initial gaussian beam distribution given in terms of Twiss parameters." << G4endl
                   << "\t If the optional part is given x,y are treated separately." << G4endl 
                   << "\t Default/current value = '" << covarianceString << "'" << G4endl << G4endl;

            G4cout << " --beamRcut <double>" << G4endl
                   << "\t Radial cutoff for the beam distribution." << G4endl
                   << "\t If given alone, generate a circular uniform distribution." << G4endl
                   << "\t If given together with --covar/-c, generate a multivariate gaussian" << G4endl
                   << "\t with all particles starting within the given radius." << G4endl
                   << "\t Default/current value = " << beam_rCut << G4endl << G4endl;

            G4cout << " --beamFile <string>" << G4endl
                   << "\t Filename to load the initial beam distribution from." << G4endl
                   << "\t Only the --zoffset flag is taken into account;" << G4endl
                   << "\t  the other beam-relevant flags are ignored." << G4endl
                   << "\t The flags --energy and --beam is only used for reference energy in Twiss etc." << G4endl
                   << "\t Note that the number of particles in the file must be at least as big as --numEvents." << G4endl
                   << "\t Expected format is determined from file ending, possibilities are:" << G4endl
                   << "\t\t .csv: Comma-separated list with 1 particle per row, fields are" << G4endl
                   << "\t\t       particle_type, x<double, mm>, x' <double,px/pz>, y, y', z, Ekin<double,MeV>" << G4endl
                   << "\t\t       Here the particle_type is specified in the same way as in --beam." << G4endl
                   << "\t\t       Please see beamFile.csv for an example." << G4endl
                   << "\t Default/current value = \"" << beam_loadFile << "\"" << G4endl << G4endl;

            G4cout << " --seed/-s <int>" << G4endl
                   << "\t Set the initial seed, 0->use the clock etc." << G4endl
                   << "\t Default/current value = " << rngSeed << G4endl << G4endl;

            G4cout << " --help/-h" << G4endl
                   << "\t Display this help-text and exit." << G4endl << G4endl;

            G4cout << " --gui/-g" << G4endl
                   << "\t Use the standard Geant4 GUI" << G4endl << G4endl;

            G4cout << " --quickmode/-q" << G4endl
                   << "\t Quickmode, skip most post-processing and plots" << G4endl
                   << "\t Default/current value = " << (quickmode?"true":"false") << G4endl << G4endl;

            G4cout << " --anaScatterTest" << G4endl
                   << "\t Compute analytical scattering angle distribution" << G4endl
                   << "\t Default/current value = " << (anaScatterTest?"true":"false") << G4endl << G4endl;

            G4cout << " --miniroot/-r" << G4endl
                   << "\t miniROOTfile, write small root file with only anlysis output, no TTrees." << G4endl
                   << "\t Default/current value = " << (miniROOTfile?"true":"false") << G4endl << G4endl;

            G4cout << " --outname/-f <string>" << G4endl
                   << "\t Output filename" << G4endl
                   << "\t Default/current value = " << filename_out << G4endl << G4endl;

            G4cout << " --outfolder/-o <string>" << G4endl
                   << "\t Output folder" << G4endl
                   << "\t Default/current value = " << foldername_out << G4endl << G4endl;

            G4cout << " --cutoffEnergyFraction <double>" << G4endl
                   << "\t Minimum of beam energy to require for 'cutoff' plots," << G4endl
                   << "\t particles with energy <= to this will not be included." << G4endl
                   << "\t Default/current value = " << cutoff_energyFraction << G4endl << G4endl;

            G4cout << " --cutoffRadius <double>" << G4endl
                   << "\t Maximum radius on target to require for 'cutoff' plots," << G4endl
                   << "\t Particles with r >= this will not be included." << G4endl
                   << "\t Default/current value = " << cutoff_radius << " [mm]" << G4endl << G4endl;

            G4cout << " --edepDZ <double>" << G4endl
                   << "\t Z bin width for energy deposit histograms" << G4endl
                   << "\t Default/current value = " << edep_dens_dz << " [mm]" << G4endl << G4endl;

            G4cout << " --engNbins <int>" << G4endl
                   << "\t Number of bins for 1D energy histograms (0 => internal default), " << G4endl
                   << "\t Default/current value = " << engNbins << G4endl << G4endl;

            G4cout << " --histPosLim <double>" << G4endl
                   << "\t Range for Position axes in ROOT histograms" << G4endl
                   << "\t Default/current value = +/-" << histPosLim << G4endl << G4endl;

            G4cout << " --histAngLim <double>" << G4endl
                   << "\t Range for Angle axes in ROOT histograms" << G4endl
                   << "\t Default/current value = +/-" << histAngLim << G4endl << G4endl;

            G4cout << " --object/--magnet (*)pos:type:length:gradient(:type=val1:specific=val2:arguments=val3)" << G4endl
                   << "\t Create an object (which may be a magnet) of the given type at the given position. " << G4endl
                   << "\t If a '*' is prepended the position (<double> [mm]), the position is the " << G4endl
                   << "\t   start of the active element relative to the end of the target;" << G4endl
                   << "\t   otherwise it is the z-position of the middle of the element." << G4endl
                   << "\t The gradient (<double> [T/m]) is the focusing gradient of the object (should be 0.0 if not a magnet)." << G4endl
                   << "\t The length (<double> [mm]) is the total length of the volumes used by the object." << G4endl
                   << "\t The rest of the arguments are given as key=value pairs, which may be type-specific." << G4endl
                   << G4endl
                   << "\t Standard key=val pairs:" << G4endl
                   << "\t     xOffset:   Center offset in X (<double> [mm]) " << G4endl
                   << "\t     yOffset:   Center offset in Y (<double> [mm]) " << G4endl
                   << "\t     xRot:      Rotation around horizontal axis (<double> [deg])" << G4endl
                   << "\t     yRot:      Rotation around vertical axis   (<double> [deg])" << G4endl
                   << "\t   Note that the offset is applied first," << G4endl
                   << "\t     then the object is rotated around the offset point." << G4endl
                   << "\t     The xRot is applied before the yRot." << G4endl
                   << G4endl
                   << "\t Accepted types and their arguments:" << G4endl
                   << G4endl
                   << "\t   'PLASMA1':" << G4endl
                   << "\t     Models a linear-field plasma lens, only vacuum, in a saphire crystal" << G4endl
                   << "\t     Specific parameters:" << G4endl
                   << "\t     radius:    Capillary radius (<double> [mm])" << G4endl
                   << "\t     totalAmps: Flag (<True/False>) to interpret the gradient parameter" << G4endl
                   << "\t                as the total current [A] instead of in [T/m]." << G4endl
                   << "\t     width:     Capillay crystal width (<double> [mm])" << G4endl
                   << "\t     height:    Capillay crystal height (<double> [mm])" << G4endl
                   << G4endl
                   << "\t   'COLLIMATOR1':" << G4endl
                   << "\t     Models a rectangular collimator with a circular hole in the middle along the z-axis, no field." << G4endl
		   << "\t     Specific parameters:" << G4endl
                   << "\t     radius:    Channel radius   (<double> [mm])" << G4endl
                   << "\t     width:     Absorber width   (<double> [mm])" << G4endl
                   << "\t     height:    Absorber height  (<double> [mm])" << G4endl
                   << "\t     material:  Absorber material (similar to -m)" << G4endl
                   << G4endl
                   << "\t   'COLLIMATORRECT':" << G4endl
                   << "\t     Models a rectangular collimator with a rectangular hole in the middle along the z-axis, no field." << G4endl
		   << "\t     Specific parameters:" << G4endl
                   << "\t     apertureWidth:    Channel width   (<double> [mm]), default: 200 [mm]" << G4endl
                   << "\t     apertureHeight:   Channel height  (<double> [mm]), default:  80 [mm]" << G4endl
                   << "\t     absorberWidth:    Absorber width  (<double> [mm]), default: 250 [mm]" << G4endl
                   << "\t     absorberHeight:   Absorber height (<double> [mm]), default: 100 [mm]" << G4endl
                   << "\t     material:         Absorber material (similar to -m), default: G4_Al" << G4endl
                   << G4endl
                   << "\t   'TARGET':" << G4endl
                   << "\t     Models a rectangular target, no field." << G4endl
                   << "\t     Specific parameters:" << G4endl
                   << "\t     width:     Target width (<double> [mm])" << G4endl
                   << "\t     height:    Target height (<double> [mm])" << G4endl
                   << "\t     material:  Target material (similar to -m)" << G4endl
                   << G4endl
                   << "\t   'TARGETR':" << G4endl
                   << "\t     Models a cylindrical target, no field." << G4endl
                   << "\t     Specific parameters:" << G4endl
                   << "\t     radius:    Target radius(<double> [mm])" << G4endl
                   << "\t     material:  Target material (similar to -m)" << G4endl
                   << G4endl
                   << "\t   'COLLIMATORHV':" << G4endl
                   << "\t     Models two collimating jaws, no field." << G4endl
                   << "\t     Specific parameters:" << G4endl
                   << "\t     gap:       Gap between collimator jaws (<double> [mm])" << G4endl
                   << "\t     HV:        Horizontal (H) or vertical (V) collimator ('H'/'V')" << G4endl
                   << "\t     jawThick:  Collimator jaw thickness (i.e. in in gap plane) (<double> [mm])" << G4endl
                   << "\t     jawHeight: Collimator jaw height (i.e. out of gap plane) (<double> [mm])" << G4endl
                   << "\t     material:  Jaw material (similar to -m)" << G4endl
                   << G4endl
                   << "\t   'SHIELDEDSCINTILLATOR':" << G4endl
                   << "\t     Models a cylindrical scintilling crystal (sensitive)" << G4endl
                   << "\t     inside a cylindrical shield (not sensitive), no field." << G4endl
                   << "\t     Specific parameters:" << G4endl
                   << "\t     scintMat:  Scintillator material (similar to -m)," << G4endl
                   << "\t                defaults to 'G4_SODIUM_IODIDE'" << G4endl
                   << "\t     shieldMat: Shielding material (similar to -m), defaults to 'G4_Pb'" << G4endl
                   << "\t     rScint:    Scintillator radius   (<double> [mm])" << G4endl
                   << "\t     lScint:    Scintillator length   (<double> [mm])" << G4endl
                   << "\t     zScint:    Scintillator position (<double> [mm])" << G4endl
                   << "\t     riShield:  Shield inner radius   (<double> [mm])" << G4endl
                   << "\t     roShield:  Shield outer radius   (<double> [mm])" << G4endl
                   << G4endl
                   << "\t   'PBW':" << G4endl
                   << "\t     Models the ESS Proton Beam Window, no field." << G4endl
                   << "\t     Please remember that 'pos' is the position of the enveloping volume," << G4endl
                   << "\t       which is not the center of the actual window." << G4endl
                   << "\t     Also please note that the standar parameter 'length' should be set to 0," << G4endl
                   << "\t       since it is auto-calculated based on the radius etc." << G4endl
                   << "\t     Specific parameters:" << G4endl
                   << "\t     radius:     Inner radius of cylinder, >0         (<double> [mm]),    default: 88.0 [mm]" << G4endl
                   << "\t     material:   Target material (similar to -m),                         default: G4_Al" << G4endl
                   << "\t     al1Thick:   Outer thickness of metal window, >0  (<double> [mm]),    default: 1.0  [mm]" << G4endl
                   << "\t     waterThick: Thickness of water channel, >0       (<double> [mm]),    default: 2.0  [mm]" << G4endl
                   << "\t     al2Thick:   Inner thickness of metal window, >0  (<double> [mm]),    default: 1.25 [mm]" << G4endl
                   << "\t     width:      Width of cylinder as seen by PBW, >0 (<double> [mm]),    default: 60.0 [mm]" << G4endl
                   << "\t     arcPhi:     Arc angle of window section          (<double> [deg]),   default: 120 [deg]" << G4endl
                   << "\t                 arcPhi should be within: 0 < arcPhi <= 180 [deg]" << G4endl
                   << G4endl;

            G4cout << "\t Currently the following magnet setups are specified:" << G4endl;
            int i = 1;
            for (auto mag : magnetDefinitions) {
                G4cout << "\t   " << i++ << ": '" << mag << "'" << G4endl;
            }

            G4cout << G4endl
                   << G4endl;

            G4cout << " Note that if both -g and -n is used, the events are ran before the GUI is opened." << G4endl;
            G4cout << " One may also use one or more arguments which does not include a '-n' -- these are forwarded untouched to Geant4." << G4endl;
            G4cout << " The first argument not in the form '-char' is interpreted as a macro to run. Don't use vis.mac, it will crash." << G4endl;
            G4cout << " Extra arguments are not compatible with -g" << G4endl;
            G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
