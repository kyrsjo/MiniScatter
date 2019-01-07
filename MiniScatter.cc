#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "DetectorConstruction.hh"
#include "ParallelWorldConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"

#include "G4PhysListFactory.hh"
#include "G4ParallelWorldPhysics.hh"

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
               G4double detector_distance,
               G4double detector_angle,
               G4double world_size,
               G4String physListName,
               G4double beam_energy,
               G4String beam_type,
               G4double beam_offset,
               G4double beam_zpos,
               G4double beam_rCut,
               G4bool   doBacktrack,
               G4int    rngSeed,
               G4String filename_out,
               G4String foldername_out,
               G4bool   quickmode,
               G4bool   miniROOTfile,
               G4double cutoff_energyFraction,
               G4double cutoff_radius,
               G4double edep_dens_dz,
               std::vector<G4String> &magnetDefinitions);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

    // Avoid problem (segfault)
    // "Error in <UnknownClass::InitInterpreter()>: LLVM SYMBOLS ARE EXPOSED TO CLING! This will cause problems; please hide them or dlopen() them after the call to TROOT::InitInterpreter()!"
    // when running with GUI.
    gROOT->Reset();

    //Parse command line arguments
    int getopt_char;
    int getopt_idx;

    G4double target_thick        = 1.0;       // Target thickness  [mm]
    G4String target_material     = "G4_Al";   // Name of target material to use

    G4double detector_distance   = 50.0;      // Detector distance at x=y=0  [mm]
    G4double detector_angle      = 0.0;       // Detectector angle around y-axis [deg]
    G4bool   detector_rotate     = false;

    G4double world_size          = 0.0;       // World size X/Y [mm]

    G4double beam_energy = 200;               // Beam energy [MeV]
    G4String beam_type   = "e-";              // Beam particle type
    G4double beam_offset = 0.0;               // Beam offset (x) [mm]
    G4double beam_zpos   = 0.0;               // Beam offset (z) [mm]
    G4bool   doBacktrack = false;             // Backtrack to the z-position?
    G4String covarianceString = "";           // Beam covariance matrix parameters
    G4double beam_rCut = 0.0;                 // Beam distribution radial cutoff

    G4String physListName = "QGSP_FTFP_BERT"; // Name of physics list to use

    G4int    numEvents    = 0;                // Number of events to generate

    G4bool   useGUI       = false;            // GUI on/off
    G4bool   quickmode    = false;            // Don't make slow plots
    G4String filename_out   = "output";       // Output filename
    G4String foldername_out = "plots";        // Output foldername
    G4bool   miniROOTfile = false;            // Write small root file
                                              // (only analysis output, no TTrees)

    G4int    rngSeed      = 123;              // RNG seed

    G4double cutoff_energyFraction = 0.95;    // [fraction]
    G4double cutoff_radius         = 1.0;     // [mm]

    G4double edep_dens_dz          = 0.0;     // Z bin width for energy deposit histograms [mm]

    std::vector<G4String> magnetDefinitions;

    static struct option long_options[] = {
                                           {"thick",                 required_argument, NULL, 't' },
                                           {"mat",                   required_argument, NULL, 'm' },
                                           {"dist",                  required_argument, NULL, 'd' },
                                           {"ang",                   required_argument, NULL, 'a' },
                                           {"worldsize",             required_argument, NULL, 'w' },
                                           {"dist",                  required_argument, NULL, 'd' },
                                           {"ang",                   required_argument, NULL, 'a' },
                                           {"phys",                  required_argument, NULL, 'p' },
                                           // -n is only short
                                           {"energy",                required_argument, NULL, 'e' },
                                           {"beam",                  required_argument, NULL, 'b' },
                                           {"xoffset",               required_argument, NULL, 'x' },
                                           {"zoffset",               required_argument, NULL, 'z' },
                                           {"covar",                 required_argument, NULL, 'c' },
                                           {"beamRcut",              required_argument, NULL, 1200},
                                           {"outname",               required_argument, NULL, 'f' },
                                           {"outfolder",             required_argument, NULL, 'o' },
                                           {"seed",                  required_argument, NULL, 's' },
                                           {"help",                  no_argument,       NULL, 'h' },
                                           {"gui",                   no_argument,       NULL, 'g' },
                                           {"quickmode",             no_argument,       NULL, 'q' },
                                           {"miniroot",              no_argument,       NULL, 'r' },
                                           {"cutoffEnergyFraction",  required_argument, NULL, 1000 },
                                           {"cutoffRadius",          required_argument, NULL, 1001 },
                                           {"edepDZ",                required_argument, NULL, 1002 },
                                           {"magnet",                required_argument, NULL, 1100 },
                                           {0,0,0,0}
    };

    while ( (getopt_char = getopt_long(argc,argv, "t:m:d:a:w:p:n:e:b:x:z:c:f:o:s:hgqr", long_options, &getopt_idx)) != -1) {
        switch(getopt_char) {
        case 'h': //Help
            printHelp(target_thick,
                      target_material,
                      detector_distance,
                      detector_angle,
                      world_size,
                      physListName,
                      beam_energy,
                      beam_type,
                      beam_offset,
                      beam_zpos,
                      beam_rCut,
                      doBacktrack,
                      rngSeed,
                      filename_out,
                      foldername_out,
                      quickmode,
                      miniROOTfile,
                      cutoff_energyFraction,
                      cutoff_radius,
                      edep_dens_dz,
                      magnetDefinitions);
            exit(1);
            break;

        case 'g': //Use GUI
            useGUI = true;
            break;

        case 't': //Target thickness
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

        case 'd': //Detector distance
            try {
                detector_distance = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading detector distance" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
            break;

        case 'a': //Detector angle
            try {
                detector_angle = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading detector angle" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
            detector_rotate = true;
            break;

        case 'w': //World size
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

        case 'p': //Named physics list
            physListName = G4String(optarg);
            break;

        case 'm': //Target material
            target_material = G4String(optarg);
            break;

        case 'n': //Number of events
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

        case 'e': //beam energy
            try {
                beam_energy = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading beam energy" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
            break;

        case 'b': //Beam type
            beam_type = G4String(optarg);
            break;

        case 'x': //beam offset (x)
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

        case 'z': //beam offset (z)
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

        case 'c': //Beam covariance matrix from Twiss parameters
            //-c epsN[um]:beta[m]:alpha(::epsN_Y[um]:betaY[m]:alphaY)
            covarianceString = G4String(optarg);
            break;

        case 1200: //Beam radial cutoff [mm]
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

        case 'f': //Output filename
            filename_out = G4String(optarg);
            break;

        case 'o': //Output foldername
            foldername_out = G4String(optarg);
            break;

        case 'r': //Mini ROOT file
            miniROOTfile = true;
            break;

        case 's': //RNG seed
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

        case 'q': // Quick mode (skip most plots)
            quickmode = true;
            break;

        case 1000: // Cutoff energy fraction
            try {
                cutoff_energyFraction = std::stod(string(optarg));
            }
            catch (const std::invalid_argument& ia) {
                G4cout << "Invalid argument when reading cutoff_energyFraction" << G4endl
                       << "Got: '" << optarg << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
            break;

        case 1001: // Cutoff radius [mm]
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

        case 1002: // Z bin width for energy deposit histograms [mm]
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

        case 1100: //Magnet definition
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
              detector_distance,
              detector_angle,
              world_size,
              physListName,
              beam_energy,
              beam_type,
              beam_offset,
              beam_zpos,
              beam_rCut,
              doBacktrack,
              rngSeed,
              filename_out,
              foldername_out,
              quickmode,
              miniROOTfile,
              cutoff_energyFraction,
              cutoff_radius,
              edep_dens_dz,
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

    G4RunManager * runManager = new G4RunManager;

    //Set the initial seed
    G4Random::setTheSeed(rngSeed);

    // Set mandatory initialization classes

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
    physlist->SetDefaultCutValue( 0.1*mm);

    DetectorConstruction* physWorld = new DetectorConstruction(target_thick,
                                                               target_material,
                                                               detector_distance,
                                                               detector_angle,
                                                               detector_rotate,
                                                               world_size,
                                                               magnetDefinitions);

    ParallelWorldConstruction* magnetSensorWorld =
        new ParallelWorldConstruction("MagnetSensorWorld",physWorld);
    physWorld->RegisterParallelWorld(magnetSensorWorld);
    physlist->RegisterPhysics(new G4ParallelWorldPhysics("MagnetSensorWorld"));

    runManager->SetUserInitialization(physWorld);

    // Set user action classes:
    //
    PrimaryGeneratorAction* gen_action = new PrimaryGeneratorAction(physWorld,
                                                                    beam_energy,
                                                                    beam_type,
                                                                    beam_offset,
                                                                    beam_zpos,
                                                                    doBacktrack,
                                                                    covarianceString,
                                                                    beam_rCut);
    runManager->SetUserAction(gen_action);
    //
    RunAction* run_action = new RunAction;
    runManager->SetUserAction(run_action);
    //
    EventAction* event_action = new EventAction(run_action);
    runManager->SetUserAction(event_action);

    // Initialize G4 kernel
    runManager->Initialize();

    physWorld->PostInitialize();

    //Set root file output filename
    RootFileWriter::GetInstance()->setFilename(filename_out);
    RootFileWriter::GetInstance()->setFoldername(foldername_out);
    RootFileWriter::GetInstance()->setQuickmode(quickmode);
    RootFileWriter::GetInstance()->setMiniFile(miniROOTfile);
    RootFileWriter::GetInstance()->setBeamEnergyCutoff(cutoff_energyFraction);
    RootFileWriter::GetInstance()->setPositionCutoffR(cutoff_radius);
    RootFileWriter::GetInstance()->setEdepDensDZ(edep_dens_dz);
    RootFileWriter::GetInstance()->setNumEvents(numEvents); // May be 0

#ifdef G4VIS_USE
    // Initialize visualization
    G4VisManager* visManager = new G4VisExecutive;
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();
#endif

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

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
#endif
    }

    //Run given number of events
    if (useGUI==false and numEvents > 0) {
        G4cout << G4String("'/run/beamOn ") + std::to_string(numEvents) << "'" << G4endl;
        UImanager->ApplyCommand(G4String("/run/beamOn ") + std::to_string(numEvents));
    }

    G4cout <<"Done." << G4endl;

    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    //                 owned and deleted by the run manager, so they should not
    //                 be deleted in the main() program !
#ifdef G4VIS_USE
    delete visManager;
#endif
    //End comment
    // Delete runManager
    delete runManager;
    G4cout << "runManager deleted" << G4endl;

    return 0;
}

//--------------------------------------------------------------------------------

void printHelp(G4double target_thick,
               G4String target_material,
               G4double detector_distance,
               G4double detector_angle,
               G4double world_size,
               G4String physListName,
               G4double beam_energy,
               G4String beam_type,
               G4double beam_offset,
               G4double beam_zpos,
               G4double beam_rCut,
               G4bool   doBacktrack,
               G4int    rngSeed,
               G4String filename_out,
               G4String foldername_out,
               G4bool   quickmode,
               G4bool   miniROOTfile,
               G4double cutoff_energyFraction,
               G4double cutoff_radius,
               G4double edep_dens_dz,
               std::vector<G4String> &magnetDefinitions) {
            G4cout << "Welcome to MiniScatter!" << G4endl
                   << G4endl
                   << "Usage/options:" << G4endl;

            G4cout << "-t <double> : Target thickness [mm],  default/current value = "
                   << target_thick << G4endl;

            G4cout << "-m <string> : Target material name,   default/current       = '"
                   << target_material << "'" << G4endl
                   << " Valid choices: 'G4_Al', 'G4_C', 'G4_Cu', 'G4_Pb', 'G4_Ti', 'G4_Si', 'G4_W', 'G4_U', "
                   << "'G4_MYLAR', 'G4_KAPTON', 'G4_STAINLESS-STEEL', 'G4_WATER', 'G4_Galactic', 'Sapphire'" << G4endl
                   << " Also possible: 'gas::pressure' "
                   << " where 'gas' is 'H_2', 'He', 'N_2', 'Ne', or 'Ar',"
                   << " and pressure is given in mbar (T=300K is assumed)." << G4endl;

            G4cout << "-d <double> : Detector distance [mm], default/current value = "
                   << detector_distance << G4endl;

            G4cout << "-a <double> : Detector angle [deg],   default/current value = "
                   << detector_angle << G4endl;

            G4cout << "-w <double> : World size X/Y [mm],    default/current value = "
                   << world_size << G4endl;

            G4cout << "-p <string> : Physics list name,      default/current       = '"
                   << physListName << G4endl;

            G4cout << "-n <int>    : Run a given number of events automatically"
                   << G4endl;

            G4cout << "-e <double> : Beam energy [MeV],      default/current value = "
                   << beam_energy << G4endl;

            G4cout << "-b <string> : Particle type,          default/current value = "
                   << beam_type << G4endl
                   << " This accepts standard Geant4 particle types (see /gun/List for all of them)," << G4endl
                   << " typcial ones are 'e-' and 'proton'." << G4endl
                   << " Ions can also be specified as 'ion::Z,A' where Z and A are the nucleus charge and mass number." << G4endl;

            G4cout << "-x <double> : Beam offset (x) [mm],   default/current value = "
                   << beam_offset << G4endl;

            G4cout << "-z (*)<double> : Beam offset (z) [mm],   default/current value = "
                   << beam_zpos
                   << ", doBacktrack = " << (doBacktrack?"true":"false") << G4endl
                   << " If set to 0.0, start at half the buffer distance. Note that target always at z=0." << G4endl
                   << " If a '*' is prepended, the distribution is to be generated at z=0," << G4endl
                   << " then backtracked to the given z value (which may be 0.0)" << G4endl;
            
            G4cout << "-c epsN[um]:beta[m]:alpha(::epsN_Y[um]:betaY[m]:alphaY) : " << G4endl
                   << " Set realistic beam distribution (on target surface); " << G4endl
                   << " if optional part given then x,y are treated separately" << G4endl;

            G4cout << "--beamRcut <double> : Radial cutoff for the beam distribution." << G4endl
                   << " If given alone, generate a circular uniform distribution." << G4endl
                   << " If given together with -c, generate a multivariate gaussian with all particles starting within the given radius." << G4endl
                   << " Default/current value = " << beam_rCut << G4endl;

            G4cout << "-s <int>    : Set the initial seed,   default/current value = "
                   << rngSeed << G4endl;

            G4cout << "-g : Use a GUI" << G4endl;

            G4cout << "-q : Quickmode, skip most post-processing and plots, default/current value = "
                   << (quickmode?"true":"false") << G4endl;

            G4cout << "-r : miniROOTfile, write small root file with only anlysis output, no TTrees, default/current value = "
                   << (miniROOTfile?"true":"false") << G4endl;

            G4cout << "-f <string> : Output filename,        default/current value = "
                   << filename_out << G4endl;

            G4cout << "-o <string : Output folder,           default/current value = "
                   << foldername_out << G4endl;

            G4cout << "--cutoffEnergyfraction : Minimum of beam energy to require for 'cutoff' plots, "
                   << "default/current value = " << cutoff_energyFraction << G4endl;

            G4cout << "--cutoffRadius         : Maximum radius on target to require for 'cutoff' plots, "
                   << "default/current value = " << cutoff_radius << " [mm]" << G4endl;

            G4cout << "--edepDZ               : Z bin width for energy deposit histograms " 
                   << "default/current value = " << edep_dens_dz << " [mm]" << G4endl;

            G4cout << "--magnet (*)pos:type:length:gradient(:type=val1:specific=val2:arguments=val3) : "
                   << " Create a magnet of the given type at the given position. " << G4endl
                   << " If a '*' is prepended the position (<double> [mm]), the position is the start of the active element "
                   << "relative to the end of the target; otherwize it is the z-position of the middle of the element." << G4endl
                   << " The gradient (<double> [T/m]) is the focusing gradient of the device." << G4endl
                   << " The length <double> [mm] is the total length of the volumes used by the device." << G4endl
                   << " The type-specific arguments are given as key=value pairs." << G4endl
                   << " Accepted types and their arguments:" << G4endl
                   << "  'PLASMA1':" << G4endl
                   << "     radius:    Capillary radius (<double> [mm])" << G4endl
                   << "     totalAmps: Flag (<True/False>) to interpret the gradient parameter"
                   << " as the total current [A] instead of in [T/m]." << G4endl
                   << "     width:     Capillay crystal width (<double> [mm])" << G4endl
                   << "     height:    Capillay crystal height (<double> [mm])" << G4endl
                   << "  'COLLIMATOR1':" << G4endl
                   << "     radius:    Channel radius (<double> [mm])" << G4endl
                   << "     width:     Absorber width (<double> [mm])" << G4endl
                   << "     height:    Absorber height (<double> [mm])" << G4endl
                   << "Currently have the following magnet setups:" << G4endl;
            for (auto mag : magnetDefinitions) {
                G4cout << mag << G4endl;
            }

            G4cout << G4endl
                   << G4endl;

            G4cout << "Note that if both -g and -n is used, the events are ran before the GUI is opened." << G4endl;
            G4cout << "One may also use one or more arguments which does not include a '-n' -- these are forwarded untouched to Geant4" << G4endl;
            G4cout << "The first argument not in the form '-char' is interpreted as a macro to run. Don't use vis.mac, it will crash." << G4endl;
            G4cout << "Extra arguments are not compatible with -g" << G4endl;
            G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
