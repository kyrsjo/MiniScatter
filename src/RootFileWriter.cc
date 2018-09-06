#include "RootFileWriter.hh"

#include "TGraph.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TBranch.h"

#include "MyEdepHit.hh"
#include "MyTrackerHit.hh"

#include "G4SDManager.hh"

#include "G4Track.hh"

#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4SystemOfUnits.hh"
#include "CLHEP/Units/SystemOfUnits.h"

#include <iostream>
#include <iomanip>
#include <experimental/filesystem> //Mainstreamed from C++17

using namespace std;
RootFileWriter* RootFileWriter::singleton = 0;
const G4String  RootFileWriter::foldername_out = "plots";

void RootFileWriter::initializeRootFile(){
    G4RunManager*           run    = G4RunManager::GetRunManager();
    DetectorConstruction*   detCon = (DetectorConstruction*)run->GetUserDetectorConstruction();
    PrimaryGeneratorAction* genAct = (PrimaryGeneratorAction*)run->GetUserPrimaryGeneratorAction();
    this->beamEnergy = genAct->get_beam_energy();

    if (not has_filename_out) {
        G4cerr << "Error: filename_out not set." << G4endl;
        exit(1);
    }
    G4String rootFileName = foldername_out + "/" + filename_out + ".root";

    //Create folder if it does not exist
    if (not experimental::filesystem::exists(foldername_out.data())) {
        G4cout << "Creating folder '" << foldername_out << "'" << G4endl;
        experimental::filesystem::create_directory(foldername_out.data());
    }
    G4cout << "Opening ROOT file '" + rootFileName +"'"<<G4endl;
    histFile = new TFile(rootFileName,"RECREATE");

    eventCounter = 0;

    // TTrees for external analysis
    targetExit = new TTree("TargetExit","TargetExit tree");
    targetExit->Branch("TargetExitBranch", &targetExitBuffer,
                       "x/D:y:z:px:py:pz:E:PDG/I:charge:eventID");

    trackerHits = new TTree("TrackerHits","TrackerHits tree");
    trackerHits->Branch("TrackerHitsBranch", &trackerHitsBuffer,
                        "x/D:y:z:px:py:pz:E:PDG/I:charge:eventID");

    // Target energy deposition
    targetEdep = new TH1D("targetEdep","targetEdep",1000,0,6);
    targetEdep->GetXaxis()->SetTitle("Total energy deposit/event [MeV]");
    targetEdep_NIEL = new TH1D("targetEdep_NIEL","targetEdep_NIEL",1000,0,1);
    targetEdep_NIEL->GetXaxis()->SetTitle("Total NIEL/event [keV]");
    targetEdep_IEL = new TH1D("targetEdep_IEL","targetEdep_IEL",1000,0,6);
    targetEdep_IEL->GetXaxis()->SetTitle("Total ionizing energy deposit/event [MeV]");

    // Target exit angle histogram
    target_exitangle_hist        = new TH1D("exitangle",
                                            "Exit angle from tracker",
                                            5001, -90, 90);
    target_exitangle_hist_cutoff = new TH1D("exitangle_cutoff",
                                            "Exit angle from tracker (charged, energy > cutoff)",
                                            50001, -10, 10);

    // Tracker histograms
    tracker_numParticles = new TH1D("numParticles","numParticles",1001,-0.5,1000.5);
    tracker_numParticles->GetXaxis()->SetTitle("Number of particles / event");

    tracker_energy       = new TH1D("energy","energy",10000,0,beamEnergy);
    tracker_energy->GetXaxis()->SetTitle("Energy of particles on the tracker [MeV]");

    tracker_hitPos        = new TH2D("trackerHitpos", "Tracker Hit position",
               1000,detCon->getDetectorSizeX()/2.0/mm,detCon->getDetectorSizeX()/2.0/mm,
               1000,detCon->getDetectorSizeY()/2.0/mm,detCon->getDetectorSizeY()/2.0/mm);
    tracker_hitPos_cutoff = new TH2D("trackerHitpos_cutoff", "Tracker Hit position (charged, energy > cutoff)",
               1000,detCon->getDetectorSizeX()/2.0/mm,detCon->getDetectorSizeX()/2.0/mm,
               1000,detCon->getDetectorSizeY()/2.0/mm,detCon->getDetectorSizeY()/2.0/mm);

    tracker_phasespaceX   =
        new TH2D("tracker_x",
                 "Tracker phase space (x)",
                 1000,detCon->getDetectorSizeX()/2.0/mm, detCon->getDetectorSizeY()/2.0/mm,
                 50001, -M_PI, M_PI);
    //tracker_phasespaceX->Sumw2();
    tracker_phasespaceY   =
        new TH2D("tracker_y",
                 "Tracker phase space (y)",
                 1000, detCon->getDetectorSizeX()/2.0/mm, detCon->getDetectorSizeY()/2.0/mm,
                 50001, -M_PI, M_PI);
    //tracker_phasespaceY->Sumw2();

    tracker_phasespaceX_cutoff   =
        new TH2D("tracker_cutoff_x",
                 "Tracker phase space (x) (charged, energy > cutoff)",
                 1000,detCon->getDetectorSizeX()/2.0/mm, detCon->getDetectorSizeY()/2.0/mm,
                 50001, -M_PI, M_PI);
    //tracker_phasespaceX_cutoff->Sumw2();
    tracker_phasespaceY_cutoff   =
        new TH2D("tracker_cutoff_y",
                 "Tracker phase space (y) (charged, energy > cutoff)",
                 1000, detCon->getDetectorSizeX()/2.0/mm, detCon->getDetectorSizeY()/2.0/mm,
                 50001, -M_PI, M_PI);
    //tracker_phasespaceY_cutoff->Sumw2();

    init_phasespaceX   =
        new TH2D("init_x",
                 "Initial phase space (x)",
                 1000,detCon->getWorldSizeX()/2.0/mm, detCon->getWorldSizeY()/2.0/mm,
                 50001, -M_PI, M_PI);
    //init_phasespaceX->Sumw2();
    init_phasespaceY   =
        new TH2D("init_y",
                 "Initial phase space (y)",
                 1000, detCon->getWorldSizeX()/2.0/mm, detCon->getWorldSizeY()/2.0/mm,
                 50001, -M_PI, M_PI);
    //init_phasespaceY->Sumw2();

    //For counting the types of particles hitting the trackers
    typeCounter["tracker"]       = particleTypesCounter();
    typeCounter["tracker_cutoff"] = particleTypesCounter();
    typeCounter["target"]        = particleTypesCounter();
    typeCounter["target_cutoff"] = particleTypesCounter();

    //Compute means and RMS of where they hit
    tracker_particleHit_x  = 0.0;
    tracker_particleHit_xx = 0.0;
    tracker_particleHit_y  = 0.0;
    tracker_particleHit_yy = 0.0;

    //Compute means and RMS of where they hit (above cutoff particles only)
    tracker_particleHit_x_cutoff  = 0.0;
    tracker_particleHit_xx_cutoff = 0.0;
    tracker_particleHit_y_cutoff  = 0.0;
    tracker_particleHit_yy_cutoff = 0.0;
    numParticles_cutoff = 0;

    //Compute RMS of target exit angle
    target_exitangle              = 0.0;
    target_exitangle2             = 0.0;
    target_exitangle_numparticles = 0;
    target_exitangle_cutoff              = 0.0;
    target_exitangle2_cutoff             = 0.0;
    target_exitangle_cutoff_numparticles = 0;

}

void RootFileWriter::doEvent(const G4Event* event){

    eventCounter++;

    G4HCofThisEvent* HCE=event->GetHCofThisEvent();
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    //**Data from TargetSD**
    G4int myTargetEdep_CollID = SDman->GetCollectionID("EdepCollection");
    if (myTargetEdep_CollID>=0){
        MyEdepHitsCollection* targetEdepHitsCollection = NULL;
        targetEdepHitsCollection = (MyEdepHitsCollection*) (HCE->GetHC(myTargetEdep_CollID));
        if (targetEdepHitsCollection != NULL) {
            G4int nEntries = targetEdepHitsCollection->entries();
            G4double edep      = 0.0;
            G4double edep_NIEL = 0.0;
            G4double edep_IEL  = 0.0;
            for (G4int i = 0; i < nEntries; i++){
                edep      += (*targetEdepHitsCollection)[i]->GetDepositedEnergy();
                edep_NIEL += (*targetEdepHitsCollection)[i]->GetDepositedEnergy_NIEL();
                edep_IEL  += (*targetEdepHitsCollection)[i]->GetDepositedEnergy() -
                    (*targetEdepHitsCollection)[i]->GetDepositedEnergy_NIEL();
            }
            targetEdep->Fill(edep/MeV);
            targetEdep_NIEL->Fill(edep_NIEL/keV);
            targetEdep_IEL->Fill(edep_IEL/MeV);
        }
        else {
            G4cout << "targetEdepHitsCollection was NULL!"<<G4endl;
        }
    }
    else {
        G4cout << "myTargetEdep_CollID was " << myTargetEdep_CollID << " < 0!"<<G4endl;
    }

    G4int myTargetExitpos_CollID = SDman->GetCollectionID("ExitposCollection");
    if (myTargetExitpos_CollID>=0) {
        MyTrackerHitsCollection* targetExitposHitsCollection = NULL;
        targetExitposHitsCollection = (MyTrackerHitsCollection*) (HCE->GetHC(myTargetExitpos_CollID));
        if (targetExitposHitsCollection != NULL) {
            G4int nEntries = targetExitposHitsCollection->entries();

            for (G4int i = 0; i < nEntries; i++) {
                //Get the data from the event
                const G4double       energy      = (*targetExitposHitsCollection)[i]->GetTrackEnergy();
                const G4int          charge      = (*targetExitposHitsCollection)[i]->GetCharge();
                const G4ThreeVector& momentum    = (*targetExitposHitsCollection)[i]->GetMomentum();
                G4double             exitangle   = atan(momentum.x()/momentum.z())/deg;
                const G4ThreeVector& hitPos      = (*targetExitposHitsCollection)[i]->GetPosition();
                const G4int          PDG         = (*targetExitposHitsCollection)[i]->GetPDG();
                const G4String&      type        = (*targetExitposHitsCollection)[i]->GetType();

                //Particle type counting
                FillParticleTypes(typeCounter["target"], PDG, type);
                if (energy >= beamEnergy*beamEnergy_cutoff) {
                    FillParticleTypes(typeCounter["target_cutoff"], PDG, type);
                }

                //Exit angle
                target_exitangle_hist->Fill(exitangle);
                if (charge != 0 and energy >= beamEnergy*beamEnergy_cutoff) {
                    target_exitangle_hist_cutoff->Fill(exitangle);
                }

                target_exitangle              += exitangle;
                target_exitangle2             += exitangle*exitangle;
                target_exitangle_numparticles += 1;

                if (charge != 0 and energy >= beamEnergy*beamEnergy_cutoff) {
                    target_exitangle_cutoff              += exitangle;
                    target_exitangle2_cutoff             += exitangle*exitangle;
                    target_exitangle_cutoff_numparticles += 1;
                }

                //Fill the TTree
                if (not miniFile) {
                    targetExitBuffer.x = hitPos.x()/mm;
                    targetExitBuffer.y = hitPos.x()/mm;
                    targetExitBuffer.z = hitPos.z()/mm;

                    targetExitBuffer.px = momentum.x()/MeV;
                    targetExitBuffer.py = momentum.y()/MeV;
                    targetExitBuffer.pz = momentum.z()/MeV;

                    targetExitBuffer.E = beamEnergy / MeV;

                    targetExitBuffer.PDG = PDG;
                    targetExitBuffer.charge = charge;

                    targetExitBuffer.eventID = eventCounter;

                    targetExit->Fill();
                }
            }

        }
        else {
            G4cout << "targetExitposHitsCollection was NULL!"<<G4endl;
        }
    }
    else {
        G4cout << "myTargetExitpos_CollID was " << myTargetExitpos_CollID << " < 0!"<<G4endl;
    }

    //**Data from detectorTrackerSD**
    G4int myTrackerSD_CollID = SDman->GetCollectionID("TrackerCollection");
    if (myTrackerSD_CollID>=0) {
        MyTrackerHitsCollection* trackerHitsCollection = NULL;
        trackerHitsCollection = (MyTrackerHitsCollection*) (HCE->GetHC(myTrackerSD_CollID));
        if (trackerHitsCollection != NULL) {
            G4int nEntries = trackerHitsCollection->entries();

            for (G4int i = 0; i < nEntries; i++) {
                //Get the data from the event
                const G4double  energy = (*trackerHitsCollection)[i]->GetTrackEnergy();
                const G4int     PDG    = (*trackerHitsCollection)[i]->GetPDG();
                const G4int     charge = (*trackerHitsCollection)[i]->GetCharge();
                const G4String& type   = (*trackerHitsCollection)[i]->GetType();
                const G4ThreeVector& hitPos   = (*trackerHitsCollection)[i]->GetPosition();
                const G4ThreeVector& momentum = (*trackerHitsCollection)[i]->GetMomentum();

                //Overall histograms
                tracker_energy->Fill(energy/MeV);

                //Hit position
                tracker_hitPos->Fill(hitPos.x()/mm, hitPos.y()/mm);
                if (charge != 0 and energy >= beamEnergy*beamEnergy_cutoff) {
                    tracker_hitPos_cutoff->Fill(hitPos.x()/mm, hitPos.y()/mm);
                }

                //Phase space
                tracker_phasespaceX->Fill(hitPos.x()/mm, momentum.x()/momentum.z());
                tracker_phasespaceY->Fill(hitPos.y()/mm, momentum.y()/momentum.z());

                if (charge != 0 and energy >= beamEnergy*beamEnergy_cutoff) {
                    tracker_phasespaceX_cutoff->Fill(hitPos.x()/mm, momentum.x()/momentum.z());
                    tracker_phasespaceY_cutoff->Fill(hitPos.y()/mm, momentum.y()/momentum.z());
                }

                //Particle type counting
                FillParticleTypes(typeCounter["tracker"], PDG, type);
                if (energy >= beamEnergy*beamEnergy_cutoff) {
                    FillParticleTypes(typeCounter["tracker_cutoff"], PDG, type);
                }

                //Hit positions
                tracker_particleHit_x  +=  hitPos.x()/mm;
                tracker_particleHit_xx += (hitPos.x()/mm)*(hitPos.x()/mm);
                tracker_particleHit_y  +=  hitPos.y()/mm;
                tracker_particleHit_yy += (hitPos.y()/mm)*(hitPos.y()/mm);

                if (charge != 0 and energy >= beamEnergy*beamEnergy_cutoff) {
                    tracker_particleHit_x_cutoff  +=  hitPos.x()/mm;
                    tracker_particleHit_xx_cutoff += (hitPos.x()/mm)*(hitPos.x()/mm);
                    tracker_particleHit_y_cutoff  +=  hitPos.y()/mm;
                    tracker_particleHit_yy_cutoff += (hitPos.y()/mm)*(hitPos.y()/mm);
                    numParticles_cutoff += 1;
                }

                //Fill the TTree
                if (not miniFile) {
                    trackerHitsBuffer.x = hitPos.x()/mm;
                    trackerHitsBuffer.y = hitPos.x()/mm;
                    trackerHitsBuffer.z = hitPos.z()/mm;

                    trackerHitsBuffer.px = momentum.x()/MeV;
                    trackerHitsBuffer.py = momentum.y()/MeV;
                    trackerHitsBuffer.pz = momentum.z()/MeV;

                    trackerHitsBuffer.E = beamEnergy / MeV;

                    trackerHitsBuffer.PDG = PDG;
                    trackerHitsBuffer.charge = charge;

                    trackerHitsBuffer.eventID = eventCounter;

                    trackerHits->Fill();
                }
            }

            tracker_numParticles->Fill(nEntries);
        }
        else{
            G4cout << "trackerHitsCollection was NULL!"<<G4endl;
        }
    }
    else{
        G4cout << "myTrackerSD_CollID was " << myTrackerSD_CollID << "<0!"<<G4endl;
    }

    // Initial particle distribution
    G4RunManager*           run    = G4RunManager::GetRunManager();
    PrimaryGeneratorAction* genAct = (PrimaryGeneratorAction*)run->GetUserPrimaryGeneratorAction();

    init_phasespaceX->Fill(genAct->x/mm,genAct->xp/rad);
    init_phasespaceY->Fill(genAct->y/mm,genAct->yp/rad);
}
void RootFileWriter::finalizeRootFile() {

    //Print out the particle types on all detector planes
    for (auto it : typeCounter) {
        PrintParticleTypes(it.second, it.first);
    }

    // ** Below cutoff **

    //Average position and RMS
    double xave  = tracker_particleHit_x / ((double)typeCounter["tracker"].numParticles);
    double yave  = tracker_particleHit_y / ((double)typeCounter["tracker"].numParticles);
    double xrms  = ( tracker_particleHit_xx - (tracker_particleHit_x*tracker_particleHit_x / ((double)typeCounter["tracker"].numParticles)) ) /
        (((double)typeCounter["tracker"].numParticles)-1.0);
    xrms = sqrt(xrms);
    double yrms  = ( tracker_particleHit_yy - (tracker_particleHit_y*tracker_particleHit_y / ((double)typeCounter["tracker"].numParticles)) ) /
        (((double)typeCounter["tracker"].numParticles)-1.0);
    yrms = sqrt(yrms);
    //Exitangle
    G4double exitangle_avg = target_exitangle /
        ((double)target_exitangle_numparticles);
    G4double exitangle_rms = ( target_exitangle2 -
                               (target_exitangle*target_exitangle /
                                ((double)target_exitangle_numparticles)) ) /
        ((double)target_exitangle_numparticles - 1.0 ) ;
    exitangle_rms = sqrt(exitangle_rms);

    G4cout << G4endl
           << "All particles (n=" << typeCounter["tracker"].numParticles << "):" << G4endl;

    G4cout << "Average x = " << xave << " [mm], RMS = " << xrms << " [mm]" << G4endl
           << "Average y = " << yave << " [mm], RMS = " << yrms << " [mm]" << G4endl;

    G4cout << G4endl
           << "Exit angle average (x) = " << exitangle_avg << " [deg]" << G4endl
           << "Exit angle RMS (x)     = " << exitangle_rms << " [deg]" << G4endl;

    // ** Above cutoff **

    // Average position and RMS
    double xave_cutoff  = tracker_particleHit_x_cutoff / ((double)numParticles_cutoff);
    double yave_cutoff  = tracker_particleHit_y_cutoff / ((double)numParticles_cutoff);
    double xrms_cutoff  = ( tracker_particleHit_xx_cutoff -
                            (tracker_particleHit_x_cutoff*tracker_particleHit_x_cutoff /
                             ((double)numParticles_cutoff)) ) /
        (((double)numParticles_cutoff)-1.0);
    xrms_cutoff = sqrt(xrms_cutoff);
    double yrms_cutoff  = ( tracker_particleHit_yy_cutoff -
                            (tracker_particleHit_y_cutoff*tracker_particleHit_y_cutoff /
                             ((double)numParticles_cutoff)) ) /
        (((double)numParticles_cutoff)-1.0);
    yrms_cutoff = sqrt(yrms_cutoff);

    //Exitangle
    G4double exitangle_avg_cutoff = target_exitangle_cutoff /
        ((double)target_exitangle_cutoff_numparticles);
    G4double exitangle_rms_cutoff = ( target_exitangle2_cutoff -
                                      (target_exitangle_cutoff*target_exitangle_cutoff /
                                       ((double)target_exitangle_cutoff_numparticles)) ) /
        ((double)target_exitangle_cutoff_numparticles - 1.0 ) ;
    exitangle_rms_cutoff = sqrt(exitangle_rms_cutoff);

    G4cout << G4endl
           << "Above cutoff (charged, energy >= "
           << beamEnergy*beamEnergy_cutoff <<" [MeV], n=" << numParticles_cutoff
           << ") only:" << G4endl << G4endl;

    G4cout << "Average x = " << xave_cutoff << " [mm], RMS = " << xrms_cutoff << " [mm]" << G4endl
           << "Average y = " << yave_cutoff << " [mm], RMS = " << yrms_cutoff << " [mm]" << G4endl;

    G4cout << G4endl
           << "Exit angle average (x) = " << exitangle_avg_cutoff << " [deg]" << G4endl
           << "Exit angle RMS (x)     = " << exitangle_rms_cutoff << " [deg]" << G4endl;
    G4cout << G4endl;

    //Compute Twiss parameters
    PrintTwissParameters(init_phasespaceX);
    PrintTwissParameters(init_phasespaceY);
    PrintTwissParameters(tracker_phasespaceX);
    PrintTwissParameters(tracker_phasespaceY);
    PrintTwissParameters(tracker_phasespaceX_cutoff);
    PrintTwissParameters(tracker_phasespaceY_cutoff);

    if (not quickmode) {
        // Compute the analytical multiple scattering angle distribution
        // Formulas from various sources:
        //
        // http://www-glast.slac.stanford.edu/software/AnaGroup/rev.pdf
        // Document describing the validation of Geant4 multiple scattering
        //
        // http://cdsweb.cern.ch/record/1279627/files/PH-EP-Tech-Note-2010-013.pdf
        // (note: Missing factor A in radiation length formula)
        //
        // https://en.wikipedia.org/wiki/Radiation_length
        // Checked: April 15th 2018
        //
        // http://pdg.lbl.gov/2004/reviews/passagerpp.pdf
        // This seems to be the base reference for the "compact fit to the data"
        // used for the radiation length, which is the formula found
        // in Wikipedia and in PH-EP-Tech-Note-2010-013.
        //
        // http://pdg.lbl.gov/2017/reviews/rpp2017-rev-passage-particles-matter.pdf
        // Multiple scattering formulas
        //
        // http://pdg.lbl.gov/2017/AtomicNuclearProperties/index.html
        // Tables used to crosscheck computed radiation lenghts

        G4cout << G4endl
               << "Computing analytical scattering..." << G4endl;
        G4RunManager*           run    = G4RunManager::GetRunManager();

        DetectorConstruction* detCon = (DetectorConstruction*)run->GetUserDetectorConstruction();
        G4double targetThickness = detCon->getTargetThickness(); //Geant units
        G4cout << "targetThickness = " << targetThickness/mm << " [mm]" << G4endl;
        G4int targetZ = detCon->GetTargetMaterialZ();
        G4cout << "targetZ         = " << targetZ << " [e]" << G4endl;
        G4double targetA = detCon->GetTargetMaterialA();
        G4cout << "targetA         = " << targetA << " [amu]" << G4endl;

        G4double targetDensity = detCon->GetTargetMaterialDensity();
        G4cout << "targetDensity   = " << targetDensity*cm3/g << " [g/cm]" << G4endl;

        G4double radiationLength = 716.4 * targetA /
            (targetZ*(targetZ+1)*log(287.0/sqrt(targetZ))); //[g/cm^2]
        G4cout << "RadiationLength = " << radiationLength << " [g/cm^2]" << G4endl;
        radiationLength /= targetDensity*cm3/g; //[cm]
        radiationLength *= cm; //Geant4 units
        G4cout << "                = " << radiationLength/cm << " [cm]" << G4endl;

        PrimaryGeneratorAction* genAct = (PrimaryGeneratorAction*)run->GetUserPrimaryGeneratorAction();
        G4double beamMass   = genAct->get_beam_particlemass(); // Geant4 units
        G4cout << "beamMass        = " << beamMass / MeV << " [MeV/c^2]" <<G4endl;
        G4double beamCharge = genAct->get_beam_particlecharge(); // [e]
        G4cout << "beamCharge      = " << beamCharge << " [e]" << G4endl;
        G4double beamBeta2   = 1 - pow(beamMass/beamEnergy*MeV, 2);
        G4cout << "beamBeta2       = " << beamBeta2 << G4endl;

        G4double theta0 = 13.6*MeV / (beamEnergy*MeV - pow(beamMass,2)/beamEnergy*MeV) *
            fabs(beamCharge) * sqrt(targetThickness/radiationLength) *
            (1 + 0.038 * log(targetThickness*pow(beamCharge,2)/(radiationLength*beamBeta2)) );
        //TF1* highland = new TF1("highland", "Highland scattering distribution",0,90);
        G4cout << "theta0          = " << theta0 << " [rad]" << G4endl
               << "                = " << theta0/deg << " [deg]" << G4endl;
        G4cout << G4endl;

        target_exitangle_hist_cutoff->Write(); // Write the unaltered histogram

        //Plot the analytical scattering distribution and the histogram together
        TCanvas* c1 = new TCanvas("scatterPlot");
        target_exitangle_hist_cutoff->GetXaxis()->SetRangeUser(-3*theta0/deg,3*theta0/deg);
        target_exitangle_hist_cutoff->GetXaxis()->SetTitle("Angle [deg]");
        target_exitangle_hist_cutoff->Scale(1.0/target_exitangle_hist_cutoff->Integral(), "width");
        target_exitangle_hist_cutoff->Draw();

        G4double exitangle_analytic_xmin = target_exitangle_hist_cutoff->GetXaxis()->GetXmin();
        G4double exitangle_analytic_xmax = target_exitangle_hist_cutoff->GetXaxis()->GetXmax();
        G4double exitangle_analytic_xrange = exitangle_analytic_xmax - exitangle_analytic_xmin;
        int exitangle_analytic_npoints = target_exitangle_hist_cutoff->GetNbinsX();
        G4double* exitangle_analytic_x = new G4double[exitangle_analytic_npoints];
        G4double* exitangle_analytic_y = new G4double[exitangle_analytic_npoints];

        for (int i = 0; i < exitangle_analytic_npoints; i++) {
            exitangle_analytic_x[i] =
                exitangle_analytic_xrange/((G4double)exitangle_analytic_npoints-1) * i
                + exitangle_analytic_xmin;

            exitangle_analytic_y[i] = 1.0/sqrt(2.0*M_PI*pow(theta0/deg,2)) *
                exp( - pow(exitangle_analytic_x[i]/(theta0/deg),2) / 2.0);
        }
        TGraph* exitangle_analytic     = new TGraph(exitangle_analytic_npoints,
                                                    exitangle_analytic_x,
                                                    exitangle_analytic_y       );
        exitangle_analytic->Draw("same");
        /*
        // Debug
        G4cout << "integral = " << exitangle_analytic->Integral() << G4endl;
        G4cout << "integral = " << target_exitangle_hist_cutoff->Integral("width") << G4endl;
        */

        G4String plot_filename = foldername_out + "/" + filename_out + "_angles.png";
        c1->SaveAs(plot_filename);
        c1->Write();

        // Write the ROOT file.
        exitangle_analytic->Write();
        delete exitangle_analytic; exitangle_analytic = NULL;

        init_phasespaceX->Write();
        init_phasespaceY->Write();

        targetEdep->Write();
        targetEdep_NIEL->Write();
        targetEdep_IEL->Write();

        target_exitangle_hist->Write();

        tracker_phasespaceX->Write();
        tracker_phasespaceY->Write();

        tracker_phasespaceX_cutoff->Write();
        tracker_phasespaceY_cutoff->Write();
    }

    //Now we have plotted, delete stuff

    delete init_phasespaceX; init_phasespaceX = NULL;
    delete init_phasespaceY; init_phasespaceY = NULL;

    delete targetEdep; targetEdep = NULL;
    delete targetEdep_NIEL; targetEdep_NIEL = NULL;
    delete targetEdep_IEL; targetEdep_IEL = NULL;

    delete target_exitangle_hist; target_exitangle_hist = NULL;

    delete tracker_phasespaceX; tracker_phasespaceX = NULL;
    delete tracker_phasespaceY; tracker_phasespaceY = NULL;

    delete tracker_phasespaceX_cutoff; tracker_phasespaceX_cutoff = NULL;
    delete tracker_phasespaceY_cutoff; tracker_phasespaceY_cutoff = NULL;

    delete target_exitangle_hist_cutoff; target_exitangle_hist_cutoff = NULL;

    delete tracker_numParticles; tracker_numParticles = NULL;
    delete tracker_energy; tracker_energy = NULL;
    delete tracker_hitPos; tracker_hitPos = NULL;
    delete tracker_hitPos_cutoff; tracker_hitPos_cutoff = NULL;

    histFile->Write();
    histFile->Close();
    delete histFile; histFile = NULL;
}

void RootFileWriter::PrintTwissParameters(TH2D* phaseSpaceHist) {
    G4cout << "Stats for '" << phaseSpaceHist->GetTitle() << "':"  << G4endl;
    double stats[7];
    phaseSpaceHist->GetStats(stats);

    // Fill used [mm] and [rad]
    double posAve   = stats[2]/stats[0];
    double angAve   = stats[4]/stats[0];
    double posVar   = (stats[3] - stats[2]*stats[2]/stats[0]) / (stats[0]-1.0) ;
    double angVar   = (stats[5] - stats[4]*stats[4]/stats[0]) / (stats[0]-1.0) ;
    double coVar    = (stats[6] - stats[2]*stats[4]/stats[0]) / (stats[0]-1.0);

    G4cout << "numHits = "  << stats[0]
           << ", posAve = " << posAve     << " [mm]"
           << ", angAve = " << angAve     << " [rad]"
           << ", posVar = " << posVar     << " [mm^2]"
           << ", angVar = " << angVar     << " [rad^2]"
           << ", coVar  = " << coVar      << " [rad*mm]"
           << G4endl;

    G4RunManager*           run    = G4RunManager::GetRunManager();
    PrimaryGeneratorAction* genAct = (PrimaryGeneratorAction*)run->GetUserPrimaryGeneratorAction();
    double gamma_rel = (genAct->get_beam_energy() * MeV) / genAct->get_beam_particlemass();
    double beta_rel = sqrt(gamma_rel*gamma_rel - 1.0) / gamma_rel;

    double det = posVar*angVar - coVar*coVar;
    double epsG = sqrt(det);
    double epsN = epsG*beta_rel*gamma_rel;  //[mm*rad]
    double beta = posVar/epsG; // [mm]
    double alpha = -coVar/epsG;

    G4cout << "Geometrical emittance  = " << epsG*1e3 << " [um]" << G4endl;
    G4cout << "Normalized emittance   = " << epsN*1e3 << " [um]"
           << ", assuming beam energy = " << genAct->get_beam_energy() << " [MeV]"
           << ", and mass = " << genAct->get_beam_particlemass()/MeV << " [MeV/c]"
           << G4endl;
    G4cout << "Twiss beta  = " << beta*1e-3  << " [m]" << G4endl
           << "Twiss alpha = " << alpha << " [-]"  << G4endl;

    G4cout << G4endl;

    // Write to root file
    TVectorD twissVector (3);
    twissVector[0] = epsN*1e3;
    twissVector[1] = beta*1e-3;
    twissVector[2] = alpha;
    twissVector.Write((G4String(phaseSpaceHist->GetName())+"_TWISS").c_str());
}

void RootFileWriter::PrintParticleTypes(particleTypesCounter& pt, G4String name) {
    //Print out the particle types hitting the tracker
    G4cout << endl;
    G4cout << "Got types at " << name << ":" << G4endl;
    if (pt.numParticles == 0){
        G4cout << "No particles hit!" << G4endl;
        return;
    }

    TVectorD particleTypes_PDG    (pt.particleTypes.size());
    TVectorD particleTypes_numpart(pt.particleTypes.size());

    size_t particleTypes_i = 0;
    for(std::map<G4int,G4int>::iterator it = pt.particleTypes.begin(); it != pt.particleTypes.end(); it++){
        G4cout << std::setw(15) << it->first << " = "
               << std::setw(15) << pt.particleNames[it->first] << ": "
               << std::setw(15) << it->second << " = ";// << G4endl;
        G4int numHashes = (G4int) ((it->second / ((double)pt.numParticles)) * 100);
        for (int i = 0; i < numHashes; i++) G4cout << "#";
        G4cout << endl;

        //Also put them in the ROOT file
        // Unfortunately, there is no TObject array type for ints (?!?),
        // and I don't want  to depend on a ROOT dictionary file.
        particleTypes_PDG     [particleTypes_i] = int(it->first);
        particleTypes_numpart [particleTypes_i] = int(it->second);

        particleTypes_i++;
    }
    particleTypes_PDG.Write((name + "_ParticleTypes_PDG").c_str());
    particleTypes_numpart.Write((name + "_ParticleTypes_numpart").c_str());

}

void RootFileWriter::FillParticleTypes(particleTypesCounter& pt, G4int PDG, G4String type) {
    if (pt.particleTypes.count(PDG) == 0) {
        pt.particleTypes[PDG] = 0;
        pt.particleNames[PDG] = type;
    }
    pt.particleTypes[PDG] += 1;
    pt.numParticles       += 1;
}
