#include "RootFileWriter.hh"

#include "MyEdepHit.hh"
#include "MyTrackerHit.hh"

#include "G4SDManager.hh"

#include "G4Track.hh"

#include "G4RunManager.hh"
#include "DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "CLHEP/Units/SystemOfUnits.h"

#include <iostream>
#include <iomanip>
#include <experimental/filesystem> //Mainstreamed from C++17

using namespace std;
RootFileWriter* RootFileWriter::singleton = 0;
const G4String  RootFileWriter::foldername_out = "plots";

void RootFileWriter::initializeRootFile(){
    //G4RunManager*         run    = G4RunManager::GetRunManager();
    //DetectorConstruction* detCon = (DetectorConstruction*)run->GetUserDetectorConstruction();

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
    histFile= new TFile(rootFileName,"RECREATE");

    targetEdep = new TH1D("targetEdep","targetEdep",1000,0,6);
    targetEdep->GetXaxis()->SetTitle("Total energy deposit/event [MeV]");
    targetEdep_NIEL = new TH1D("targetEdep_NIEL","targetEdep_NIEL",1000,0,1);
    targetEdep_NIEL->GetXaxis()->SetTitle("Total NIEL/event [keV]");
    targetEdep_IEL = new TH1D("targetEdep_IEL","targetEdep_IEL",1000,0,6);
    targetEdep_IEL->GetXaxis()->SetTitle("Total ionizing energy deposit/event [MeV]");

    // Tracker histograms
    tracker_numParticles = new TH1D("numParticles","numParticles",1001,-0.5,1000.5);
    tracker_numParticles->GetXaxis()->SetTitle("Number of particles / event");

    tracker_energy       = new TH1D("energy","energy",10000,0,10.0);
    tracker_energy->GetXaxis()->SetTitle("Energy of particles on the tracker [MeV]");

    //For counting the types of particles hitting the tracker
    tracker_particleTypes.clear();
    tracker_particleNames.clear();
    numParticles_total = 0;
}

void RootFileWriter::doEvent(const G4Event* event){

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
        else{
            G4cout << "targetEdepHitsCollection was NULL!"<<G4endl;
        }
    }
    else{
        G4cout << "myTargetEdep_CollID was " << myTargetEdep_CollID << " < 0!"<<G4endl;
    }

    //**Data from detectorTrackerSD**
    G4int myTrackerSD_CollID = SDman->GetCollectionID("TrackerCollection");
    if (myTrackerSD_CollID>=0){
        MyTrackerHitsCollection* trackerHitsCollection = NULL;
        trackerHitsCollection = (MyTrackerHitsCollection*) (HCE->GetHC(myTrackerSD_CollID));
        if (trackerHitsCollection != NULL) {
            G4int nEntries = trackerHitsCollection->entries();

            for (G4int i = 0; i < nEntries; i++){
                const G4double  energy = (*trackerHitsCollection)[i]->GetTrackEnergy();
                //const G4double  angle  = (*trackerHitsCollection)[i]->GetTrackAngle(); //UNUSED
                const G4int     PDG    = (*trackerHitsCollection)[i]->GetPDG();
                //const G4int     charge = (*trackerHitsCollection)[i]->GetCharge(); // UNUSED
                const G4String& type  = (*trackerHitsCollection)[i]->GetType();

                //Overall histograms
                tracker_energy->Fill(energy/TeV);

                //Particle type counting
                if (tracker_particleTypes.count(PDG) == 0){
                    tracker_particleTypes[PDG] = 0;
                    tracker_particleNames[PDG] = type;
                }
                tracker_particleTypes[PDG] += 1;
                numParticles_total += 1;
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


}
void RootFileWriter::finalizeRootFile(){
    targetEdep->Write();
    delete targetEdep; targetEdep = NULL;
    targetEdep_NIEL->Write();
    delete targetEdep_NIEL; targetEdep_NIEL = NULL;
    targetEdep_IEL->Write();
    delete targetEdep_IEL; targetEdep_IEL = NULL;

    tracker_numParticles->Write();
    delete tracker_numParticles; tracker_numParticles = NULL;
    tracker_energy->Write();
    delete tracker_energy; tracker_energy = NULL;

    //Print out the particle types hitting the tracker
    G4cout << "Got types at tracker:" << G4endl;
    for(std::map<G4int,G4int>::iterator it=tracker_particleTypes.begin(); it !=tracker_particleTypes.end(); it++){
        G4cout << std::setw(15) << it->first << " = "
               << std::setw(15) << tracker_particleNames[it->first] << ": "
               << std::setw(15) << it->second << " = ";// << G4endl;
        G4int numHashes = (G4int) ((it->second / ((double)numParticles_total)) * 100);
        for (int i = 0; i < numHashes; i++) G4cout << "#";
        G4cout << endl;
    }
    tracker_particleTypes.clear();

    histFile->Write();
    histFile->Close();
    delete histFile; histFile = NULL;
}
