#include "analysis.hh"
//#include "AntiPHit.hh"
#include "MyEdepHit.hh"
#include "MyTrackerHit.hh"
//#include "DetectorConstruction.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include <iostream>

//#include "TDatabasePDG.h"

using namespace std;
analysis* analysis::singleton = 0;

void analysis::makeHistograms(){
  histFile= new TFile("plots/histo.root","RECREATE");
  
  targetEdep = new TH1D("targetEdep","targetEdep",1000,0,3);
  targetEdep->GetXaxis()->SetTitle("Total energy deposit/event [MeV]");
  targetEdep_NIEL = new TH1D("targetEdep_NIEL","targetEdep_NIEL",1000,0,1);
  targetEdep_NIEL->GetXaxis()->SetTitle("Total NIEL/event [MeV]");
  targetEdep_IEL = new TH1D("targetEdep_IEL","targetEdep_IEL",1000,0,3);
  targetEdep_IEL->GetXaxis()->SetTitle("Total ionizing energy deposit/event [MeV]");

  tracker_numParticles = new TH1D("numParticles","numParticles",1000,0,1000);
  tracker_angle        = new TH1D("angle","angle",1000,0,3.14/3.0);
}

void analysis::writePerEvent(const G4Event* event){
  
  G4HCofThisEvent* HCE=event->GetHCofThisEvent();
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  //**Data from targetEdepSD**
  G4int myTargetEdepSD_CollID = SDman->GetCollectionID("EdepCollection");
  if (myTargetEdepSD_CollID>=0){
    MyEdepHitsCollection* targetEdepHitsCollection = NULL;
    targetEdepHitsCollection = (MyEdepHitsCollection*) (HCE->GetHC(myTargetEdepSD_CollID));
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
      targetEdep_NIEL->Fill(edep_NIEL/MeV);
      targetEdep_IEL->Fill(edep_IEL/MeV);
    }
    else{
      cout << "targetEdepHitsCollection was NULL!"<<endl;
    }
  }
  else{
    cout << "myTargetEdepSD_CollID was " << myTargetEdepSD_CollID << "<0!"<<endl;
  }

  //**Data from detectorTrackerSD**
  G4int myTrackerSD_CollID = SDman->GetCollectionID("TrackerCollection");
  if (myTargetEdepSD_CollID>=0){
    MyTrackerHitsCollection* trackerHitsCollection = NULL;
    trackerHitsCollection = (MyTrackerHitsCollection*) (HCE->GetHC(myTrackerSD_CollID));
    if (trackerHitsCollection != NULL) {
      G4int nEntries = trackerHitsCollection->entries();
      // G4double edep      = 0.0;
      // G4double edep_NIEL = 0.0;
      // G4double edep_IEL  = 0.0;
      for (G4int i = 0; i < nEntries; i++){
      	
	tracker_angle->Fill((*trackerHitsCollection)[i]->GetTrackAngle());

	G4int PDG = (*trackerHitsCollection)[i]->GetPDG();
	//cout << (*trackerHitsCollection)[i]->GetPDG() << endl;
	if (tracker_particleTypes.count(PDG) == 0){
	  tracker_particleTypes[PDG] = 0;
	}
	tracker_particleTypes[PDG] += 1;
      }
      tracker_numParticles->Fill(nEntries);
    }
    else{
      cout << "trackerHitsCollection was NULL!"<<endl;
    }
  }
  else{
    cout << "myTrackerSD_CollID was " << myTrackerSD_CollID << "<0!"<<endl;
  }

  
}
void analysis::writeHistograms(){
  targetEdep->Write();
  delete targetEdep;
  targetEdep_NIEL->Write();
  delete targetEdep_NIEL;
  targetEdep_IEL->Write();
  delete targetEdep_IEL;

  tracker_numParticles->Write();
  delete tracker_numParticles;
  
  tracker_angle->Write();
  delete tracker_angle;

  //TDatabasePDG PDG_DB();
  cout << "Got types at tracker:" << endl;
  for(std::map<G4int,G4int>::iterator it=tracker_particleTypes.begin(); it !=tracker_particleTypes.end(); it++){
    cout << it->first << " " << it->second << endl;
  }
  tracker_particleTypes.clear();

  histFile->Write();
  histFile->Close();
}

