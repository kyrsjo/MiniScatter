#include "analysis.hh"
#include "AntiPHit.hh"
#include "MyEdepHit.hh"
#include "DetectorConstruction.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4DigiManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
//#include "TStyle.h"
#include <iostream>
#include <sstream>
//#include "TCanvas.h"
//#include "TPaveText.h"
using namespace std;
analysis* analysis::singleton = 0;

//Constructor: it initializes parameters
analysis::analysis(){;
}
void analysis::makeHistograms(){
  // G4RunManager*run= G4RunManager::GetRunManager();
  // DetectorConstruction*DetC= (DetectorConstruction*)run->GetUserDetectorConstruction();

  histFile= new TFile("plots/histo.root","RECREATE");
  
  targetEdep = new TH1D("targetEdep","targetEdep",1000,0,3);
  targetEdep->GetXaxis()->SetTitle("Total energy deposit/event [MeV]");
  targetEdep_NIEL = new TH1D("targetEdep_NIEL","targetEdep_NIEL",1000,0,1);
  targetEdep_NIEL->GetXaxis()->SetTitle("Total NIEL/event [MeV]");
  targetEdep_IEL = new TH1D("targetEdep_IEL","targetEdep_IEL",1000,0,3);
  targetEdep_IEL->GetXaxis()->SetTitle("Total ionizing energy deposit/event [MeV]");
}

void analysis::writePerEvent(const G4Event* event){
  
  G4HCofThisEvent* HCE=event->GetHCofThisEvent();
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  //Data from the AntiPCollection
  
  G4int AntiPCollectionID = SDman->GetCollectionID("AntiPCollection");
  AntiPHitsCollection* hc=NULL;
  hc=(AntiPHitsCollection*)(HCE->GetHC(AntiPCollectionID));
  if(hc==NULL){
    G4cout << "hc was NULL!"<<endl;
  }
  else{
    G4int nEntries = hc->entries();
    energyHisto=new TH2D("EnergyGeant","EnergyGeant",256,-0.703999,0.703999,256,-0.7039999,0.703999);
    for(G4int itr  = 0 ; itr < nEntries ; itr++) {
      energyHisto->Fill((*hc)[itr]->GetHitPosition().z()/cm,(*hc)[itr]->GetHitPosition().y()/cm,(*hc)[itr]->GetDepositedEnergy()/cm);
    }
    energyHisto->Write();
    delete energyHisto;
  }

  //Data from targetEdepSD
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
  
}
void analysis::writeHistograms(){
  targetEdep->Write();
  delete targetEdep;
  targetEdep_NIEL->Write();
  delete targetEdep_NIEL;
  targetEdep_IEL->Write();
  delete targetEdep_IEL;
  
  histFile->Write();
  histFile->Close();
}

