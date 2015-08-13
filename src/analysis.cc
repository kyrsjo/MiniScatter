#include "analysis.hh"
#include "AntiPHit.hh"
#include "DetectorConstruction.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4DigiManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "TStyle.h"
#include <iostream>
#include <sstream>
#include "TCanvas.h"
#include "TrackingAction.hh"
#include "TPaveText.h"
using namespace std;
analysis* analysis::singleton = 0;

//Constructor: it initializes parameters
analysis::analysis(){;
}
void analysis::makeHistograms(){
  // G4RunManager*run= G4RunManager::GetRunManager();
  // DetectorConstruction*DetC= (DetectorConstruction*)run->GetUserDetectorConstruction();

  histFile= new TFile("plots/histo.root","RECREATE");
}

void analysis::writePerEvent(const G4Event* event){

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4int AntiPCollectionID = SDman->GetCollectionID("AntiPCollection");
  G4HCofThisEvent * HCE=event->GetHCofThisEvent();
  AntiPHitsCollection* hc=0;
  hc=(AntiPHitsCollection*)(HCE->GetHC(AntiPCollectionID));
  if(hc==NULL){
    G4cout<<"NULL"<<G4endl;
  }
  G4int nEntries = hc->entries();
  energyHisto=new TH2D("EnergyGeant","EnergyGeant",256,-0.703999,0.703999,256,-0.7039999,0.703999);
  for(G4int itr  = 0 ; itr < nEntries ; itr++) {
    energyHisto->Fill((*hc)[itr]->GetHitPosition().z()/cm,(*hc)[itr]->GetHitPosition().y()/cm,(*hc)[itr]->GetDepositedEnergy()/cm);
    G4cout<<"position "<<(*hc)[itr]->GetHitPosition().x()/cm<<"  "<<(*hc)[itr]->GetHitPosition().y()/cm<<"  "<<(*hc)[itr]->GetHitPosition().z()/cm<<G4endl;
  }
  energyHisto->Write();
  delete energyHisto;
  
}
void analysis::writeHistograms(){
  histFile->Write();
  histFile->Close();
}

