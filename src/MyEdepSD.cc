#include "MyEdepSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"    
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>

//using namespace std;


MyEdepSD::MyEdepSD(const G4String& name) :
 G4VSensitiveDetector(name) {
  
  collectionName.insert("EdepCollection");
  fHitsCollectionID = -1;
}

MyEdepSD::~MyEdepSD() {}

// Called at the beginning of each event
void MyEdepSD::Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent) { 
  
  //For every event, make a new hits collection
  fHitsCollection = new MyEdepHitsCollection(SensitiveDetectorName, collectionName[0]);
  if (fHitsCollectionID < 0) {
    fHitsCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }

  // Add collection to the event
  hitsCollectionOfThisEvent->AddHitsCollection(fHitsCollectionID, fHitsCollection);

}

// Called each step in the scoring logical volume
G4bool MyEdepSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
  
  MyEdepHit* aHit = new MyEdepHit();
  aHit->AddDepositedEnergy(aStep->GetTotalEnergyDeposit());
  aHit->AddDepositedEnergy_NIEL(aStep->GetNonIonizingEnergyDeposit());
  //std::cout<< "Creating hit! edep="<< aStep->GetTotalEnergyDeposit()/GeV << "[GeV]" << std::endl;

  fHitsCollection->insert(aHit);
  
  return true;
}
