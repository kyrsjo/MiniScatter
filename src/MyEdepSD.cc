#include "MyEdepSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"    
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include <iostream>
#include "G4SystemOfUnits.hh"


using namespace std;


MyEdepSD::MyEdepSD(const G4String& name) :
 G4VSensitiveDetector(name) {
  
  //collectionName.insert("EdepCollection");
  //fHitsCollectionID = -1;
}

MyEdepSD::~MyEdepSD() {}


// Called at the beginning of each event
void MyEdepSD::Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent)
{
  //For every event, make a new hits collection
  fHitsCollection = new MyEdepHitsCollection(SensitiveDetectorName, collectionName[0]);
  if(fHitsCollectionID < -1){
    fHitsCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }
  
  // Add collection to the event
  hitsCollectionOfThisEvent->AddHitsCollection(fHitsCollectionID, fHitsCollection);
}

// Called each step in the scoring logical volume
G4bool MyEdepSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
  
  // Get energy deposited in this step
  G4double depositedEnergy = aStep->GetTotalEnergyDeposit()-aStep->GetNonIonizingEnergyDeposit();

  MyEdepHit* aHit = new MyEdepHit();
  aHit->AddDepositedEnergy(depositedEnergy);                  // Set energy deposit to the Hit object.

  fHitsCollection->insert(aHit);// Store (insert) the Hit object to hit collection container.
  return true;
}

void MyEdepSD::EndOfEvent(G4HCofThisEvent*) {}
