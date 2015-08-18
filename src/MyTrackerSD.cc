#include "MyTrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"    
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

#include <iostream>

//using namespace std;


MyTrackerSD::MyTrackerSD(const G4String& name) :
  G4VSensitiveDetector(name) {
  
  collectionName.insert("TrackerCollection");
  fHitsCollectionID = -1;
  
  G4RunManager*     run= G4RunManager::GetRunManager();
  detectorConstruction = (DetectorConstruction*)run->GetUserDetectorConstruction();
}

MyTrackerSD::~MyTrackerSD() {}

// Called at the beginning of each event
void MyTrackerSD::Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent) { 
  
  //For every event, make a new hits collection
  fHitsCollection = new MyTrackerHitsCollection(SensitiveDetectorName, collectionName[0]);
  if (fHitsCollectionID < 0) {
    fHitsCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }

  // Add collection to the event
  hitsCollectionOfThisEvent->AddHitsCollection(fHitsCollectionID, fHitsCollection);

}

// Called each step in the scoring logical volume
G4bool MyTrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
  G4Track* theTrack = aStep->GetTrack();
  G4double energy = theTrack->GetTotalEnergy();
  const G4ThreeVector& hitPos = aStep->GetPreStepPoint()->GetPosition();
  
  G4ParticleDefinition* particleType = theTrack->GetDefinition();
  G4int particleID = particleType->GetPDGEncoding();
  G4double angle = atan2(sqrt(hitPos.x()*hitPos.x() + hitPos.y()*hitPos.y()),detectorConstruction->GetDetectorDistance());

  MyTrackerHit* aHit = new MyTrackerHit(energy,angle,particleID);
  fHitsCollection->insert(aHit);
  
  return true;
}
