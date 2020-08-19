/*
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

#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

//#include <iostream>
//using namespace std;


TrackerSD::TrackerSD(const G4String& name) :
  G4VSensitiveDetector(name) {

  fHitsCollection = NULL; //Note: The HitsCollection objects are deleted by Geant4.
  collectionName.insert(name+"_exitpos");
  fHitsCollectionID = -1;

  G4RunManager*     run= G4RunManager::GetRunManager();
  detectorConstruction = (DetectorConstruction*)run->GetUserDetectorConstruction();
}

TrackerSD::~TrackerSD() {}

// Called at the beginning of each event
void TrackerSD::Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent) {

  //For every event, make a new hits collection
  fHitsCollection = new TrackerHitsCollection(SensitiveDetectorName, collectionName[0]);
  if (fHitsCollectionID < 0) {
    fHitsCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }

  // Add collection to the event
  hitsCollectionOfThisEvent->AddHitsCollection(fHitsCollectionID, fHitsCollection);

}

// Called each step in the scoring logical volume
G4bool TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
  //Only use incoming tracks
  if (aStep->GetPreStepPoint()->GetStepStatus()!=fGeomBoundary) {
    //G4cout << "TrackerSD::ProcessHits(): SKIP!" << G4endl;
    return false;
  }

  G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();
  const G4ThreeVector momentum = aStep->GetPreStepPoint()->GetMomentum();
  const G4ThreeVector& hitPos = aStep->GetPreStepPoint()->GetPosition();

  G4Track* theTrack = aStep->GetTrack();
  G4ParticleDefinition* particleType = theTrack->GetDefinition();
  G4int particleID = particleType->GetPDGEncoding();
  G4int particleCharge = particleType->GetPDGCharge();

  TrackerHit* aHit = new TrackerHit(hitPos, momentum, energy, particleID, particleCharge);
  aHit->SetType(particleType->GetParticleSubType());
  fHitsCollection->insert(aHit);

  return true;
}
