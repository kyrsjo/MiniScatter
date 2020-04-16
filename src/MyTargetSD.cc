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

#include "MyTargetSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

#include <iostream>

//using namespace std;


MyTargetSD::MyTargetSD(const G4String& name) :
    G4VSensitiveDetector(name) {

    fHitsCollection_edep = NULL; //Note: The HitsCollection objects are deleted by Geant4.
    collectionName.insert(name+"_edep");
    fHitsCollectionID_edep = -1;

    fHitsCollection_exitpos = NULL; //Note: The HitsCollection objects are deleted by Geant4.
    collectionName.insert(name+"_exitpos");
    fHitsCollectionID_exitpos = -1;
}

MyTargetSD::~MyTargetSD() {}

// Called at the beginning of each event, making the hits collections
void MyTargetSD::Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent) {

    //Energy deposits
    fHitsCollection_edep = new MyEdepHitsCollection(SensitiveDetectorName, collectionName[0]);
    if (fHitsCollectionID_edep < 0) {
        fHitsCollectionID_edep = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection_edep);
    }
    hitsCollectionOfThisEvent->AddHitsCollection(fHitsCollectionID_edep, fHitsCollection_edep);

    //Exit positions
    fHitsCollection_exitpos = new MyTrackerHitsCollection(SensitiveDetectorName, collectionName[1]);
    if (fHitsCollectionID_exitpos < 0) {
        fHitsCollectionID_exitpos = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection_exitpos);
    }
    hitsCollectionOfThisEvent->AddHitsCollection(fHitsCollectionID_exitpos, fHitsCollection_exitpos);
}

// Called each step in the scoring logical volume
G4bool MyTargetSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {

    //Always do the energy deposit
    MyEdepHit* aHit_edep = new MyEdepHit(aStep->GetTotalEnergyDeposit(),
                                         aStep->GetNonIonizingEnergyDeposit(),
                                         aStep->GetPreStepPoint()->GetPosition(),
                                         aStep->GetPostStepPoint()->GetPosition());
    fHitsCollection_edep->insert(aHit_edep);

    //Only use outgoing tracks
    if (aStep->GetPostStepPoint()->GetStepStatus()==fGeomBoundary) {
        G4double energy = aStep->GetPostStepPoint()->GetKineticEnergy();
        const G4ThreeVector momentum = aStep->GetPostStepPoint()->GetMomentum();
        const G4ThreeVector& hitPos = aStep->GetPostStepPoint()->GetPosition();

        G4Track* theTrack = aStep->GetTrack();
        G4ParticleDefinition* particleType = theTrack->GetDefinition();
        G4int particleID = particleType->GetPDGEncoding();
        G4int particleCharge = particleType->GetPDGCharge();

        MyTrackerHit* aHit = new MyTrackerHit(hitPos, momentum, energy, particleID, particleCharge);
        aHit->SetType(particleType->GetParticleSubType());
        fHitsCollection_exitpos->insert(aHit);
    }

    return true;

}
