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

    collectionName.insert("EdepCollection");
    fHitsCollectionID_edep = -1;
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

}

// Called each step in the scoring logical volume
G4bool MyTargetSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {

    //Always do the energy deposit
    MyEdepHit* aHit_edep = new MyEdepHit();
    aHit_edep->SetDepositedEnergy(aStep->GetTotalEnergyDeposit());
    aHit_edep->SetDepositedEnergy_NIEL(aStep->GetNonIonizingEnergyDeposit());
    fHitsCollection_edep->insert(aHit_edep);

    return true;

}
