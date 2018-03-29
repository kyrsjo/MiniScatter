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

    collectionName.insert("MomentumCollection_in");
    fHitsCollectionID_in = -1;

    collectionName.insert("MomentumCollection_out");
    fHitsCollectionID_out = -1;

    collectionName.insert("EdepCollection");
    fHitsCollectionID_edep = -1;
}

MyTargetSD::~MyTargetSD() {}

// Called at the beginning of each event, making the hits collections
void MyTargetSD::Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent) {

    //Inbound momentum
    fHitsCollection_in = new MyMomentumHitsCollection(SensitiveDetectorName, collectionName[0]);
    if (fHitsCollectionID_in < 0) {
        fHitsCollectionID_in = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection_in);
    }
    hitsCollectionOfThisEvent->AddHitsCollection(fHitsCollectionID_in, fHitsCollection_in);

    //Outbound momentum
    fHitsCollection_out = new MyMomentumHitsCollection(SensitiveDetectorName, collectionName[1]);
    if (fHitsCollectionID_out < 0) {
        fHitsCollectionID_out = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection_out);
    }
    hitsCollectionOfThisEvent->AddHitsCollection(fHitsCollectionID_out, fHitsCollection_out);

    //Energy deposits
    fHitsCollection_edep = new MyEdepHitsCollection(SensitiveDetectorName, collectionName[2]);
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

    //Momentum: Inbound?
    if (aStep->GetPreStepPoint()->GetStepStatus()==fGeomBoundary) {
        const G4ThreeVector momentum = aStep->GetPreStepPoint()->GetMomentum();
        MyMomentumHit* aHit = new MyMomentumHit(momentum);
        fHitsCollection_in->insert(aHit);
        //G4cout << "Pre" << G4endl;
    }

    //Momentum: Outbound?
    if (aStep->GetPostStepPoint()->GetStepStatus()==fGeomBoundary) {
        const G4ThreeVector momentum = aStep->GetPostStepPoint()->GetMomentum();
        MyMomentumHit* aHit = new MyMomentumHit(momentum);
        fHitsCollection_out->insert(aHit);
        //G4cout << "Post" << G4endl;
    }

    return true;

}
