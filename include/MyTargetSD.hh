#ifndef TargetSD
#define TargetSD

#include "G4VSensitiveDetector.hh"
#include "MyEdepHit.hh"
#include "MyTrackerHit.hh"

class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;

class MyTargetSD : public G4VSensitiveDetector {

public:

    MyTargetSD (const G4String& name);
    virtual ~MyTargetSD();

    // Methods
    virtual void Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent);

    virtual G4bool ProcessHits(G4Step* aStep,G4TouchableHistory* history);

    virtual void EndOfEvent(G4HCofThisEvent*) {};
private:

    // Data members
    MyEdepHitsCollection* fHitsCollection_edep;
    G4int fHitsCollectionID_edep;

    MyTrackerHitsCollection* fHitsCollection_exitpos;
    G4int fHitsCollectionID_exitpos;
};

#endif
