#ifndef TrackerSD
#define TrackerSD

#include "DetectorConstruction.hh"
#include "G4VSensitiveDetector.hh"
#include "MyTrackerHit.hh"

class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;

class MyTrackerSD : public G4VSensitiveDetector {

public:

    MyTrackerSD (const G4String& name);
    virtual ~MyTrackerSD();

    // Methods
    virtual void Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent);

    virtual G4bool ProcessHits(G4Step* aStep,G4TouchableHistory* history);

    virtual void EndOfEvent(G4HCofThisEvent*) {};
private:

    DetectorConstruction* detectorConstruction;

    // Data members
    MyTrackerHitsCollection* fHitsCollection;
    G4int fHitsCollectionID;
};

#endif
