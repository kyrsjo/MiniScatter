#ifndef TargetSD
#define TargetSD

#include "G4VSensitiveDetector.hh"
#include "MyMomentumHit.hh"
#include "MyEdepHit.hh"

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
  MyMomentumHitsCollection* fHitsCollection_in;
  G4int fHitsCollectionID_in;

  MyMomentumHitsCollection* fHitsCollection_out;
  G4int fHitsCollectionID_out;

  MyEdepHitsCollection* fHitsCollection_edep;
  G4int fHitsCollectionID_edep;
};

#endif

