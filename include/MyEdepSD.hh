#ifndef EdepSD
#define EdepSD

#include "G4VSensitiveDetector.hh"
#include "MyEdepHit.hh"

class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;

class MyEdepSD : public G4VSensitiveDetector {

public:
  
  MyEdepSD (const G4String& name);
  virtual ~MyEdepSD();

  
  // Methods
  virtual void Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent);
  
  virtual G4bool ProcessHits(G4Step* aStep,G4TouchableHistory* history);
  
  virtual void EndOfEvent(G4HCofThisEvent* hitsCollectionOfThisEvent);
  
 private:
  
  // Data members
  MyEdepHitsCollection* fHitsCollection;
  
  G4int fHitsCollectionID;
};

#endif

