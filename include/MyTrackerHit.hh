#ifndef MyTrackerHit_HH
#define MyTrackerHit_HH

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4VHit.hh"
#include "G4ThreeVector.hh"

class G4AttDef;
class G4AttValue;

class MyTrackerHit : public G4VHit {

public:
  
  // Constructors
  MyTrackerHit();
  MyTrackerHit(G4double energy, G4double angle, G4int id, G4ThreeVector momentum);

  // Destructor
  virtual ~MyTrackerHit();  
  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  // Methods
  virtual void Print();
  
  inline G4double GetTrackEnergy() const {return trackEnergy;}
  inline G4double GetTrackAngle()  const {return trackAngle;}
  inline G4int    GetPDG()         const {return PDG;}
  inline const G4ThreeVector& GetMomentum() const {return trackMomentum;}
  
private:
  
  G4double trackEnergy;
  G4double trackAngle;
  G4int    PDG;
  
  G4ThreeVector trackMomentum;
};

typedef G4THitsCollection<MyTrackerHit> MyTrackerHitsCollection;

extern G4Allocator<MyTrackerHit> MyTrackerHitAllocator;

inline void* MyTrackerHit::operator new(size_t) {
  void* aHit;
  aHit = (void*)MyTrackerHitAllocator.MallocSingle();
  return aHit;
}

inline void MyTrackerHit::operator delete(void* aHit) {
  MyTrackerHitAllocator.FreeSingle((MyTrackerHit*) aHit);
}

#endif


