#ifndef MyMomentumHit_HH
#define MyMomentumHit_HH

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4VHit.hh"
#include "G4ThreeVector.hh"

class G4AttDef;
class G4AttValue;

class MyMomentumHit : public G4VHit {

public:
  
  // Constructors
  MyMomentumHit(G4ThreeVector momentum);

  // Destructor
  virtual ~MyMomentumHit();  
  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  // Methods
  virtual void Print();
  
  inline const G4ThreeVector& GetMomentum() const {return trackMomentum;}
  
private:
  
  G4ThreeVector trackMomentum;
};

typedef G4THitsCollection<MyMomentumHit> MyMomentumHitsCollection;

extern G4Allocator<MyMomentumHit> MyMomentumHitAllocator;

inline void* MyMomentumHit::operator new(size_t) {
  void* aHit;
  aHit = (void*)MyMomentumHitAllocator.MallocSingle();
  return aHit;
}

inline void MyMomentumHit::operator delete(void* aHit) {
  MyMomentumHitAllocator.FreeSingle((MyMomentumHit*) aHit);
}

#endif


