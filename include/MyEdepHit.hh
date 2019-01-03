#ifndef MyEdepHit_HH
#define MyEdepHit_HH

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4VHit.hh"
#include "G4ThreeVector.hh"

class G4AttDef;
class G4AttValue;

class MyEdepHit : public G4VHit {

public:

    // Constructors
    MyEdepHit() :
        fDepositedEnergy(0.0), fDepositedEnergy_NIEL(0.0),
        preStepPoint(0,0,0), postStepPoint(0,0,0) {};
    MyEdepHit(G4double Edep, G4double Edep_NIEL,
              G4ThreeVector preStepPoint_in,G4ThreeVector postStepPoint_in) :
        fDepositedEnergy(Edep), fDepositedEnergy_NIEL(Edep_NIEL),
        preStepPoint(preStepPoint_in), postStepPoint(postStepPoint_in) {};

    // Destructor
    virtual ~MyEdepHit() {};
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);

    // Methods
    virtual void Print() {};

    // Deposited energy
    inline void SetDepositedEnergy(G4double energy) {fDepositedEnergy = energy;}
    inline G4double GetDepositedEnergy() const {return fDepositedEnergy;}

    // Deposited energy
    inline void SetDepositedEnergy_NIEL(G4double energy) {fDepositedEnergy_NIEL = energy;}
    inline G4double GetDepositedEnergy_NIEL() const {return fDepositedEnergy_NIEL;}

    //Pre- and post-step points
    inline const G4ThreeVector& GetPreStepPoint()  const {return preStepPoint;}
    inline const G4ThreeVector& GetPostStepPoint() const {return postStepPoint;}

private:

    G4double fDepositedEnergy;      // Energy deposit
    G4double fDepositedEnergy_NIEL; // Energy deposit (NIEL)

    G4ThreeVector preStepPoint;
    G4ThreeVector postStepPoint;
};

typedef G4THitsCollection<MyEdepHit> MyEdepHitsCollection;

extern G4Allocator<MyEdepHit> MyEdepHitAllocator;

inline void* MyEdepHit::operator new(size_t) {
    void* aHit;
    aHit = (void*)MyEdepHitAllocator.MallocSingle();
    return aHit;
}

inline void MyEdepHit::operator delete(void* aHit) {
    MyEdepHitAllocator.FreeSingle((MyEdepHit*) aHit);
}

#endif
