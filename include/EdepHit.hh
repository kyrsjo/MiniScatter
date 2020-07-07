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
#ifndef EdepHit_HH
#define EdepHit_HH

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4VHit.hh"
#include "G4ThreeVector.hh"

class G4AttDef;
class G4AttValue;

class EdepHit : public G4VHit {

public:

    // Constructors
    EdepHit() :
        fDepositedEnergy(0.0), fDepositedEnergy_NIEL(0.0),
        preStepPoint(0,0,0), postStepPoint(0,0,0) {};
    EdepHit(G4double Edep, G4double Edep_NIEL,
              G4ThreeVector preStepPoint_in,G4ThreeVector postStepPoint_in) :
        fDepositedEnergy(Edep), fDepositedEnergy_NIEL(Edep_NIEL),
        preStepPoint(preStepPoint_in), postStepPoint(postStepPoint_in) {};

    // Destructor
    virtual ~EdepHit() {};
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

typedef G4THitsCollection<EdepHit> EdepHitsCollection;

extern G4Allocator<EdepHit> EdepHitAllocator;

inline void* EdepHit::operator new(size_t) {
    void* aHit;
    aHit = (void*)EdepHitAllocator.MallocSingle();
    return aHit;
}

inline void EdepHit::operator delete(void* aHit) {
    EdepHitAllocator.FreeSingle((EdepHit*) aHit);
}

#endif
