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
    MyTrackerHit() :
        trackPosition(0,0,0), trackMomentum(0,0,0), trackEnergy(0.0), PDG(0), particleCharge(0.0) {};
    MyTrackerHit(G4ThreeVector position, G4ThreeVector momentum, G4double energy, G4int id, G4int charge) :
        trackPosition(position), trackMomentum(momentum), trackEnergy(energy), PDG(id), particleCharge(charge) {};

    // Destructor
    virtual ~MyTrackerHit() {};
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);

    // Methods
    virtual void Print() {};

    inline const G4ThreeVector& GetPosition() const {return trackPosition;}
    inline const G4ThreeVector& GetMomentum() const {return trackMomentum;}

    inline G4double GetTrackEnergy()    const {return trackEnergy;}
    inline G4int    GetPDG()            const {return PDG;}
    inline G4int    GetCharge()         const {return particleCharge;}

    inline void SetType(G4String type) {particleType = type;}
    inline const G4String& GetType() const {return particleType;}
private:

    G4ThreeVector trackPosition; //Global coordinates [G4 units]
    G4ThreeVector trackMomentum;

    G4double trackEnergy;
    G4int    PDG;
    G4int    particleCharge;

    G4String particleType;
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
