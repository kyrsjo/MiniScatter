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
#ifndef TargetSD_hh
#define TargetSD_hh 1

#include "G4VSensitiveDetector.hh"
#include "EdepHit.hh"
#include "TrackerHit.hh"

class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;

class TargetSD : public G4VSensitiveDetector {

public:

    TargetSD (const G4String& name);
    virtual ~TargetSD();

    // Methods
    virtual void Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent);

    virtual G4bool ProcessHits(G4Step* aStep,G4TouchableHistory* history);

    virtual void EndOfEvent(G4HCofThisEvent*) {};
private:

    // Data members
    EdepHitsCollection* fHitsCollection_edep;
    G4int fHitsCollectionID_edep;

    TrackerHitsCollection* fHitsCollection_exitpos;
    G4int fHitsCollectionID_exitpos;
};

#endif
