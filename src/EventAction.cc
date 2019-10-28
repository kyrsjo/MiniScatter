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

#include "EventAction.hh"
#include "G4DigiManager.hh"
#include "RootFileWriter.hh"
#include "RunAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

#include <iostream>
#include <fstream>

using namespace std;

//------------------------------------------------------------------------------

EventAction::EventAction(RunAction* run) : runAct(run) {}


//------------------------------------------------------------------------------

void EventAction::EndOfEventAction(const G4Event* event) {
    RootFileWriter::GetInstance()->doEvent(event);

    G4int eventID = event->GetEventID();
    if (eventID % 500 == 0) {
        G4cout << "Event# "<<event->GetEventID() << G4endl;
    }
}

//------------------------------------------------------------------------------
