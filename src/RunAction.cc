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

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include "RootFileWriter.hh"
#include "PrimaryGeneratorAction.hh"



//--------------------------------------------------------------------------------

RunAction::RunAction(){
    // initialize root
}

//--------------------------------------------------------------------------------

RunAction::~RunAction() {
    // Delete histogram object
}

//--------------------------------------------------------------------------------

void RunAction::BeginOfRunAction(const G4Run*) {
    RootFileWriter::GetInstance()->initializeRootFile();
}


//--------------------------------------------------------------------------------

void RunAction::EndOfRunAction(const G4Run*) {
    RootFileWriter::GetInstance()->finalizeRootFile();

    G4RunManager*           run    = G4RunManager::GetRunManager();
    PrimaryGeneratorAction* genAct = (PrimaryGeneratorAction*)run->GetUserPrimaryGeneratorAction();
    genAct->endOfRun();
}

//--------------------------------------------------------------------------------
