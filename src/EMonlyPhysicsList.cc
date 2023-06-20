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

#include "EMonlyPhysicsList.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include "G4EmStandardPhysicsSS.hh"

EMonlyPhysicsList::EMonlyPhysicsList(G4String EMlistName) :
    G4VModularPhysicsList() {

  if (EMlistName == "EM") {
      RegisterPhysics( new G4EmStandardPhysics() );
  }
  else if(EMlistName == "EMV") {
      RegisterPhysics( new G4EmStandardPhysics_option1() );
  }
  else if(EMlistName == "EMX") {
      RegisterPhysics( new G4EmStandardPhysics_option2() );
  }
  else if(EMlistName == "EMY") {
      RegisterPhysics( new G4EmStandardPhysics_option3() );
  }
  else if(EMlistName == "EMZ") {
      RegisterPhysics( new G4EmStandardPhysics_option4() );
  }
  else if(EMlistName == "LIV") {
      RegisterPhysics( new G4EmLivermorePhysics() );
  }
  else if(EMlistName == "PEN") {
      RegisterPhysics( new G4EmPenelopePhysics() );
  }
  else if(EMlistName == "GS") {
      RegisterPhysics( new G4EmStandardPhysicsGS() );
  }
  else if(EMlistName == "LE") {
      RegisterPhysics( new G4EmLowEPPhysics() );
  }
  else if(EMlistName == "WVI") {
      RegisterPhysics( new G4EmStandardPhysicsWVI() );
  }
  else if(EMlistName == "SS") {
      RegisterPhysics( new G4EmStandardPhysicsSS() );
  }
  else {
      G4String errormessage = "No EM physics list named '" + EMlistName + "' was found";
      G4Exception("EMonlyPhysicsList::EMonlyPhysicsList()", "MSPhysList2000", FatalException,
		  errormessage);
  }

  //

}
