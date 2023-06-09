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
#ifndef EMonlyPhysicsList_hh
#define EMonlyPhysicsList_hh 1

#include "G4VModularPhysicsList.hh"

/* 
 * The purpose of this class is to make it possible to build a physics list
 * which only contains the electromagnetic part, no hadronics.
 * This is not really correct, but will be used for testing multiple scattering.
 */

class EMonlyPhysicsList : public G4VModularPhysicsList {
public:
     EMonlyPhysicsList(G4String EMlistName);
};

#endif
