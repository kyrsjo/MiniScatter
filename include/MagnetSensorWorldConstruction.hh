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

#ifndef MagnetSensorWorldConstruction_hh
#define ParallelWolrdConstruction_hh 1

#include "G4VUserParallelWorld.hh"

#include "DetectorConstruction.hh"

// Parallel geometry in which the magnet sensitive detectors are defined;
// this gets around the problem of how to attach the SDs correctly onto complex magnet geometries.

class MagnetSensorWorldConstruction : public G4VUserParallelWorld {
public:
    MagnetSensorWorldConstruction(G4String worldName, DetectorConstruction* mainGeometryConstruction_in);
    virtual ~MagnetSensorWorldConstruction();

    virtual void Construct();
    virtual void ConstructSD();

private:
    DetectorConstruction* mainGeometryConstruction = NULL;
    std::vector <G4VPhysicalVolume*> magnetDetectorPVs;
};

#endif
