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

#include "MagnetSensorWorldConstruction.hh"
#include "MagnetClasses.hh"

#include "G4PVPlacement.hh"

MagnetSensorWorldConstruction::MagnetSensorWorldConstruction(G4String worldName, DetectorConstruction* mainGeometryConstruction_in)  :
    G4VUserParallelWorld(worldName), mainGeometryConstruction(mainGeometryConstruction_in) {

}

MagnetSensorWorldConstruction::~MagnetSensorWorldConstruction() {}

void MagnetSensorWorldConstruction::Construct() {
    //Define the parallel geometries of the magnets

    G4VPhysicalVolume* ghostWorld        = GetWorld();
    G4LogicalVolume*   ghostWorldLogical = ghostWorld->GetLogicalVolume();

    for (auto magnet : mainGeometryConstruction->magnets) {
        //Pretty much a copy of what goes on in DetectorConstruction::Construct()
        if (magnet->GetUseGhostDetector()) {
            G4VPhysicalVolume* magnetDetectorPV   =
                new G4PVPlacement(magnet->GetMainPV_transform(),
                                  magnet->GetDetectorLV(),
                                  magnet->magnetName + "_detectorPV",
                                  ghostWorldLogical,
                                  false,
                                  0,
                                  true);
            magnetDetectorPVs.push_back(magnetDetectorPV);
        }
        else {
            magnetDetectorPVs.push_back(NULL);
        }

    }
}

void MagnetSensorWorldConstruction::ConstructSD() {

    for (auto magnet : mainGeometryConstruction->magnets) {
        magnet->AddSD();
    }
}
