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

#include "VirtualTrackerWorldConstruction.hh"
//#include "MagnetClasses.hh"

#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"

#include "TrackerSD.hh"

VirtualTrackerWorldConstruction* VirtualTrackerWorldConstruction::singleton = NULL;

VirtualTrackerWorldConstruction::VirtualTrackerWorldConstruction(G4String worldName, DetectorConstruction* mainGeometryConstruction_in,
                                                                   std::vector <G4double>* trackerDistances_in, G4double trackerAngle_in)  :
    G4VUserParallelWorld(worldName), mainGeometryConstruction(mainGeometryConstruction_in) {

    trackerAngle = trackerAngle_in*deg;
    for (auto d : *trackerDistances_in) {
        trackerDistances.push_back(d*mm);
    }

    TrackerThickness = 1.0*um;
    G4double trackerAngle_pos = fabs(trackerAngle/rad); // Helper variable -- we assume positive angle in the calculations,
                                                        // but if negative the same restrictions would just come from the other side
    TrackerSizeX = ( mainGeometryConstruction->getWorldSizeX()/2.0 / cos(trackerAngle_pos) - trackerThickness * tan(trackerAngle_pos) ) * 2.0;

    TrackerSizeY = mainGeometryConstruction->getWorldSizeY();

    this->singleton = this;
}

VirtualTrackerWorldConstruction::~VirtualTrackerWorldConstruction() {

}

void VirtualTrackerWorldConstruction::Construct() {

    G4VPhysicalVolume* ghostWorld        = GetWorld();
    G4LogicalVolume*   ghostWorldLogical = ghostWorld->GetLogicalVolume();


    //Define the parallel geometries of the trackers

    G4Material* vacuumMaterial = G4Material::GetMaterial("G4_Galactic");
    if (not vacuumMaterial) {
        G4cerr << "Internal error -- material G4_Galactic not found in VirtualTrackerWorldConstruction::Construct()!" << G4endl;
        exit(1);
    }

    int idx = 0;
    for (auto dist : trackerDistances) {
        idx++;
        G4Box* solidTracker = new G4Box(G4String("TrackerS_")+std::to_string(idx), TrackerSizeX/2,TrackerSizeY/2,TrackerThickness/2);
        G4LogicalVolume* logicTracker = new G4LogicalVolume(solidTracker, vacuumMaterial, G4String("TrackerLV")+std::to_string(idx));
        G4RotationMatrix* trackerRot = new G4RotationMatrix();
        trackerRot->rotateY(trackerAngle);
        G4ThreeVector zTrans(0.0-sin(trackerAngle/rad)*TrackerThickness/2,
                             0.0, 
                             dist-cos(trackerAngle/rad)*TrackerThickness/2); //Center of tracker plane should be at the specified position
        G4VPhysicalVolume* physiTracker = 
                       new G4PVPlacement(G4Transform3D(*trackerRot,zTrans),          //Translate then rotate
                                         logicTracker,                               //its logical volume
                                         G4String("TrackerPV")+std::to_string(idx),  //its name
                                         ghostWorldLogical,                          //its mother
                                         false,                                      //pMany not used
                                         0,                                          //copy number
                                         true);                                      //Check for overlaps

        virtualTrackerLVs.push_back(logicTracker);
        virtualTrackerPVs.push_back(physiTracker);
    }

}

void VirtualTrackerWorldConstruction::ConstructSD() {
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    int idx = 0;
    for (auto tLV : virtualTrackerLVs) {
        idx++;
        G4VSensitiveDetector* trackerSD = new TrackerSD(G4String("tracker_")+std::to_string(idx));
        SDman->AddNewDetector(trackerSD);
        tLV->SetSensitiveDetector(trackerSD);
    }

}
