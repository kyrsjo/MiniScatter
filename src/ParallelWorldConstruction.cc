#include "ParallelWorldConstruction.hh"
#include "MagnetClasses.hh"

#include "G4PVPlacement.hh"

ParallelWorldConstruction::ParallelWorldConstruction(G4String worldName, DetectorConstruction* mainGeometryConstruction_in)  :
    G4VUserParallelWorld(worldName), mainGeometryConstruction(mainGeometryConstruction_in) {

}

ParallelWorldConstruction::~ParallelWorldConstruction() {}

void ParallelWorldConstruction::Construct() {
    //Define the parallel geometries of the magnets

    G4VPhysicalVolume* ghostWorld        = GetWorld();
    G4LogicalVolume*   ghostWorldLogical = ghostWorld->GetLogicalVolume();

    for (auto magnet : mainGeometryConstruction->magnets) {
        //Pretty much a copy of what goes on in DetectorConstruction::Construct()
        G4VPhysicalVolume* magnetDetectorPV   =
            new G4PVPlacement(NULL,
                              G4ThreeVector(magnet->GetXOffset(),
                                            magnet->GetYOffset(),
                                            magnet->getZ0()      ),
                              magnet->GetDetectorLV(),
                              magnet->magnetName + "_detectorPV",
                              ghostWorldLogical,
                              false,
                              0,
                              true);
        magnetDetectorPVs.push_back(magnetDetectorPV);

    }
}

void ParallelWorldConstruction::ConstructSD() {

    for (auto magnet : mainGeometryConstruction->magnets) {
        magnet->AddSD();
    }
}
