#include "DetectorConstruction.hh"
#include "MyTargetSD.hh"
#include "MyTrackerSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4UniformMagField.hh"
#include "G4Element.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

#include <cmath>

//------------------------------------------------------------------------------

DetectorConstruction::DetectorConstruction(G4double TargetThickness_in,
                                           G4String TargetMaterial_in,
                                           G4double DetectorDistance_in,
                                           G4double DetectorAngle_in,
                                           G4bool   DetectorRotated_in) :
    AlMaterial(0), TargetMaterial(0),
    solidWorld(0),logicWorld(0),physiWorld(0),
    solidTarget(0),logicTarget(0),physiTarget(0),
    magField(0) {

    TargetThickness = TargetThickness_in*mm;
    DetectorThickness = 1*um;
    DetectorDistance = DetectorDistance_in*mm;
    DetectorRotated = DetectorRotated_in;
    DetectorAngle = DetectorAngle_in*pi/180.0;

    if (not DetectorRotated) { // No detector angle, make a simple 200x200cm world
        WorldSizeZ   = 200*cm;
        if (DetectorDistance > WorldSizeZ / 2.0) {
            //If neccessary, expand the z-size of the volume to fit the detector
            WorldSizeZ   = (DetectorDistance + DetectorThickness + 10*cm)*2.0;
        }

        WorldSizeX = 200*cm;
        WorldSizeY = 200*cm;

        DetectorSizeX = WorldSizeX;
        DetectorSizeY = WorldSizeY;

        TargetSizeX = WorldSizeX;
        TargetSizeY = WorldSizeY;
    }
    else { //Detector has angle -- make sure that everything fits together
        double theta = abs(DetectorAngle);
        // Offset of detector center vs. end of target
        double dz = DetectorDistance-TargetThickness/2.0;
        // Transverse size of volume needed to contain the detector
        double dx = dz/tan(theta);
        // Total length of detector volume (including part to be "chopped off" to fit the thickness)
        double rp = dz/sin(theta);

        // How much to cut off the end so that it doesn't crash with the wall when thickness > 0?
        double dr = 0.0;
        if (theta < pi/4.0) {
            dr = DetectorThickness/(2.0*tan(theta));
        }
        else if (theta > pi/4.0) {
            dr = (DetectorThickness*tan(theta))/2.0;
        }
        else if (theta == pi/4.0) {
            dr = DetectorThickness/2.0;
        }
        else if (theta > pi/2.0) {
            G4cerr << "Error: Detector angle  should be within +/- pi/2.0" << G4endl;
            exit(1);
        }

        DetectorSizeX = (rp-dr)*2.0;
        DetectorSizeY = DetectorSizeX;

        WorldSizeX = dx*2.0;
        WorldSizeY = DetectorSizeY;

        WorldSizeZ = 2*(TargetThickness/2+2*dz);

        TargetSizeX = WorldSizeX;
        TargetSizeY = WorldSizeY;
    }

    // materials
    DefineMaterials();
    SetTargetMaterial(TargetMaterial_in);
    DetectorMaterial = vacuumMaterial;
}

//------------------------------------------------------------------------------

G4VPhysicalVolume* DetectorConstruction::Construct() {
    // Clean old geometry, if any
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    // World volume
    solidWorld = new G4Box("WorldS", WorldSizeX/2.0, WorldSizeY/2.0, WorldSizeZ/2.0);
    logicWorld = new G4LogicalVolume(solidWorld, vacuumMaterial, "WorldLV");

    physiWorld = new G4PVPlacement(0,               //no rotation
                                   G4ThreeVector(), //World volume must be centered at the origin
                                   logicWorld,      //its logical volume
                                   "World",         //its name
                                   0,               //its mother  volume
                                   false,           //pMany not used
                                   0,               //copy number
                                   true);           //Check for overlaps

    //constructing the target
    solidTarget = new G4Box("TargetS", TargetSizeX/2,TargetSizeY/2, TargetThickness/2);
    logicTarget = new G4LogicalVolume(solidTarget, TargetMaterial,"TargetLV");
    physiTarget = new G4PVPlacement(NULL,                        //no rotation
                                    G4ThreeVector(0.0,0.0,0.0),  //its position
                                    logicTarget,                 //its logical volume
                                    "TargetPV",                  //its name
                                    logicWorld,                  //its mother
                                    false,                       //pMany not used
                                    0,                           //copy number
                                    true);                       //Check for overlaps

    //The "detector"
    solidDetector = new G4Box("DetectorS", DetectorSizeX/2,DetectorSizeY/2,DetectorThickness/2);
    logicDetector = new G4LogicalVolume(solidDetector, DetectorMaterial, "DetectorLV");

    if (DetectorRotated) {

        G4RotationMatrix* detectorRot = new G4RotationMatrix();
        detectorRot->rotateY(DetectorAngle*rad);
        G4ThreeVector zTrans(0.0,0.0,DetectorDistance);
        physiDetector = new G4PVPlacement(G4Transform3D(*detectorRot,zTrans),        //Translate then rotate
                                          logicDetector,                            //its logical volume
                                          "DetectorPV",                             //its name
                                          logicWorld,                               //its mother
                                          false,                                    //pMany not used
                                          0,                                        //copy number
                                          true);                                    //Check for overlaps

    }
    else{
        physiDetector = new G4PVPlacement(NULL,                                     //No rotation
                                          G4ThreeVector(0.0,0.0,DetectorDistance),  //its position
                                          logicDetector,                            //its logical volume
                                          "DetectorPV",                             //its name
                                          logicWorld,                               //its mother
                                          false,                                    //pMany not used
                                          0,                                        //copy number
                                          true);                                    //Check for overlaps
    }

    // Get pointer to detector manager
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    G4VSensitiveDetector* targetSD = new MyTargetSD("TargetSD_target");
    SDman->AddNewDetector(targetSD);
    logicTarget->SetSensitiveDetector(targetSD);
    G4VSensitiveDetector* detectorSD = new MyTrackerSD("TrackerSD_tracker");
    SDman->AddNewDetector(detectorSD);
    logicDetector->SetSensitiveDetector(detectorSD);

    return physiWorld;
}

//------------------------------------------------------------------------------

void DetectorConstruction::DefineMaterials() {
    G4NistManager* man = G4NistManager::Instance();
    man->SetVerbose(1);

    G4Material* Al = man->FindOrBuildMaterial("G4_Al");
    G4Material* C  = man->FindOrBuildMaterial("G4_C");
    G4Material* Cu = man->FindOrBuildMaterial("G4_Cu");
    G4Material* Pb = man->FindOrBuildMaterial("G4_Pb");
    G4Material* Ti = man->FindOrBuildMaterial("G4_Ti");

    G4Material* Vacuum = man->FindOrBuildMaterial("G4_Galactic");

    //default materials
    vacuumMaterial   = Vacuum;
    AlMaterial       = Al;
    CMaterial        = C;
    CuMaterial       = Cu;
    PbMaterial       = Pb;
    TiMaterial       = Ti;
}

//------------------------------------------------------------------------------

void DetectorConstruction::SetTargetMaterial(G4String materialChoice) {
    // search the material by its name
    G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
    if (pttoMaterial) TargetMaterial = pttoMaterial;
    else {
        G4cerr << "Error when setting material '"
               << materialChoice << "' -- not found!" << G4endl;
        exit(1);
    }
}

//------------------------------------------------------------------------------
