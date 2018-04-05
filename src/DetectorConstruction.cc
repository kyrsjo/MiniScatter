#include "DetectorConstruction.hh"
#include "MyTargetSD.hh"
#include "MyTrackerSD.hh"

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4UniformMagField.hh"

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
#include <string>

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

    G4cout << G4endl;

    // Compute geometry

    TargetThickness = TargetThickness_in*mm;
    DetectorThickness = 1*um;
    DetectorDistance = DetectorDistance_in*mm;
    DetectorRotated = DetectorRotated_in;
    DetectorAngle = DetectorAngle_in*pi/180.0;

    //Set the z-size of the world volume to fit the detector + buffer,
    // in case of no rotation
    G4double WorldSizeZ_minimum = (DetectorDistance + DetectorThickness + WorldSizeZ_buffer)*2.0;

    if (not DetectorRotated) { // No detector angle, make a simple 200x200cm world

        WorldSizeZ = WorldSizeZ_minimum;

        G4double DetectorTargetDistance = (DetectorDistance - TargetThickness/2.0 - DetectorThickness/2.0);
        if (DetectorTargetDistance < 0.0) {
            G4cerr << "DetectorTargetDistance < 0.0 => Detector is inside target :(" << G4endl;
            exit(1);
        }
        G4cout << "Creating an unrotated detector; "
               << "distance target end to detector start = "
               << DetectorTargetDistance/mm << " [mm]" << G4endl;

        WorldSizeX = 5*cm;
        WorldSizeY = 5*cm;

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
        if (WorldSizeZ < WorldSizeZ_minimum) {
            WorldSizeZ = WorldSizeZ_minimum;
        }

        TargetSizeX = WorldSizeX;
        TargetSizeY = WorldSizeY;
    }

    // materials
    DefineMaterials(); // The standard ones

    if (TargetMaterial_in.contains("::")) {
        // Gas target
        DefineGas(TargetMaterial_in);
    }
    else {
        // Solid target
        SetTargetMaterial(TargetMaterial_in);
    }
    DetectorMaterial = vacuumMaterial;

    G4cout << G4endl;
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
    G4Material* Si = man->FindOrBuildMaterial("G4_Si");

    G4Material* Vacuum = man->FindOrBuildMaterial("G4_Galactic");

    //default materials
    vacuumMaterial   = Vacuum;
    AlMaterial       = Al;
    CMaterial        = C;
    CuMaterial       = Cu;
    PbMaterial       = Pb;
    TiMaterial       = Ti;
    SiMaterial       = Si;
}

//------------------------------------------------------------------------------

void DetectorConstruction::DefineGas(G4String TargetMaterial_in) {
    G4cout << G4endl;

    if (not TargetMaterial_in.contains("::")) {
        G4cerr << "Error in DetectorConstruction::DefineGas() -- no '::' was found"
               << "in material name = '" << TargetMaterial_in << "'"
               << G4endl;
        exit(1);
    }

    str_size colonPos = TargetMaterial_in.index("::");
    str_size pressurePos = colonPos+2;

    if (pressurePos >= TargetMaterial_in.length()) {
        G4cerr << "Error in DetectorConstruction::DefineGas() -- no pressure was found"
               << " after '::' in material name = '" << TargetMaterial_in << "'"
               << G4endl;
        exit(1);
    }
    G4String material_in = TargetMaterial_in(0, colonPos);
    G4String pressure_in = TargetMaterial_in(pressurePos, TargetMaterial_in.length());

    G4double pressure = 0.0;
    try {
        pressure = std::stod(std::string(pressure_in));
    }
    catch (const std::invalid_argument& ia) {
        G4cerr << "Invalid argument when reading pressure" << G4endl
               << "Got: '" << pressure_in << "'" << G4endl
               << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
        exit(1);
    }

    // ** Define the gas **
    // Compute properties
    constexpr G4double temperature = 300*kelvin;

    //Define materials
    G4double aHe       = 4.002602*g/mole;
    G4double densityHe = aHe*(pressure*bar*1e-3)/(temperature*Avogadro*CLHEP::k_Boltzmann);
    G4Isotope* isHe4 = new G4Isotope("He4", //Name
                                     2,     //iz
                                     4,     //n
                                     aHe);  //a
    G4Element* elHe  = new G4Element("Helium", //name
                                     "He",     //symbol
                                     1);       //ncomponents
    elHe->AddIsotope(isHe4, 1.0);
    this->gasHe = new G4Material("HeliumGas",        //name
                                 densityHe,          //density
                                 1,                  //ncomponents
                                 kStateGas,          //state
                                 temperature*kelvin, //temp
                                 pressure*bar*1e-3 );//pressure
    gasHe->AddElement(elHe, 1.0);
    G4cout << "Built He gas, pressure = "<< pressure
           << " [mbar], temperature = " << temperature
           << " [K], density = " << densityHe / g * meter3 << " [g/m3]"
           << G4endl;

    //Define materials
    G4double aAr       = 39.948*g/mole;
    G4double densityAr = aAr*(pressure*bar*1e-3)/(temperature*Avogadro*CLHEP::k_Boltzmann);
    G4Isotope* isAr40 = new G4Isotope("Ar40", //Name
                                     18,      //iz
                                     40,      //n
                                     aAr);    //a
    G4Element* elAr  = new G4Element("Argon", //name
                                     "Ar",    //symbol
                                     1);      //ncomponents
    elAr->AddIsotope(isAr40, 1.0);
    this->gasAr = new G4Material("ArgonGas",         //name
                                 densityAr,          //density
                                 1,                  //ncomponents
                                 kStateGas,          //state
                                 temperature*kelvin, //temp
                                 pressure*bar*1e-3 );//pressure
    gasAr->AddElement(elAr, 1.0);
    G4cout << "Built Ar gas, pressure = "<< pressure
           << " [mbar], temperature = " << temperature
           << " [K], density = " << densityAr / g * meter3 << " [g/m3]"
           << G4endl;

    if (material_in == "He") {
        TargetMaterial = this->gasHe;
    }
    else if (material_in == "Ar") {
        TargetMaterial = this->gasAr;
    }
    else {
        G4cerr << "Error in DetectorConstricton::DefineGas()" << G4endl;
        G4cerr << "Gas type '" << material_in << "' unknown." << G4endl;
        exit(1);
    }

    G4cout << G4endl;
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
