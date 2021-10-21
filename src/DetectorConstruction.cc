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

#include "DetectorConstruction.hh"
#include "TargetSD.hh"
#include "TrackerSD.hh"

#include "MagnetClasses.hh"

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

#include "G4VisAttributes.hh"

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
                                           G4double TargetAngle_in,
                                           G4String BackgroundMaterial_in,
                                           G4double WorldSize_in,
                                           G4double WorldMinLength_in,
                                           std::vector <G4String> &magnetDefinitions_in) :
    solidWorld(0),logicWorld(0),physiWorld(0),
    solidTarget(0),logicTarget(0),physiTarget(0),
    magnetDefinitions(magnetDefinitions_in) {

    // Process input

    TargetThickness = TargetThickness_in*mm;
    TargetAngle = TargetAngle_in*deg;

    // Helper variable for target  -- we assume positive angle in the calculations,
    // but if negative the same restrictions would just come from the other side
    G4double TargetAngle_pos = fabs(TargetAngle/rad); 

    if ( fabs(TargetAngle_in) > 45.0) {
        G4cerr << "abs(TargetAngle) too close to or above 45 degrees, current setting = " << TargetAngle_in << " [deg]" << G4endl;
        G4cerr << "This would result in a very or infinitely long volume or inverted left/right axis." << G4endl;
        exit(1);
    }

    WorldSizeX = WorldSize_in*mm;
    WorldSizeY = WorldSize_in*mm;

    if (TargetThickness == 0.0) {
        /*
        // Makes sense to allow no targets since we can have non-vacuum background materials...
        if (magnetDefinitions.size() == 0) {
            G4cerr << "Error: Object/magnet definitions must be used if TargetThickness=0." << G4endl;
            exit(1);
        }
        */
        HasTarget = false;
    }
    else if (TargetThickness < 0.0) {
        G4cerr << "Error in DetectorConstruction::DetectorConstruction():" << G4endl
               << "TargetThickness = " << TargetThickness/mm << " < 0.0; this is not allowed." << G4endl;
        exit(1);
    }
    else {
        HasTarget = true;
    }
    //Find the size of the target
    TargetSizeX = ( WorldSizeX/2.0/cos(TargetAngle_pos) - tan(TargetAngle_pos) * TargetThickness/2 ) * 2.0;
    if (TargetSizeX < 0.0) {
        G4cerr << "Error: Target thickness too large for World Size and angle" << G4endl;
        exit(1);
    }

    TargetSizeY = WorldSizeY;

    //Find the minimum length of the world
    WorldSizeZ = WorldMinLength_in*mm;
    G4double maxZ_target = sqrt(TargetSizeX*TargetSizeX + TargetThickness*TargetThickness)/2.0; // Usually an overestimate, but safe
    if (maxZ_target > WorldSizeZ ) { WorldSizeZ = maxZ_target; }
    // Times 2 and round up to nearest 10 mm
    WorldSizeZ = ceil((2*WorldSizeZ/mm)/10.0)*10.0*mm;


    G4cout << "World construction details:" << G4endl
           << "World size  X/Y/Z [mm] : " << WorldSizeX/mm  << ", " << WorldSizeY/mm  << ", " << WorldSizeZ/mm      << ";" << G4endl;
    G4cout << "Target size X/Y/Z [mm] : " << TargetSizeX/mm << ", " << TargetSizeY/mm << ", " << TargetThickness/mm << ";" << G4endl;

    // Note: No explicit check for magnet overlaps etc. here.
    //       This must be taken care of by the user!
    int i = 1;
    for (auto mds : magnetDefinitions) {
        G4String magnetName = "magnet_" + std::to_string(i++);
        magnets.push_back( MagnetBase::MagnetFactory(mds, this, magnetName) );
    }

    // materials
    DefineMaterials(); // The standard ones
    if (TargetThickness == 0.0) {
        // No target, only magnets!
        TargetMaterial = NULL;
    }
    else if (TargetMaterial_in.contains("::")) {
        // Gas target
        TargetMaterial = DefineGas(TargetMaterial_in);
    }
    else {
        // Solid target
        SetTargetMaterial(TargetMaterial_in);
    }

    SetBackgroundMaterial(BackgroundMaterial_in);

    DetectorMaterial = BackgroundMaterial;

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
    logicWorld = new G4LogicalVolume(solidWorld, BackgroundMaterial, "WorldLV");

    physiWorld = new G4PVPlacement(0,               //no rotation
                                   G4ThreeVector(), //World volume must be centered at the origin
                                   logicWorld,      //its logical volume
                                   "World",         //its name
                                   0,               //its mother  volume
                                   false,           //pMany not used
                                   0,               //copy number
                                   true);           //Check for overlaps

    //logicWorld->SetVisAttributes(G4VisAttributes(false));
    logicWorld->SetVisAttributes(G4Colour(0.5,0.5,0.5,0.1));

    //constructing the target
    if (TargetThickness > 0.0) {
        solidTarget = new G4Box("TargetS", TargetSizeX/2,TargetSizeY/2, TargetThickness/2);
        logicTarget = new G4LogicalVolume(solidTarget, TargetMaterial,"TargetLV");

        G4RotationMatrix* targetRot = new G4RotationMatrix();
        targetRot->rotateY(TargetAngle);
        G4ThreeVector zTrans(0.0,0.0,0.0);
        physiTarget = new G4PVPlacement(G4Transform3D(*targetRot,zTrans), //Translate (by zero)
                                                                          // then rotate
                                          logicTarget,                    //its logical volume
                                          "TargetPV",                     //its name
                                          logicWorld,                     //its mother
                                          false,                          //pMany not used
                                          0,                              //copy number
                                          true);                          //Check for overlaps
    }
    else {
        solidTarget = NULL;
        logicTarget = NULL;
        physiTarget = NULL;
    }

    // Get pointer to detector manager
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    if (logicTarget != NULL) {
        G4VSensitiveDetector* targetSD = new TargetSD("target");
        SDman->AddNewDetector(targetSD);
        logicTarget->SetSensitiveDetector(targetSD);
    }

    // Build magnets
    for (auto magnet : magnets) {
        // More or less repeated in ParallelWorldConstruction::Construct()
        magnet->Construct();

        G4VPhysicalVolume* magnetPV   = new G4PVPlacement(magnet->GetMainPV_transform(),
                                                          magnet->GetMainLV(),
                                                          magnet->magnetName + "_mainPV",
                                                          logicWorld,
                                                          false,
                                                          0,
                                                          true);

        magnetPVs.push_back(magnetPV);
    }

    return physiWorld;
}

//------------------------------------------------------------------------------

void DetectorConstruction::DefineMaterials() {
    // List of available materials:
    // http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/Appendix/materialNames.html
    G4NistManager* man = G4NistManager::Instance();
    man->SetVerbose(1);

    AlMaterial = man->FindOrBuildMaterial("G4_Al");
    CuMaterial = man->FindOrBuildMaterial("G4_C");
    CuMaterial = man->FindOrBuildMaterial("G4_Cu");
    PbMaterial = man->FindOrBuildMaterial("G4_Pb");
    TiMaterial = man->FindOrBuildMaterial("G4_Ti");
    SiMaterial = man->FindOrBuildMaterial("G4_Si");
    WMaterial  = man->FindOrBuildMaterial("G4_W");
    UMaterial  = man->FindOrBuildMaterial("G4_U");
    FeMaterial = man->FindOrBuildMaterial("G4_Fe");

    MylarMaterial          = man->FindOrBuildMaterial("G4_MYLAR");
    KaptonMaterial         = man->FindOrBuildMaterial("G4_KAPTON");
    StainlessSteelMaterial = man->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    WaterMaterial          = man->FindOrBuildMaterial("G4_WATER");
    NaIMaterial            = man->FindOrBuildMaterial("G4_SODIUM_IODIDE");

    vacuumMaterial = man->FindOrBuildMaterial("G4_Galactic");
    airMaterial    = man->FindOrBuildMaterial("G4_AIR");

    G4Element* elAl = new G4Element("Aluminium", "Al", 13.0, 26.9815385*g/mole);
    G4Element* elO  = new G4Element("Oxygen",    "O",  8.0,  15.999*g/mole);
    SapphireMaterial = new G4Material("Sapphire", 4.0*g/cm3, 2);
    SapphireMaterial->AddElement(elAl, 2);
    SapphireMaterial->AddElement(elO,  3);

    G4Element* elCr = new G4Element("Chromium", "Cr", 24.0, 51.9961*g/mole);
    ChromoxMaterial = new G4Material("ChromoxPure",5.22*g/cm3,2);
    ChromoxMaterial->AddElement(elCr,2);
    ChromoxMaterial->AddElement(elO,3);

    // Chromox parameters from https://www.advatech-uk.co.uk/chromox.html
    ChromoxScreenMaterial = new G4Material("ChromoxScreen", 3.85*g/cm3,2);
    ChromoxScreenMaterial->AddMaterial(SapphireMaterial,0.995);
    ChromoxScreenMaterial->AddMaterial(ChromoxMaterial,0.005);
}

//------------------------------------------------------------------------------

G4Material* DetectorConstruction::DefineGas(G4String gasMaterialName) {
    G4cout << G4endl;

    if (not gasMaterialName.contains("::")) {
        G4cerr << "Error in DetectorConstruction::DefineGas() -- no '::' was found"
               << "in material name = '" << gasMaterialName << "'"
               << G4endl;
        exit(1);
    }

    str_size colonPos = gasMaterialName.index("::");
    str_size pressurePos = colonPos+2;

    if (pressurePos >= gasMaterialName.length()) {
        G4cerr << "Error in DetectorConstruction::DefineGas() -- no pressure was found"
               << " after '::' in material name = '" << gasMaterialName << "'"
               << G4endl;
        exit(1);
    }
    G4String material_in = gasMaterialName(0, colonPos);
    G4String pressure_in = gasMaterialName(pressurePos, gasMaterialName.length()); // Bug, 2nd argument is length not position

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

    //Define materials (Hydrogen-1 / H_2)
    G4double aH         = 1.008*g/mole;
    G4double densityH_2 = 2*aH*(pressure*bar*1e-3)/(temperature*Avogadro*CLHEP::k_Boltzmann);
    G4Isotope* isH1    = new G4Isotope("H1",    //Name
                                       1,       //iz
                                       1,       //n
                                       aH);     //a
    G4Element* elH  = new G4Element("Hydrogen",   //name
                                    "H",          //symbol
                                    1);           //ncomponents
    elH->AddIsotope(isH1, 1.0);
    this->gasH_2 = new G4Material("HydrogenGas",      //name
                                densityH_2,         //density
                                1,                  //ncomponents
                                kStateGas,          //state
                                temperature,        //temp
                                pressure*bar*1e-3 );//pressure
    gasH_2->AddElement(elH, 2);
    G4cout << "Built H_2 gas, pressure = "<< pressure
           << " [mbar], temperature = " << temperature/kelvin
           << " [K], density = " << densityH_2 / g * meter3 << " [g/m3]"
           << G4endl;

    //Define materials (Helium-4)
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
                                 temperature,        //temp
                                 pressure*bar*1e-3 );//pressure
    gasHe->AddElement(elHe, 1);
    G4cout << "Built He gas, pressure = "<< pressure
           << " [mbar], temperature = " << temperature/kelvin
           << " [K], density = " << densityHe / g * meter3 << " [g/m3]"
           << G4endl;

    //Define materials (Nitrogen-14 / N_2)
    G4double aN        = 14.007*g/mole;
    G4double densityN_2 = 2*aN*(pressure*bar*1e-3)/(temperature*Avogadro*CLHEP::k_Boltzmann);
    G4Isotope* isN14   = new G4Isotope("N14",    //Name
                                         7,       //iz
                                         14,      //n
                                         aN);     //a
    G4Element* elN  = new G4Element("Nitrogen",   //name
                                    "N",          //symbol
                                    1);           //ncomponents
    elN->AddIsotope(isN14, 1.0);
    this->gasN_2 = new G4Material("NitrogenGas",      //name
                                  densityN_2,         //density
                                  1,                  //ncomponents
                                  kStateGas,          //state
                                  temperature,        //temp
                                  pressure*bar*1e-3 );//pressure
    gasN_2->AddElement(elN, 2);
    G4cout << "Built N_2 gas, pressure = "<< pressure
           << " [mbar], temperature = " << temperature/kelvin
           << " [K], density = " << densityN_2 / g * meter3 << " [g/m3]"
           << G4endl;

    //Define materials (Neon-20)
    G4double aNe       = 19.9924401754*g/mole;
    G4double densityNe = aNe*(pressure*bar*1e-3)/(temperature*Avogadro*CLHEP::k_Boltzmann);
    G4Isotope* isNe20 = new G4Isotope("Ne20", //Name
                                     10,      //iz
                                     20,      //n
                                     aNe);    //a
    G4Element* elNe  = new G4Element("Neon", //name
                                     "Ne",    //symbol
                                     1);      //ncomponents
    elNe->AddIsotope(isNe20, 1.0);
    this->gasNe = new G4Material("NeonGas",          //name
                                 densityNe,          //density
                                 1,                  //ncomponents
                                 kStateGas,          //state
                                 temperature,        //temp
                                 pressure*bar*1e-3 );//pressure
    gasNe->AddElement(elNe, 1);
    G4cout << "Built Ne gas, pressure = "<< pressure
           << " [mbar], temperature = " << temperature/kelvin
           << " [K], density = " << densityNe / g * meter3 << " [g/m3]"
           << G4endl;

    //Define materials (Argon-40)
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
                                 temperature,        //temp
                                 pressure*bar*1e-3 );//pressure
    gasAr->AddElement(elAr, 1);
    G4cout << "Built Ar gas, pressure = "<< pressure
           << " [mbar], temperature = " << temperature/kelvin
           << " [K], density = " << densityAr / g * meter3 << " [g/m3]"
           << G4endl;

    G4Material* returnMaterial = NULL;

    if (material_in == "H_2") {
        returnMaterial = this->gasH_2;
    }
    else if (material_in == "He") {
        returnMaterial = this->gasHe;
    }
    else if (material_in == "N_2") {
        returnMaterial = this->gasN_2;
    }
    else if (material_in == "Ne") {
        returnMaterial = this->gasNe;
    }
    else if (material_in == "Ar") {
        returnMaterial = this->gasAr;
    }
    else {
        G4cerr << "Error in DetectorConstruction::DefineGas()" << G4endl;
        G4cerr << "Gas type '" << material_in << "' unknown." << G4endl;
        exit(1);
    }

    G4cout << G4endl;

    return returnMaterial;
}

//------------------------------------------------------------------------------

void DetectorConstruction::SetTargetMaterial(G4String materialChoice) {
    // TODO: Merge with SetBackgroundMaterial()

    // search the material by its name
    if(!GetHasTarget() || TargetThickness == 0.0){
        G4cerr << "Error in DetectorConstruction::SetTargetMaterial():" << G4endl
               << " No target material is actually defined; probably target thickness is 0.0." << G4endl;
        exit(1);
    }

    G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
    if (pttoMaterial) {
        TargetMaterial = pttoMaterial;
    }
    else {
        G4cerr << "Error when setting material '"
               << materialChoice << "' -- not found!" << G4endl;
        G4MaterialTable* materialTable = G4Material::GetMaterialTable();
        G4cerr << "Valid choices:" << G4endl;
        for (auto mat : *materialTable) {
            G4cerr << mat->GetName() << G4endl;
        }
        exit(1);
    }
}

//------------------------------------------------------------------------------

void DetectorConstruction::SetBackgroundMaterial(G4String materialChoice) {
    // TODO: Merge with SetTargetMaterial()

    // search the material by its name
    G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
    if (pttoMaterial) {
        BackgroundMaterial = pttoMaterial;
    }
    else {
        G4cerr << "Error when setting material '"
               << materialChoice << "' -- not found!" << G4endl;
        G4MaterialTable* materialTable = G4Material::GetMaterialTable();
        G4cerr << "Valid choices:" << G4endl;
        for (auto mat : *materialTable) {
            G4cerr << mat->GetName() << G4endl;
        }
        exit(1);
    }
}
//------------------------------------------------------------------------------

G4int DetectorConstruction::GetTargetMaterialZ() {
    //Return the nuclear charge of the most common species in the target

    if (!GetHasTarget() || TargetMaterial == NULL){
        G4cerr << "Error in DetectorConstruction::GetTargetMaterialZ():" << G4endl
               << " No target material is actually defined; probably target thickness is 0.0." << G4endl;
        exit(1);
    }

    const size_t numElements             = TargetMaterial->GetNumberOfElements();
    const G4ElementVector* elementVector = TargetMaterial->GetElementVector();
    const G4int* atomsVector             = TargetMaterial->GetAtomsVector();

    G4int maxAtoms = 0;
    size_t maxAtomsIndex = -1;

    for (size_t i = 0; i<numElements; i++) {
        if (atomsVector[i] > maxAtoms) {
            maxAtoms = atomsVector[i];
            maxAtomsIndex = i;
        }
    }

    return (*elementVector)[maxAtomsIndex]->GetZ();
}

//------------------------------------------------------------------------------

G4double DetectorConstruction::GetTargetMaterialA() {
    //Return the average mass number of the most common species in the target

    if (!GetHasTarget() || TargetMaterial == NULL){
        G4cerr << "Error in DetectorConstruction::GetTargetMAterialZ():" << G4endl
               << " No target material is actually defined; probably target thickness is 0.0." << G4endl;
        exit(1);
    }

    const size_t numElements             = TargetMaterial->GetNumberOfElements();
    const G4ElementVector* elementVector = TargetMaterial->GetElementVector();
    const G4int* atomsVector             = TargetMaterial->GetAtomsVector();

    G4int maxAtoms = 0;
    size_t maxAtomsIndex = -1;

    for (size_t i = 0; i<numElements; i++) {
        if (atomsVector[i] > maxAtoms) {
            maxAtoms = atomsVector[i];
            maxAtomsIndex = i;
        }
    }

    return (*elementVector)[maxAtomsIndex]->GetAtomicMassAmu();

}

//------------------------------------------------------------------------------

G4double DetectorConstruction::GetTargetMaterialDensity() {
    if (!GetHasTarget() || TargetMaterial == NULL){
        G4cerr << "Error in DetectorConstruction::GetTargetMAterialZ():" << G4endl
               << " No target material is actually defined; probably target thickness is 0.0." << G4endl;
        exit(1);
    }
    return TargetMaterial->GetDensity();
}

//------------------------------------------------------------------------------

G4Material* DetectorConstruction::GetTargetMaterial() {
    if (!GetHasTarget() || TargetMaterial == NULL) {
        G4cerr << "Error in DetectorConstruction::GetTargetMaterial():" << G4endl
               << " No target material is actually defined; probably target thickness is 0.0." << G4endl;
        exit(1);
    }
    return TargetMaterial;
}

//------------------------------------------------------------------------------

G4Material* DetectorConstruction::GetBackgroundMaterial() {
    if (BackgroundMaterial == NULL) {
        G4cerr << "Error in DetectorConstruction::GetBackgroundMaterial():" << G4endl
               << " No background material is actually defined?" << G4endl;
        exit(1);
    }
    return BackgroundMaterial;
}

//------------------------------------------------------------------------------

void DetectorConstruction::PostInitialize() {
    // Setup the magnet fields.
    for (auto mag : magnets) {
        mag->PostInitialize();
    }
}
