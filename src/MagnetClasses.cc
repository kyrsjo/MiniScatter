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

#include "MagnetClasses.hh"

#include "G4String.hh"
#include <string>

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"

#include "G4PVPlacement.hh"

#include "MyTargetSD.hh"
#include "G4SDManager.hh"

#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4Mag_UsualEqRhs.hh"

#include "G4PhysicalConstants.hh"

MagnetBase* MagnetBase:: MagnetFactory(G4String inputString, DetectorConstruction* detCon, G4String magnetName) {

    //Split by '::'
    std::vector<G4String> argList;
    str_size startPos = 0;
    str_size endPos = 0;
    do {
        endPos   = inputString.index(":",startPos);
        argList.push_back(inputString(startPos,endPos-startPos));
        startPos = endPos+1;
    } while (endPos != std::string::npos);

    //Debug
    for (auto arg : argList) {
        G4cout << "'" << arg << "'" << G4endl;
    }

    if (argList.size() < 4) {
        G4cerr << "Error when parsing magnet input string '"<<inputString<<"'" << G4endl;
        G4cerr << "Expected at least 4 arguments, as in 'pos:type:length:gradient'." << G4endl;
        exit(1);
    }

    //Type-specific key=val pairs
    std::map<G4String,G4String> keyValPairs;
    for (size_t i = 4; i < argList.size(); i++) {
        str_size eqPos = argList[i].index("=",0);
        if (eqPos == std::string::npos) {
            G4cerr << "Error when parsing key=val pair '" << argList[i] << "', no '=' found!" << G4endl;
            exit(1);
        }
        keyValPairs[argList[i](0,eqPos)] = argList[i](eqPos+1,std::string::npos);
    }

    //Debug
    for (auto arg : keyValPairs) {
        G4cout << "'" << arg.first << "' = '" << arg.second << "'" << G4endl;
    }

    //Parse the common part
    G4bool doRelPos = false;
    G4double magnetPos = false;
    if (strlen(argList[0]) > 0 && argList[0][0] == '*') {
        doRelPos = true;
    }
    try {
        if (doRelPos) {
            magnetPos = std::stod(std::string(argList[0]).substr(1));
        }
        else {
            magnetPos = std::stod(std::string(argList[0]));
        }
        magnetPos *= mm;
    }
    catch (const std::invalid_argument& ia) {
        G4cout << "Invalid argument when reading magnet position" << G4endl
               << "Got: '" << argList[0] << "'" << G4endl
               << "Expected a floating point number! (exponential notation and prepended '*' is accepted)" << G4endl;
        exit(1);
    }

    G4String magnetType = argList[1];

    G4double magnetLength = 0.0;
    try {
        magnetLength = std::stod(std::string(argList[2]));
        magnetLength *= mm;
    }
    catch (const std::invalid_argument& ia) {
        G4cerr << "Invalid argument when reading magnet length" << G4endl
               << "Got: '" << argList[2] << "'" << G4endl
               << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
        exit(1);
    }

    G4double magnetGradient = 0.0;
    try {
        magnetGradient = std::stod(std::string(argList[3]));
    }
    catch (const std::invalid_argument& ia) {
        G4cerr << "Invalid argument when reading magnet gradient" << G4endl
               << "Got: '" << argList[3] << "'" << G4endl
               << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
        exit(1);
    }

    // Build new Magnets
    MagnetBase* theMagnet = NULL;
    if (magnetType == "PLASMA1") {
        theMagnet = new MagnetPLASMA1(magnetPos, doRelPos, magnetLength, magnetGradient,
                                      keyValPairs, detCon, magnetName);
    }
    else if(magnetType == "COLLIMATOR1") {
        theMagnet = new MagnetCOLLIMATOR1(magnetPos, doRelPos, magnetLength, magnetGradient,
                                          keyValPairs, detCon, magnetName);
    }
    else {
        G4cerr << "Uknown magnet type '" << magnetType << "'" << G4endl;
        exit(1);
    }

    //exit(0);
    return theMagnet;
}

G4LogicalVolume* MagnetBase::MakeNewMainLV(G4String name_postfix){
    // Builds an outer volume
    // Note: The mainLV's physical volume(s) are created in the DetectorConstruction classes

    G4Material* vacuumMaterial = G4Material::GetMaterial("G4_Galactic");
    if (not vacuumMaterial) {
        G4cerr << "Internal error -- material G4_Galactic not found in MagnetBase::MakeMainLV()!" << G4endl;
        exit(1);
    }

    mainLV_w = detCon->getWorldSizeX()-2*xOffset-2*(length/2.0)*sin(xRot/rad);
    mainLV_h = detCon->getWorldSizeX()-2*xOffset-2*(length/2.0)*sin(yRot/rad);

    if (mainLV_w < 0.0 or mainLV_h < 0.0) {
        G4cerr << "Error in MagnetBase::MakeNewMainLV():" << G4endl
               << " It is not possible to fit the rotated and translated" << G4endl
               << " outer volume for magnet '" << magnetName << "' inside the world volume." << G4endl;
        exit(1);
    }

    G4VSolid* mainBox = new G4Box(magnetName+"_"+name_postfix+"S",
                                  mainLV_w/2.0,
                                  mainLV_h/2.0,
                                  length/2.0);

    G4LogicalVolume* newMainLV = new G4LogicalVolume(mainBox,vacuumMaterial,
                                                     magnetName+"_"+name_postfix+"LV");

    return newMainLV;
}

G4LogicalVolume* MagnetBase::GetMainLV() const {
    //Returns the mainLV.

    if (this->mainLV != NULL) {
        return this->mainLV;
    }
    else {
        G4cout << "Error in MagnetClasses::GetMainLV(): MainLV has not been constructed!" << G4endl;
        exit(1);
    }
}

void MagnetBase::BuildMainPV_transform() {
    G4RotationMatrix mainPV_rotate;
    mainPV_rotate.rotateX(xRot);
    mainPV_rotate.rotateY(yRot);
    G4ThreeVector    mainPV_translate;
    mainPV_translate.setX(xOffset);
    mainPV_translate.setY(yOffset);
    mainPV_transform = G4Transform3D(mainPV_rotate, mainPV_translate);
}

void MagnetBase::ConstructDetectorLV() {
    if (this->detectorLV != NULL) {
        G4cerr << "Error in MagnetBase::ConstructDetectorLV(): The detectorLV has already been constructed?" << G4endl;
        exit(1);
    }

    this->detectorLV = MakeNewMainLV("detector");
}
void MagnetBase::AddSD(){
    // Add the TargetSD to the virtual logical volume of the magnet.
    // This records the outgoing position and energy deposit in the magnet.

    // Get pointer to detector manager
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    magnetSD = new MyTargetSD(magnetName);
    SDman->AddNewDetector(magnetSD);
    this->detectorLV->SetSensitiveDetector(magnetSD);
}
G4LogicalVolume* MagnetBase::GetDetectorLV() const {
    //Returns the detectorLV.

    if (this->detectorLV != NULL) {
        return this->detectorLV;
    }
    else {
        G4cout << "Error in MagnetClasses::GetDetectorLV(): DetectorLV has not been constructed!" << G4endl;
        exit(1);
    }
}

MagnetPLASMA1::MagnetPLASMA1(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                             std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                             G4String magnetName_in) :
    MagnetBase(zPos_in, doRelPos_in, length_in, gradient_in, keyValPairs_in, detCon_in, magnetName_in) {

    // Default values
    G4bool inputIsTotalAmps = false;
    capRadius= 1.0*mm;
    cryWidth = 10.0*mm;
    cryHeight = 20.0*mm;

    for (auto it : keyValPairs) {
        if (it.first == "radius") {
            try {
                capRadius = std::stod(std::string(it.second)) * mm;
            }
            catch (const std::invalid_argument& ia) {
                G4cerr << "Invalid argument when reading capillary radius" << G4endl
                       << "Got: '" << it.second << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
        }
        else if (it.first == "totalAmps") {
            if (it.second == "True") {
                inputIsTotalAmps = true;
            }
            else if (it.second == "False") {
                inputIsTotalAmps = false;
            }
            else {
                G4cerr << "Invalid argument when reading capillary totalAmps flag" << G4endl
                       << "Got: " << it.second << "'" << G4endl
                       << "Expected 'True' or 'False'" << G4endl;
                exit(1);
            }
        }
        else if (it.first == "width") {
            try {
                cryWidth = std::stod(std::string(it.second)) * mm;
            }
            catch (const std::invalid_argument& ia) {
                G4cerr << "Invalid argument when reading crystal width" << G4endl
                       << "Got: '" << it.second << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
        }
        else if (it.first == "height") {
            try {
                cryHeight = std::stod(std::string(it.second)) * mm;
            }
            catch (const std::invalid_argument& ia) {
                G4cerr << "Invalid argument when reading crystal height" << G4endl
                       << "Got: '" << it.second << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
        }
        else if (it.first == "xOffset") {
            try {
                xOffset = std::stod(std::string(it.second)) * mm;
            }
            catch (const std::invalid_argument& ia) {
                G4cerr << "Invalid argument when reading xOffset" << G4endl
                       << "Got: '" << it.second << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
        }
        else if (it.first == "yOffset") {
            try {
                yOffset = std::stod(std::string(it.second)) * mm;
            }
            catch (const std::invalid_argument& ia) {
                G4cerr << "Invalid argument when reading xOffset" << G4endl
                       << "Got: '" << it.second << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
        }
        else if (it.first == "xRot") {
            try {
                xRot = std::stod(std::string(it.second)) * deg;
            }
            catch (const std::invalid_argument& ia) {
                G4cerr << "Invalid argument when reading xRot" << G4endl
                       << "Got: '" << it.second << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
        }
        else if (it.first == "yRot") {
            try {
                yRot = std::stod(std::string(it.second)) * deg;
            }
            catch (const std::invalid_argument& ia) {
                G4cerr << "Invalid argument when reading yRot" << G4endl
                       << "Got: '" << it.second << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
        }
        else {
            G4cerr << "MagnetPLASMA1 did not understand key=value pair '"
                   << it.first << "'='" << it.second << "'." << G4endl;
            exit(1);
        }
    }

    //if the current is given, compute and set the gradient; otherwise compute and set the current
    if (inputIsTotalAmps) {
        plasmaTotalCurrent = gradient; //[A]
        // Gradient is always stored as [T/m]
        gradient = ( mu0*(plasmaTotalCurrent*ampere) / (twopi*capRadius*capRadius) ) / (tesla/meter);
    }
    else {
        // [A]
        plasmaTotalCurrent = ( (gradient*tesla/meter) * twopi*capRadius*capRadius / mu0 ) / ampere;
    }

    G4cout << "Initialized a MagnetPLASMA1, parameters:" << G4endl;
    G4cout << "\t magnetName         = " << magnetName << G4endl;
    G4cout << "\t Z0                 = " << getZ0()/mm         << " [mm]"  << G4endl;
    G4cout << "\t length             = " << length/mm          << " [mm]"  << G4endl;
    G4cout << "\t gradient           = " << gradient           << " [T/m]" << G4endl;
    G4cout << "\t plasmaTotalcurrent = " << plasmaTotalCurrent << " [A]"   << G4endl;
    G4cout << "\t capRadius          = " << capRadius/mm       << " [mm]"  << G4endl;
    G4cout << "\t cryWidth           = " << cryWidth/mm        << " [mm]"  << G4endl;
    G4cout << "\t cryHeight          = " << cryHeight/mm       << " [mm]"  << G4endl;
    G4cout << "\t xOffset            = " << xOffset/mm         << " [mm]"  << G4endl;
    G4cout << "\t yOffset            = " << yOffset/mm         << " [mm]"  << G4endl;
    G4cout << "\t xRot               = " << xRot/deg           << " [deg]" << G4endl;
    G4cout << "\t yRot               = " << yRot/deg           << " [deg]" << G4endl;

}

void MagnetPLASMA1::Construct() {

    if (this->mainLV != NULL) {
        G4cerr << "Error in MagnetPLASMA1::Construct(): The mainLV has already been constructed?" << G4endl;
        exit(1);
    }

    this->mainLV = MakeNewMainLV("main");

    /*
    G4Material* vacuumMaterial = G4Material::GetMaterial("G4_Galactic");
    if (not vacuumMaterial) {
        G4cerr << "Internal error -- material G4_Galactic not found in MagnetPLASMA1::Construct()!" << G4endl;
        exit(1);
    }

    //Field box
    // TODO: Not really neccessary, use the mainLV directly
    // (unless the goal is to make it smaller and contain the gas...)
    G4VSolid* fieldBox            = new G4Box(magnetName+"_fieldBoxS",
                                              detCon->getWorldSizeX()/2.0,
                                              detCon->getWorldSizeX()/2.0,
                                              length/2.0);
    G4Material* fieldBoxMaterial  = vacuumMaterial; // TODO
    G4LogicalVolume* fieldBoxLV   = new G4LogicalVolume(fieldBox,
                                                        fieldBoxMaterial,
                                                        magnetName+"_fieldBoxLV");
    G4VPhysicalVolume* fieldBoxPV = new G4PVPlacement(NULL,
                                                      G4ThreeVector(0.0,0.0,0.0),
                                                      fieldBoxLV,
                                                      magnetName + "_fieldBoxPV",
                                                      this->mainLV,
                                                      false,
                                                      0,
                                                      true); */
    field = new FieldPLASMA1(plasmaTotalCurrent, capRadius,
                             G4ThreeVector(xOffset, yOffset, getZ0()),mainLV);
    G4FieldManager* fieldMgr = new G4FieldManager(field);
    G4Mag_UsualEqRhs* fieldEquation = new G4Mag_UsualEqRhs(field);
    G4MagIntegratorStepper* fieldStepper = new G4ClassicalRK4(fieldEquation);
    G4ChordFinder* fieldChordFinder = new G4ChordFinder(field, capRadius/4.0, fieldStepper);
    fieldMgr->SetChordFinder(fieldChordFinder);
    //fieldBoxLV->SetFieldManager(fieldMgr,true);
    mainLV->SetFieldManager(fieldMgr,true);

    if (cryWidth > mainLV_w || cryHeight > mainLV_h) {
        G4cerr << "Error in MagnetPLASMA1::Construct():" << G4endl
               << " The crystal is wider than it's allowed envelope "
               << " including offsets and rotations."  << G4endl;
        exit(1);
    }

    if (capRadius > cryWidth/2.0 or capRadius > cryHeight/2.0) {
        G4cerr << "Error in MagnetPLASMA1::Construct():" << G4endl
               << " The capillary doesn't fit in the crystal!" << G4endl;
        exit(1);
    }

    //TODO: Insert here a "gas box" that is exactly the same size as the crystal
    // and made of gas material, but has no hole, OR is exactly the size of the hole.

    // The crystal
    G4VSolid* crystalBox      = new G4Box(magnetName+"_crystalBoxS",
                                          cryWidth/2.0, cryHeight/2.0, length/2.0);
    G4VSolid* crystalCylinder = new G4Tubs(magnetName+"_crystalCylinderS",
                                           0.0, capRadius, length,
                                           0.0, 360.0*deg);
    G4VSolid* crystalSolid    = new G4SubtractionSolid(magnetName+"_crystalS",
                                                       crystalBox, crystalCylinder);

    G4Material* sapphireMaterial = G4Material::GetMaterial("Sapphire");
    if (not sapphireMaterial) {
        G4cerr << "Internal error -- material Sapphire not found in MagnetPLASMA1::Construct()!" << G4endl;
        exit(1);
    }
    G4LogicalVolume*   crystalLV = new G4LogicalVolume(crystalSolid,sapphireMaterial, magnetName+"_crystalLV");
    G4VPhysicalVolume* crystalPV = new G4PVPlacement(NULL,
                                                     G4ThreeVector(0.0,0.0,0.0),
                                                     crystalLV,
                                                     magnetName + "_crystalPV",
                                                     mainLV,
                                                     false,
                                                     0,
                                                     true);

    ConstructDetectorLV();
    BuildMainPV_transform();
}


G4Navigator* FieldBase::fNavigator = NULL;

void FieldBase::SetupTransform() {
    // Initialization of global->local transform based on the Geant4 example  "extended/field/field04"
    // Ran the first time the GetFieldValue is called, since it needs to know the whole volume tree.

    //If neccessary, create our own navigator
    G4Navigator* theNavigator =
        G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    if (!fNavigator) {
        // fNaviator is a shared static object. If it does not exist, create it.
        FieldBase::fNavigator = new G4Navigator();
        if( theNavigator->GetWorldVolume() )
            fNavigator->SetWorldVolume(theNavigator->GetWorldVolume());
    }

    //·set·fGlobalToLocal·transform
    fNavigator->LocateGlobalPointAndSetup(centerPoint,0,false);
    G4TouchableHistoryHandle touchable = fNavigator->CreateTouchableHistoryHandle();
    G4int depth = touchable->GetHistoryDepth();
    G4bool foundVolume = false;
    for (G4int i = 0; i<depth; ++i) {
        if(touchable->GetVolume()->GetLogicalVolume() == fieldLV) {
            foundVolume = true;
            break;
        }
        touchable->MoveUpHistory();
    }
    if (not foundVolume) {
        G4cerr << "Internal error in FieldBase::InitializeTransform(): Could not find the volume!" << G4endl;
        exit(1);
    }
    fGlobalToLocal = touchable->GetHistory()->GetTopTransform();
}

FieldPLASMA1::FieldPLASMA1(G4double current_in, G4double radius_in,
                           G4ThreeVector centerPoint_in, G4LogicalVolume* fieldLV_in) :
    FieldBase(centerPoint_in, fieldLV_in), plasmaTotalCurrent(current_in), capRadius(radius_in) {
    //Gradient is always stored as [T/m]. Current stored as [A]. Distances stored in G4 units.
    gradient = ( mu0*(plasmaTotalCurrent*ampere) / (twopi*capRadius*capRadius) ) / (tesla/meter);
}

void FieldPLASMA1::GetFieldValue(const G4double point[4], G4double field[6]) const {

    G4ThreeVector global(point[0],point[1],point[2]);
    G4ThreeVector local = fGlobalToLocal.TransformPoint(global);

    //G4cout << "GetFieldValue, global = " << global << " local = " << local << G4endl;

    G4ThreeVector B(0.0,0.0,0.0);

    G4double r = sqrt(local[0]*local[0] + local[1]*local[1]); // Position in [G4 units],
                                                              // centered on the magnet center

    G4double Btheta;
    if (r < capRadius) { //Inside the capillary
        Btheta = gradient*(r/meter) * tesla; //Field strength in [G4 units]; convert r to [m]
    }
    else {               //Outside the capillary
        // Current stored as A
        Btheta = mu0*(plasmaTotalCurrent*ampere) / (twopi*r); //Field strength in [G4 units]
    }

    G4double theta = atan2(local[1],local[0]);
    B[0] = -sin(theta)*Btheta;
    B[1] =  cos(theta)*Btheta;
    B[2] = 0.0;

    B = fGlobalToLocal.Inverse().TransformAxis(B);
    field[0] = B[0];
    field[1] = B[1];
    field[2] = B[2];
}

MagnetCOLLIMATOR1::MagnetCOLLIMATOR1(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                                     std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                                     G4String magnetName_in) :
    MagnetBase(zPos_in, doRelPos_in, length_in, gradient_in, keyValPairs_in, detCon_in, magnetName_in){

    for (auto it : keyValPairs) {
        if (it.first == "radius") {
            try {
                radius = std::stod(std::string(it.second)) * mm;
            }
            catch (const std::invalid_argument& ia) {
                G4cerr << "Invalid argument when reading collimator radius" << G4endl
                       << "Got: '" << it.second << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
        }
        else if (it.first == "width") {
            try {
                width = std::stod(std::string(it.second)) * mm;
            }
            catch (const std::invalid_argument& ia) {
                G4cerr << "Invalid argument when reading absorber width" << G4endl
                       << "Got: '" << it.second << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
        }
        else if (it.first == "height") {
            try {
                height = std::stod(std::string(it.second)) * mm;
            }
            catch (const std::invalid_argument& ia) {
                G4cerr << "Invalid argument when reading absorber height" << G4endl
                       << "Got: '" << it.second << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
        }
        else if (it.first == "material") {
            absorberMaterialName = it.second;
        }
        else if (it.first == "xOffset") {
            try {
                xOffset = std::stod(std::string(it.second)) * mm;
            }
            catch (const std::invalid_argument& ia) {
                G4cerr << "Invalid argument when reading xOffset" << G4endl
                       << "Got: '" << it.second << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
        }
        else if (it.first == "yOffset") {
            try {
                yOffset = std::stod(std::string(it.second)) * mm;
            }
            catch (const std::invalid_argument& ia) {
                G4cerr << "Invalid argument when reading yOffset" << G4endl
                       << "Got: '" << it.second << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
        }
        else if (it.first == "xRot") {
            try {
                xRot = std::stod(std::string(it.second)) * deg;
            }
            catch (const std::invalid_argument& ia) {
                G4cerr << "Invalid argument when reading xRot" << G4endl
                       << "Got: '" << it.second << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
        }
        else if (it.first == "yRot") {
            try {
                yRot = std::stod(std::string(it.second)) * deg;
            }
            catch (const std::invalid_argument& ia) {
                G4cerr << "Invalid argument when reading yRot" << G4endl
                       << "Got: '" << it.second << "'" << G4endl
                       << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
                exit(1);
            }
        }
        else {
            G4cerr << "MagnetCOLLIMATOR1 did not understand key=value pair '"
                   << it.first << "'='" << it.second << "'." << G4endl;
            exit(1);
        }
    }

    if (gradient != 0.0) {
        G4cerr << "Invalid gradient for COLLIMATOR1: Gradient must be 0.0, but was "
               << gradient << "[T/m]" << G4endl;
        exit(1);
    }

    G4cout << "Initialized a MagnetCOLLIMATOR1, parameters:" << G4endl;
    G4cout << "\t magnetName         = " << magnetName << G4endl;
    G4cout << "\t Z0                 = " << getZ0()/mm         << " [mm]"  << G4endl;
    G4cout << "\t length             = " << length/mm          << " [mm]"  << G4endl;
    G4cout << "\t gradient           = " << gradient           << " [T/m]" << G4endl;
    G4cout << "\t radius             = " << radius/mm          << " [mm]"  << G4endl;
    G4cout << "\t width              = " << width/mm           << " [mm]"  << G4endl;
    G4cout << "\t height             = " << height/mm          << " [mm]"  << G4endl;
    G4cout << "\t xOffset            = " << xOffset/mm         << " [mm]"  << G4endl;
    G4cout << "\t yOffset            = " << yOffset/mm         << " [mm]"  << G4endl;
    G4cout << "\t xRot               = " << xRot/deg           << " [deg]" << G4endl;
    G4cout << "\t yRot               = " << yRot/deg           << " [deg]" << G4endl;
}

void MagnetCOLLIMATOR1::Construct() {
    if (this->mainLV != NULL) {
        G4cerr << "Error in MagnetCOLLIMATOR1::Construct(): The mainLV has already been constructed?" << G4endl;
        exit(1);
    }

    this->mainLV = MakeNewMainLV("main");

    //Sanity checks on dimensions
    if (width > mainLV_w || height > mainLV_h) {
        G4cerr << "Error in MagnetCOLLIMATOR1::Construct():" << G4endl
               << " The absorber is wider than it's allowed envelope "
               << " including offsets and rotations."  << G4endl;
        exit(1);
    }

    if (radius > width/2.0 or radius > height/2.0) {
        G4cerr << "Error in MagnetCOLLIMATOR1::Construct():" << G4endl
               << " The channel doesn't fit in the absorber!" << G4endl;
        exit(1);
    }

    // Build the absorber
    G4VSolid* absorberBox      = new G4Box(magnetName+"_absorberBoxS",
                                           width/2.0, height/2.0, length/2.0);
    G4VSolid* absorberCylinder = new G4Tubs(magnetName+"_absorberCylinderS",
                                            0.0, radius, length,
                                            0.0, 360.0*deg);
    G4VSolid* absorberSolid    = new G4SubtractionSolid(magnetName+"_absorberS",
                                                       absorberBox, absorberCylinder);

    absorberMaterial = G4Material::GetMaterial(absorberMaterialName);
    if (not absorberMaterial){
        G4cerr << "Error when setting material '"
               << absorberMaterialName << "' for MagnetCollimator '"
               << magnetName << "' -- not found!" << G4endl;
        G4MaterialTable* materialTable = G4Material::GetMaterialTable();
        G4cerr << "Valid choices:" << G4endl;
        for (auto mat : *materialTable) {
            G4cerr << mat->GetName() << G4endl;
        }
        exit(1);
    }

    G4LogicalVolume*   absorberLV = new G4LogicalVolume(absorberSolid,absorberMaterial, magnetName+"_absorberLV");
    G4VPhysicalVolume* absorberPV = new G4PVPlacement(NULL,
                                                      G4ThreeVector(0.0,0.0,0.0),
                                                      absorberLV,
                                                      magnetName + "_absorberPV",
                                                      mainLV,
                                                      false,
                                                      0,
                                                      true);

    ConstructDetectorLV();
    BuildMainPV_transform();
}
