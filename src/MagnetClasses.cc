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

#include "TargetSD.hh"
#include "G4SDManager.hh"

#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4Mag_UsualEqRhs.hh"

#include "G4PhysicalConstants.hh"

#include "G4VisAttributes.hh"

#include "G4Exception.hh"

MagnetBase* MagnetBase::MagnetFactory(G4String inputString, DetectorConstruction* detCon, G4String magnetName) {

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
    //for (auto arg : argList) {
    //    G4cout << "'" << arg << "'" << G4endl;
    //}

    if (argList.size() < 4) {
        G4cerr << "Error when parsing object/magnet input string '"<<inputString<<"'" << G4endl;
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
    //for (auto arg : keyValPairs) {
    //    G4cout << "'" << arg.first << "' = '" << arg.second << "'" << G4endl;
    //}

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
        theMagnet = new MagnetPLASMA1             (magnetPos, doRelPos, magnetLength, magnetGradient,
                                                   keyValPairs, detCon, magnetName);
    }
    else if(magnetType == "COLLIMATOR1") {
        theMagnet = new MagnetCOLLIMATOR1         (magnetPos, doRelPos, magnetLength, magnetGradient,
                                                   keyValPairs, detCon, magnetName);
    }
    else if(magnetType == "COLLIMATORRECT") {
        theMagnet = new MagnetCOLLIMATORRECT      (magnetPos, doRelPos, magnetLength, magnetGradient,
                                                   keyValPairs, detCon, magnetName);
    }
    else if(magnetType == "TARGET") {
        theMagnet = new MagnetTARGET              (magnetPos, doRelPos, magnetLength, magnetGradient,
                                                   keyValPairs, detCon, magnetName);
    }
    else if(magnetType == "TARGETR") {
        theMagnet = new MagnetTARGETR             (magnetPos, doRelPos, magnetLength, magnetGradient,
                                                   keyValPairs, detCon, magnetName);
    }
    else if (magnetType == "COLLIMATORHV") {
        theMagnet = new MagnetCOLLIMATORHV        (magnetPos, doRelPos, magnetLength, magnetGradient,
                                                   keyValPairs, detCon, magnetName);
    }
    else if (magnetType == "SHIELDEDSCINTILLATOR") {
        theMagnet = new MagnetSHIELDEDSCINTILLATOR(magnetPos, doRelPos, magnetLength, magnetGradient,
                                                   keyValPairs, detCon, magnetName);
    }
    else if(magnetType == "PBW") {
        theMagnet = new MagnetPBW                 (magnetPos, doRelPos, magnetLength, magnetGradient,
                                                   keyValPairs, detCon, magnetName);
    }
    else {
        G4cerr << "Unknown magnet type '" << magnetType << "'" << G4endl;
        exit(1);
    }

    //exit(0);
    return theMagnet;
}

G4LogicalVolume* MagnetBase::MakeNewMainLV(G4String name_postfix, G4double width, G4double height){
    // Builds an outer volume
    // Note: The mainLV's physical volume(s) are created in the DetectorConstruction class
    // Set width and height both <= 0.0 to autogenerate the wrapping volume width and height

    G4Material* BackgroundMaterial = detCon->GetBackgroundMaterial();

    if (width <= 0.0 and height <= 0) {
        //Autogenerate a maximalistic "slice",
        // which will later be rotated when it is placed in the DetectorConstruction
        mainLV_w = detCon->getWorldSizeX()-2*xOffset-2*(length/2.0)*sin(xRot/rad);
        mainLV_h = detCon->getWorldSizeX()-2*xOffset-2*(length/2.0)*sin(yRot/rad);
    }
    else {
        //Use input width/height
        mainLV_w = width;
        mainLV_h = height;
    }

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

    G4LogicalVolume* newMainLV = new G4LogicalVolume(mainBox,BackgroundMaterial,
                                                     magnetName+"_"+name_postfix+"LV");

    // Make the bounding volume invisible by default
    newMainLV->SetVisAttributes(G4VisAttributes(false));

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
    mainPV_translate.setZ(getZ0());
    mainPV_transform = G4Transform3D(mainPV_rotate, mainPV_translate);
}

void MagnetBase::ConstructDetectorLV() {
    if (this->detectorLV != NULL) {
        G4cerr << "Error in MagnetBase::ConstructDetectorLV(): The detectorLV has already been constructed?" << G4endl;
        exit(1);
    }
    if(this->mainLV == NULL) {
        G4cerr << "Error in MagnetBase::ConstructDetectorLV(): The mainLV is not yet constructed?" << G4endl;
        exit(1);
    }

    this->detectorLV = MakeNewMainLV("detector",mainLV_w,mainLV_h);
}
void MagnetBase::AddSD(){
    // Add the TargetSD to the virtual logical volume of the magnet.
    // This records the outgoing position and energy deposit in the magnet.

    // Get pointer to detector manager
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    magnetSD = new TargetSD(magnetName);
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

// Input parsing helpers
G4double MagnetBase::ParseDouble(G4String inStr, G4String readWhat) {
    try {
        return std::stod(std::string(inStr));
    }
    catch (const std::invalid_argument& ia) {
        G4cerr << "Invalid argument when reading " << readWhat << G4endl
               << "Got: '" << inStr << "'" << G4endl
               << "Expected a floating point number! (exponential notation is accepted)" << G4endl;
        exit(1);
    }
}

G4bool MagnetBase::ParseBool(G4String inStr, G4String readWhat) {
    if (inStr == "True") {
        return true;
    }
    else if (inStr == "False") {
        return false;
    }
    else {
        G4cerr << "Invalid argument when reading " << readWhat << " flag" << G4endl
               << "Got: " << inStr << "'" << G4endl
               << "Expected 'True' or 'False'" << G4endl;
        exit(1);
    }
}

void MagnetBase::ParseOffsetRot(G4String k, G4String v) {
    if (k == "xOffset") {
        xOffset = ParseDouble(v, "xOffset") * mm;
    }
    else if (k == "yOffset") {
        yOffset = ParseDouble(v, "yOffset") * mm;
    }
    else if (k == "xRot") {
        xRot = ParseDouble(v, "xRot") * deg;
    }
    else if (k == "yRot") {
        yRot = ParseDouble(v, "yRot") * deg;
    }
    else {
        G4cerr << "MagnetBase::ParseOffsetRot() cannot parse key '" << k << "' "
               << "(value = '" << v << "')" << G4endl;
        exit(1);
    }
}

void MagnetBase::PrintCommonParameters() {
    G4cout << " Initialized a '" << magnetType << "', parameters:"  <<             G4endl;
    G4cout << "\t magnetName              = " << magnetName         <<             G4endl;
    G4cout << "\t Z0                      = " << getZ0()/mm         << " [mm]"  << G4endl;
    G4cout << "\t length                  = " << length/mm          << " [mm]"  << G4endl;
    G4cout << "\t gradient                = " << gradient           << " [T/m]" << G4endl;
    G4cout << "\t xOffset                 = " << xOffset/mm         << " [mm]"  << G4endl;
    G4cout << "\t yOffset                 = " << yOffset/mm         << " [mm]"  << G4endl;
    G4cout << "\t xRot                    = " << xRot/deg           << " [deg]" << G4endl;
    G4cout << "\t yRot                    = " << yRot/deg           << " [deg]" << G4endl;
}

/** FIELD PATTERN BASE CLASS **/

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