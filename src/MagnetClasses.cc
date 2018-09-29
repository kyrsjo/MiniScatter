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
    else {
        G4cerr << "Uknown magnet type '" << magnetType << "'" << G4endl;
        exit(1);
    }

    //exit(0);
    return theMagnet;
}

void MagnetBase::AddSD(G4LogicalVolume* mainLV) {
    // Adds a TargetSD to the main logical volume of the magnet.
    // This records the outgoing position and energy deposit in the magnet.

    // Get pointer to detector manager
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    magnetSD = new MyTargetSD(magnetName);
    SDman->AddNewDetector(magnetSD);
    mainLV->SetSensitiveDetector(magnetSD);
}

MagnetPLASMA1::MagnetPLASMA1(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                             std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                             G4String magnetName_in) :
    MagnetBase(zPos_in, doRelPos_in, length_in, gradient_in, keyValPairs_in, detCon_in, magnetName_in) {

    G4bool inputIsTotalAmps = false;

    capRadius= 1.0*mm;

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
        else {
            G4cerr << "Capillary did not understand key=value pair '"
                   << it.first << "'='" << it.second << "'." << G4endl;
            exit(1);
        }
    }

    if (inputIsTotalAmps) {
        plasmaTotalCurrent = gradient;
        gradient = ( mu0*(plasmaTotalCurrent*ampere) / (twopi*capRadius*capRadius) ) *meter/tesla;
    }
    else {
        plasmaTotalCurrent = ( (gradient*tesla/meter) * twopi*capRadius*capRadius / mu0 ) / ampere;
    }

    G4cout << "Initialized a MagnetPLASMA1, parameters:" << G4endl;
    G4cout << "\t magnetName         = " << magnetName << G4endl;
    G4cout << "\t Z0                 = " << getZ0()/mm         << " [mm]"  << G4endl;
    G4cout << "\t length             = " << length/mm          << " [mm]"  << G4endl;
    G4cout << "\t gradient           = " << gradient           << " [T/m]" << G4endl;
    G4cout << "\t plasmaTotalcurrent = " << plasmaTotalCurrent << " [A]"   << G4endl;
    G4cout << "\t capRadius          = " << capRadius/mm       << " [mm]"  << G4endl;
    //if the current is given, compute and set the gradient; otherwise compute and set the current
}

G4LogicalVolume* MagnetPLASMA1::Construct() {

    //Build the outer volume (TODO: Factorize this into a common function)
    G4VSolid* mainBox = new G4Box(magnetName+"_mainS",
                                detCon->getWorldSizeX()/2.0, detCon->getWorldSizeX()/2.0, length/2.0);
    // search the material by its name
    G4Material* vacuumMaterial = G4Material::GetMaterial("G4_Galactic");
    if (not vacuumMaterial) {
        G4cerr << "Internal error -- material G4_Galactic not found in MagnetPLASMA1::Construct()!" << G4endl;
        exit(1);
    }
    G4LogicalVolume* mainLV = new G4LogicalVolume(mainBox,vacuumMaterial, magnetName+"_mainLV");

    //Field and gas box
    G4VSolid* fieldBox            = new G4Box(magnetName+"_fieldBoxS",
                                              20.0*mm/2.0, 20.0*mm/2.0, length/2.0);
    G4Material* gasMaterial       = vacuumMaterial; // TODO
    G4LogicalVolume* fieldBoxLV   = new G4LogicalVolume(fieldBox,gasMaterial,magnetName+"_fieldBoxLV");
    G4VPhysicalVolume* fieldBoxPV = new G4PVPlacement(NULL,
                                                      G4ThreeVector(0.0,0.0,0.0),
                                                      fieldBoxLV,
                                                      magnetName + "_fieldBoxPV",
                                                      mainLV,
                                                      false,
                                                      0,
                                                      true);
    field = new FieldPLASMA1(plasmaTotalCurrent, capRadius,
                             G4ThreeVector(0.0, 1.0*mm, getZ0()),fieldBoxLV);
    G4FieldManager* fieldMgr = new G4FieldManager(field);
    G4Mag_UsualEqRhs* fieldEquation = new G4Mag_UsualEqRhs(field);
    G4MagIntegratorStepper* fieldStepper = new G4ClassicalRK4(fieldEquation);
    G4ChordFinder* fieldChordFinder = new G4ChordFinder(field, capRadius/4.0, fieldStepper);
    fieldMgr->SetChordFinder(fieldChordFinder);
    fieldBoxLV->SetFieldManager(fieldMgr,true);

    // The crystal
    G4VSolid* crystalBox      = new G4Box(magnetName+"_crystalBoxS",
                                          10.0*mm/2.0, 20.0*mm/2.0, length/2.0);
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
                                                     fieldBoxLV,
                                                     false,
                                                     0,
                                                     true);

    AddSD(mainLV);

    return mainLV;
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
    gradient = ( mu0*(plasmaTotalCurrent*ampere) / (twopi*capRadius*capRadius) ) *meter/tesla;
}

void FieldPLASMA1::GetFieldValue(const G4double point[4], G4double field[6]) const {

    G4ThreeVector global(point[0],point[1],point[2]);
    G4ThreeVector local = fGlobalToLocal.TransformPoint(global);

    //G4cout << "GetFieldValue, global = " << global << " local = " << local << G4endl;

    G4ThreeVector B(0.0,0.0,0.0); // TODO: Compute in a reasonable manner

    G4double r = sqrt(local[0]*local[0] + local[1]*local[1]);
    G4double Btheta;
    if (r < capRadius) {
        Btheta = gradient*(r/meter) * tesla;
    }
    else {
        Btheta = mu0*plasmaTotalCurrent*ampere / (twopi*capRadius);
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
