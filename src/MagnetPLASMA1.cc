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

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"

#include "G4PVPlacement.hh"

#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4Mag_UsualEqRhs.hh"

#include "G4PhysicalConstants.hh"


MagnetPLASMA1::MagnetPLASMA1(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                             std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                             G4String magnetName_in) :
    MagnetBase(zPos_in, doRelPos_in, length_in, gradient_in, keyValPairs_in, detCon_in, magnetName_in, "PLASMA1") {

    // Default values
    G4bool inputIsTotalAmps = false;
    capRadius = 1.0*mm;
    cryWidth  = 10.0*mm;
    cryHeight = 20.0*mm;

    for (auto it : keyValPairs) {
        if (it.first == "radius") {
            capRadius = ParseDouble(it.second, "capillary radius") * mm;
        }
        else if (it.first == "totalAmps") {
            inputIsTotalAmps = ParseBool(it.second, "capillary totalAmps");
        }
        else if (it.first == "width") {
            cryWidth = ParseDouble(it.second, "crystal width") * mm;
        }
        else if (it.first == "height") {
            cryWidth = ParseDouble(it.second, "crystal height") * mm;
        }
        else if (it.first == "xOffset" || it.first == "yOffset" || it.first == "xRot" || it.first == "yRot") {
            ParseOffsetRot(it.first, it.second);
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

    PrintCommonParameters();
    G4cout << "\t plasmaTotalcurrent      = " << plasmaTotalCurrent << " [A]"   << G4endl;
    G4cout << "\t capRadius               = " << capRadius/mm       << " [mm]"  << G4endl;
    G4cout << "\t cryWidth                = " << cryWidth/mm        << " [mm]"  << G4endl;
    G4cout << "\t cryHeight               = " << cryHeight/mm       << " [mm]"  << G4endl;


}

void MagnetPLASMA1::Construct() {

    if (this->mainLV != NULL) {
        G4cerr << "Error in MagnetPLASMA1::Construct(): The mainLV has already been constructed?" << G4endl;
        exit(1);
    }

    //Sanity checks on dimensions
    if (cryWidth > detCon->getWorldSizeX() || cryHeight > detCon->getWorldSizeY()) {
        G4cerr << "Error in MagnetPLASMA1::Construct():" << G4endl
               << " The absorber is wider than the world volume."  << G4endl;
        exit(1);
    }

    this->mainLV = MakeNewMainLV("main", cryWidth,cryHeight);

    //Field
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

    sapphireMaterial = G4Material::GetMaterial("Sapphire");
    if (not sapphireMaterial) {
        G4cerr << "Internal error -- material Sapphire not found in MagnetPLASMA1::Construct()!" << G4endl;
        exit(1);
    }
    G4LogicalVolume*   crystalLV = new G4LogicalVolume(crystalSolid,sapphireMaterial, magnetName+"_crystalLV");
    //G4VPhysicalVolume* crystalPV =
                                   new G4PVPlacement  (NULL,
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


/** FIELD PATTERN CLASS **/

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

