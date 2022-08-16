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
 *  You should have received a cop//180,180 + 240,60 is greaty of the GNU General Public License
 *  along with MiniScatter.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "MagnetClasses.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"

#include "G4PVPlacement.hh"

MagnetPBW::MagnetPBW(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                            std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                            G4String magnetName_in) :
    MagnetBase(zPos_in, doRelPos_in, length_in, gradient_in, keyValPairs_in, detCon_in, magnetName_in, "PBW"){

    for (auto it : keyValPairs) {
        if (it.first == "radius") {
            radius  = ParseDouble(it.second, "target radius") * mm;
        }
        else if (it.first == "material") {
            targetMaterialName = it.second;
        }
        else if (it.first == "al1Thick") {
            al1Thick = ParseDouble(it.second, "Al 1 Thickness") * mm;
        }
        else if (it.first == "waterThick") {
            waterThick = ParseDouble(it.second, "water Thickness") * mm;
        }
        else if (it.first == "al2Thick") {
            al2Thick = ParseDouble(it.second, "Al 2 Thickness") * mm;
        }
        else if (it.first == "xOffset" || it.first == "yOffset" || it.first == "xRot" || it.first == "yRot") {
            ParseOffsetRot(it.first, it.second);
        }
        else {
            G4cerr << "MagnetPBW did not understand key=value pair '"
                   << it.first << "'='" << it.second << "'." << G4endl;
            exit(1);
        }
    }

    if (gradient != 0.0) {
        G4cerr << "Invalid gradient for PBW: Gradient must be 0.0, but was "
               << gradient << " [T/m]" << G4endl;
        exit(1);
    }

    thickness = al1Thick + waterThick + al2Thick;

    //!!!!!Overriding length - if not 0 raise error  -overwrites.
    //Because PBW will be rotated later!
    width = 60; //assumes "width" < radius
    //Assumes a PBW angle range from 30 deg to 150 deg!
    length = (radius + thickness) * 0.5 + (radius + thickness + 5.0); //decimal is sin30 and z translation is needed
    height = 2* (radius + thickness) * 0.86602540378; //decimal is cos30 (sqrt(3)/2)

    PrintCommonParameters();
    G4cout << "\t targetMaterialName      = " << targetMaterialName <<             G4endl;
    G4cout << "\t inner radius            = " << radius/mm          << " [mm]"  << G4endl;
    G4cout << "\t al1Thick                = " << al1Thick/mm        << " [mm]"  << G4endl;
    G4cout << "\t waterThick              = " << waterThick/mm      << " [mm]"  << G4endl;
    G4cout << "\t al2Thick                = " << al2Thick/mm        << " [mm]"  << G4endl;
    G4cout << "\t PBW thickness           = " << thickness/mm       << " [mm]"  << G4endl;
    G4cout << "\t MainLV Width            = " << width/mm           << " [mm]"  << G4endl;
    G4cout << "\t MainLV Height           = " << height/mm          << " [mm]"  << G4endl;
    G4cout << "\t MainLV Length           = " << length/mm          << " [mm]"  << G4endl;

}

void MagnetPBW::Construct() {
    if (this->mainLV != NULL) {
        G4cerr << "Error in MagnetPBW::Construct(): The mainLV has already been constructed?" << G4endl;
        exit(1);
    }

    //Sanity checks on dimensions
    if (radius > detCon->getWorldSizeX()/2 || radius > detCon->getWorldSizeY()/2) {
        G4cerr << "Error in MagnetPBW::Construct():" << G4endl
               << " The absorber is bigger than the world volume."  << G4endl;
        exit(1);
    }

    this->mainLV = MakeNewMainLV("main",length,height);

    // Build the target (PBW)
    G4VSolid* targetSolid      = new G4Tubs(magnetName+"_targetS",
                                            radius, radius + thickness, width*0.5,
                                            30.0*deg, 120.0*deg);

    targetMaterial = G4Material::GetMaterial(targetMaterialName);
    if (not targetMaterial){
        G4cerr << "Error when setting material '"
               << targetMaterialName << "' for MagnetPBW '"
               << magnetName << "' -- not found!" << G4endl;
        G4MaterialTable* materialTable = G4Material::GetMaterialTable();
        G4cerr << "Valid choices:" << G4endl;
        for (auto mat : *materialTable) {
            G4cerr << mat->GetName() << G4endl;
        }
        exit(1);
    }

    //Define rotation so PBW is oriented correct with no user modification
    G4RotationMatrix* pRot = new G4RotationMatrix();
    pRot->rotateX(270.0*deg);

    G4LogicalVolume*  targetLV  = new G4LogicalVolume (targetSolid,targetMaterial, magnetName+"_targetLV");
    G4PVPlacement*    targetPV  = new G4PVPlacement   (pRot,
                                                      G4ThreeVector(0.0,0.0,(radius + thickness) + 5.0), //beam origin is 5mm before PBW
                                                      targetLV,
                                                      magnetName + "_targetPV",
                                                      mainLV,
                                                      false,
                                                      0,
                                                      true);

    G4VSolid*        waterSolid = new G4Tubs          (magnetName+"_waterS",
                                                      radius + al2Thick, radius + (al2Thick + waterThick), width*0.5,
                                                      60.0*deg,60.0*deg);

    G4LogicalVolume*    waterLV = new G4LogicalVolume (waterSolid,G4Material::GetMaterial("G4_WATER"),
                                                      magnetName + "_waterLV");
                                  new G4PVPlacement   (NULL,
                                                      G4ThreeVector(0.0,0.0,0.0),
                                                      waterLV,
                                                      magnetName + "_waterPV",
                                                      targetLV,               //mother volume is targetSolidLV!
                                                      false,
                                                      0,
                                                      true);

    ConstructDetectorLV();
    BuildMainPV_transform();
}