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
#include "G4SubtractionSolid.hh"

#include "G4PVPlacement.hh"

MagnetTARGET::MagnetTARGET(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                                     std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                                     G4String magnetName_in) :
    MagnetBase(zPos_in, doRelPos_in, length_in, gradient_in, keyValPairs_in, detCon_in, magnetName_in, "TARGET"){

    for (auto it : keyValPairs) {
        if (it.first == "width") {
            width  = ParseDouble(it.second, "target width") * mm;
        }
        else if (it.first == "height") {
            height = ParseDouble(it.second, "target height") * mm;
        }
        else if (it.first == "material") {
            targetMaterialName = it.second;
        }
        else if (it.first == "xOffset" || it.first == "yOffset" || it.first == "xRot" || it.first == "yRot") {
            ParseOffsetRot(it.first, it.second);
        }
        else {
            G4cerr << "MagnetTARGET did not understand key=value pair '"
                   << it.first << "'='" << it.second << "'." << G4endl;
            exit(1);
        }
    }

    if (gradient != 0.0) {
        G4cerr << "Invalid gradient for TARGET: Gradient must be 0.0, but was "
               << gradient << " [T/m]" << G4endl;
        exit(1);
    }

    PrintCommonParameters();
    G4cout << "\t targetMaterialName      = " << targetMaterialName <<             G4endl;
    G4cout << "\t width                   = " << width/mm           << " [mm]"  << G4endl;
    G4cout << "\t height                  = " << height/mm          << " [mm]"  << G4endl;

}

void MagnetTARGET::Construct() {
    if (this->mainLV != NULL) {
        G4cerr << "Error in MagnetTARGET::Construct(): The mainLV has already been constructed?" << G4endl;
        exit(1);
    }

    //Sanity checks on dimensions
    if (width > detCon->getWorldSizeX() || height > detCon->getWorldSizeY()) {
        G4cerr << "Error in MagnetTARGET::Construct():" << G4endl
               << " The absorber is bigger than the world volume."  << G4endl;
        exit(1);
    }

    this->mainLV = MakeNewMainLV("main",width,height);

    // Build the target
    G4VSolid* targetSolid      = new G4Box(magnetName+"_targetS",
                                           width/2.0, height/2.0, length/2.0);

    targetMaterial = G4Material::GetMaterial(targetMaterialName);
    if (not targetMaterial){
        G4cerr << "Error when setting material '"
               << targetMaterialName << "' for MagnetTARGET '"
               << magnetName << "' -- not found!" << G4endl;
        G4MaterialTable* materialTable = G4Material::GetMaterialTable();
        G4cerr << "Valid choices:" << G4endl;
        for (auto mat : *materialTable) {
            G4cerr << mat->GetName() << G4endl;
        }
        exit(1);
    }

    G4LogicalVolume*   targetLV = new G4LogicalVolume(targetSolid,targetMaterial, magnetName+"_targetLV");
    G4VPhysicalVolume* targetPV = new G4PVPlacement  (NULL,
                                                      G4ThreeVector(0.0,0.0,0.0),
                                                      targetLV,
                                                      magnetName + "_targetPV",
                                                      mainLV,
                                                      false,
                                                      0,
                                                      false);
    if(targetPV->CheckOverlaps()) {
        G4String errormessage = "Overlap detected when placing targetPV for magnet \n"
            "\t'" + magnetName + "' of type '" + magnetType + "'\n"
            "\t, see error message above for more info.";
        G4Exception("MagnetTARGET::Construct()", "MSDetConMagnet1001",FatalException,errormessage);
    }

    ConstructDetectorLV();
    BuildMainPV_transform();
}

