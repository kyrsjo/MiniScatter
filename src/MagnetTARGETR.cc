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

MagnetTARGETR::MagnetTARGETR(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                            std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                            G4String magnetName_in) :
    MagnetBase(zPos_in, doRelPos_in, length_in, gradient_in, keyValPairs_in, detCon_in, magnetName_in, "TARGETR"){

    for (auto it : keyValPairs) {
        if (it.first == "radius") {
            radius  = ParseDouble(it.second, "target radius") * mm;
        }
        else if (it.first == "material") {
            targetMaterialName = it.second;
        }
        else if (it.first == "xOffset" || it.first == "yOffset" || it.first == "xRot" || it.first == "yRot") {
            ParseOffsetRot(it.first, it.second);
        }
        else {
            G4cerr << "MagnetTARGETR did not understand key=value pair '"
                   << it.first << "'='" << it.second << "'." << G4endl;
            exit(1);
        }
    }

    if (gradient != 0.0) {
        G4cerr << "Invalid gradient for TARGETR: Gradient must be 0.0, but was "
               << gradient << " [T/m]" << G4endl;
        exit(1);
    }

    PrintCommonParameters();
    G4cout << "\t targetMaterialName      = " << targetMaterialName <<             G4endl;
    G4cout << "\t radius                  = " << radius/mm          << " [mm]"  << G4endl;

}

void MagnetTARGETR::Construct() {
    if (this->mainLV != NULL) {
        G4cerr << "Error in MagnetTARGETR::Construct(): The mainLV has already been constructed?" << G4endl;
        exit(1);
    }

    //Sanity checks on dimensions
    if (radius > detCon->getWorldSizeX()/2 || radius > detCon->getWorldSizeY()/2) {
        G4cerr << "Error in MagnetTARGETR::Construct():" << G4endl
               << " The absorber is bigger than the world volume."  << G4endl;
        exit(1);
    }

    this->mainLV = MakeNewMainLV("main",2*radius, 2*radius);

    // Build the target
    G4VSolid* targetSolid      = new G4Tubs(magnetName+"_targetS",
                                            0.0, radius, length/2.0,
                                            0.0, 360.0*deg);

    targetMaterial = G4Material::GetMaterial(targetMaterialName);
    if (not targetMaterial){
        G4cerr << "Error when setting material '"
               << targetMaterialName << "' for MagnetTARGETR '"
               << magnetName << "' -- not found!" << G4endl;
        G4MaterialTable* materialTable = G4Material::GetMaterialTable();
        G4cerr << "Valid choices:" << G4endl;
        for (auto mat : *materialTable) {
            G4cerr << mat->GetName() << G4endl;
        }
        exit(1);
    }

    G4LogicalVolume*   targetLV = new G4LogicalVolume(targetSolid,targetMaterial, magnetName+"_targetLV");
    G4PVPlacement*     targetPV = new G4PVPlacement  (NULL,
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
        G4Exception("MagnetTARGETR::Construct()", "MSDetConMagnet1001",FatalException,errormessage);
    }

    ConstructDetectorLV();
    BuildMainPV_transform();
}

