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

MagnetCOLLIMATORHV::MagnetCOLLIMATORHV(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                                       std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                                       G4String magnetName_in) :
    MagnetBase(zPos_in, doRelPos_in, length_in, gradient_in, keyValPairs_in, detCon_in, magnetName_in, "TARGETR"){

    for (auto it : keyValPairs) {
        if (it.first == "gap") {
            gap  = ParseDouble(it.second, "collimator gap") * mm;
        }
        else if (it.first == "material") {
            targetMaterialName = it.second;
        }
        else if (it.first == "HV") {
            if (it.second == "H") {
                isH = true;
            }
            else if (it.second == "V") {
                isH = false;
            }
            else {
                G4cerr << "Invalid argument when reading horizontal/vertical flag" << G4endl
                       << "Got: '" << it.second << "'" << G4endl
                       << "Expected 'H' or 'V'."   << G4endl;
                exit(1);
            }
        }
        else if (it.first == "jawThick") {
            jawThick  = ParseDouble(it.second, "collimator jaw thickness") * mm;
        }
        else if (it.first == "jawHeight") {
            jawHeight = ParseDouble(it.second, "collimator jaw height") * mm;
        }
        else if (it.first == "xOffset" || it.first == "yOffset" || it.first == "xRot" || it.first == "yRot") {
            ParseOffsetRot(it.first, it.second);
        }
        else {
            G4cerr << "MagnetCOLLIMATORHV did not understand key=value pair '"
                   << it.first << "'='" << it.second << "'." << G4endl;
            exit(1);
        }
    }

    if (gradient != 0.0) {
        G4cerr << "Invalid gradient for MagnetCOLLIMATORHV: Gradient must be 0.0, but was "
               << gradient << " [T/m]" << G4endl;
        exit(1);
    }

    PrintCommonParameters();
    G4cout << "\t HV                      = " << (isH ? 'H' : 'V')  <<             G4endl;
    G4cout << "\t targetMaterialName      = " << targetMaterialName <<             G4endl;
    G4cout << "\t gap                     = " << gap/mm             << " [mm]"  << G4endl;
    G4cout << "\t jaw thickness           = " << jawThick/mm        << " [mm]"  << G4endl;
    G4cout << "\t jaw height              = " << jawHeight/mm       << " [mm]"  << G4endl;

}

void MagnetCOLLIMATORHV::Construct() {
    if (this->mainLV != NULL) {
        G4cerr << "Error in MagnetCOLLIMATORHV::Construct(): The mainLV has already been constructed?" << G4endl;
        exit(1);
    }

    this->mainLV = MakeNewMainLV("main");

    //Sanity checks on dimensions
    if (gap+2*jawThick > mainLV_w || mainLV_h > mainLV_h) {
        G4cerr << "Error in MagnetCOLLIMATORHV::Construct():" << G4endl
               << " The absorber is bigger than it's allowed envelope "
               << " including offsets and rotations."  << G4endl;
        G4cerr << "mainLW_w = " << mainLV_w/mm << " [mm]" << G4endl;
        G4cerr << "mainLW_h = " << mainLV_h/mm << " [mm]" << G4endl;
        exit(1);
    }

    // Build the target
    targetMaterial = G4Material::GetMaterial(targetMaterialName);
    if (not targetMaterial){
        G4cerr << "Error when setting material '"
               << targetMaterialName << "' for MagnetCOLLIMATORHV '"
               << magnetName << "' -- not found!" << G4endl;
        G4MaterialTable* materialTable = G4Material::GetMaterialTable();
        G4cerr << "Valid choices:" << G4endl;
        for (auto mat : *materialTable) {
            G4cerr << mat->GetName() << G4endl;
        }
        exit(1);
    }

    G4VSolid* jawPlu  = new G4Box(magnetName+"_jawPluS", (isH?jawThick:jawHeight)/2,(isH?jawHeight:jawThick)/2, length/2);
    G4VSolid* jawMin  = new G4Box(magnetName+"_jawMinS", (isH?jawThick:jawHeight)/2,(isH?jawHeight:jawThick)/2, length/2);

    G4LogicalVolume*   targetPluLV = new G4LogicalVolume(jawPlu,targetMaterial, magnetName+"_targetPluLV");
    G4LogicalVolume*   targetMinLV = new G4LogicalVolume(jawMin,targetMaterial, magnetName+"_targetMinLV");
    
    G4double offset = gap/2+jawThick/2;
    new G4PVPlacement  (NULL,
                        G4ThreeVector(isH?offset:0.0,isH?0.0:offset,0.0),
                        targetPluLV,
                        magnetName + "_targetPluPV",
                        mainLV,
                        false,
                        0,
                        true);
    new G4PVPlacement  (NULL,
                        G4ThreeVector(isH?-offset:0.0,isH?0.0:-offset,0.0),
                        targetMinLV,
                        magnetName + "_targetMinPV",
                        mainLV,
                        false,
                        0,
                        true);

    ConstructDetectorLV();
    BuildMainPV_transform();
}

