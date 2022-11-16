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

MagnetCOLLIMATORRECT::MagnetCOLLIMATORRECT(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                                     std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                                     G4String magnetName_in) :
    MagnetBase(zPos_in, doRelPos_in, length_in, gradient_in, keyValPairs_in, detCon_in, magnetName_in, "COLLIMATORRECT"){

    for (auto it : keyValPairs) {
        if (it.first == "holeWidth") {
            holeWidth = ParseDouble(it.second, "collimator hole width") * mm;
        }
        else if (it.first == "holeHeight") {
            holeHeight = ParseDouble(it.second, "collimator hole height") * mm;
        }
        else if (it.first == "absorberWidth") {
            absorberWidth  = ParseDouble(it.second, "absorber width") * mm;
        }
        else if (it.first == "absorberHeight") {
            absorberHeight = ParseDouble(it.second, "absorber height") * mm;
        }
        else if (it.first == "material") {
            absorberMaterialName = it.second;
        }
        else if (it.first == "xOffset" || it.first == "yOffset" || it.first == "xRot" || it.first == "yRot") {
            ParseOffsetRot(it.first, it.second);
        }
        else {
            G4cerr << "MagnetCOLLIMATORRECT did not understand key=value pair '"
                   << it.first << "'='" << it.second << "'." << G4endl;
            exit(1);
        }
    }

    if (gradient != 0.0) {
        G4cerr << "Invalid gradient for COLLIMATORRECT: Gradient must be 0.0, but was "
               << gradient << " [T/m]" << G4endl;
        exit(1);
    }

    PrintCommonParameters();
    G4cout << "\t absorberMaterialName    = " << absorberMaterialName <<             G4endl;
    G4cout << "\t hole width              = " << holeWidth/mm         << " [mm]"  << G4endl;
    G4cout << "\t hole height             = " << holeHeight/mm        << " [mm]"  << G4endl;
    G4cout << "\t absorberWidth           = " << absorberWidth/mm     << " [mm]"  << G4endl;
    G4cout << "\t absorberHeight          = " << absorberHeight/mm    << " [mm]"  << G4endl;

}

void MagnetCOLLIMATORRECT::Construct() {
    if (this->mainLV != NULL) {
        G4cerr << "Error in MagnetCOLLIMATORRECT::Construct(): The mainLV has already been constructed?" << G4endl;
        exit(1);
    }

    //Sanity checks on dimensions
    if (absorberWidth > detCon->getWorldSizeX() || absorberHeight > detCon->getWorldSizeY()) {
        G4cerr << "Error in MagnetCOLLIMATORRECT::Construct():" << G4endl
               << " The absorber is wider than the world volume."  << G4endl;
        exit(1);
    }

    this->mainLV = MakeNewMainLV("main", absorberWidth,absorberHeight);

    if (holeWidth > absorberWidth/2.0 or holeHeight > absorberHeight/2.0) {
        G4cerr << "Error in MagnetCOLLIMATORRECT::Construct():" << G4endl
               << " The channel doesn't fit in the absorber!" << G4endl;
        exit(1);
    }

    // Build the absorber
    G4VSolid* absorberBox      = new G4Box(magnetName+"_absorberBoxS",
                                          absorberWidth/2.0, absorberHeight/2.0, length/2.0);
    G4VSolid* holeBox          = new G4Box(magnetName+"_holeBox",
                                          holeWidth/2.0, holeHeight/2.0, (length+1)/2.0); //increase to not have rounding error thick film on entrance
    G4VSolid* absorberSolid    = new G4SubtractionSolid(magnetName+"_absorberS",
                                          absorberBox, holeBox);

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
    G4VPhysicalVolume* absorberPV = new G4PVPlacement  (NULL,
                                                        G4ThreeVector(0.0,0.0,0.0),
                                                        absorberLV,
                                                        magnetName + "_absorberPV",
                                                        mainLV,
                                                        false,
                                                        0,
                                                        false);

    if(absorberPV->CheckOverlaps()) {
        G4String errormessage = "Overlap detected when placing absorberPV for magnet \n"
            "\t'" + magnetName + "' of type '" + magnetType + "'\n"
            "\t, see error message above for more info.";
        G4Exception("MagnetCOLLIMATORRECT::Construct()", "MSDetConMagnet1001",FatalException,errormessage);
    }

    ConstructDetectorLV();
    BuildMainPV_transform();
}