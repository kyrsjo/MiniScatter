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
        if (it.first == "apertureWidth") {
            apertureWidth = ParseDouble(it.second, "collimator aperture width") * mm;
        }
        else if (it.first == "apertureHeight") {
            apertureHeight = ParseDouble(it.second, "collimator aperture height") * mm;
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
            G4ExceptionDescription errormessage;
            errormessage << "MagnetCOLLIMATORRECT did not understand key=value pair '"
                   << it.first << "'='" << it.second << "'." << G4endl;
            G4Exception("MagnetCOLLIMATORRECT::Construct()", "MSDetConMagnet1001",FatalException,errormessage);
        }
    }

    if (gradient != 0.0) {
        G4ExceptionDescription errormessage;
        errormessage << "Invalid gradient for COLLIMATORRECT: Gradient must be 0.0, but was "
               << gradient << " [T/m]" << G4endl;
        G4Exception("MagnetCOLLIMATORRECT::Construct()", "MSDetConMagnet1002",FatalException,errormessage);
    }

    PrintCommonParameters();
    G4cout << "\t absorberMaterialName    = " << absorberMaterialName <<             G4endl;
    G4cout << "\t aperture width          = " << apertureWidth/mm     << " [mm]"  << G4endl;
    G4cout << "\t aperture height         = " << apertureHeight/mm    << " [mm]"  << G4endl;
    G4cout << "\t absorberWidth           = " << absorberWidth/mm     << " [mm]"  << G4endl;
    G4cout << "\t absorberHeight          = " << absorberHeight/mm    << " [mm]"  << G4endl;

}

void MagnetCOLLIMATORRECT::Construct() {
    if (this->mainLV != NULL) {
        G4ExceptionDescription errormessage;
        errormessage << "Error in MagnetCOLLIMATORRECT::Construct(): The mainLV has already been constructed?" << G4endl;
        G4Exception("MagnetCOLLIMATORRECT::Construct()", "MSDetConMagnet1003",FatalException,errormessage);
    }

    //Sanity checks on dimensions
    if (absorberWidth > detCon->getWorldSizeX() || absorberHeight > detCon->getWorldSizeY()) {
        G4ExceptionDescription errormessage;
        errormessage << "Error in MagnetCOLLIMATORRECT::Construct():" << G4endl
               << " The absorber is wider than the world volume."  << G4endl;
        G4Exception("MagnetCOLLIMATORRECT::Construct()", "MSDetConMagnet1004",FatalException,errormessage);
    }

    this->mainLV = MakeNewMainLV("main", absorberWidth,absorberHeight);

    if (apertureWidth > absorberWidth or apertureHeight > absorberHeight) {
        G4ExceptionDescription errormessage;
        errormessage << "Error in MagnetCOLLIMATORRECT::Construct():" << G4endl
               << " The channel doesn't fit in the absorber!" << G4endl;
        G4Exception("MagnetCOLLIMATORRECT::Construct()", "MSDetConMagnet1005",FatalException,errormessage);
    }

    // Build the absorber
    G4VSolid* absorberBox      = new G4Box(magnetName+"_absorberBoxS",
                                          absorberWidth/2.0, absorberHeight/2.0, length/2.0);
    //Note: Length of the cutout is increased by 1 mm to not have rounding error "film" on entrance
    G4VSolid* apertureBox      = new G4Box(magnetName+"_apertureBox",
                                          apertureWidth/2.0, apertureHeight/2.0, (length+1*mm)/2.0);
    G4VSolid* absorberSolid    = new G4SubtractionSolid(magnetName+"_absorberS",
                                          absorberBox, apertureBox);

    absorberMaterial = G4Material::GetMaterial(absorberMaterialName);
    if (not absorberMaterial){
        G4ExceptionDescription errormessage;
        errormessage << "Error when setting material '"
               << absorberMaterialName << "' for MagnetCollimator '"
               << magnetName << "' -- not found!" << G4endl;
        G4MaterialTable* materialTable = G4Material::GetMaterialTable();
        errormessage << "Valid choices:" << G4endl;
        for (auto mat : *materialTable) {
            errormessage << mat->GetName() << G4endl;
        }
        G4Exception("MagnetCOLLIMATORRECT::Construct()", "MSDetConMagnet1006",FatalException,errormessage);
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
        G4ExceptionDescription errormessage;
        errormessage << "Overlap detected when placing absorberPV for magnet \n"
            << "\t'" << magnetName << "' of type '" << magnetType << "'\n"
            << "\t, see error message above for more info.";
        G4Exception("MagnetCOLLIMATORRECT::Construct()", "MSDetConMagnet1007",FatalException,errormessage);
    }

    ConstructDetectorLV();
    BuildMainPV_transform();
}