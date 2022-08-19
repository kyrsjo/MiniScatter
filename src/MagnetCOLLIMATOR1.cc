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

MagnetCOLLIMATOR1::MagnetCOLLIMATOR1(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                                     std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                                     G4String magnetName_in) :
    MagnetBase(zPos_in, doRelPos_in, length_in, gradient_in, keyValPairs_in, detCon_in, magnetName_in, "COLLIMATOR1"){

    for (auto it : keyValPairs) {
        if (it.first == "radius") {
            radius = ParseDouble(it.second, "collimator radius") * mm;
        }
        else if (it.first == "width") {
            width  = ParseDouble(it.second, "absorber width") * mm;
        }
        else if (it.first == "height") {
            height = ParseDouble(it.second, "absorber height") * mm;
        }
        else if (it.first == "material") {
            absorberMaterialName = it.second;
        }
        else if (it.first == "xOffset" || it.first == "yOffset" || it.first == "xRot" || it.first == "yRot") {
            ParseOffsetRot(it.first, it.second);
        }
        else {
            G4cerr << "MagnetCOLLIMATOR1 did not understand key=value pair '"
                   << it.first << "'='" << it.second << "'." << G4endl;
            exit(1);
        }
    }

    if (gradient != 0.0) {
        G4cerr << "Invalid gradient for COLLIMATOR1: Gradient must be 0.0, but was "
               << gradient << " [T/m]" << G4endl;
        exit(1);
    }

    PrintCommonParameters();
    G4cout << "\t absorberMaterialName    = " << absorberMaterialName <<             G4endl;
    G4cout << "\t radius                  = " << radius/mm            << " [mm]"  << G4endl;
    G4cout << "\t width                   = " << width/mm             << " [mm]"  << G4endl;
    G4cout << "\t height                  = " << height/mm            << " [mm]"  << G4endl;

}

void MagnetCOLLIMATOR1::Construct() {
    if (this->mainLV != NULL) {
        G4Exception("MagnetCOLLIMATOR1::Construct()", "MSDetConMagnetCollimator1000",FatalException,"Internal error -- The mainLV has already been constructed?");
    }

    //Sanity checks on dimensions
    if (width > detCon->getWorldSizeX() || height > detCon->getWorldSizeY()) {
        G4Exception("MagnetCOLLIMATOR1::Construct()", "MSDetConMagnetCollimator1001",FatalException,"The absorber is wider than the world volume.");

    }

    this->mainLV = MakeNewMainLV("main", width,height);

    if (radius > width/2.0 or radius > height/2.0) {
        G4Exception("MagnetCOLLIMATOR1::Construct()", "MSDetConMagnetCollimator1002",FatalException,"The channel doesn't fit in the absorber!");
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
        G4String errormessage = "Error when setting material '" + absorberMaterialName + "', it was not found.\n";
        errormessage += "Valid choices:\n";
        G4MaterialTable* materialTable = G4Material::GetMaterialTable();
        for (auto mat : *materialTable) {
            errormessage += mat->GetName() + "\n";
        }
        G4Exception("MagnetCOLLIMATOR1::Construct()", "MSDetConMagnetCollimator1003",FatalException,errormessage);
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
        G4Exception("MagnetCOLLIMATOR1::Construct()", "MSDetConMagnet1001",FatalException,errormessage);
    }

    ConstructDetectorLV();
    BuildMainPV_transform();
}

