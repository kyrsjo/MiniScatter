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
        G4cerr << "Error in MagnetCOLLIMATOR1::Construct(): The mainLV has already been constructed?" << G4endl;
        exit(1);
    }

    //Sanity checks on dimensions
    if (width > detCon->getWorldSizeX() || height > detCon->getWorldSizeY()) {
        G4cerr << "Error in MagnetCOLLIMATOR1::Construct():" << G4endl
               << " The absorber is wider than the world volume."  << G4endl;
        exit(1);
    }

    this->mainLV = MakeNewMainLV("main", width,height);

    if (radius > width/2.0 or radius > height/2.0) {
        G4cerr << "Error in MagnetCOLLIMATOR1::Construct():" << G4endl
               << " The channel doesn't fit in the absorber!" << G4endl;
        exit(1);
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
    //G4VPhysicalVolume* absorberPV =
                                    new G4PVPlacement  (NULL,
                                                        G4ThreeVector(0.0,0.0,0.0),
                                                        absorberLV,
                                                        magnetName + "_absorberPV",
                                                        mainLV,
                                                        false,
                                                        0,
                                                        true);

    ConstructDetectorLV();
    BuildMainPV_transform();
}

