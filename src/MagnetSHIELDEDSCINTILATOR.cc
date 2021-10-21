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

MagnetSHIELDEDSCINTILLATOR::MagnetSHIELDEDSCINTILLATOR(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                                                       std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                                                       G4String magnetName_in) :
    MagnetBase(zPos_in, doRelPos_in, length_in, gradient_in, keyValPairs_in, detCon_in, magnetName_in, "SHIELDEDSCINTILLATOR") {

    useGhostDetector = false;

    for (auto it : keyValPairs) {
        if (it.first == "scintMat") {
            scintillatorMaterialName = it.second;
        }
        else if (it.first == "shieldMat") {
            shieldingMaterialName = it.second;
        }
        else if (it.first == "rScint") {
            r_scint  = ParseDouble(it.second, "scintillator radius") * mm;
        }
        else if (it.first == "lScint") {
            l_scint  = ParseDouble(it.second, "scintillator length") * mm;
        }
        else if (it.first == "zScint") {
            z_scint  = ParseDouble(it.second, "scintillator position") * mm;
        }
        else if (it.first == "riShield") {
            ri_shield  = ParseDouble(it.second, "shield inner radius") * mm;
        }
        else if (it.first == "roShield") {
            ro_shield  = ParseDouble(it.second, "shield outer radius") * mm;
        }
        else if (it.first == "xOffset" || it.first == "yOffset" || it.first == "xRot" || it.first == "yRot") {
            ParseOffsetRot(it.first, it.second);
        }
        else {
            G4cerr << "MagnetSHIELDEDSCINTILLATOR did not understand key=value pair '"
                   << it.first << "'='" << it.second << "'." << G4endl;
            exit(1);
        }
    }

    if (gradient != 0.0) {
        G4cerr << "Invalid gradient for MagnetSHIELDEDSCINTILLATOR: Gradient must be 0.0, but was "
               << gradient << " [T/m]" << G4endl;
        exit(1);
    }

    PrintCommonParameters();
    G4cout << "\t scintillatorMaterialName= " << scintillatorMaterialName <<             G4endl;
    G4cout << "\t shieldingMaterialName   = " << shieldingMaterialName    <<             G4endl;
    G4cout << "\t r_scint                 = " << r_scint/mm               << " [mm]"  << G4endl;
    G4cout << "\t l_scint                 = " << l_scint/mm               << " [mm]"  << G4endl;
    G4cout << "\t z_scint                 = " << z_scint/mm               << " [mm]"  << G4endl;
    G4cout << "\t ri_shield               = " << ri_shield/mm             << " [mm]"  << G4endl;
    G4cout << "\t ro_shield               = " << ro_shield/mm             << " [mm]"  << G4endl;

}

void MagnetSHIELDEDSCINTILLATOR::Construct() {
    if (this->mainLV != NULL) {
        G4cerr << "Error in MagnetSHIELDEDSCINTILLATOR::Construct(): The mainLV has already been constructed?" << G4endl;
        exit(1);
    }

    //Sanity checks on dimensions; it could still fit with rotations etc.
    if (ro_shield > detCon->getWorldSizeX()/2 || ro_shield > detCon->getWorldSizeY()/2) {
        G4cerr << "Error in MagnetSHIELDEDSCINTILLATOR::Construct():" << G4endl
               << " The shield radius is too big for the  world volume." << G4endl;
        G4cerr << "mainLW_w = " << mainLV_w/mm << " [mm]" << G4endl;
        G4cerr << "mainLW_h = " << mainLV_h/mm << " [mm]" << G4endl;
        exit(1);
    }
    if (ri_shield >= ro_shield || r_scint >= ri_shield) {
        G4cerr << "Error in MagnetSHIELDEDSCINTILLATOR::Construct():" << G4endl
               << " Nesting of radii is wrong." << G4endl
               << " Expected ro_shield >= ri_shield >= r_scint." << G4endl;
        G4cerr << "r_scint   = " << r_scint/mm   << " [mm]" << G4endl;
        G4cerr << "ri_shield = " << ri_shield/mm << " [mm]" << G4endl;
        G4cerr << "ro_shield = " << ro_shield/mm << " [mm]" << G4endl;
        exit(1);
    }
    if (l_scint > length || z_scint+l_scint/2 > length/2)  {
        G4cerr << "Error in MagnetSHIELDEDSCINTILLATOR::Construct():" << G4endl
               << " Nesting of lengths or scintillator position is wrong." << G4endl
               << " Expected l_scint <= length and l_scint+z_scint/2 <= length/2." << G4endl;
        G4cerr << "length  = "  << length/mm   << " [mm]" << G4endl;
        G4cerr << "l_scint = "  << l_scint/mm  << " [mm]" << G4endl;
        G4cerr << "z_scint = "  << z_scint/mm  << " [mm]" << G4endl;
        exit(1);
    }

    this->mainLV = MakeNewMainLV("main", 2*ro_shield, 2*ro_shield);

    // Build the object
    scintillatorMaterial = G4Material::GetMaterial(scintillatorMaterialName);
    if (not scintillatorMaterial){
        G4cerr << "Error when setting scintillatorMaterial '"
               << scintillatorMaterial << "' for MagnetSHIELDEDSCINTILLATOR '"
               << magnetName << "' -- not found!" << G4endl;
        G4MaterialTable* materialTable = G4Material::GetMaterialTable();
        G4cerr << "Valid choices:" << G4endl;
        for (auto mat : *materialTable) {
            G4cerr << mat->GetName() << G4endl;
        }
        exit(1);
    }

    shieldingMaterial = G4Material::GetMaterial(shieldingMaterialName);
    if (not shieldingMaterial){
        G4cerr << "Error when setting shieldingMaterial '"
               << shieldingMaterial << "' for MagnetSHIELDEDSCINTILLATOR '"
               << magnetName << "' -- not found!" << G4endl;
        G4MaterialTable* materialTable = G4Material::GetMaterialTable();
        G4cerr << "Valid choices:" << G4endl;
        for (auto mat : *materialTable) {
            G4cerr << mat->GetName() << G4endl;
        }
        exit(1);
    }

    // Scintillator
    G4VSolid*         scintillatorSolid = new G4Tubs         (magnetName+"_scintS",
                                                              0.0, r_scint, l_scint/2.0,
                                                              0.0, 360.0*deg);

                      scintillatorLV    = new G4LogicalVolume(scintillatorSolid,scintillatorMaterial, magnetName+"_scintLV");

                                          new G4PVPlacement  (NULL,
                                                              G4ThreeVector(0.0,0.0,z_scint),
                                                              scintillatorLV,
                                                              magnetName + "_scintPV",
                                                              mainLV,
                                                              false,
                                                              0,
                                                              true);

    //Shielding
    G4VSolid* shieldingCylOuter = new G4Tubs(magnetName+"_shieldingCylOuterS",
                                             0.0, ro_shield,length/2,
                                             0.0, 360.0*deg);
    G4VSolid* shieldingCylInner = new G4Tubs(magnetName+"_shieldingCylInnerS",
                                             0.0, ri_shield,length/2,
                                             0.0, 360.0*deg);
    G4VSolid* shieldingSolid    = new G4SubtractionSolid(magnetName+"_shieldingS",
                                                         shieldingCylOuter, shieldingCylInner);

                    shieldingLV = new G4LogicalVolume(shieldingSolid, shieldingMaterial, magnetName+"_shieldingLV");

                                  new G4PVPlacement(NULL,
                                                    G4ThreeVector(0.0,0.0,0.0),
                                                    shieldingLV,
                                                    magnetName+"shieldingPV",
                                                    mainLV,
                                                    false,
                                                    0,
                                                    true);

    //Other
    ConstructDetectorLV();
    BuildMainPV_transform();
}

void MagnetSHIELDEDSCINTILLATOR::ConstructDetectorLV() {
    if (this->detectorLV != NULL) {
        G4cerr << "Error in MagnetSHIELDEDSCINTILLATOR::ConstructDetectorLV(): The detectorLV has already been constructed?" << G4endl;
        exit(1);
    }
    if(this->mainLV == NULL) {
        G4cerr << "Error in MagnetSHIELDEDSCINTILLATOR::ConstructDetectorLV(): The mainLV is not yet constructed?" << G4endl;
        exit(1);
    }

    //In this case, we just want the scintillator to be detecting things (not the lead shielding), so there is no "overall" detectorLV
    this->detectorLV = this->scintillatorLV;
    //MagnetBase::ConstructDetectorLV();
}