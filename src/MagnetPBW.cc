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

#include "G4Tubs.hh"

#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"

#include <cmath>

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
        else if (it.first == "arcPhi") {
            arcPhi = ParseDouble(it.second, "Window Arc Angle") * deg;
        }
        else if (it.first == "width") {
            arcPhi = ParseDouble(it.second, "Window Width") * mm;
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

    //Error catching
    if (length != 0.0) {
        G4cerr << "Invalid length for PBW: Length must be 0.0, but was "
               << length / mm << " [mm]" << G4endl; 
        exit(1);
    }
    if (radius <= 0.0) {
        G4cerr << "Invalid radius for PBW: Radius must be > 0.0, but was "
               << radius / mm << " [mm]" << G4endl; 
        exit(1);
    }
    if (al1Thick <= 0.0) {
        G4cerr << "Invalid al1Thick for PBW: Al1Thick must be > 0.0, but was "
               << al1Thick / mm << " [mm]" << G4endl; 
        exit(1);
    }
    if (waterThick <= 0.0) {
        G4cerr << "Invalid waterThick for PBW: WaterThick must be > 0.0, but was "
               << waterThick / mm << " [mm]" << G4endl; 
        exit(1);
    }
    if (al2Thick <= 0.0) {
        G4cerr << "Invalid al2Thick for PBW: Al2Thick must be > 0.0, but was "
               << al2Thick / mm << " [mm]" << G4endl; 
        exit(1);
    }
    if (width <= 0.0) {
        G4cerr << "Invalid width for PBW: Width must be > 0.0, but was "
               << width / mm << " [mm]" << G4endl; 
        exit(1);
    }
    if (arcPhi / deg < 0.0 || arcPhi / deg > 180.0) {
        G4cerr << "Invalid arc angle for PBW: ArcPhi must be within: 0 <= arcPhi <= 180, but was "
               << arcPhi / deg << " [deg]" << G4endl; 
        exit(1);
    }

    //Calculate dimensions for mainLV box and positioning
    thickness = al1Thick + waterThick + al2Thick;
    startPhi = rightAng/rad - (arcPhi/rad * 0.5);
    waterStartPhi = rightAng/rad - (arcPhi/rad * 0.25);
    length = radius * (1 - cos(arcPhi/rad * 0.5)) + thickness;
    height = 2 * sin(arcPhi/rad * 0.5) * (radius + thickness);
    boxCenter = radius * cos(arcPhi/rad * 0.5) + length * 0.5;

    PrintCommonParameters();
    G4cout << "\t targetMaterialName      = " << targetMaterialName <<             G4endl;
    G4cout << "\t inner radius            = " << radius/mm          << " [mm]"  << G4endl;
    G4cout << "\t al1Thick                = " << al1Thick/mm        << " [mm]"  << G4endl;
    G4cout << "\t waterThick              = " << waterThick/mm      << " [mm]"  << G4endl;
    G4cout << "\t al2Thick                = " << al2Thick/mm        << " [mm]"  << G4endl;
    G4cout << "\t PBW thickness           = " << thickness/mm       << " [mm]"  << G4endl;
    G4cout << "\t arcPhi                  = " << arcPhi/deg         << " [deg]" << G4endl;
    G4cout << "\t PBW Width               = " << width/mm           << " [mm]"  << G4endl;
    G4cout << "\t Calculated Values:" << G4endl;
    G4cout << "\t MainLV Height           = " << height/mm          << " [mm]"  << G4endl;
    G4cout << "\t MainLV Length           = " << length/mm          << " [mm]"  << G4endl;
    G4cout << "\t Box Center              = " << boxCenter/mm       << " [mm]"  << G4endl;
    G4cout << "\t The PBW first surface " << G4endl;
    G4cout << "\t    position is shifted by " << -length/mm * 0.5   <<  " [mm]"  << G4endl;

}

void MagnetPBW::Construct() {
    if (this->mainLV != NULL) {
        G4cerr << "Error in MagnetPBW::Construct(): The mainLV has already been constructed?" << G4endl;
        exit(1);
    }

    //Sanity checks on dimensions
    if (radius > detCon->getWorldSizeX()/2 || radius > detCon->getWorldSizeY()/2) {
        G4cerr << "Error in MagnetPBW::Construct():" << G4endl
               << "  The absorber radius is bigger than the world volume (--worldsize)."  << G4endl;
        exit(1);
    }

    this->mainLV = MakeNewMainLV("main",width,height);

    // Build the target (PBW)
    G4VSolid* targetSolid = new G4Tubs(magnetName+"_targetS",
                                      radius,
                                      radius + thickness,
                                      width * 0.5,
                                      startPhi/rad,
                                      arcPhi/rad);

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

    G4RotationMatrix* pRot = new G4RotationMatrix();
    pRot->rotateX(90.0*deg);
    pRot->rotateY(90.0*deg);

    G4LogicalVolume*  targetLV  = new G4LogicalVolume (targetSolid,
                                                      targetMaterial,
                                                      magnetName+"_targetLV");

    G4PVPlacement*    targetPV  = new G4PVPlacement   (pRot,
                                                      G4ThreeVector(0.0, 0.0, boxCenter),
                                                      targetLV,
                                                      magnetName + "_targetPV",
                                                      mainLV,
                                                      false,
                                                      0,
                                                      true);
    if(targetPV->CheckOverlaps()) {
        G4String errormessage = "Overlap detected when placing targetPV for magnet \n"
            "\t'" + magnetName + "' of type '" + magnetType + "'\n"
            "\t, see error message above for more info.";
        G4Exception("MagnetPBW::Construct()", "MSDetConMagnet1001",FatalException,errormessage);
    }

    G4VSolid*        waterSolid = new G4Tubs          (magnetName+"_waterS",
                                                      radius + al2Thick,
                                                      radius + (al2Thick + waterThick),
                                                      width * 0.5,
                                                      waterStartPhi/rad,
                                                      arcPhi/rad * 0.5);

    G4LogicalVolume*    waterLV = new G4LogicalVolume (waterSolid,G4Material::GetMaterial("G4_WATER"),
                                                      magnetName + "_waterLV");

    G4PVPlacement*      waterPV = new G4PVPlacement   (NULL,
                                                      G4ThreeVector(0.0, 0.0, 0.0),
                                                      waterLV,
                                                      magnetName + "_waterPV",
                                                      targetLV,           //mother volume is targetLV!
                                                      false,
                                                      0,
                                                      true);
    if(waterPV->CheckOverlaps()) {
        G4String errormessage = "Overlap detected when placing waterPV for magnet \n"
            "\t'" + magnetName + "' of type '" + magnetType + "'\n"
            "\t, see error message above for more info.";
        G4Exception("MagnetPBW::Construct()", "MSDetConMagnet1001",FatalException,errormessage);
    }

    //Set color and line segments per circle for Visualization
    G4VisAttributes* aluminum = new G4VisAttributes(G4Colour(0.66,0.67,0.71));
    aluminum->SetForceLineSegmentsPerCircle(100);
    targetLV->SetVisAttributes(aluminum);
    
    G4VisAttributes* water = new G4VisAttributes(G4Colour(0,1,1));
    water->SetForceLineSegmentsPerCircle(100);
    waterLV->SetVisAttributes(water);

    ConstructDetectorLV();
    BuildMainPV_transform();
}