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
#ifndef MagnetClasses_h$
#define MagnetClasses_h$

#include "G4SystemOfUnits.hh"

#include "DetectorConstruction.hh"
#include "FieldClasses.hh"

class MagnetBase {
public:
    static MagnetBase* MagnetFactory(G4String inputString, DetectorConstruction* detCon, G4String magnetName);

    MagnetBase(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
               std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
               G4String magnetName_in) :
        zPos(zPos_in), doRelPos(doRelPos_in), length(length_in),gradient(gradient_in),
        keyValPairs(keyValPairs_in), detCon(detCon_in), magnetName(magnetName_in) {};

    virtual G4double getZ0() {
        //Get the global z position of the center of the magnet in G4 units.

        if (doRelPos) {
            return detCon->getTargetThickness()/2.0 + zPos + length/2.0;
        }
        else {
            return zPos;
        }
    };

    virtual void PostInitialize() {
        if (field != NULL) {
            // There are "magnets" with no field
            field->PostInitialize();
        }
    }

    G4double GetLength()    const { return length;  };
    G4double GetXOffset()   const { return xOffset; };
    G4double GetYOffset()   const { return yOffset; };

protected:
    G4double zPos;     // [G4 units]
    G4bool   doRelPos;
    G4double length;   // [G4units]
    G4double gradient; // [T/m]

    G4double xOffset = 0.0; // [G4 length units]
    G4double yOffset = 0.0; // [G4 length units]

    //Rotations are applied after moving, around the new position
    G4double xRot = 0.0; // Rotation around horizontal axis [G4 angle units]
    G4double yRot = 0.0; // Rotation around vertical axis   [G4 angle units]

    std::map<G4String,G4String> keyValPairs;
    DetectorConstruction* detCon;

    FieldBase* field = NULL;

    G4LogicalVolume* mainLV = NULL;
    G4LogicalVolume* MakeNewMainLV(G4String name_postfix);

    G4double mainLV_w = 0.0; // Width  of mainLV after removing what is needed for trans and rot [G4 units]
    G4double mainLV_h = 0.0; // Height of mainLV after removing what is needed for trans and rot [G4 units]

    G4Transform3D mainPV_transform;
    void BuildMainPV_transform();

    virtual void ConstructDetectorLV();
    G4LogicalVolume* detectorLV = NULL;
    G4VSensitiveDetector* magnetSD = NULL;

public:
    const G4String magnetName;

    const G4Transform3D GetMainPV_transform() const {return mainPV_transform;};

    virtual void Construct() = 0;
    G4LogicalVolume* GetMainLV() const;
    G4LogicalVolume* GetDetectorLV() const;
    void AddSD(); // Adds an SD to the detectorLV
};

class MagnetPLASMA1 : public MagnetBase {
public:
    MagnetPLASMA1(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                  std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                  G4String magnetName_in);

    virtual void Construct();
private:
    G4double plasmaTotalCurrent; // [A]
    G4double capRadius;          // [G4 length units]
    G4double cryWidth;           // [G4 length units]
    G4double cryHeight;          // [G4 length units]
};

class MagnetCOLLIMATOR1 : public MagnetBase {
public:
    MagnetCOLLIMATOR1(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                      std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                      G4String magnetName_in);
    virtual void Construct();
private:
    G4String absorberMaterialName;
    G4Material* absorberMaterial = NULL;
    G4double width  = 10.0*mm;  //[G4 length units]
    G4double height = 50.0*mm;  //[G4 length units]
    G4double radius = 50.0*mm;  //[G4 length units]
};

class MagnetTARGET : public MagnetBase {
public:
    MagnetTARGET(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                      std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                      G4String magnetName_in);
    virtual void Construct();
private:
    G4String targetMaterialName;
    G4Material* targetMaterial = NULL;
    G4double width  = 10.0*mm;  //[G4 length units]
    G4double height = 10.0*mm;  //[G4 length units]
};

#endif
