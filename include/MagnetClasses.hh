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
#ifndef MagnetClasses_h
#define MagnetClasses_h

#include "G4SystemOfUnits.hh"

#include "DetectorConstruction.hh"
//#include "FieldClasses.hh"

// For the field classes
#include "G4MagneticField.hh"
#include "G4Navigator.hh"

/** Magnet field base classes
 *  (must be on top, as it is used in the MagnetBase)
 * **/

class FieldBase : public G4MagneticField {
    // Class that has the navigator and transform
public:
    FieldBase(G4ThreeVector centerPoint_in, G4LogicalVolume* fieldLV_in) :
        centerPoint(centerPoint_in), fieldLV(fieldLV_in) {};
    virtual void PostInitialize() {
        SetupTransform();
    }
private:
    G4ThreeVector centerPoint;
    G4LogicalVolume* fieldLV;
    static G4Navigator* fNavigator;
protected:
    void SetupTransform();
    G4AffineTransform fGlobalToLocal;

    G4double gradient; // [T/m]
};

/** Magnet/Object geometry base class **/

class MagnetBase {
public:
    static MagnetBase* MagnetFactory(G4String inputString, DetectorConstruction* detCon, G4String magnetName);

    MagnetBase(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
               std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
               G4String magnetName_in, G4String magnetType_in) :
        zPos(zPos_in), doRelPos(doRelPos_in), length(length_in),gradient(gradient_in),
        keyValPairs(keyValPairs_in), detCon(detCon_in),
        magnetName(magnetName_in), magnetType(magnetType_in) {};

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

    virtual G4double GetTypicalDensity() const = 0;

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
    G4LogicalVolume* MakeNewMainLV(G4String name_postfix, G4double width, G4double height);

    G4double mainLV_w = -1.0; // Width  of mainLV, in local (pre-rotation) coordinate axis [G4 units]
    G4double mainLV_h = -1.0; // Height of mainLV, in local (pre-rotation) coordinate axis [G4 units]

    G4Transform3D mainPV_transform;
    void BuildMainPV_transform();

    virtual void ConstructDetectorLV();
    G4LogicalVolume* detectorLV = NULL;
    G4VSensitiveDetector* magnetSD = NULL;
    //By default, the detectors will be attached to dedicated volumes
    // in the MagnetSensorWorld, which are copies of the mainLV.
    // This means all steps in the mainLV are recorded.
    // Set this flag to false in the constructor if one of the physically existing
    // logical volumes are used instead.
    G4bool useGhostDetector = true;

public:
    const G4String magnetName;
    const G4String magnetType;

    const G4Transform3D GetMainPV_transform() const {return mainPV_transform;};

    virtual void Construct() = 0;
    G4LogicalVolume* GetMainLV() const;
    G4LogicalVolume* GetDetectorLV() const;
    virtual void AddSD(); // Adds an SD to the detectorLV
    G4bool GetUseGhostDetector() {return useGhostDetector;}

public:
    //Parsing helpers
    void ParseOffsetRot(G4String k, G4String v);

    G4bool   ParseBool  (G4String inStr, G4String readWhat);
    G4double ParseDouble(G4String inStr, G4String readWhat);

    void PrintCommonParameters();
};

/** Various magnets or objects **/

// PLASMA LENS
class MagnetPLASMA1 : public MagnetBase {
public:
    MagnetPLASMA1(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                  std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                  G4String magnetName_in);

    virtual void Construct();

    virtual G4double GetTypicalDensity() const { return sapphireMaterial->GetDensity(); };
private:
    G4double plasmaTotalCurrent; // [A]
    G4double capRadius;          // [G4 length units]
    G4double cryWidth;           // [G4 length units]
    G4double cryHeight;          // [G4 length units]

    G4Material* sapphireMaterial = NULL;
};

class FieldPLASMA1 : public FieldBase {
public:
    FieldPLASMA1(G4double current_in, G4double radius_in,
                 G4ThreeVector centerPoint_in, G4LogicalVolume* fieldLV_in);
    virtual void GetFieldValue(const G4double point[4], G4double field[6]) const;

private:
    G4double plasmaTotalCurrent = 0.0;    // [A]
    G4double capRadius          = 0.5*mm; // [G4 units]

};

// CIRCULAR OPENING COLLIMATOR
class MagnetCOLLIMATOR1 : public MagnetBase {
public:
    MagnetCOLLIMATOR1(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                      std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                      G4String magnetName_in);
    virtual void Construct();

    virtual G4double GetTypicalDensity() const { return absorberMaterial->GetDensity(); };

private:
    G4String absorberMaterialName;
    G4Material* absorberMaterial = NULL;
    G4double width  = 10.0*mm;  //[G4 length units]
    G4double height = 50.0*mm;  //[G4 length units]
    G4double radius = 50.0*mm;  //[G4 length units]
};

// TARGET -- JUST A SIMPLE SLAB OF MATERIAL
class MagnetTARGET : public MagnetBase {
public:
    MagnetTARGET(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                      std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                      G4String magnetName_in);
    virtual void Construct();

    virtual G4double GetTypicalDensity() const { return targetMaterial->GetDensity(); };

private:
    G4String targetMaterialName;
    G4Material* targetMaterial = NULL;
    G4double width  = 10.0*mm;  //[G4 length units]
    G4double height = 10.0*mm;  //[G4 length units]
};

// TARGETR -- JUST A SIMPLE SLAB OF MATERIAL, BUT ROUND
class MagnetTARGETR : public MagnetBase {
public:
    MagnetTARGETR(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                    std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                    G4String magnetName_in);
    virtual void Construct();

    virtual G4double GetTypicalDensity() const { return targetMaterial->GetDensity(); };

private:
    G4String targetMaterialName;
    G4Material* targetMaterial = NULL;
    G4double radius = 10.0*mm;  //[G4 length units]
};

// COLLIMATORHV -- TWO SIMPLE SLABS OF MATERIALS, WITH EITHER
// A HORISONAL (HCOLL) OR VERTICAL (VCOLL) GAP
class MagnetCOLLIMATORHV : public MagnetBase {
public:
    MagnetCOLLIMATORHV(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                       std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                       G4String magnetName_in);
    virtual void Construct();

    virtual G4double GetTypicalDensity() const { return targetMaterial->GetDensity(); };

private:
    G4String    targetMaterialName;
    G4Material* targetMaterial = NULL;
    G4double    gap       = 10.0*mm; //[G4 length units]
    G4double    jawThick  = 50.0*mm; //[G4 length units]
    G4double    jawHeight = 50.0*mm; //[G4 length units]
    G4bool      isH       = true;    //Horizontal collimator -> vertical gap opening? If not, then vertical gap
};

// SHIELDEDSCINTILATOR
// A cylindrical scintillator with a hollow cylindrical shell
// This defines the detector so it is only sensitive in the scintillator
class MagnetSHIELDEDSCINTILLATOR : public MagnetBase {
public:
    MagnetSHIELDEDSCINTILLATOR(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                               std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                               G4String magnetName_in);
    virtual void Construct();

    virtual G4double GetTypicalDensity() const { return scintillatorMaterial->GetDensity(); };

protected:
    virtual void ConstructDetectorLV(); //Overide the general ConstructDetectorLV()

private:
    G4String    scintillatorMaterialName = "G4_SODIUM_IODIDE";
    G4Material* scintillatorMaterial = NULL;
    G4String    shieldingMaterialName = "G4_Pb";
    G4Material* shieldingMaterial = NULL;
    G4double    r_scint   = 25.0*mm; //[G4 length units] Radius of scintillating cyrstal
    G4double    l_scint   = 50.0*mm; //[G4 length units] Length of scintillating crystal
    G4double    z_scint   = 0.0*mm;  //[G4 length units] Longitudinal position of center of scintillating crystal, relative to center of shield
    G4double    ri_shield = 30.0*mm; //[G4 length units] Inner radius of shield
    G4double    ro_shield = 50.0*mm; //[G4 length units] Outer radius of shield

    G4LogicalVolume* scintillatorLV = NULL;
    G4LogicalVolume* shieldingLV    = NULL;
};

// PBW -- Replica of the ESS PBW
// A section of a cylindrical shell from startPhi to (startPhi + arcPhi) of a thickness
// Includes a water channel in the middle from (startPhi * 2) to (startPhi * 2 + arcPhi * 0.5)
class MagnetPBW : public MagnetBase {
public:
    MagnetPBW(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                    std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                    G4String magnetName_in);
    virtual void Construct();

    virtual G4double GetTypicalDensity() const { return targetMaterial->GetDensity(); };

private:
    G4String targetMaterialName = "G4_Al";
    G4Material* targetMaterial;
    G4double radius        = 88.0*mm;   //[G4 length units] Inner radius of cylinder, >0
    G4double al1Thick      = 1.0*mm;    //[G4 length units] Outer thickness of metal window, >0
    G4double waterThick    = 2.0*mm;    //[G4 length units] Thickness of water channel, >0
    G4double al2Thick      = 1.25*mm;   //[G4 length units] Inner thickness of metal window, >0
    G4double thickness     = 4.25*mm;   //[G4 length units] Total thicnkess
    G4double width         = 200.0*mm;  //[G4 length units] Width of cylinder as seen by PBW, >0
    G4double height        = 160.0*mm;  //[G4 length units] Height of window section, >0
    G4double arcPhi        = 120.0*deg; //[deg] Arc angle of window section, within: 0 <= arcPhi <= 180
    G4double startPhi      = 30*deg;    //[deg] Start angle of window section, = 90 - arcPhi/2 for geometry
    G4double waterStartPhi = 60*deg;    //[deg] Start angle of water channel, = 90 - arcPhi/4 for geometry
    G4double boxCenter     = 68.125*mm; //[G4 length units] Distance to translate PBW so center lines up with pos
    G4double rightAng      = 90*deg;    //[deg] Reference for calculations
};

// Rectangular Opening Collimator
// Default dimensions are for ESS PBIP with minimum absorber volume
class MagnetCOLLIMATORRECT : public MagnetBase {
public:
    MagnetCOLLIMATORRECT(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                      std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                      G4String magnetName_in);
    virtual void Construct();

    virtual G4double GetTypicalDensity() const { return absorberMaterial->GetDensity(); };

private:
    G4String absorberMaterialName;
    G4Material* absorberMaterial = NULL;
    G4double absorberWidth  = 250.0*mm; //[G4 length units]
    G4double absorberHeight = 100.0*mm; //[G4 length units]
    G4double holeWidth      = 200.0*mm; //[G4 length units]
    G4double holeHeight     = 80.0*mm;  //[G4 length units]
};

#endif
