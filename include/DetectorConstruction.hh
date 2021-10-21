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
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"

#include <vector>

// Forward declarations
class MagnetBase; 
//--------------------------------------------------------------------------------

class DetectorConstruction : public G4VUserDetectorConstruction {

public:
    DetectorConstruction(G4double TargetThickness_in,
                         G4String TargetMaterial_in,
                         G4double TargetAngle_in,
                         G4String BackgroundMaterial_in,
                         G4double WorldSize_in,
                         G4double WorldMinLength_in,
                         std::vector <G4String> &magnetDefinitions_in);
    ~DetectorConstruction(){};

private:
    void SetTargetMaterial (G4String);
    void SetBackgroundMaterial (G4String);
    //  void SetDetectorMaterial (G4String);
public:

    G4bool   GetHasTarget() {return HasTarget;};
    G4int    GetTargetMaterialZ();
    G4double GetTargetMaterialA();
    G4double GetTargetMaterialDensity();

    G4Material* GetTargetMaterial();
    G4Material* GetBackgroundMaterial();

    //    void SetMagField(G4double);
    G4VPhysicalVolume* Construct();
    void PostInitialize(); // To be called after construct, but before tracking starts
public:

    const G4VPhysicalVolume* getphysiWorld() {return physiWorld;};
    const G4VPhysicalVolume* getTargetPV()   {return physiTarget;};

    //All in G4 units
    inline G4double getTargetThickness() const {return TargetThickness;};
    inline G4double getTargetSizeX()     const {return TargetSizeX;};
    inline G4double getTargetSizeY()     const {return TargetSizeY;};

    inline G4double getWorldSizeZ()       const {return WorldSizeZ;};
    inline G4double getWorldSizeX()       const {return WorldSizeX;};
    inline G4double getWorldSizeY()       const {return WorldSizeY;};

private:
    G4Material*        vacuumMaterial = NULL;
    G4Material*        airMaterial    = NULL;

    G4Material*        AlMaterial     = NULL;
    G4Material*        CMaterial      = NULL;
    G4Material*        CuMaterial     = NULL;
    G4Material*        PbMaterial     = NULL;
    G4Material*        TiMaterial     = NULL;
    G4Material*        SiMaterial     = NULL;
    G4Material*        WMaterial      = NULL;
    G4Material*        UMaterial      = NULL;
    G4Material*        FeMaterial     = NULL;

    G4Material*        MylarMaterial          = NULL;
    G4Material*        KaptonMaterial         = NULL;
    G4Material*        StainlessSteelMaterial = NULL;
    G4Material*        WaterMaterial          = NULL;

    G4Material*        SapphireMaterial = NULL;

    G4Material*        ChromoxMaterial       = NULL;
    G4Material*        ChromoxScreenMaterial = NULL;

    G4Material*        gasH_2         = NULL;
    G4Material*        gasHe          = NULL;
    G4Material*        gasN_2         = NULL;
    G4Material*        gasNe          = NULL;
    G4Material*        gasAr          = NULL;

    //These are all in G4 units
    G4double           WorldSizeX;
    G4double           WorldSizeY;
    G4double           WorldSizeZ;

    G4double           TargetSizeX;
    G4double           TargetSizeY;
    G4double           TargetThickness;

    G4double           TargetAngle;

    G4bool             HasTarget      = false;
    G4Material*        TargetMaterial = NULL;

    G4Material*        BackgroundMaterial = NULL;

    G4double           DetectorDistance;
    G4Material*        DetectorMaterial;

    G4Box*             solidWorld;    //pointer to the solid World
    G4LogicalVolume*   logicWorld;    //pointer to the logical World
    G4VPhysicalVolume* physiWorld;    //pointer to the physical World

    G4Box*             solidTarget; //pointer to the solid target
    G4LogicalVolume*   logicTarget; //pointer to the logical target
    G4VPhysicalVolume* physiTarget; //pointer to the physical target
private:
    void DefineMaterials();
    G4Material* DefineGas(G4String TargetMaterial_in);

private:
    std::vector <G4String> &magnetDefinitions;
public:
    // This one needs to be accessed by e.g. the rootFileWriter
    std::vector <MagnetBase*> magnets;
private:
    std::vector <G4VPhysicalVolume*> magnetPVs;

};

//--------------------------------------------------------------------------------

#endif
