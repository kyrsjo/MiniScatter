#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4UniformMagField.hh"

//--------------------------------------------------------------------------------

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    DetectorConstruction(G4double TargetThickness_in  = 1.0,
                         G4String TargetMaterial_in   = "G4_Cu",
                         G4double DetectorDistance_in = 500.0,
                         G4double DetectorAngle_in    = 0.0,
                         G4bool   DetectorRotated_in  = false
                         );
    ~DetectorConstruction(){};

public:
    void SetTargetMaterial (G4String);
    //  void SetDetectorMaterial (G4String);

    void SetMagField(G4double);
    G4VPhysicalVolume* Construct();
public:

    const G4VPhysicalVolume* getphysiWorld() {return physiWorld;};
    const G4VPhysicalVolume* getTargetPV()   {return physiTarget;};

    inline G4double getTargetThickness()   {return TargetThickness;};
    inline G4double getTargetSizeX() const {return TargetSizeX;};
    inline G4double getTargetSizeY() const {return TargetSizeY;};

    inline G4double getDetectorDistance()  {return DetectorDistance;};
    inline G4double getDetectorSizeX() const {return TargetSizeX;};
    inline G4double getDetectorSizeY() const {return TargetSizeY;};

private:
    G4Material*        vacuumMaterial;
    G4Material*        AlMaterial;
    G4Material*        CMaterial;
    G4Material*        CuMaterial;
    G4Material*        PbMaterial;
    G4Material*        TiMaterial;
    G4Material*        StainlessSteel;

    G4double           WorldSizeX;
    G4double           WorldSizeY;
    G4double           WorldSizeZ;

    G4double           TargetSizeX;
    G4double           TargetSizeY;
    G4double           TargetThickness;

    G4Material*        TargetMaterial;

    G4double           DetectorSizeX;
    G4double           DetectorSizeY;
    // Distance from center of target to center of detector:
    G4double           DetectorThickness;
    G4double           DetectorAngle; //[rad]
    G4bool             DetectorRotated;

    G4double           DetectorDistance;
    G4Material*        DetectorMaterial;

    G4Box*             solidWorld;    //pointer to the solid World
    G4LogicalVolume*   logicWorld;    //pointer to the logical World
    G4VPhysicalVolume* physiWorld;    //pointer to the physical World

    G4Box*             solidTarget; //pointer to the solid target
    G4LogicalVolume*   logicTarget; //pointer to the logical target
    G4VPhysicalVolume* physiTarget; //pointer to the physical target

    G4Box*             solidDetector; //pointer to the solid detector
    G4LogicalVolume*   logicDetector; //pointer to the logical detector
    G4VPhysicalVolume* physiDetector; //pointer to the physical detector

    G4UniformMagField* magField;      //pointer to the magnetic field

private:
    void DefineMaterials();
};

//--------------------------------------------------------------------------------

#endif
