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

    G4double GetLength()    const { return length; };

protected:
    G4double zPos;     // [G4 units]
    G4bool   doRelPos;
    G4double length;   // [G4units]
    G4double gradient; // [T/m]

    std::map<G4String,G4String> keyValPairs;
    DetectorConstruction* detCon;

    FieldBase* field = NULL;

    G4LogicalVolume* mainLV = NULL;
    G4LogicalVolume* MakeNewMainLV(G4String name_postfix);

    virtual void ConstructDetectorLV();
    G4LogicalVolume* detectorLV = NULL;
    G4VSensitiveDetector* magnetSD = NULL;

public:
    const G4String magnetName;

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
    G4double height = 50*mm;    //[G4 length units]
    G4double radius = 50*mm;    //[G4 length units]
};

#endif
