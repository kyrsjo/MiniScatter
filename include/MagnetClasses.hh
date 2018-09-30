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
            return detCon->getTargetThickness()/2.0 + zPos + length/2.0 + lengthPad/2.0;
        }
        else {
            return zPos;
        }
    };

    virtual void PostInitialize() {
        if (field == NULL) {
            G4cerr << "Internal error: field==NULL!" << G4endl;
            exit(1);
        }
        field->PostInitialize();
    }

    G4double GetLength()    const { return length; };
    G4double GetLengthPad() const { return lengthPad; };

protected:
    G4double zPos;     // [G4 units]
    G4bool   doRelPos;
    G4double length;   // [G4units]
    G4double lengthPad = 1e-7*mm; // [G4units]
    G4double gradient; // [T/m]

    std::map<G4String,G4String> keyValPairs;
    DetectorConstruction* detCon;

    FieldBase* field = NULL;

    void AddSD(G4LogicalVolume* mainLV);
    G4VSensitiveDetector* magnetSD = NULL;

public:
    const G4String magnetName;

    virtual G4LogicalVolume* Construct() = 0;
};

class MagnetPLASMA1 : public MagnetBase {
public:
    MagnetPLASMA1(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                  std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in,
                  G4String magnetName_in);

    virtual G4LogicalVolume* Construct();
private:
    G4double plasmaTotalCurrent; // [A]
    G4double capRadius;          // [G4 lenght units]
    G4double cryWidth;           // [G4 length units]
    G4double cryHeight;          // [G4 length units]
};

#endif
