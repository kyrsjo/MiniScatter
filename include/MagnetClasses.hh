#ifndef MagnetClasses_h$
#define MagnetClasses_h$

#include "G4SystemOfUnits.hh"

#include "DetectorConstruction.hh"
#include "FieldClasses.hh"

class MagnetBase {
public:
    static MagnetBase* MagnetFactory(G4String inputString, DetectorConstruction* detCon);

    MagnetBase(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
               std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in) :
        zPos(zPos_in), doRelPos(doRelPos_in), length(length_in),gradient(gradient_in),
        keyValPairs(keyValPairs_in), detCon(detCon_in) {};

    virtual G4double getZ0();
    virtual void PostInitialize() {
        if (field == NULL) {
            G4cerr << "Internal error: field==NULL!" << G4endl;
            exit(1);
        }
        field->PostInitialize();
    }
protected:
    G4double zPos;     // [G4 units]
    G4bool   doRelPos;
    G4double length;   // [G4units]
    G4double gradient; // [T/m]
    std::map<G4String,G4String> keyValPairs;
    DetectorConstruction* detCon;
    FieldBase* field = NULL;
public:
    virtual G4LogicalVolume* Construct(G4String magnetName) = 0;
};

class MagnetPLASMA1 : public MagnetBase {
public:
    MagnetPLASMA1(G4double zPos_in, G4bool doRelPos_in, G4double length_in, G4double gradient_in,
                  std::map<G4String,G4String> &keyValPairs_in, DetectorConstruction* detCon_in);

    virtual G4LogicalVolume* Construct(G4String magnetName);
private:
    G4double plasmaTotalCurrent; // [A]
    G4double capRadius;          // [G4 units]
};

#endif