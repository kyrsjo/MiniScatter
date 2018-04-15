#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class DetectorConstruction;

// -----------------------------------------------------------------------------------------

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction(DetectorConstruction* DC,
                           G4double beam_energy_in,
                           G4String beam_type_in,
                           G4double beam_offset_in);
    virtual ~PrimaryGeneratorAction();
    void GeneratePrimaries(G4Event*);

    G4double get_beam_energy() const { return beam_energy; };
    G4double get_beam_particlemass() const { return particle->GetPDGMass(); };
    G4double get_beam_particlecharge() const { return particle->GetPDGCharge(); };
private:
    G4ParticleGun*           particleGun;  //pointer a to G4  class
    DetectorConstruction*    Detector;     //pointer to the geometry

    G4double beam_energy;    // Beam energy [MeV]
    G4String beam_type;      // Beam particle type
    G4double beam_offset;    // Beam offset (x) [mm]

    G4ParticleDefinition* particle; // Particle type

};

// -----------------------------------------------------------------------------------------

#endif
