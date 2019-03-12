#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

//System of units defines variables like "s" and "m" in the global scope,
// which are then shadowed inside functions in the header. Let's ignore it.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "TMatrixDfwd.h"
#include "TDecompChol.h"
#pragma GCC diagnostic pop

#include "TRandom.h"

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
                           G4double beam_offset_in,
                           G4double beam_zpos_in,
                           G4bool   doBacktrack_in,
                           G4String covarianceString_in,
                           G4double Rcut_in );
    virtual ~PrimaryGeneratorAction();
    void GeneratePrimaries(G4Event*);

    G4double get_beam_energy()         const { return beam_energy; };
    G4double get_beam_particlemass()   const { return particle->GetPDGMass(); };
    G4double get_beam_particlecharge() const { return particle->GetPDGCharge(); };
private:
    G4ParticleGun*           particleGun;  //pointer a to G4  class
    DetectorConstruction*    Detector;     //pointer to the geometry

    G4double beam_energy;    // Beam energy [MeV]
    G4String beam_type;      // Beam particle type
    G4double beam_offset;    // Beam offset (x) [mm]
    G4double beam_zpos;      // Beam initial z position [converted to G4 units in constructor]
    G4bool   doBacktrack;    // Generate at z=0 then backtrack to injection position?

    G4ParticleDefinition* particle; // Particle type

    // Setup for covariance
    G4bool hasCovariance = false;
    void setupCovariance();
    G4double convertColons(str_size startPos, str_size endPos, G4String paramName);

    G4String covarianceString; // String defining the covariance matrix via Twiss parameters
    G4double epsN_x;  // Normalized emittance  (x) [um]
    G4double epsG_x;  // Geometrical emittance (x) [um]
    G4double beta_x;  // Beta function         (x) [m]
    G4double alpha_x; // Alpha function        (x) [-]

    G4double epsN_y;  // Normalized emittance  (y) [um]
    G4double epsG_y;  // Geometrical emittance (y) [um]
    G4double beta_y;  // Beta function         (y) [m]
    G4double alpha_y; // Alpha function        (y) [-]

    // Beam covariance matrices [m,rad]
    TMatrixD covarX;
    TMatrixD covarY;

    // covar = covar_L*covar_L^T, which is used for generation.
    // Note that ROOT's cholesky routine return the upper-triangular covar_LU^T,
    // so a transpose is neccessary.
    TMatrixD covarX_L;
    TMatrixD covarY_L;

    //Setup for circular uniform distribution / Rcut
    G4double Rcut; // [mm]

    TRandom* RNG;

public:
    //Leave the generated positions where RootFileWriter can pick it up [G4 units]
    G4double x,xp,y,yp;
};

// -----------------------------------------------------------------------------------------

#endif
