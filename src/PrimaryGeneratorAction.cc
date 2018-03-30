#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLatex.h"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <cmath>
//#include "CLHEP/Random/RandGauss.h"
//#include "/Applications/geant4.9.6.p02/source/externals/clhep/include/CLHEP/Random/RandLandau.h" //(include this before method definition)

// -----------------------------------------------------------------------------------------

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* DC,
                                               G4double beam_energy_in,
                                               G4String beam_type_in,
                                               G4double beam_offset_in) :
    Detector(DC), beam_energy(beam_energy_in), beam_type(beam_type_in), beam_offset(beam_offset_in){
    G4int n_particle = 1;
    particleGun  = new G4ParticleGun(n_particle);
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle(beam_type);
    if (particle == NULL) {
        G4cerr << "Error - particle named '" << beam_type << "'not found" << G4endl;
        //particleTable->DumpTable();
        exit(1);
    }
    particleGun->SetParticleDefinition(particle);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    particleGun->SetParticlePosition(G4ThreeVector(beam_offset*mm,0,-30*cm));
    particleGun->SetParticleEnergy(beam_energy*TeV);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
    particleGun->GeneratePrimaryVertex(anEvent);
}
