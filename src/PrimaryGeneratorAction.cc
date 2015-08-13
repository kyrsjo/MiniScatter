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

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* DC)
  :Detector(DC)
{
  G4String particleName = "anti_proton";
  //G4String particleName = "geantino";
  
  //G4String particleName = "pi-";
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName);		
  particleGun->SetParticleDefinition(particle);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{ 
  delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
  //set energy and stuff
  particleGun->SetParticlePosition(G4ThreeVector(-2*cm,0,0));	
  particleGun->SetParticleEnergy(0.5*MeV);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1,0,0));
  particleGun->GeneratePrimaryVertex(anEvent);
}

