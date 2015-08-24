#ifndef ANALYSIS_HH_
#define ANALYSIS_HH_
#include "G4Event.hh"
#include "G4Run.hh"
#include "DetectorConstruction.hh"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <map>

class analysis {
public:
  //! Singleton pattern
  static analysis* GetInstance() {
    if ( analysis::singleton == NULL ) analysis::singleton = new analysis();
    return analysis::singleton;
  }

  void makeHistograms();
  void fillHistograms(const G4Event *anEvent);
  void writeHistograms();
  void writePerEvent(const G4Event* event);

private: 
  analysis(){};
  //! Singleton static instance
  static analysis* singleton;

  //The histogram file
  TFile *histFile;
  
  //Histograms
  TH1D* targetEdep;
  TH1D* targetEdep_NIEL;
  TH1D* targetEdep_IEL;

  TH1D* target_sumMomentum_z;

  TH1D* tracker_numParticles;
  TH1D* tracker_angle;
  TH1D* tracker_energy;
  
  std::map<G4int,G4int> tracker_particleTypes;
  std::map<G4int,G4String> tracker_particleNames;
  G4int numParticles_total;
  
  TH1D* tracker_protonAngle;
  TH1D* tracker_protonEnergy;
  
  TH1D* tracker_sumMomentum_z;
  
  TH1D* tracker_energyAngle;
  TH1D* tracker_energyAngle_charged;
  TH1D* tracker_energyAngle_neutral;
  
  //Used to generate filenames
  G4String physListName;
  // (Target thickness extracted directly from the DetectorConstruction)
public:
  
  void SetMetadata(const G4String physListName_in);
};

  
#endif
