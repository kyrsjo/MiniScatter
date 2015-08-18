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

  TFile *histFile;

  TH2D* energyHisto;
  TH1D* pixelHisto;
  

private: 
  analysis(){};
  //! Singleton static instance
  static analysis* singleton;

  TH1D* targetEdep;
  TH1D* targetEdep_NIEL;
  TH1D* targetEdep_IEL;
  
  TH1D* tracker_numParticles;
  TH1D* tracker_angle;
  TH1D* tracker_energy;
  std::map<G4int,G4int> tracker_particleTypes;
  
  //TH1D* tracker_energyflux_neutrals;
  
  TH1D* tracker_protonAngle;
  TH1D* tracker_protonEnergy;
};

  
#endif
