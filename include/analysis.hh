#ifndef ANALYSIS_HH_
#define ANALYSIS_HH_
#include "G4Event.hh"
#include "G4Run.hh"
#include "DetectorConstruction.hh"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;
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
  TFile *file;
  TFile *histFile;
  vector<TH2D*> pixelMap;
  vector<TH2D*> particleMap;
  int numberOfFrames;
  ofstream outfile;
  map <string,int> data;
  map <vector<int>,double> largestEnergy;
  int protons;
  int alpha; 
  int gamma; 
  int electrons;
  int pions;
  int kaons;
  int others;
  //TH2D* energyHisto;
  TH2D* energyHisto;
  TH1D* pixelHisto;
  

private: 
  analysis();
  //! Singleton static instance
  static analysis* singleton;
  double halfDetX;
  double halfDetY;
  double halfDetZ;
  int xpixels,ypixels;
  double xpitch_pixel;
  double ypitch_pixel;
};

  
#endif
