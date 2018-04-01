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
    void writeHistograms();
    void writePerEvent(const G4Event* event);

    void setFilename(G4String filename_out_arg) {
        this->filename_out = filename_out_arg;
        has_filename_out = true;
    };
private:
    analysis(){
        has_filename_out = false;
    };

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

    //Output file naming
    G4String filename_out;
    G4bool has_filename_out;
    static const G4String foldername_out;

};

#endif
