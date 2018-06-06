#ifndef ROOTFILEWRITER_HH_
#define ROOTFILEWRITERANALYSIS_HH_
#include "G4Event.hh"
#include "G4Run.hh"
#include "DetectorConstruction.hh"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <map>

class RootFileWriter {
public:
    //! Singleton pattern
    static RootFileWriter* GetInstance() {
        if ( RootFileWriter::singleton == NULL ) RootFileWriter::singleton = new RootFileWriter();
        return RootFileWriter::singleton;
    }

    void initializeRootFile();
    void finalizeRootFile();
    void doEvent(const G4Event* event);

    void setFilename(G4String filename_out_arg) {
        this->filename_out = filename_out_arg;
        has_filename_out = true;
    };
    void setQuickmode(G4bool quickmode_arg) {
        this->quickmode = quickmode_arg;
    };

private:
    RootFileWriter(){
        has_filename_out = false;
    };

    //! Singleton static instance
    static RootFileWriter* singleton;

    //The ROOT file
    TFile *histFile;

    //Histograms

    // Target histograms
    TH1D* targetEdep;
    TH1D* targetEdep_NIEL;
    TH1D* targetEdep_IEL;

    TH1D* target_exitangle_hist;
    TH1D* target_exitangle_hist_cutoff;

    //Tracker histograms
    TH1D* tracker_numParticles;
    TH1D* tracker_energy;
    TH2D* tracker_hitPos;
    TH2D* tracker_hitPos_cutoff;
    TH2D* tracker_phasespaceX_cutoff;
    TH2D* tracker_phasespaceY_cutoff;

    //Initial distribution
    TH2D* init_phasespaceX;
    TH2D* init_phasespaceY;

    // End-of-run statistics

    // Count the number of each particle type that hits the tracker
    std::map<G4int,G4int> tracker_particleTypes;
    std::map<G4int,G4String> tracker_particleNames;
    G4int numParticles_total;

    // Compute means and standard deviations of where the particles hit the tracker
    G4double tracker_particleHit_x;
    G4double tracker_particleHit_xx;
    G4double tracker_particleHit_y;
    G4double tracker_particleHit_yy;

    G4double tracker_particleHit_x_cutoff;
    G4double tracker_particleHit_xx_cutoff;
    G4double tracker_particleHit_y_cutoff;
    G4double tracker_particleHit_yy_cutoff;

    G4int numParticles_cutoff;

    //Target exit angle RMS
    G4double target_exitangle;
    G4double target_exitangle2;
    G4int target_exitangle_numparticles;
    G4double target_exitangle_cutoff;
    G4double target_exitangle2_cutoff;
    G4int target_exitangle_cutoff_numparticles;

    // Internal stuff
    //Output file naming
    G4String filename_out;
    G4bool has_filename_out;
    static const G4String foldername_out;

    G4bool quickmode = false;

    G4double beamEnergy; // [MeV]
    // Compute statistics for charged particles with energy > this cutoff
    static constexpr G4double beamEnergy_cutoff = 0.95;

    void PrintTwissParameters(TH2D* phaseSpaceHist);
};

#endif
