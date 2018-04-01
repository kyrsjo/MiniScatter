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

    //Tracker histograms
    TH1D* tracker_numParticles;
    TH1D* tracker_energy;

    // Count the number of each particle type that hits the tracker
    std::map<G4int,G4int> tracker_particleTypes;
    std::map<G4int,G4String> tracker_particleNames;
    G4int numParticles_total;

    //Output file naming
    G4String filename_out;
    G4bool has_filename_out;
    static const G4String foldername_out;

};

#endif
