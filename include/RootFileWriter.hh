/*
 * This file is part of MiniScatter.
 *
 *  MiniScatter is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MiniScatter is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MiniScatter.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef ROOTFILEWRITER_HH
#define ROOTFILEWRITER_HH 1
#include "G4Event.hh"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include <map>

class TRandom;

// Use a simple struct for writing to ROOT file,
// since this requires no dictionary to read.
struct trackerHitStruct {
    Double_t x; // [mm]
    Double_t y; // [mm]
    Double_t z; // [mm]

    Double_t px; // [MeV/c]
    Double_t py; // [MeV/c]
    Double_t pz; // [MeV/c]

    Double_t E; // [MeV]

    Int_t PDG;
    Int_t charge;

    Int_t eventID;
};

class particleTypesCounter {
public:
    particleTypesCounter(){
        particleTypes.clear();
        particleNames.clear();
        numParticles = 0;
    }
    // In both cases, the index is the PDG id.
    std::map<G4int,G4int> particleTypes;    // The number of particles of each type
    std::map<G4int,G4String> particleNames; // The name of each particle type
    G4int numParticles;
};

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
    void setFoldername(G4String foldername_out_arg) {
        this->foldername_out = foldername_out_arg;
    };

    void setQuickmode(G4bool quickmode_arg) {
        this->quickmode = quickmode_arg;
    };
    void setanaScatterTest(G4bool anaScatterTest_arg) {
        this->anaScatterTest = anaScatterTest_arg;
    };
    void setMiniFile(G4bool miniFile_arg) {
        this->miniFile = miniFile_arg;
    };

    void setBeamEnergyCutoff(G4double cutFrac){
        this->beamEnergy_cutoff = cutFrac;
    }
    void setPositionCutoffR(G4double cutR) {
        this->position_cutoffR = cutR;
    }

    void setNumEvents(G4int numEvents_in) {
        this->numEvents = numEvents_in;
    }

    void setEdepDensDZ(G4double edep_dens_dz_in) {
        this->edep_dens_dz = edep_dens_dz_in;
    }
    void setEngNbins(G4int edepNbins_in);

    void setRNGseed(G4int rngSeed_in) {
        this->rngSeed = rngSeed_in;
    }

    void setHistPosLim(G4double histPos_in) {
        this->phasespacehist_posLim = histPos_in;
    }

    void setHistAngLim(G4double histAngLim_in) {
        this->phasespacehist_angLim = histAngLim_in;
    }

private:
    RootFileWriter(){
        has_filename_out = false;
    };

    //! Singleton static instance
    static RootFileWriter* singleton;

    //The ROOT file
    TFile *histFile                                                             = NULL;

    // TTrees //
    TTree* targetExit                                                           = NULL;
    TTree* trackerHits                                                          = NULL;
    TTree* initParts                                                            = NULL;
    trackerHitStruct targetExitBuffer;
    trackerHitStruct trackerHitsBuffer;
    trackerHitStruct initPartsBuffer;

    Double_t* magnetEdepsBuffer                                                 = NULL;
    TTree* magnetEdeps                                                          = NULL;

    // Histograms //

    // Target histograms
    TH1D* targetEdep                                                            = NULL;
    TH1D* targetEdep_NIEL                                                       = NULL;
    TH1D* targetEdep_IEL                                                        = NULL;

    std::map<G4int,TH1D*> target_exit_energy;
    std::map<G4int,TH1D*> target_exit_cutoff_energy;

    TH1D* target_exitangle_hist                                                 = NULL;
    TH1D* target_exitangle_hist_cutoff                                          = NULL;

    TH2D* target_exit_phasespaceX                                               = NULL;
    TH2D* target_exit_phasespaceY                                               = NULL;
    TH2D* target_exit_phasespaceX_cutoff                                        = NULL;
    TH2D* target_exit_phasespaceY_cutoff                                        = NULL;
    TH2D* target_exit_phasespaceXY                                              = NULL;
    TH2D* target_exit_phasespaceXY_cutoff                                       = NULL;
  
    std::map<G4int,TH1D*> target_exit_Rpos;
    std::map<G4int,TH1D*> target_exit_Rpos_cutoff;

    TH3D* target_edep_dens                                                      = NULL;
    TH2D* target_edep_rdens                                                     = NULL;

    // Magnet histograms
    std::vector<TH1D*> magnet_edep;
    std::vector<TH3D*> magnet_edep_dens;
    std::vector<TH2D*> magnet_edep_rdens;

    std::vector<std::map<G4int,TH1D*>> magnet_exit_Rpos;
    std::vector<std::map<G4int,TH1D*>> magnet_exit_Rpos_cutoff;
    std::vector<TH2D*> magnet_exit_phasespaceX;
    std::vector<TH2D*> magnet_exit_phasespaceY;
    std::vector<TH2D*> magnet_exit_phasespaceX_cutoff;
    std::vector<TH2D*> magnet_exit_phasespaceY_cutoff;
    std::vector<std::map<G4int,TH2D*>> magnet_exit_phasespaceX_cutoff_PDG;
    std::vector<std::map<G4int,TH2D*>> magnet_exit_phasespaceY_cutoff_PDG;
    std::vector<std::map<G4int,TH1D*>> magnet_exit_energy;
    std::vector<std::map<G4int,TH1D*>> magnet_exit_cutoff_energy;

    //Tracker histograms
    std::vector<TH1D*> tracker_numParticles;
    std::vector<TH1D*> tracker_energy;
    std::vector<std::map<G4int,TH1D*>> tracker_type_energy;
    std::vector<std::map<G4int,TH1D*>> tracker_type_cutoff_energy;
    std::vector<TH2D*> tracker_phasespaceX;
    std::vector<TH2D*> tracker_phasespaceY;
    std::vector<TH2D*> tracker_phasespaceX_cutoff;
    std::vector<TH2D*> tracker_phasespaceY_cutoff;
    std::vector<TH2D*> tracker_phasespaceXY;
    std::vector<TH2D*> tracker_phasespaceXY_cutoff;
    std::vector<std::map<G4int,TH2D*>> tracker_phasespaceX_cutoff_PDG;
    std::vector<std::map<G4int,TH2D*>> tracker_phasespaceY_cutoff_PDG;
    std::vector<std::map<G4int,TH2D*>> tracker_phasespaceXY_cutoff_PDG;

    std::vector<std::map<G4int,TH1D*>> tracker_Rpos;
    std::vector<std::map<G4int,TH1D*>> tracker_Rpos_cutoff;

    //Initial distribution
    TH2D* init_phasespaceX                                                      = NULL;
    TH2D* init_phasespaceY                                                      = NULL;
    TH2D* init_phasespaceXY                                                     = NULL;
    TH1D* init_E                                                                = NULL;

    // End-of-run statistics

    // Count the number of each particle type that hits the tracker
    std::map<G4String,particleTypesCounter> typeCounter;

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

    G4String foldername_out = "plots";

    G4String rootFileName = "!ROOTFILENAME_NOT_INITIALIZED!";

    G4bool quickmode      = false;
    G4bool anaScatterTest = false;
    G4bool miniFile       = false;

    G4double beamEnergy; // [MeV]

    // Compute statistics for charged particles with energy > this cutoff
    G4double beamEnergy_cutoff = 0.95;
    // Compute statistics for charged particles with position inside this radius
    G4double position_cutoffR = 1.0; // [mm]

    //Set default limits for histograms
    G4double phasespacehist_posLim = 10;
    G4double phasespacehist_angLim = 5;

    //Delta z for the energy deposition density radial TH2Ds and TH3Ds [mm, 0 => Disable; < 0 => TH2Ds only]
    G4double edep_dens_dz = 0.0;

    //Energy deposition and energy-remaining number of bins
    G4int engNbins = 1000;

    // RNG for sampling over the step
    G4int rngSeed;
    TRandom* RNG                                                                = NULL;

    Int_t eventCounter; // Used for EventID-ing and metadata
    Int_t numEvents;    // Used for comparing to eventCounter with metadata;
                        // only reflects the -n <int> command line flag
                        // so it may be 0 if this was not set.

    void PrintTwissParameters(TH2D* phaseSpaceHist);
    void PrintParticleTypes(particleTypesCounter& pt, G4String name);
    void FillParticleTypes(particleTypesCounter& pt, G4int PDG, G4String type);
};

#endif
