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

#include "RootFileWriter.hh"

#include "TGraph.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TBranch.h"

#include "TRandom1.h"

#include "MyEdepHit.hh"
#include "MyTrackerHit.hh"

#include "G4SDManager.hh"

#include "G4Track.hh"
#include "G4RunManager.hh"

#include "DetectorConstruction.hh"
#include "MagnetClasses.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4SystemOfUnits.hh"

#include <iostream>
#include <iomanip>

#include <unistd.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
struct stat stat_info;

using namespace std;
RootFileWriter* RootFileWriter::singleton = 0;

const G4double RootFileWriter::phasespacehist_posLim = 10.0*mm;
const G4double RootFileWriter::phasespacehist_angLim = 5.0*deg;

void RootFileWriter::initializeRootFile(){
    G4RunManager*           run    = G4RunManager::GetRunManager();
    DetectorConstruction*   detCon = (DetectorConstruction*)run->GetUserDetectorConstruction();
    PrimaryGeneratorAction* genAct = (PrimaryGeneratorAction*)run->GetUserPrimaryGeneratorAction();
    this->beamEnergy = genAct->get_beam_energy();

    //Count all particles that are Fill'ed for the stats used to compute the twiss parameters,
    // even if they are outside the phasespacehist_posLim / phaspacehist_angLim.
    TH2D::StatOverflows(true);

    if (not has_filename_out) {
        G4cerr << "Error: filename_out not set." << G4endl;
        exit(1);
    }
    G4String rootFileName = foldername_out + "/" + filename_out + ".root";
    G4cout << "foldername = '" << foldername_out << "'" << G4endl;

    //Create folder if it does not exist
    // Note: This code assumes UNIX filesystem conventions.
    // Linux, Mac, etc. should be OK. Windows is NOT supported, and it may eat all your cheese.
    // On a cheese-eating platform, replace the recusive folder-creation logic below
    // (everything inside the if (errno==ENOENT) { .... }) with
    // mkdir(foldername_out.c_str(), 0755);
    // This will not do recursive folder creation, but it should work.
    // Note: no guarantees that some other part of the code won't drink all your beer, this is untested...
    #if not ( defined(unix) || defined(__unix__) || defined(__unix) )
    #error Only UNIX is supported.
    #endif
    if(stat(foldername_out.c_str(), &stat_info) != 0) {
        if (errno == ENOENT) {
            G4cout << "Creating folder '" << foldername_out << "'" << G4endl;

            G4String foldername_out_full;
            if (foldername_out.c_str()[0] != '/') {

                char* cwd = get_current_dir_name();
                if (cwd == NULL) {
                    perror("Error getting the current path");
                    exit(1);
                }

                if (cwd[strlen(cwd)-1] == '/') {
                    foldername_out_full = G4String(cwd) + foldername_out;
                }
                else {
                    foldername_out_full = G4String(cwd) + G4String('/') + foldername_out;
                }

                delete[] cwd;

                G4cout << "Converted relative path '" << foldername_out << "' to absolute path '" 
                       << foldername_out_full << "'." << G4endl;
            }
            else {
                foldername_out_full = foldername_out;
            }
            const char* path_full = foldername_out_full.c_str();

            if( strlen(path_full) < 2 ) {
                G4cerr << "ERROR: Path string must be more than '/' (and why is '/' not existing?), got '"
                       << path_full << "' - aborting!" << G4endl;
                exit(1);
            }

            //Create the folders recursively
            size_t idx_stop;
            for (idx_stop = 1; idx_stop < strlen(path_full); ++idx_stop) {
                if (path_full[idx_stop] == '/') { // Got it!

                    //Extract the path up to here
                    char* path_part = new char[idx_stop+1];
                    strncpy(path_part, path_full, idx_stop);
                    path_part[idx_stop] = '\0';

                    // Create it if it does not exist
                    if(stat(path_part, &stat_info) != 0) {
                        if (errno == ENOENT) {
                            mkdir(path_part, 0755);
                        }
                    }

                    delete[] path_part;
                }
            }
            //3. Get the end if the string is not terminated by a "/"
            if (path_full[strlen(path_full)-1] != '/') {
                // We have already checked that the folder doesn't exist
                mkdir (path_full, 0755);
            }
        }
        else {
            G4cerr << "ERROR: Could not lookup folder " << foldername_out << " - aborting!" << G4endl;
            exit(1);
        }
    }
    else if(not S_ISDIR(stat_info.st_mode)) {
        G4cerr << "ERROR: An non-folder entity named " << foldername_out << " already exist- aborting!" << G4endl;
        exit(1);
    }
    else {
        G4cout << "Folder '" << foldername_out << "' already exists -- using it!" << G4endl;
    }

    G4cout << "Opening ROOT file '" + rootFileName +"'"<<G4endl;
    histFile = new TFile(rootFileName,"RECREATE");
    if ( not histFile->IsOpen() ) {
        G4cerr << "Opening TFile '" << rootFileName << "' failed; quitting." << G4endl;
        exit(1);
    }

    eventCounter = 0;

    RNG = new TRandom1((UInt_t) rngSeed);

    // TTrees for external analysis
    if (not miniFile) {
        if (detCon->GetHasTarget()) {
            targetExit = new TTree("TargetExit","TargetExit tree");
            targetExit->Branch("TargetExitBranch", &targetExitBuffer,
                               "x/D:y:z:px:py:pz:E:PDG/I:charge:eventID");
        }

        trackerHits = new TTree("TrackerHits","TrackerHits tree");
        trackerHits->Branch("TrackerHitsBranch", &trackerHitsBuffer,
                            "x/D:y:z:px:py:pz:E:PDG/I:charge:eventID");

        magnetEdeps = new TTree("magnetEdeps", "Magnet Edeps tree");
    }

    // Target energy deposition
    if (detCon->GetHasTarget()) {
        targetEdep = new TH1D("targetEdep","targetEdep",engNbins,0,beamEnergy);
        targetEdep->GetXaxis()->SetTitle("Total energy deposit/event [MeV]");
        targetEdep_NIEL = new TH1D("targetEdep_NIEL","targetEdep_NIEL",1000,0,1);
        targetEdep_NIEL->GetXaxis()->SetTitle("Total NIEL/event [keV]");
        targetEdep_IEL = new TH1D("targetEdep_IEL","targetEdep_IEL",engNbins,0,beamEnergy);
        targetEdep_IEL->GetXaxis()->SetTitle("Total ionizing energy deposit/event [MeV]");

        if(edep_dens_dz != 0.0) {
            G4int target_edep_nbins_dz = (int) ceil((detCon->getTargetThickness()/mm)/fabs(this->edep_dens_dz));
            G4cout << "NBINS_DZ for target_edep_dens = " << target_edep_nbins_dz << G4endl;

            if ( this->edep_dens_dz > 0.0 ) {
                target_edep_dens = new TH3D("target_edep_dens",
                                            "Target energy deposition density [MeV/bin]",
                                            100, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                            100, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                            target_edep_nbins_dz, 0.0, detCon->getTargetThickness()/mm
                                            );
                target_edep_dens->GetXaxis()->SetTitle("X position [mm]");
                target_edep_dens->GetYaxis()->SetTitle("Y position [mm]");
                target_edep_dens->GetZaxis()->SetTitle("Z position [mm]");
            }
            else {
                target_edep_dens = NULL;
            }

            target_edep_rdens = new TH2D("target_edep_rdens",
                                        "Target radial energy deposition density [MeV/bin]",
                                        target_edep_nbins_dz, 0.0,detCon->getTargetThickness()/mm,
                                        1000, 0.0, 2*phasespacehist_posLim / mm
                                        );
            target_edep_rdens->GetXaxis()->SetTitle("Z position [mm]");
            target_edep_rdens->GetYaxis()->SetTitle("R position [mm]");
        }
        else {
            target_edep_dens  = NULL;
            target_edep_rdens = NULL;
        }

        // Target tracking info
        target_exit_energy[11]  = new TH1D("target_exit_energy_PDG11",
                                        "Particle energy when exiting target (electrons)",
                                        engNbins,0,beamEnergy);
        target_exit_energy[-11] = new TH1D("target_exit_energy_PDG-11",
                                        "Particle energy when exiting target (positrons)",
                                        engNbins,0,beamEnergy);
        target_exit_energy[22]  = new TH1D("target_exit_energy_PDG22",
                                        "Particle energy when exiting target (photons)",
                                        engNbins,0,beamEnergy);
        target_exit_energy[2212]= new TH1D("target_exit_energy_PDG2212",
                                        "Particle energy when exiting target (protons)",
                                        engNbins,0,beamEnergy);
        target_exit_energy[0]   = new TH1D("target_exit_energy_PDGother",
                                        "Particle energy when exiting target (other)",
                                        engNbins,0,beamEnergy);
        for (auto it : target_exit_energy) {
            it.second->GetXaxis()->SetTitle("Energy [MeV]");
        }

        target_exit_cutoff_energy[11]  = new TH1D("target_exit_cutoff_energy_PDG11",
                                                "Particle energy when exiting target (electrons) (r < Rcut, E > Ecut)",
                                                engNbins,0,beamEnergy);
        target_exit_cutoff_energy[-11] = new TH1D("target_exit_cutoff_energy_PDG-11",
                                                "Particle energy when exiting target (positrons) (r < Rcut, E > Ecut)",
                                                engNbins,0,beamEnergy);
        target_exit_cutoff_energy[22]  = new TH1D("target_exit_cutoff_energy_PDG22",
                                                "Particle energy when exiting target (photons) (r < Rcut, E > Ecut)",
                                                engNbins,0,beamEnergy);
        target_exit_cutoff_energy[2212]= new TH1D("target_exit_cutoff_energy_PDG2212",
                                                "Particle energy when exiting target (protons) (r < Rcut, E > Ecut)",
                                                engNbins,0,beamEnergy);
        target_exit_cutoff_energy[0]   = new TH1D("target_exit_cutoff_energy_PDGother",
                                                "Particle energy when exiting target (other) (r < Rcut)",
                                                engNbins,0,beamEnergy);
        for (auto it : target_exit_cutoff_energy) {
            it.second->GetXaxis()->SetTitle("Energy [MeV]");
        }

        // Target exit angle histogram
        target_exitangle_hist        = new TH1D("target_exit_angle",
                                                "Exit angle from target",
                                                5001, -90, 90);
        target_exitangle_hist_cutoff = new TH1D("target_exit_angle_cutoff",
                                                "Exit angle from target (charged, energy > Ecut, r < Rcut)",
                                                5001, -90, 90);

        // Target exit phasespace histograms
        target_exit_phasespaceX        = new TH2D("target_exit_x",
                                                "Target exit phase space (x)",
                                                1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        target_exit_phasespaceX->GetXaxis()->SetTitle("Position x [mm]");
        target_exit_phasespaceX->GetYaxis()->SetTitle("Angle dx/dz [rad]");

        target_exit_phasespaceY        = new TH2D("target_exit_y",
                                                "Target exit phase space (y)",
                                                1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        target_exit_phasespaceY->GetXaxis()->SetTitle("Position y [mm]");
        target_exit_phasespaceY->GetYaxis()->SetTitle("Angle dy/dz [rad]");

        target_exit_phasespaceX_cutoff = new TH2D("target_exit_cutoff_x",
                                                "Target exit phase space (x) (charged, energy > Ecut, r < Rcut)",
                                                1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        target_exit_phasespaceX_cutoff->GetXaxis()->SetTitle("Position x [mm]");
        target_exit_phasespaceX_cutoff->GetYaxis()->SetTitle("Angle dx/dz [rad]");

        target_exit_phasespaceY_cutoff = new TH2D("target_exit_cutoff_y",
                                                "Target exit phase space (y) (charged, energy > Ecut, r < Rcut)",
                                                1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        target_exit_phasespaceY_cutoff->GetXaxis()->SetTitle("Position y [mm]");
        target_exit_phasespaceY_cutoff->GetYaxis()->SetTitle("Angle dy/dz [rad]");
    }

    // Tracker histograms
    tracker_numParticles = new TH1D("numParticles","numParticles",1001,-0.5,1000.5);
    tracker_numParticles->GetXaxis()->SetTitle("Number of particles / event");

    tracker_energy       = new TH1D("energy","Energy of all particles hitting the tracker",10000,0,beamEnergy);
    tracker_energy->GetXaxis()->SetTitle("Energy per particle [MeV]");

    tracker_type_energy[11]  = new TH1D("tracker_energy_PDG11",
                                      "Particle energy when hitting tracker (electrons)",
                                      engNbins,0,beamEnergy);
    tracker_type_energy[-11] = new TH1D("tracker_energy_PDG-11",
                                      "Particle energy when hitting tracker (positrons)",
                                      engNbins,0,beamEnergy);
    tracker_type_energy[22]  = new TH1D("tracker_energy_PDG22",
                                      "Particle energy when hitting tracker (photons)",
                                      engNbins,0,beamEnergy);
    tracker_type_energy[2212]= new TH1D("tracker_energy_PDG2212",
                                      "Particle energy when hitting tracker (protons)",
                                      engNbins,0,beamEnergy);
    tracker_type_energy[0]   = new TH1D("tracker_energy_PDGother",
                                      "Particle energy when hitting tracker (other)",
                                      engNbins,0,beamEnergy);
    for (auto it : tracker_type_energy) {
        it.second->GetXaxis()->SetTitle("Energy [MeV]");
    }

    tracker_type_cutoff_energy[11]  = new TH1D("tracker_cutoff_energy_PDG11",
                                              "Particle energy when hitting tracker (electrons) (r < Rcut)",
                                              engNbins,0,beamEnergy);
    tracker_type_cutoff_energy[-11] = new TH1D("tracker_cutoff_energy_PDG-11",
                                              "Particle energy when hitting tracker (positrons) (r < Rcut)",
                                              engNbins,0,beamEnergy);
    tracker_type_cutoff_energy[22]  = new TH1D("tracker_cutoff_energy_PDG22",
                                              "Particle energy when hitting tracker (photons) (r < Rcut)",
                                              engNbins,0,beamEnergy);
    tracker_type_cutoff_energy[2212]= new TH1D("tracker_cutoff_energy_PDG2212",
                                              "Particle energy when hitting tracker (protons) (r < Rcut)",
                                              engNbins,0,beamEnergy);
    tracker_type_cutoff_energy[0]   = new TH1D("tracker_cutoff_energy_PDGother",
                                              "Particle energy when hitting tracker (other) (r < Rcut)",
                                              engNbins,0,beamEnergy);
    for (auto it : tracker_type_cutoff_energy) {
        it.second->GetXaxis()->SetTitle("Energy [MeV]");
    }

    tracker_hitPos        = new TH2D("trackerHitpos", "Tracker Hit position",
               1000,-detCon->getDetectorSizeX()/2.0/mm,detCon->getDetectorSizeX()/2.0/mm,
               1000,-detCon->getDetectorSizeY()/2.0/mm,detCon->getDetectorSizeY()/2.0/mm);
    tracker_hitPos_cutoff = new TH2D("trackerHitpos_cutoff", "Tracker Hit position (charged, energy > Ecut)",
               1000,-position_cutoffR, position_cutoffR,
               1000,-position_cutoffR, position_cutoffR);

    tracker_phasespaceX   =
        new TH2D("tracker_x",
                 "Tracker phase space (x)",
                 1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                 1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
    //tracker_phasespaceX->Sumw2();
    tracker_phasespaceY   =
        new TH2D("tracker_y",
                 "Tracker phase space (y)",
                 1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                 1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
    //tracker_phasespaceY->Sumw2();

    tracker_phasespaceX_cutoff   =
        new TH2D("tracker_cutoff_x",
                 "Tracker phase space (x) (charged, energy > Ecut, r < Rcut)",
                 1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                 1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
    //tracker_phasespaceX_cutoff->Sumw2();
    tracker_phasespaceY_cutoff   =
        new TH2D("tracker_cutoff_y",
                 "Tracker phase space (y) (charged, energy > Ecut, r < Rcut)",
                 1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                 1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
    //tracker_phasespaceY_cutoff->Sumw2();

    init_phasespaceX   =
        new TH2D("init_x",
                 "Initial phase space (x)",
                 1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                 1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
    //init_phasespaceX->Sumw2();
    init_phasespaceY   =
        new TH2D("init_y",
                 "Initial phase space (y)",
                 1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                 1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
    //init_phasespaceY->Sumw2();
    init_phasespaceXY   =
        new TH2D("init_xy",
                 "Initial phase space (x,y)",
                 1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                 1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm);
    init_E =
        new TH1D("init_E",
                 "Initial particle energy",
                 1000, 0.0, max(beamEnergy*1.1,genAct->get_beam_energy_flatMax()));

    // Limit for radial histograms
    G4double minR = min(detCon->getWorldSizeX(),detCon->getWorldSizeY())/mm;

    // Target R position
    if (detCon->GetHasTarget()) {
        target_exit_Rpos[11]  = new TH1D("target_exit_rpos_PDG11",
                                        "target exit rpos (electrons)",
                                        1000,0,minR);
        target_exit_Rpos[-11] = new TH1D("target_exit_rpos_PDG-11",
                                        "target exit rpos (positrons)",
                                        1000,0,minR);
        target_exit_Rpos[22]  = new TH1D("target_exit_rpos_PDG22",
                                        "target exit rpos (photons)",
                                        1000,0,minR);
        target_exit_Rpos[2212]= new TH1D("target_exit_rpos_PDG2212",
                                        "target exit rpos (protons)",
                                        1000,0,minR);
        target_exit_Rpos[0]   = new TH1D("target_exit_rpos_PDGother",
                                        "target exit rpos (other)",
                                        1000,0,minR);
        for (auto hist : target_exit_Rpos) {
            hist.second->GetXaxis()->SetTitle("R [mm]");
        }
        target_exit_Rpos_cutoff[11]  = new TH1D("target_exit_rpos_cutoff_PDG11",
                                                "target exit rpos (electrons, energy > E_cut)",
                                                1000,0,minR);
        target_exit_Rpos_cutoff[-11] = new TH1D("target_exit_rpos_cutoff_PDG-11",
                                                "target exit rpos (positrons, energy > E_cut)",
                                                1000,0,minR);
        target_exit_Rpos_cutoff[22]  = new TH1D("target_exit_rpos_cutoff_PDG22",
                                                "target exit rpos (photons, energy > E_cut)",
                                                1000,0,minR);
        target_exit_Rpos_cutoff[2212]= new TH1D("target_exit_rpos_cutoff_PDG2212",
                                                "target exit rpos (protons, energy > E_cut)",
                                                1000,0,minR);
        target_exit_Rpos_cutoff[0]   = new TH1D("target_exit_rpos_cutoff_PDGother",
                                                "target exit rpos (other, energy > E_cut)",
                                                1000,0,minR);
        for (auto hist : target_exit_Rpos_cutoff) {
            hist.second->GetXaxis()->SetTitle("R [mm]");
        }
    }

    // Tracker R position
    tracker_Rpos[11]  = new TH1D("tracker_rpos_PDG11",
                                 "tracker rpos (electrons)",
                                 1000,0,minR);
    tracker_Rpos[-11] = new TH1D("tracker_rpos_PDG-11",
                                 "tracker rpos (positrons)",
                                 1000,0,minR);
    tracker_Rpos[22]  = new TH1D("tracker_rpos_PDG22",
                                 "tracker rpos (photons)",
                                 1000,0,minR);
    tracker_Rpos[2212]= new TH1D("tracker_rpos_PDG2212",
                                 "tracker rpos (protons)",
                                 1000,0,minR);
    tracker_Rpos[0]   = new TH1D("tracker_rpos_PDGother",
                                 "tracker rpos (other)",
                                 1000,0,minR);
    for (auto hist : tracker_Rpos) {
        hist.second->GetXaxis()->SetTitle("R [mm]");
    }
    tracker_Rpos_cutoff[11]  = new TH1D("tracker_rpos_cutoff_PDG11",
                                        "tracker rpos (electrons, energy > E_cut)",
                                        1000,0,minR);
    tracker_Rpos_cutoff[-11] = new TH1D("tracker_rpos_cutoff_PDG-11",
                                        "tracker rpos (positrons, energy > E_cut)",
                                        1000,0,minR);
    tracker_Rpos_cutoff[22]  = new TH1D("tracker_rpos_cutoff_PDG22",
                                        "tracker rpos (photons, energy > E_cut)",
                                        1000,0,minR);
    tracker_Rpos_cutoff[2212]= new TH1D("tracker_rpos_cutoff_PDG2212",
                                        "tracker rpos (protons, energy > E_cut)",
                                        1000,0,minR);
    tracker_Rpos_cutoff[0]   = new TH1D("tracker_rpos_cutoff_PDGother",
                                        "tracker rpos (other, energy > E_cut)",
                                        1000,0,minR);
    for (auto hist : tracker_Rpos_cutoff) {
        hist.second->GetXaxis()->SetTitle("R [mm]");
    }

    //For counting the types of particles hitting the detectors (for magnets it is defined elsewhere)
    if (detCon->GetHasTarget()){
        typeCounter["tracker"]       = particleTypesCounter();
        typeCounter["tracker_cutoff"] = particleTypesCounter();
    }
    typeCounter["target"]        = particleTypesCounter();
    typeCounter["target_cutoff"] = particleTypesCounter();

    //Compute means and RMS of where they hit
    tracker_particleHit_x  = 0.0;
    tracker_particleHit_xx = 0.0;
    tracker_particleHit_y  = 0.0;
    tracker_particleHit_yy = 0.0;

    //Compute means and RMS of where they hit (above cutoff particles only)
    tracker_particleHit_x_cutoff  = 0.0;
    tracker_particleHit_xx_cutoff = 0.0;
    tracker_particleHit_y_cutoff  = 0.0;
    tracker_particleHit_yy_cutoff = 0.0;
    numParticles_cutoff = 0;

    //Compute RMS of target exit angle
    if (detCon->GetHasTarget()) {
        target_exitangle              = 0.0;
        target_exitangle2             = 0.0;
        target_exitangle_numparticles = 0;
        target_exitangle_cutoff              = 0.0;
        target_exitangle2_cutoff             = 0.0;
        target_exitangle_cutoff_numparticles = 0;
    }

    // Magnet histograms
    for (auto mag : detCon->magnets) {
        const G4String magName = mag->magnetName;

        magnet_edep.push_back( new TH1D((magName + "_edep").c_str(),(magName + " edep").c_str(),
                                        1000,0,beamEnergy) );
        magnet_edep.back()->GetXaxis()->SetTitle("Total energy deposit/event [MeV]");

        G4int mag_edep_nbins_dz = (int) ceil((mag->GetLength()/mm) / fabs(this->edep_dens_dz));
        if(edep_dens_dz != 0.0) {
            G4cout << "NBINS_DZ for " << magName << "_edep_dens = " << mag_edep_nbins_dz << G4endl;

            if (this->edep_dens_dz > 0.0 ) {
               TH3D* mag_edep_dens = new TH3D((magName + "_edep_dens").c_str(),
                                              (magName + " energy deposition density [MeV/bin]").c_str(),
                                              100,-phasespacehist_posLim/mm, phasespacehist_posLim/mm,
                                              100,-phasespacehist_posLim/mm, phasespacehist_posLim/mm,
                                              mag_edep_nbins_dz, 0.0, mag->GetLength()/mm);
                mag_edep_dens->GetXaxis()->SetTitle("X position [mm]");
                mag_edep_dens->GetYaxis()->SetTitle("Y position [mm]");
                mag_edep_dens->GetZaxis()->SetTitle("Z position [mm]");
            
                magnet_edep_dens.push_back(mag_edep_dens);
            }
            else {
                magnet_edep_dens.push_back(NULL);
            }

            TH2D* mag_edep_rdens = new TH2D((magName + "_edep_rdens").c_str(),
                                            (magName + " radial energy deposition density [MeV/bin]").c_str(),
                                            mag_edep_nbins_dz, 0.0, mag->GetLength()/mm,
                                            1000, 0.0, 2*phasespacehist_posLim / mm );

            mag_edep_rdens->GetXaxis()->SetTitle("Z position [mm]");
            mag_edep_rdens->GetYaxis()->SetTitle("R position [mm]");

            magnet_edep_rdens.push_back(mag_edep_rdens);
        }
        else {
            magnet_edep_dens.push_back(NULL);
            magnet_edep_rdens.push_back(NULL);
        }

        //G4double minR = min(detCon->getWorldSizeX(),detCon->getWorldSizeY())/mm;
        magnet_exit_Rpos.push_back(std::map<G4int,TH1D*>());
        magnet_exit_Rpos.back()[11]  = new TH1D((magName + "_rpos_PDG11").c_str(),
                                                (magName + " rpos (electrons)").c_str(),
                                                1000,0,minR);
        magnet_exit_Rpos.back()[-11] = new TH1D((magName + "_rpos_PDG-11").c_str(),
                                                (magName + " rpos (positrons)").c_str(),
                                                1000,0,minR);
        magnet_exit_Rpos.back()[22]  = new TH1D((magName + "_rpos_PDG22").c_str(),
                                                (magName + " rpos (photons)").c_str(),
                                                1000,0,minR);
        magnet_exit_Rpos.back()[2212]= new TH1D((magName + "_rpos_PDG2212").c_str(),
                                                (magName + " rpos (protons)").c_str(),
                                                1000,0,minR);
        magnet_exit_Rpos.back()[0]   = new TH1D((magName + "_rpos_PDGother").c_str(),
                                                (magName + " rpos (other)").c_str(),
                                                1000,0,minR);
        for (auto hist : magnet_exit_Rpos.back()) {
            hist.second->GetXaxis()->SetTitle("R [mm]");
        }

        magnet_exit_Rpos_cutoff.push_back(std::map<G4int,TH1D*>());
        magnet_exit_Rpos_cutoff.back()[11]  = new TH1D((magName + "_rpos_cutoff_PDG11").c_str(),
                                                       (magName + " rpos (electrons, energy > Ecut)").c_str(),
                                                       1000,0,minR);
        magnet_exit_Rpos_cutoff.back()[-11] = new TH1D((magName + "_rpos_cutoff_PDG-11").c_str(),
                                                       (magName + " rpos (positrons, energy > Ecut)").c_str(),
                                                       1000,0,minR);
        magnet_exit_Rpos_cutoff.back()[22]  = new TH1D((magName + "_rpos_cutoff_PDG22").c_str(),
                                                       (magName + " rpos (photons, energy > Ecut)").c_str(),
                                                       1000,0,minR);
        magnet_exit_Rpos_cutoff.back()[2212]= new TH1D((magName + "_rpos_cutoff_PDG2212").c_str(),
                                                       (magName + " rpos (protons, energy > Ecut)").c_str(),
                                                       1000,0,minR);
        magnet_exit_Rpos_cutoff.back()[0]   = new TH1D((magName + "_rpos_cutoff_PDGother").c_str(),
                                                       (magName + " rpos (other, energy > Ecut)").c_str(),
                                                       1000,0,minR);
        for (auto hist : magnet_exit_Rpos_cutoff.back()) {
            hist.second->GetXaxis()->SetTitle("R [mm]");
        }

        magnet_exit_phasespaceX.push_back
            ( new TH2D((magName+"_x").c_str(),
                       (magName+" phase space (x)").c_str(),
                       1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                       1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad)
              );
        magnet_exit_phasespaceX.back()->GetXaxis()->SetTitle("X [mm]");
        magnet_exit_phasespaceX.back()->GetYaxis()->SetTitle("X'");

        magnet_exit_phasespaceY.push_back
            ( new TH2D((magName+"_y").c_str(),
                       (magName+" phase space (y)").c_str(),
                       1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                       1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad)
              );
        magnet_exit_phasespaceY.back()->GetXaxis()->SetTitle("Y [mm]");
        magnet_exit_phasespaceY.back()->GetYaxis()->SetTitle("Y'");

        magnet_exit_phasespaceX_cutoff.push_back
            ( new TH2D((magName+"_cutoff_x").c_str(),
                       (magName+" phase space (x) (charged, energy > Ecut, r < Rcut)").c_str(),
                       1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                       1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad)
              );
        magnet_exit_phasespaceX_cutoff.back()->GetXaxis()->SetTitle("X [mm]");
        magnet_exit_phasespaceX_cutoff.back()->GetYaxis()->SetTitle("X'");

        magnet_exit_phasespaceY_cutoff.push_back
            ( new TH2D((magName+"_cutoff_y").c_str(),
                       (magName+" phase space (y) (charged, energy > Ecut, r < Rcut)").c_str(),
                       1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                       1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad)
              );
        magnet_exit_phasespaceY_cutoff.back()->GetXaxis()->SetTitle("Y [mm]");
        magnet_exit_phasespaceY_cutoff.back()->GetYaxis()->SetTitle("Y'");

        typeCounter[magName]             = particleTypesCounter();
        typeCounter[magName + "_cutoff"] = particleTypesCounter();

        magnet_exit_energy.push_back(std::map<G4int,TH1D*>());
        magnet_exit_energy.back()[11]  = new TH1D((magName+"_exit_energy_PDG11").c_str(),
                                                  ("Particle energy when exiting "+magName+" (electrons)").c_str(),
                                                  engNbins,0,beamEnergy);
        magnet_exit_energy.back()[-11] = new TH1D((magName+"_exit_energy_PDG-11").c_str(),
                                                  ("Particle energy when exiting "+magName+" (positrons)").c_str(),
                                                  engNbins,0,beamEnergy);
        magnet_exit_energy.back()[22]  = new TH1D((magName+"_exit_energy_PDG22").c_str(),
                                                  ("Particle energy when exiting "+magName+" (photons)").c_str(),
                                                  engNbins,0,beamEnergy);
        magnet_exit_energy.back()[2212]= new TH1D((magName+"_exit_energy_PDG2212").c_str(),
                                                  ("Particle energy when exiting "+magName+" (protons)").c_str(),
                                                  engNbins,0,beamEnergy);
        magnet_exit_energy.back()[0]   = new TH1D((magName+"_exit_energy_PDGother").c_str(),
                                                  ("Particle energy when exiting "+magName+" (other)").c_str(),
                                                  engNbins,0,beamEnergy);
        for (auto PDG : magnet_exit_energy.back()) {
            PDG.second->GetXaxis()->SetTitle("Energy [MeV]");
        }

        magnet_exit_cutoff_energy.push_back(std::map<G4int,TH1D*>());
        magnet_exit_cutoff_energy.back()[11]  = new TH1D((magName+"_exit_cutoff_energy_PDG11").c_str(),
                                                  ("Particle energy when exiting "+magName+" (electrons, r < Rcut)").c_str(),
                                                  engNbins,0,beamEnergy);
        magnet_exit_cutoff_energy.back()[-11] = new TH1D((magName+"_exit_cutoff_energy_PDG-11").c_str(),
                                                  ("Particle energy when exiting "+magName+" (positrons, r < Rcut)").c_str(),
                                                  engNbins,0,beamEnergy);
        magnet_exit_cutoff_energy.back()[22]  = new TH1D((magName+"_exit_cutoff_energy_PDG22").c_str(),
                                                  ("Particle energy when exiting "+magName+" (photons, r < Rcut)").c_str(),
                                                  engNbins,0,beamEnergy);
        magnet_exit_cutoff_energy.back()[2212]= new TH1D((magName+"_exit_cutoff_energy_PDG2212").c_str(),
                                                  ("Particle energy when exiting "+magName+" (protons, r < Rcut)").c_str(),
                                                  engNbins,0,beamEnergy);
        magnet_exit_cutoff_energy.back()[0]   = new TH1D((magName+"_exit_cutoff_energy_PDGother").c_str(),
                                                  ("Particle energy when exiting "+magName+" (other, r < Rcut)").c_str(),
                                                  engNbins,0,beamEnergy);
        for (auto PDG : magnet_exit_cutoff_energy.back()) {
            PDG.second->GetXaxis()->SetTitle("Energy [MeV]");
        }
    }

    if (not miniFile) {
        size_t numMagnets = detCon->magnets.size();
        if (numMagnets > 0) {
            magnetEdepsBuffer = new Double_t[numMagnets];
            size_t i = 0;
            for (auto mag : detCon->magnets) {
                G4String magName = mag->magnetName;
                magnetEdeps->Branch(magName, &(magnetEdepsBuffer[i]), (magName+"/D").c_str());
                i++;
            }
        }
    }
}

void RootFileWriter::doEvent(const G4Event* event){

    G4RunManager*           run    = G4RunManager::GetRunManager();
    DetectorConstruction*   detCon = (DetectorConstruction*)run->GetUserDetectorConstruction();
    PrimaryGeneratorAction* genAct = (PrimaryGeneratorAction*)run->GetUserPrimaryGeneratorAction();

    eventCounter++;

    G4HCofThisEvent* HCE=event->GetHCofThisEvent();
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    // *** Data from TargetSD ***
    if (detCon->GetHasTarget()) {
        G4int myTargetEdep_CollID = SDman->GetCollectionID("target_edep");
        if (myTargetEdep_CollID>=0){
            MyEdepHitsCollection* targetEdepHitsCollection = NULL;
            targetEdepHitsCollection = (MyEdepHitsCollection*) (HCE->GetHC(myTargetEdep_CollID));
            if (targetEdepHitsCollection != NULL) {
                G4int nEntries = targetEdepHitsCollection->entries();
                G4double edep      = 0.0; // G4 units, normalized before Fill()
                G4double edep_NIEL = 0.0; // G4 units, normalized before Fill()
                G4double edep_IEL  = 0.0; // G4 units, normalized before Fill()
                for (G4int i = 0; i < nEntries; i++){
                    MyEdepHit* edepHit = (*targetEdepHitsCollection)[i];
                    if (edepHit->GetDepositedEnergy() < 1e-20*MeV) continue;

                    edep      += edepHit->GetDepositedEnergy();
                    edep_NIEL += edepHit->GetDepositedEnergy_NIEL();
                    edep_IEL  += edepHit->GetDepositedEnergy() - edepHit->GetDepositedEnergy_NIEL();

                    //Randomly spread the energy deposits over the step
                    if (target_edep_rdens != NULL) {
                        G4ThreeVector edepStep = edepHit->GetPostStepPoint() - edepHit->GetPreStepPoint();
                        G4double edepStepLen = edepStep.mag();
                        int numSamples = (int) ceil(2*(edepStepLen/mm)/fabs(edep_dens_dz));
                        for (int j = 0; j < numSamples; j++){
                            G4ThreeVector posSample = edepHit->GetPreStepPoint() + RNG->Uniform(edepStepLen)*edepStep;
                            G4double sample_z = posSample.z() + detCon->getTargetThickness()/2.0;
                            if ( this->edep_dens_dz > 0.0 ) {
                                target_edep_dens->Fill(posSample.x()/mm,
                                                       posSample.y()/mm,
                                                       sample_z/mm,
                                                       edepHit->GetDepositedEnergy()/numSamples/MeV);
                            }
                            G4double sample_r = sqrt(posSample.x()*posSample.x() + posSample.y()*posSample.y());
                            target_edep_rdens->Fill(sample_z/mm, sample_r/mm, edepHit->GetDepositedEnergy()/numSamples/MeV);
                        }
                    }
                }

                targetEdep->Fill(edep/MeV);
                targetEdep_NIEL->Fill(edep_NIEL/keV);
                targetEdep_IEL->Fill(edep_IEL/MeV);
            }
            else {
                G4cout << "targetEdepHitsCollection was NULL!"<<G4endl;
            }
        }
        else {
            G4cout << "myTargetEdep_CollID was " << myTargetEdep_CollID << " < 0!"<<G4endl;
        }

        G4int myTargetExitpos_CollID = SDman->GetCollectionID("target_exitpos");
        if (myTargetExitpos_CollID>=0) {
            MyTrackerHitsCollection* targetExitposHitsCollection = NULL;
            targetExitposHitsCollection = (MyTrackerHitsCollection*) (HCE->GetHC(myTargetExitpos_CollID));
            if (targetExitposHitsCollection != NULL) {
                G4int nEntries = targetExitposHitsCollection->entries();

                for (G4int i = 0; i < nEntries; i++) {
                    //Get the data from the event
                    const G4double       energy      = (*targetExitposHitsCollection)[i]->GetTrackEnergy();
                    const G4int          charge      = (*targetExitposHitsCollection)[i]->GetCharge();
                    const G4ThreeVector& momentum    = (*targetExitposHitsCollection)[i]->GetMomentum();
                    const G4double       exitangle   = atan(momentum.x()/momentum.z())/deg;
                    const G4ThreeVector& hitPos      = (*targetExitposHitsCollection)[i]->GetPosition();
                    const G4double       hitR        = sqrt(hitPos.x()*hitPos.x() + hitPos.y()*hitPos.y());
                    const G4int          PDG         = (*targetExitposHitsCollection)[i]->GetPDG();
                    const G4String&      type        = (*targetExitposHitsCollection)[i]->GetType();

                    //Particle type counting
                    FillParticleTypes(typeCounter["target"], PDG, type);
                    if (energy/MeV > beamEnergy*beamEnergy_cutoff and hitR/mm < position_cutoffR) {
                        FillParticleTypes(typeCounter["target_cutoff"], PDG, type);
                    }

                    //Exit angle
                    target_exitangle_hist->Fill(exitangle);
                    if (charge != 0 and energy/MeV > beamEnergy*beamEnergy_cutoff and hitR/mm < position_cutoffR) {
                        target_exitangle_hist_cutoff->Fill(exitangle);
                    }

                    target_exitangle              += exitangle;
                    target_exitangle2             += exitangle*exitangle;
                    target_exitangle_numparticles += 1;

                    if (charge != 0 and energy/MeV > beamEnergy*beamEnergy_cutoff  and hitR/mm < position_cutoffR) {
                        target_exitangle_cutoff              += exitangle;
                        target_exitangle2_cutoff             += exitangle*exitangle;
                        target_exitangle_cutoff_numparticles += 1;
                    }

                    //Phase space
                    target_exit_phasespaceX->Fill(hitPos.x()/mm, momentum.x()/momentum.z());
                    target_exit_phasespaceY->Fill(hitPos.y()/mm, momentum.y()/momentum.z());

                    if (charge != 0 and energy/MeV > beamEnergy*beamEnergy_cutoff and hitR/mm < position_cutoffR) {
                        target_exit_phasespaceX_cutoff->Fill(hitPos.x()/mm, momentum.x()/momentum.z());
                        target_exit_phasespaceY_cutoff->Fill(hitPos.y()/mm, momentum.y()/momentum.z());
                    }

                    //Energy
                    if (target_exit_energy.find(PDG) != target_exit_energy.end()) {
                        target_exit_energy[PDG]->Fill(energy/MeV);
                    }
                    else {
                        target_exit_energy[0]->Fill(energy/MeV);
                    }

                    if (hitR/mm < position_cutoffR and energy/MeV > beamEnergy*beamEnergy_cutoff) {
                        if (target_exit_cutoff_energy.find(PDG) != target_exit_cutoff_energy.end()) {
                            target_exit_cutoff_energy[PDG]->Fill(energy/MeV);
                        }
                        else {
                            target_exit_cutoff_energy[0]->Fill(energy/MeV);
                        }
                    }

                    //R position
                    if (target_exit_Rpos.find(PDG) != target_exit_Rpos.end()) {
                        target_exit_Rpos[PDG]->Fill(hitR/mm);
                    }
                    else {
                        target_exit_Rpos[0]->Fill(hitR/mm);
                    }
                    if (energy/MeV > beamEnergy*beamEnergy_cutoff) {
                        if (target_exit_Rpos_cutoff.find(PDG) != target_exit_Rpos_cutoff.end()) {
                            target_exit_Rpos_cutoff[PDG]->Fill(hitR/mm);
                        }
                        else {
                            target_exit_Rpos_cutoff[0]->Fill(hitR/mm);
                        }
                    }

                    //Fill the TTree
                    if (not miniFile) {
                        targetExitBuffer.x = hitPos.x()/mm;
                        targetExitBuffer.y = hitPos.x()/mm;
                        targetExitBuffer.z = hitPos.z()/mm;

                        targetExitBuffer.px = momentum.x()/MeV;
                        targetExitBuffer.py = momentum.y()/MeV;
                        targetExitBuffer.pz = momentum.z()/MeV;

                        targetExitBuffer.E = energy / MeV;

                        targetExitBuffer.PDG = PDG;
                        targetExitBuffer.charge = charge;

                        targetExitBuffer.eventID = eventCounter;

                        targetExit->Fill();
                    }
                }

            }
            else {
                G4cout << "targetExitposHitsCollection was NULL!"<<G4endl;
            }
        }
        else {
            G4cout << "myTargetExitpos_CollID was " << myTargetExitpos_CollID << " < 0!"<<G4endl;
        }
    }

    //**Data from detectorTrackerSD**
    G4int myTrackerSD_CollID = SDman->GetCollectionID("TrackerCollection");
    if (myTrackerSD_CollID>=0) {
        MyTrackerHitsCollection* trackerHitsCollection = NULL;
        trackerHitsCollection = (MyTrackerHitsCollection*) (HCE->GetHC(myTrackerSD_CollID));
        if (trackerHitsCollection != NULL) {
            G4int nEntries = trackerHitsCollection->entries();

            for (G4int i = 0; i < nEntries; i++) {
                //Get the data from the event
                const G4double  energy = (*trackerHitsCollection)[i]->GetTrackEnergy();
                const G4int     PDG    = (*trackerHitsCollection)[i]->GetPDG();
                const G4int     charge = (*trackerHitsCollection)[i]->GetCharge();
                const G4String& type   = (*trackerHitsCollection)[i]->GetType();
                const G4ThreeVector& hitPos   = (*trackerHitsCollection)[i]->GetPosition();
                const G4ThreeVector& momentum = (*trackerHitsCollection)[i]->GetMomentum();
                const G4double       hitR     = sqrt(hitPos.x()*hitPos.x() + hitPos.y()*hitPos.y());

                //Overall histograms
                tracker_energy->Fill(energy/MeV);

                if (tracker_type_energy.find(PDG) != tracker_type_energy.end()) {
                    tracker_type_energy[PDG]->Fill(energy/MeV);
                }
                else {
                    tracker_type_energy[0]->Fill(energy/MeV);
                }

                if (hitR/mm < position_cutoffR) {
                    if (tracker_type_cutoff_energy.find(PDG) != tracker_type_cutoff_energy.end()) {
                        tracker_type_cutoff_energy[PDG]->Fill(energy/MeV);
                    }
                    else {
                        tracker_type_cutoff_energy[0]->Fill(energy/MeV);
                    }
                }

                //Hit position
                tracker_hitPos->Fill(hitPos.x()/mm, hitPos.y()/mm);
                if (charge != 0 and energy/MeV > beamEnergy*beamEnergy_cutoff and hitR/mm < position_cutoffR) {
                    tracker_hitPos_cutoff->Fill(hitPos.x()/mm, hitPos.y()/mm);
                }

                //Phase space
                tracker_phasespaceX->Fill(hitPos.x()/mm, momentum.x()/momentum.z());
                tracker_phasespaceY->Fill(hitPos.y()/mm, momentum.y()/momentum.z());

                if (charge != 0 and energy/MeV > beamEnergy*beamEnergy_cutoff and hitR/mm < position_cutoffR) {
                    tracker_phasespaceX_cutoff->Fill(hitPos.x()/mm, momentum.x()/momentum.z());
                    tracker_phasespaceY_cutoff->Fill(hitPos.y()/mm, momentum.y()/momentum.z());
                }

                //Particle type counting
                FillParticleTypes(typeCounter["tracker"], PDG, type);
                if (energy/MeV > beamEnergy*beamEnergy_cutoff and hitR/mm < position_cutoffR) {
                    FillParticleTypes(typeCounter["tracker_cutoff"], PDG, type);
                }

                //Hit positions
                tracker_particleHit_x  +=  hitPos.x()/mm;
                tracker_particleHit_xx += (hitPos.x()/mm)*(hitPos.x()/mm);
                tracker_particleHit_y  +=  hitPos.y()/mm;
                tracker_particleHit_yy += (hitPos.y()/mm)*(hitPos.y()/mm);

                if (charge != 0 and energy/MeV > beamEnergy*beamEnergy_cutoff) {
                    tracker_particleHit_x_cutoff  +=  hitPos.x()/mm;
                    tracker_particleHit_xx_cutoff += (hitPos.x()/mm)*(hitPos.x()/mm);
                    tracker_particleHit_y_cutoff  +=  hitPos.y()/mm;
                    tracker_particleHit_yy_cutoff += (hitPos.y()/mm)*(hitPos.y()/mm);
                    numParticles_cutoff += 1;
                }

                //R position
                if (tracker_Rpos.find(PDG) != tracker_Rpos.end()) {
                    tracker_Rpos[PDG]->Fill(hitR/mm);
                }
                else {
                    tracker_Rpos[0]->Fill(hitR/mm);
                }
                if (energy/MeV > beamEnergy*beamEnergy_cutoff) {
                    if (tracker_Rpos_cutoff.find(PDG) != tracker_Rpos_cutoff.end()) {
                        tracker_Rpos_cutoff[PDG]->Fill(hitR/mm);
                    }
                    else {
                        tracker_Rpos_cutoff[0]->Fill(hitR/mm);
                    }
                }

                //Fill the TTree
                if (not miniFile) {
                    trackerHitsBuffer.x = hitPos.x()/mm;
                    trackerHitsBuffer.y = hitPos.x()/mm;
                    trackerHitsBuffer.z = hitPos.z()/mm;

                    trackerHitsBuffer.px = momentum.x()/MeV;
                    trackerHitsBuffer.py = momentum.y()/MeV;
                    trackerHitsBuffer.pz = momentum.z()/MeV;

                    trackerHitsBuffer.E = energy / MeV;

                    trackerHitsBuffer.PDG = PDG;
                    trackerHitsBuffer.charge = charge;

                    trackerHitsBuffer.eventID = eventCounter;

                    trackerHits->Fill();
                }
            }

            tracker_numParticles->Fill(nEntries);
        }
        else{
            G4cout << "trackerHitsCollection was NULL!"<<G4endl;
        }
    }
    else{
        G4cout << "myTrackerSD_CollID was " << myTrackerSD_CollID << "<0!"<<G4endl;
    }

    // Initial particle distribution
    init_phasespaceX->Fill(genAct->x/mm,genAct->xp/rad);
    init_phasespaceY->Fill(genAct->y/mm,genAct->yp/rad);
    init_phasespaceXY->Fill(genAct->x/mm,genAct->y/mm);
    init_E->Fill(genAct->E/MeV);

    // *** Data from Magnets, which use a TargetSD ***
    size_t magIdx = -1;
    for (auto mag : detCon->magnets) {
        const G4String magName = mag->magnetName;
        magIdx++;

        //Edep data collection
        G4int myMagnetEdep_CollID = SDman->GetCollectionID(magName+"_edep");
        if (myMagnetEdep_CollID>=0){
            MyEdepHitsCollection* magnetEdepHitsCollection = NULL;
            magnetEdepHitsCollection = (MyEdepHitsCollection*) (HCE->GetHC(myMagnetEdep_CollID));
            if (magnetEdepHitsCollection != NULL) {
                G4int nEntries = magnetEdepHitsCollection->entries();
                G4double edep      = 0.0;
                for (G4int i = 0; i < nEntries; i++){
                    MyEdepHit* edepHit = (*magnetEdepHitsCollection)[i];
                    if (edepHit->GetDepositedEnergy() < 1e-20*MeV) continue;

                    edep      += edepHit->GetDepositedEnergy();

                    //Randomly spread the energy deposits over the step
                    if (magnet_edep_rdens[magIdx] != NULL) {
                        G4ThreeVector edepStep = edepHit->GetPostStepPoint() - edepHit->GetPreStepPoint();
                        G4double edepStepLen = edepStep.mag();
                        int numSamples = (int) ceil(2*(edepStepLen/mm)/fabs(edep_dens_dz));
                        for (int j = 0; j < numSamples; j++){
                            G4ThreeVector posSample = edepHit->GetPreStepPoint() + RNG->Uniform(edepStepLen)*edepStep;
                            G4double sample_z = posSample.z() + mag->GetLength()/2.0;
                            if (this->edep_dens_dz > 0.0 ) {
                                magnet_edep_dens[magIdx]->Fill(posSample.x()/mm,
                                                              posSample.y()/mm,
                                                              sample_z/mm,
                                                              edepHit->GetDepositedEnergy()/numSamples/MeV);
                            }

                            G4double sample_r = sqrt(posSample.x()*posSample.x() + posSample.y()*posSample.y());
                            magnet_edep_rdens[magIdx]->Fill(sample_z/mm,
                                                            sample_r/mm,
                                                            edepHit->GetDepositedEnergy()/numSamples/MeV);

                            //Debug
                            /*
                            G4cout << eventCounter << " : " << j << "/" << numSamples << ", " << mag->magnetName << " : "
                                   << posSample.x()/mm << ", " << posSample.y()/mm << ", " << sample_z/mm << ", "
                                   << edepHit->GetDepositedEnergy()/numSamples/MeV << G4endl;
                            */
                        }

                        //DEBUG
                        if (edepHit->GetPreStepPoint().z() < -mag->GetLength()/2.0 || edepHit->GetPostStepPoint().z() < -mag->GetLength()/2.0) {
                            G4cout << eventCounter << ", " << mag->magnetName << " : "
                                   << edepHit->GetPreStepPoint() << " -> " << edepHit->GetPreStepPoint() 
                                   << ", " << edepHit->GetDepositedEnergy()/MeV << G4endl;
                        }
                    }
                }

                magnet_edep[magIdx]->Fill(edep/MeV);

                //TTree, for event-by-event analysis
                if (not miniFile){
                    magnetEdepsBuffer[magIdx] = edep/MeV;
                }
            }
            else {
                G4cout << "magnetEdepHitsCollection was NULL for '" << magName << "'!" <<G4endl;
            }
        }
        else {
            G4cout << "myMagnetEdep_CollID was " << myMagnetEdep_CollID << " < 0 for '" << magName << "'!"<<G4endl;
        }


        // Exitpos data collection
        G4int myMagnetExitpos_CollID = SDman->GetCollectionID(magName + "_exitpos");
        if (myMagnetExitpos_CollID>=0) {
            MyTrackerHitsCollection* magnetExitposHitsCollection = NULL;
            magnetExitposHitsCollection = (MyTrackerHitsCollection*) (HCE->GetHC(myMagnetExitpos_CollID));
            if (magnetExitposHitsCollection != NULL) {
                G4int nEntries = magnetExitposHitsCollection->entries();

                for (G4int i = 0; i < nEntries; i++) {
                    //Get the data from the event
                    const G4double       energy      = (*magnetExitposHitsCollection)[i]->GetTrackEnergy();
                    const G4int          charge      = (*magnetExitposHitsCollection)[i]->GetCharge();
                    const G4ThreeVector& momentum    = (*magnetExitposHitsCollection)[i]->GetMomentum();
                    //const G4double       exitangle   = atan(momentum.x()/momentum.z())/deg;
                    const G4ThreeVector& hitPos      = (*magnetExitposHitsCollection)[i]->GetPosition();
                    const G4double       hitR        = sqrt(hitPos.x()*hitPos.x() + hitPos.y()*hitPos.y());
                    const G4int          PDG         = (*magnetExitposHitsCollection)[i]->GetPDG();
                    const G4String&      type        = (*magnetExitposHitsCollection)[i]->GetType();

                    if ( abs( hitPos.z() -
                              (detCon->magnets[magIdx]->GetLength()/2.0 +
                               detCon->magnets[magIdx]->getZ0()          )
                              ) < 1e-7 ) {
                        // We are on the downstream exit face.
                        // Note: Coordinates in global coordinates.

                        //Particle type counting
                        FillParticleTypes(typeCounter[magName], PDG, type);
                        if (energy/MeV > beamEnergy*beamEnergy_cutoff and hitR/mm < position_cutoffR) {
                            FillParticleTypes(typeCounter[magName + "_cutoff"], PDG, type);
                        }

                        //Phase space
                        magnet_exit_phasespaceX[magIdx]->
                            Fill(hitPos.x()/mm, momentum.x()/momentum.z());
                        magnet_exit_phasespaceY[magIdx]->
                            Fill(hitPos.y()/mm, momentum.y()/momentum.z());

                        if ( charge != 0 and
                             energy/MeV > beamEnergy*beamEnergy_cutoff and
                             hitR/mm < position_cutoffR
                             ) {
                            magnet_exit_phasespaceX_cutoff[magIdx]->
                                Fill(hitPos.x()/mm, momentum.x()/momentum.z());
                            magnet_exit_phasespaceY_cutoff[magIdx]->
                                Fill(hitPos.y()/mm, momentum.y()/momentum.z());
                        }

                        //R position
                        if (magnet_exit_Rpos[magIdx].find(PDG) != magnet_exit_Rpos[magIdx].end()) {
                            magnet_exit_Rpos[magIdx][PDG]->Fill(hitR/mm);
                        }
                        else {
                            magnet_exit_Rpos[magIdx][0]->Fill(hitR/mm);
                        }
                        if (energy/MeV > beamEnergy*beamEnergy_cutoff) {
                            if (magnet_exit_Rpos_cutoff[magIdx].find(PDG) !=
                                magnet_exit_Rpos_cutoff[magIdx].end()) {
                                magnet_exit_Rpos_cutoff[magIdx][PDG]->Fill(hitR/mm);
                            }
                            else {
                                magnet_exit_Rpos_cutoff[magIdx][0]->Fill(hitR/mm);
                            }
                        }

                        //Energy
                        if (magnet_exit_energy[magIdx].find(PDG) !=
                            magnet_exit_energy[magIdx].end()) {
                            magnet_exit_energy[magIdx][PDG]->Fill(energy/MeV);
                        }
                        else {
                            magnet_exit_energy[magIdx][0]->Fill(energy/MeV);
                        }

                        if (hitR/mm < position_cutoffR) {
                            if (magnet_exit_cutoff_energy[magIdx].find(PDG) !=
                                magnet_exit_cutoff_energy[magIdx].end()) {
                                magnet_exit_cutoff_energy[magIdx][PDG]->Fill(energy/MeV);
                            }
                            else {
                                magnet_exit_cutoff_energy[magIdx][0]->Fill(energy/MeV);
                            }
                        }
                    }

                    /*
                    //Fill the TTree (TODO: Use correct buffer etc.)
                    if (not miniFile) {
                        targetExitBuffer.x = hitPos.x()/mm;
                        targetExitBuffer.y = hitPos.x()/mm;
                        targetExitBuffer.z = hitPos.z()/mm;

                        targetExitBuffer.px = momentum.x()/MeV;
                        targetExitBuffer.py = momentum.y()/MeV;
                        targetExitBuffer.pz = momentum.z()/MeV;

                        targetExitBuffer.E = energy / MeV;

                        targetExitBuffer.PDG = PDG;
                        targetExitBuffer.charge = charge;

                        targetExitBuffer.eventID = eventCounter;

                        targetExit->Fill();
                    }
                    */
                }
            }
            else {
                G4cout << "magnetExitposHitsCollection was NULL! for '" << magName << "'"<<G4endl;
            }
        }
        else {
            G4cout << "myMagnetExitpos_CollID was " << myMagnetExitpos_CollID << " < 0 for '" << magName << "'!"<<G4endl;
        }
    } // END loop over magnets
    if (not miniFile) {
        magnetEdeps->Fill(); // Outside loop over magnets
    }
}
void RootFileWriter::finalizeRootFile() {

    //Needed for some of the processing
    G4RunManager*           run  = G4RunManager::GetRunManager();
    DetectorConstruction* detCon = (DetectorConstruction*)run->GetUserDetectorConstruction();

    //Print out the particle types on all detector planes
    for (auto it : typeCounter) {
        PrintParticleTypes(it.second, it.first);
    }

    // ** Below cutoff **

    //Tracker average position and RMS
    double xave  = tracker_particleHit_x / ((double)typeCounter["tracker"].numParticles);
    double yave  = tracker_particleHit_y / ((double)typeCounter["tracker"].numParticles);
    double xrms  = ( tracker_particleHit_xx - (tracker_particleHit_x*tracker_particleHit_x / ((double)typeCounter["tracker"].numParticles)) ) /
        (((double)typeCounter["tracker"].numParticles)-1.0);
    xrms = sqrt(xrms);
    double yrms  = ( tracker_particleHit_yy - (tracker_particleHit_y*tracker_particleHit_y / ((double)typeCounter["tracker"].numParticles)) ) /
        (((double)typeCounter["tracker"].numParticles)-1.0);
    yrms = sqrt(yrms);

    //Exitangle
    G4double exitangle_avg=NAN;
    G4double exitangle_rms=NAN;
    if(detCon->GetHasTarget()) {
        exitangle_avg = target_exitangle / ((double)target_exitangle_numparticles);
        exitangle_rms = ( target_exitangle2 -
                            (target_exitangle*target_exitangle /
                            ((double)target_exitangle_numparticles)) ) /
                        ( (double)target_exitangle_numparticles - 1.0 ) ;
        exitangle_rms = sqrt(exitangle_rms);
    }

    G4cout << G4endl
           << "All particles (n=" << typeCounter["tracker"].numParticles << "):" << G4endl;

    G4cout << "Average x = " << xave << " [mm], RMS = " << xrms << " [mm]" << G4endl
           << "Average y = " << yave << " [mm], RMS = " << yrms << " [mm]" << G4endl;

    if (detCon->GetHasTarget()) {
        G4cout << G4endl
                << "Exit angle average (x) = " << exitangle_avg << " [deg]" << G4endl
                << "Exit angle RMS (x)     = " << exitangle_rms << " [deg]" << G4endl;
    }

    // ** Above cutoff **

    // Average position and RMS
    double xave_cutoff  = tracker_particleHit_x_cutoff / ((double)numParticles_cutoff);
    double yave_cutoff  = tracker_particleHit_y_cutoff / ((double)numParticles_cutoff);
    double xrms_cutoff  = ( tracker_particleHit_xx_cutoff -
                            (tracker_particleHit_x_cutoff*tracker_particleHit_x_cutoff /
                             ((double)numParticles_cutoff)) ) /
        (((double)numParticles_cutoff)-1.0);
    xrms_cutoff = sqrt(xrms_cutoff);
    double yrms_cutoff  = ( tracker_particleHit_yy_cutoff -
                            (tracker_particleHit_y_cutoff*tracker_particleHit_y_cutoff /
                             ((double)numParticles_cutoff)) ) /
        (((double)numParticles_cutoff)-1.0);
    yrms_cutoff = sqrt(yrms_cutoff);

    //Exitangle
    G4double exitangle_avg_cutoff = NAN;
    G4double exitangle_rms_cutoff = NAN;
    if(detCon->GetHasTarget()) {
        exitangle_avg_cutoff = target_exitangle_cutoff /
            ( (double)target_exitangle_cutoff_numparticles );
        exitangle_rms_cutoff =  ( target_exitangle2_cutoff -
                                    (target_exitangle_cutoff*target_exitangle_cutoff /
                                    ((double)target_exitangle_cutoff_numparticles)) ) /
                                ( (double)target_exitangle_cutoff_numparticles - 1.0 ) ;
        exitangle_rms_cutoff = sqrt(exitangle_rms_cutoff);
    }

    G4cout << G4endl
           << "Above cutoff (charged, energy > "
           << beamEnergy*beamEnergy_cutoff <<" [MeV], n=" << numParticles_cutoff
           << ") only:" << G4endl << G4endl;

    G4cout << "Average x = " << xave_cutoff << " [mm], RMS = " << xrms_cutoff << " [mm]" << G4endl
           << "Average y = " << yave_cutoff << " [mm], RMS = " << yrms_cutoff << " [mm]" << G4endl;

    if (detCon->GetHasTarget()) {
        G4cout << G4endl
            << "Exit angle average (x) = " << exitangle_avg_cutoff << " [deg]" << G4endl
            << "Exit angle RMS (x)     = " << exitangle_rms_cutoff << " [deg]" << G4endl;
        G4cout << G4endl;
    }

    //General metadata
    // Ugly hack: Use a double to store an int,
    // since there are no streamable int arrays without making a dict.
    G4cout << "** Metadata **" << G4endl;
    G4cout << "eventCounter  = " << eventCounter << G4endl;
    G4cout << "numEvents     = " << numEvents    << G4endl;
    if (detCon->GetHasTarget()) {
        G4cout << "targetDensity = " << detCon->GetTargetMaterialDensity()*cm3/g
                                     << " [g/cm^3]" << G4endl;
    }
    else {
        G4cout << "targetDensity = " << 0.0
                                     << " [g/cm^3]" << G4endl;
    }

    TVectorD metadataVector (3);
    metadataVector[0] = double(eventCounter);
    metadataVector[1] = double(numEvents);
    if (detCon->GetHasTarget()) {
        metadataVector[2] = detCon->GetTargetMaterialDensity()*cm3/g;
    }
    else {
        metadataVector[2] = 0.0;
    }
    metadataVector.Write("metadata");
    G4cout << G4endl;

    // Magnet metadata
    for (auto mag : detCon->magnets) {
        TVectorD magnetMetadataVector(1);
        magnetMetadataVector[0] = double(mag->GetTypicalDensity()*cm3/g);
        magnetMetadataVector.Write((mag->magnetName + "_metadata").c_str());
    }
    

    //Compute Twiss parameters
    PrintTwissParameters(init_phasespaceX);
    PrintTwissParameters(init_phasespaceY);
    if (detCon->GetHasTarget()) {
        PrintTwissParameters(target_exit_phasespaceX);
        PrintTwissParameters(target_exit_phasespaceY);
        PrintTwissParameters(target_exit_phasespaceX_cutoff);
        PrintTwissParameters(target_exit_phasespaceY_cutoff);
    }
    for (size_t magIdx = 0; magIdx < magnet_exit_phasespaceX.size(); magIdx++) {
        PrintTwissParameters(magnet_exit_phasespaceX[magIdx]);
        PrintTwissParameters(magnet_exit_phasespaceY[magIdx]);
        PrintTwissParameters(magnet_exit_phasespaceX_cutoff[magIdx]);
        PrintTwissParameters(magnet_exit_phasespaceY_cutoff[magIdx]);
    }
    PrintTwissParameters(tracker_phasespaceX);
    PrintTwissParameters(tracker_phasespaceY);
    PrintTwissParameters(tracker_phasespaceX_cutoff);
    PrintTwissParameters(tracker_phasespaceY_cutoff);

    if (anaScatterTest and detCon->GetHasTarget()) {
        // Compute the analytical multiple scattering angle distribution
        // Formulas from various sources:
        //
        // http://www-glast.slac.stanford.edu/software/AnaGroup/rev.pdf
        // Document describing the validation of Geant4 multiple scattering
        //
        // http://cdsweb.cern.ch/record/1279627/files/PH-EP-Tech-Note-2010-013.pdf
        // (note: Missing factor A in radiation length formula)
        //
        // https://en.wikipedia.org/wiki/Radiation_length
        // Checked: April 15th 2018
        //
        // http://pdg.lbl.gov/2004/reviews/passagerpp.pdf
        // This seems to be the base reference for the "compact fit to the data"
        // used for the radiation length, which is the formula found
        // in Wikipedia and in PH-EP-Tech-Note-2010-013.
        //
        // http://pdg.lbl.gov/2017/reviews/rpp2017-rev-passage-particles-matter.pdf
        // Multiple scattering formulas
        //
        // http://pdg.lbl.gov/2017/AtomicNuclearProperties/index.html
        // Tables used to crosscheck computed radiation lenghts

        G4cout << G4endl
               << "Computing analytical scattering..." << G4endl;

        G4double targetThickness = detCon->getTargetThickness(); //Geant units
        G4cout << "targetThickness = " << targetThickness/mm << " [mm]" << G4endl;
        G4int targetZ = detCon->GetTargetMaterialZ();
        G4cout << "targetZ         = " << targetZ << " [e]" << G4endl;
        G4double targetA = detCon->GetTargetMaterialA();
        G4cout << "targetA         = " << targetA << " [amu]" << G4endl;

        G4double targetDensity = detCon->GetTargetMaterialDensity();
        G4cout << "targetDensity   = " << targetDensity*cm3/g << " [g/cm]" << G4endl;

        G4double radiationLength = 716.4 * targetA /
            (targetZ*(targetZ+1)*log(287.0/sqrt(targetZ))); //[g/cm^2]
        G4cout << "RadiationLength = " << radiationLength << " [g/cm^2]" << G4endl;
        radiationLength /= targetDensity*cm3/g; //[cm]
        radiationLength *= cm; //Geant4 units
        G4cout << "                = " << radiationLength/cm << " [cm]" << G4endl;

        PrimaryGeneratorAction* genAct = (PrimaryGeneratorAction*)run->GetUserPrimaryGeneratorAction();
        G4double beamMass   = genAct->get_beam_particlemass(); // Geant4 units
        G4cout << "beamMass        = " << beamMass / MeV << " [MeV/c^2]" <<G4endl;
        G4double beamCharge = genAct->get_beam_particlecharge(); // [e]
        G4cout << "beamCharge      = " << beamCharge << " [e]" << G4endl;
        G4double beamBeta2   = 1 - pow(beamMass/beamEnergy*MeV, 2);
        G4cout << "beamBeta2       = " << beamBeta2 << G4endl;

        G4double theta0 = 13.6*MeV / (beamEnergy*MeV - pow(beamMass,2)/beamEnergy*MeV) *
            fabs(beamCharge) * sqrt(targetThickness/radiationLength) *
            (1 + 0.038 * log(targetThickness*pow(beamCharge,2)/(radiationLength*beamBeta2)) );
        //TF1* highland = new TF1("highland", "Highland scattering distribution",0,90);
        G4cout << "theta0          = " << theta0 << " [rad]" << G4endl
               << "                = " << theta0/deg << " [deg]" << G4endl;
        G4cout << G4endl;

        target_exitangle_hist_cutoff->Write(); // Write the unaltered histogram

        //Plot the analytical scattering distribution and the histogram together
        TCanvas* c1 = new TCanvas("scatterPlot");
        target_exitangle_hist_cutoff->GetXaxis()->SetRangeUser(-3*theta0/deg,3*theta0/deg);
        target_exitangle_hist_cutoff->GetXaxis()->SetTitle("Angle [deg]");
        target_exitangle_hist_cutoff->Scale(1.0/target_exitangle_hist_cutoff->Integral(), "width");
        target_exitangle_hist_cutoff->Draw();

        G4double exitangle_analytic_xmin = target_exitangle_hist_cutoff->GetXaxis()->GetXmin();
        G4double exitangle_analytic_xmax = target_exitangle_hist_cutoff->GetXaxis()->GetXmax();
        G4double exitangle_analytic_xrange = exitangle_analytic_xmax - exitangle_analytic_xmin;
        int exitangle_analytic_npoints = target_exitangle_hist_cutoff->GetNbinsX();
        G4double* exitangle_analytic_x = new G4double[exitangle_analytic_npoints];
        G4double* exitangle_analytic_y = new G4double[exitangle_analytic_npoints];

        for (int i = 0; i < exitangle_analytic_npoints; i++) {
            exitangle_analytic_x[i] =
                exitangle_analytic_xrange/((G4double)exitangle_analytic_npoints-1) * i
                + exitangle_analytic_xmin;

            exitangle_analytic_y[i] = 1.0/sqrt(2.0*M_PI*pow(theta0/deg,2)) *
                exp( - pow(exitangle_analytic_x[i]/(theta0/deg),2) / 2.0);
        }
        TGraph* exitangle_analytic     = new TGraph(exitangle_analytic_npoints,
                                                    exitangle_analytic_x,
                                                    exitangle_analytic_y       );
        exitangle_analytic->Draw("same");
        /*
        // Debug
        G4cout << "integral = " << exitangle_analytic->Integral() << G4endl;
        G4cout << "integral = " << target_exitangle_hist_cutoff->Integral("width") << G4endl;
        */

        G4String plot_filename = foldername_out + "/" + filename_out + "_angles.png";
        //c1->SaveAs(plot_filename);
        c1->Write();

        // Write the ROOT file.
        exitangle_analytic->Write();
        delete exitangle_analytic; exitangle_analytic = NULL;
    }


    if (not quickmode) {
        //Write the 2D histograms to the ROOT file (slow)
        G4cout << "Writing 2D histograms..." << G4endl;
        init_phasespaceX->Write();
        init_phasespaceY->Write();
        init_phasespaceXY->Write();

        if (detCon->GetHasTarget()) {
            target_exit_phasespaceX->Write();
            target_exit_phasespaceY->Write();

            target_exit_phasespaceX_cutoff->Write();
            target_exit_phasespaceY_cutoff->Write();

            if (target_edep_rdens != NULL) {
                target_edep_rdens->Write();
            }
        }

        tracker_hitPos->Write();
        tracker_hitPos_cutoff->Write();

        tracker_phasespaceX->Write();
        tracker_phasespaceY->Write();

        tracker_phasespaceX_cutoff->Write();
        tracker_phasespaceY_cutoff->Write();

        for (auto it : magnet_edep_rdens) {
            it->Write();
        }
        for (auto it : magnet_exit_phasespaceX) {
            it->Write();
        }
        for (auto it : magnet_exit_phasespaceY) {
            it->Write();
        }
        for (auto it : magnet_exit_phasespaceX_cutoff) {
            it->Write();
        }
        for (auto it : magnet_exit_phasespaceY_cutoff) {
            it->Write();
        }

        // Write the 3D histograms to the root file (slower)
        G4cout << "Writing 3D histograms..." << G4endl;
        if (detCon->GetHasTarget()) {
            if (target_edep_dens != NULL) {
                target_edep_dens->Write();
            }
        }
        if (this->edep_dens_dz > 0.0 ) {
            for (auto it : magnet_edep_dens) {
                it->Write();
            }
        }

    }

    if (not miniFile) {
        G4cout << "Writing TTrees..." << G4endl;

        if (detCon->GetHasTarget()) {
            targetExit->Write();
        }
        trackerHits->Write();

        magnetEdeps->Write();
        delete magnetEdeps;
        magnetEdeps=NULL;

        if (magnetEdepsBuffer != NULL) {
            delete magnetEdepsBuffer;
            magnetEdepsBuffer = NULL;
        }
    }

    G4cout << "Writing 1D histograms..." << G4endl;

    init_E->Write();

    if (detCon->GetHasTarget()) {
        targetEdep->Write();
        targetEdep_NIEL->Write();
        targetEdep_IEL->Write();

        // (Loops over particle types)
        for (auto it : target_exit_energy) {
            it.second->Write();
            delete it.second;
        }
        target_exit_energy.clear();
        for (auto it : target_exit_cutoff_energy) {
            it.second->Write();
            delete it.second;
        }
        target_exit_cutoff_energy.clear();
        for (auto it: target_exit_Rpos) {
            it.second->Write();
            delete it.second;
        }
        target_exit_Rpos.clear();
        for (auto it: target_exit_Rpos_cutoff) {
            it.second->Write();
            delete it.second;
        }
        target_exit_Rpos_cutoff.clear();
    }

    for (auto it : tracker_type_energy) {
        it.second->Write();
        delete it.second;
    }
    tracker_type_energy.clear();
    for (auto it : tracker_type_cutoff_energy) {
        it.second->Write();
        delete it.second;
    }
    tracker_type_cutoff_energy.clear();

    for (auto it: tracker_Rpos) {
        it.second->Write();
        delete it.second;
    }
    tracker_Rpos.clear();
    for (auto it: tracker_Rpos_cutoff) {
        it.second->Write();
        delete it.second;
    }
    tracker_Rpos_cutoff.clear();

    if (detCon->GetHasTarget()) {
        target_exitangle_hist->Write();
    }

    // Clear magnet 3D hists
    if (this->edep_dens_dz > 0.0 ) {
        for (auto it : magnet_edep_dens) {
            delete it;
        }
    }
    magnet_edep_dens.clear();
    // Clear magnet 2D hists
    for (auto it : magnet_edep_rdens) {
        delete it;
    }
    magnet_edep_rdens.clear();
    for (auto it : magnet_exit_phasespaceX) {
        delete it;
    }
    magnet_exit_phasespaceX.clear();
    for (auto it : magnet_exit_phasespaceY) {
        delete it;
    }
    magnet_exit_phasespaceY.clear();
    for (auto it : magnet_exit_phasespaceX_cutoff) {
        delete it;
    }
    magnet_exit_phasespaceX_cutoff.clear();
    for (auto it : magnet_exit_phasespaceY_cutoff) {
        delete it;
    }
    magnet_exit_phasespaceY_cutoff.clear();

    // Write and clear magnet 1D hists
    for (auto it : magnet_edep) {
        it->Write();
        delete it;
    }
    magnet_edep.clear();

    for (auto mag : magnet_exit_Rpos) {
        for (auto PDG : mag) {
            PDG.second->Write();
            delete PDG.second;
        }
        mag.clear();
    }
    magnet_exit_Rpos.clear();

    for (auto mag : magnet_exit_Rpos_cutoff) {
        for (auto PDG : mag) {
            PDG.second->Write();
            delete PDG.second;
        }
        mag.clear();
    }
    magnet_exit_Rpos_cutoff.clear();

    for (auto mag : magnet_exit_energy) {
        for (auto PDG : mag) {
            PDG.second->Write();
            delete PDG.second;
        }
        mag.clear();
    }
    magnet_exit_energy.clear();

    for (auto mag : magnet_exit_cutoff_energy) {
        for (auto PDG : mag) {
            PDG.second->Write();
            delete PDG.second;
        }
        mag.clear();
    }
    magnet_exit_cutoff_energy.clear();

    //Now that we have plotted, delete stuff

    delete init_phasespaceX; init_phasespaceX = NULL;
    delete init_phasespaceY; init_phasespaceY = NULL;
    delete init_phasespaceXY; init_phasespaceXY = NULL;

    if (detCon->GetHasTarget()) {
        delete targetEdep; targetEdep = NULL;
        delete targetEdep_NIEL; targetEdep_NIEL = NULL;
        delete targetEdep_IEL; targetEdep_IEL = NULL;

        delete target_exitangle_hist; target_exitangle_hist = NULL;

        delete target_exit_phasespaceX; target_exit_phasespaceX = NULL;
        delete target_exit_phasespaceY; target_exit_phasespaceY = NULL;

        delete target_exit_phasespaceX_cutoff; target_exit_phasespaceX_cutoff = NULL;
        delete target_exit_phasespaceY_cutoff; target_exit_phasespaceY_cutoff = NULL;
    }

    delete tracker_phasespaceX; tracker_phasespaceX = NULL;
    delete tracker_phasespaceY; tracker_phasespaceY = NULL;

    delete tracker_phasespaceX_cutoff; tracker_phasespaceX_cutoff = NULL;
    delete tracker_phasespaceY_cutoff; tracker_phasespaceY_cutoff = NULL;

    if (detCon->GetHasTarget()) {
        if (target_edep_dens != NULL) {
            delete target_edep_dens; target_edep_dens = NULL;
        }
        if (target_edep_rdens != NULL) {
            delete target_edep_rdens; target_edep_rdens = NULL;
        }

        delete target_exitangle_hist_cutoff; target_exitangle_hist_cutoff = NULL;
    }

    delete tracker_numParticles; tracker_numParticles = NULL;
    delete tracker_energy; tracker_energy = NULL;
    delete tracker_hitPos; tracker_hitPos = NULL;
    delete tracker_hitPos_cutoff; tracker_hitPos_cutoff = NULL;

    if(not miniFile) {
        delete trackerHits; trackerHits = NULL;
        if (detCon->GetHasTarget()) {
            delete targetExit; targetExit = NULL;
        }
    }

    histFile->Write();
    histFile->Close();
    delete histFile; histFile = NULL;
}

void RootFileWriter::PrintTwissParameters(TH2D* phaseSpaceHist) {
    G4cout << "Stats for '" << phaseSpaceHist->GetTitle() << "':"  << G4endl;
    double stats[7];
    phaseSpaceHist->GetStats(stats);

    // Fill used [mm] and [rad]
    double posAve   = stats[2]/stats[0];
    double angAve   = stats[4]/stats[0];
    double posVar   = (stats[3] - stats[2]*stats[2]/stats[0]) / (stats[0]-1.0) ;
    double angVar   = (stats[5] - stats[4]*stats[4]/stats[0]) / (stats[0]-1.0) ;
    double coVar    = (stats[6] - stats[2]*stats[4]/stats[0]) / (stats[0]-1.0);

    G4cout << "numHits = "  << stats[0]
           << ", posAve = " << posAve     << " [mm]"
           << ", angAve = " << angAve     << " [rad]"
           << ", posVar = " << posVar     << " [mm^2]"
           << ", angVar = " << angVar     << " [rad^2]"
           << ", coVar  = " << coVar      << " [rad*mm]"
           << G4endl;

    G4RunManager*           run    = G4RunManager::GetRunManager();
    PrimaryGeneratorAction* genAct = (PrimaryGeneratorAction*)run->GetUserPrimaryGeneratorAction();
    double gamma_rel = 1 + ( genAct->get_beam_energy() * MeV / genAct->get_beam_particlemass() );
    double beta_rel = sqrt(gamma_rel*gamma_rel - 1.0) / gamma_rel;

    double det = posVar*angVar - coVar*coVar;
    double epsG = sqrt(det);
    double epsN = epsG*beta_rel*gamma_rel;  //[mm*rad]
    double beta = posVar/epsG; // [mm]
    double alpha = -coVar/epsG;

    G4cout << "Geometrical emittance          = " << epsG*1e3 << " [um]" << G4endl;
    G4cout << "Normalized emittance           = " << epsN*1e3 << " [um]"
           << ", assuming beam kinetic energy = " << genAct->get_beam_energy() << " [MeV]"
           << ", and mass = " << genAct->get_beam_particlemass()/MeV << " [MeV/c^2]"
           << G4endl;
    G4cout << "Twiss beta  = " << beta*1e-3  << " [m]" << G4endl
           << "Twiss alpha = " << alpha << " [-]"  << G4endl;

    G4cout << G4endl;

    // Write to root file
    TVectorD twissVector (8);
    twissVector[0] = epsN*1e3;
    twissVector[1] = beta*1e-3;
    twissVector[2] = alpha;
    twissVector[3] = posAve;
    twissVector[4] = angAve;
    twissVector[5] = posVar;
    twissVector[6] = angVar;
    twissVector[7] = coVar;
    twissVector.Write((G4String(phaseSpaceHist->GetName())+"_TWISS").c_str());
}

void RootFileWriter::PrintParticleTypes(particleTypesCounter& pt, G4String name) {
    //Print out the particle types hitting the tracker
    G4cout << endl;
    G4cout << "Got types at " << name << ":" << G4endl;
    if (pt.numParticles == 0){
        G4cout << "No particles hit!" << G4endl;
        return;
    }

    TVectorD particleTypes_PDG    (pt.particleTypes.size());
    TVectorD particleTypes_numpart(pt.particleTypes.size());

    size_t particleTypes_i = 0;
    for(std::map<G4int,G4int>::iterator it = pt.particleTypes.begin(); it != pt.particleTypes.end(); it++){
        G4cout << std::setw(15) << it->first << " = "
               << std::setw(15) << pt.particleNames[it->first] << ": "
               << std::setw(15) << it->second << " = ";// << G4endl;
        G4int numHashes = (G4int) ((it->second / ((double)pt.numParticles)) * 100);
        for (int i = 0; i < numHashes; i++) G4cout << "#";
        G4cout << endl;

        //Also put them in the ROOT file
        // Unfortunately, there is no TObject array type for ints (?!?),
        // and I don't want  to depend on a ROOT dictionary file.
        particleTypes_PDG     [particleTypes_i] = int(it->first);
        particleTypes_numpart [particleTypes_i] = int(it->second);

        particleTypes_i++;
    }
    particleTypes_PDG.Write((name + "_ParticleTypes_PDG").c_str());
    particleTypes_numpart.Write((name + "_ParticleTypes_numpart").c_str());

}

void RootFileWriter::FillParticleTypes(particleTypesCounter& pt, G4int PDG, G4String type) {
    if (pt.particleTypes.count(PDG) == 0) {
        pt.particleTypes[PDG] = 0;
        pt.particleNames[PDG] = type;
    }
    pt.particleTypes[PDG] += 1;
    pt.numParticles       += 1;
}

void RootFileWriter::setEngNbins(G4int edepNbins_in) {
    if (edepNbins_in > 0) {
        this->engNbins = edepNbins_in;
    }
    else if (edepNbins_in < 0) {
        G4cerr << "Error: edepNbins must be > 0 (or 0 for auto)" << G4endl;
        exit(1);
    }
}
