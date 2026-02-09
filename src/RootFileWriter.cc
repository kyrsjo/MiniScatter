
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
#include "MiniScatterVersion.hh"

#include "TNamed.h"

#include "TGraph.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TBranch.h"

#include "TRandom1.h"

#include "EdepHit.hh"
#include "TrackerHit.hh"

#include "G4SDManager.hh"

#include "G4Track.hh"
#include "G4RunManager.hh"

#include "DetectorConstruction.hh"
#include "MagnetClasses.hh"
#include "VirtualTrackerWorldConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4SystemOfUnits.hh"

#include "G4Exception.hh"

#include <iostream>
#include <iomanip>

#include <unistd.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
struct stat stat_info;

using namespace std;
RootFileWriter* RootFileWriter::singleton = NULL;

void RootFileWriter::initializeRootFile(){
    G4RunManager*           run    = G4RunManager::GetRunManager();
    DetectorConstruction*   detCon = (DetectorConstruction*)run->GetUserDetectorConstruction();
    PrimaryGeneratorAction* genAct = (PrimaryGeneratorAction*)run->GetUserPrimaryGeneratorAction();
    this->beamEnergy = genAct->get_beam_energy();

    //Count all particles that are Fill'ed for the stats used to compute the twiss parameters,
    // even if they are outside the phasespacehist_posLim / phaspacehist_angLim.
    TH2D::StatOverflows(true);

    if (not has_filename_out) {
        G4String errormessage = "filename_out not set";
        G4Exception("RootFileWriter::initializeRootFile()", "MSRootFile1000",FatalException,errormessage);
    }
    rootFileName = foldername_out + "/" + filename_out + ".root";
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
    // #error Only UNIX is supported.
    #endif
    if(stat(foldername_out.c_str(), &stat_info) != 0) {
        if (errno == ENOENT) {
            G4cout << "Creating folder '" << foldername_out << "'" << G4endl;

            G4String foldername_out_full;
            if (foldername_out.c_str()[0] != '/') {

		#ifdef __APPLE__
		#include <unistd.h>
		#include <limits.h>

		char cwd[PATH_MAX];
		if (getcwd(cwd, sizeof(cwd)) == nullptr) {
   		 perror("getcwd failed");
		}

		#else
		char* cwd = get_current_dir_name();
		//#endif

                //char* cwd = get_current_dir_name();
                if (cwd == NULL) {
                    G4String errormessage = "Error getting the current path";
                    G4Exception("RootFileWriter::initializeRootFile()", "MSRootFile1001",FatalException,errormessage);
                }
		#endif

                if (cwd[strlen(cwd)-1] == '/') {
                    foldername_out_full = G4String(cwd) + foldername_out;
                }
                else {
                    foldername_out_full = G4String(cwd) + G4String("/") + foldername_out;
                }
		
		#ifndef __APPLE__
                delete[] cwd;
		#endif

                G4cout << "Converted relative path '" << foldername_out << "' to absolute path '" 
                       << foldername_out_full << "'." << G4endl;
            }
            else {
                foldername_out_full = foldername_out;
            }
            const char* path_full = foldername_out_full.c_str();

            if( strlen(path_full) < 2 ) {
                G4String errormessage = "Path string must be more than '/' (and why is '/' not existing?), got '" + G4String(path_full) + "'";
                G4Exception("RootFileWriter::initializeRootFile()", "MSRootFile1002",FatalException,errormessage);
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
            G4String errormessage = "Could not lookup folder '" + foldername_out + "'";
            G4Exception("RootFileWriter::initializeRootFile()", "MSRootFile1003",FatalException,errormessage);
        }
    }
    else if(not S_ISDIR(stat_info.st_mode)) {
        G4String errormessage = "A non-folder entity named '" + foldername_out + "' already exist";
        G4Exception("RootFileWriter::initializeRootFile()", "MSRootFile1004",FatalException,errormessage);
    }
    else {
        G4cout << "Folder '" << foldername_out << "' already exists -- using it!" << G4endl;
    }

    G4cout << "Opening ROOT file '" + rootFileName +"'"<<G4endl;
    histFile = new TFile(rootFileName,"RECREATE");
    if ( not histFile->IsOpen() ) {
        G4String errormessage = "Opening TFile '" + rootFileName + "' failed.";
        G4Exception("RootFileWriter::initializeRootFile()", "MSRootFile1005",FatalException,errormessage);
    }
    G4cout << G4endl;

    // Add some metadata to the newly opened file
    TNamed* version_number = new TNamed("MiniScatter_version_number", miniscatter_version);
    version_number->Write();
    delete version_number;

    eventCounter = 0;

    // Limit for radial histograms
    G4double minR = min(detCon->getWorldSizeX(),detCon->getWorldSizeY())/mm;

    RNG = new TRandom1((UInt_t) rngSeed);

    // TTrees for external analysis
    if (not miniFile) {
        if (detCon->GetHasTarget()) {
            targetExit = new TTree("TargetExit","TargetExit tree");
            targetExit->Branch("TargetExitBranch", &targetExitBuffer,
                               "x/D:y:z:px:py:pz:E:PDG/I:charge:eventID");
        }

        initParts = new TTree("InitParts","InitParts tree");
        initParts->Branch("InitPartsBranch",&initPartsBuffer,
                            "x/D:y:z:px:py:pz:E:PDG/I:charge:eventID");

        trackerHits = new TTree("TrackerHits","TrackerHits tree");
        trackerHits->Branch("TrackerHitsBranch", &trackerHitsBuffer,
                            "x/D:y:z:px:py:pz:E:PDG/I:charge:eventID");

        // For magnets, branches are created elsewhere, one per magnet.
        magnetEdeps = new TTree("magnetEdeps", "Magnet Edeps tree");

        // For magnetExit phase space data Tree and Branches are created elsewhere
    }

    // Target energy deposition
    if (detCon->GetHasTarget()) {
        targetEdep = new TH1D("targetEdep","targetEdep",engNbins+1,0,beamEnergy*(1+1/engNbins));
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
        target_exitangle_hist           = new TH1D("target_exit_angle",
						   "Exit angle from target",
						   5001, -90, 90);
        target_exitangle_hist_cutoff    = new TH1D("target_exit_angle_cutoff",
						   "Exit angle from target (charged, energy > Ecut, r < Rcut)",
						   5001, -90, 90);

        // Target exit phasespace histograms
        target_exit_phasespaceX         = new TH2D("target_exit_x",
						   "Target exit phase space (x)",
						   1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
						   1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        target_exit_phasespaceX->GetXaxis()->SetTitle("Position x [mm]");
        target_exit_phasespaceX->GetYaxis()->SetTitle("Angle dx/dz [rad]");

        target_exit_phasespaceY         = new TH2D("target_exit_y",
						   "Target exit phase space (y)",
						   1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
						   1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        target_exit_phasespaceY->GetXaxis()->SetTitle("Position y [mm]");
        target_exit_phasespaceY->GetYaxis()->SetTitle("Angle dy/dz [rad]");

        target_exit_phasespaceX_cutoff  = new TH2D("target_exit_cutoff_x",
						   "Target exit phase space (x) (charged, energy > Ecut, r < Rcut)",
						   1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
						   1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        target_exit_phasespaceX_cutoff->GetXaxis()->SetTitle("Position x [mm]");
        target_exit_phasespaceX_cutoff->GetYaxis()->SetTitle("Angle dx/dz [rad]");

        target_exit_phasespaceY_cutoff  = new TH2D("target_exit_cutoff_y",
						   "Target exit phase space (y) (charged, energy > Ecut, r < Rcut)",
						   1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
						   1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        target_exit_phasespaceY_cutoff->GetXaxis()->SetTitle("Position y [mm]");
        target_exit_phasespaceY_cutoff->GetYaxis()->SetTitle("Angle dy/dz [rad]");

        target_exit_phasespaceXY        = new TH2D("target_exit_xy",
						   "Target exit phase space (x,y)",
						   1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
						   1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm);
        target_exit_phasespaceXY->GetXaxis()->SetTitle("Position x [mm]");
        target_exit_phasespaceXY->GetYaxis()->SetTitle("Position y [mm]");
        target_exit_phasespaceXY_cutoff = new TH2D("target_exit_cutoff_xy",
						   "Target exit phase space (x,y) (charged, energy > Ecut, r < Rcut)",
						   1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
						   1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm);
        target_exit_phasespaceXY_cutoff->GetXaxis()->SetTitle("Position x [mm]");
        target_exit_phasespaceXY_cutoff->GetYaxis()->SetTitle("Position y [mm]");

    }

    // Tracker histograms
    VirtualTrackerWorldConstruction* traCon = VirtualTrackerWorldConstruction::getInstance();
    for (int idx = 0; idx < traCon->getNumTrackers(); idx++) {

        std::string trackerName;
        if (traCon->getNumTrackers() == 1) { trackerName = "tracker"; }
        else                               { trackerName = std::string("tracker_") + std::to_string(idx+1); }

        typeCounter[trackerName]             = particleTypesCounter();
        typeCounter[trackerName + "_cutoff"] = particleTypesCounter();

        tracker_numParticles.push_back(new TH1D((trackerName+"_numParticles").c_str(),(trackerName+" numParticles").c_str(),1001,-0.5,1000.5));
        tracker_numParticles.back()->GetXaxis()->SetTitle("Number of particles / event");

        tracker_energy.push_back( new TH1D((trackerName+"_energy").c_str(),("Energy of all particles hitting "+trackerName).c_str(),10000,0,beamEnergy));
        tracker_energy.back()->GetXaxis()->SetTitle("Energy per particle [MeV]");

        tracker_type_energy.push_back(std::map<G4int,TH1D*>());
        tracker_type_energy.back()[11]  = new TH1D((trackerName+"_energy_PDG11").c_str(),
                                                   ("Particle energy when hitting "+trackerName+" (electrons)").c_str(),
                                                   engNbins,0,beamEnergy);
        tracker_type_energy.back()[-11] = new TH1D((trackerName+"_energy_PDG-11").c_str(),
                                                   ("Particle energy when hitting "+trackerName+" (positrons)").c_str(),
                                                   engNbins,0,beamEnergy);
        tracker_type_energy.back()[22]  = new TH1D((trackerName+"_energy_PDG22").c_str(),
                                                   ("Particle energy when hitting "+trackerName+" (photons)").c_str(),
                                                   engNbins,0,beamEnergy);
        tracker_type_energy.back()[2212]= new TH1D((trackerName+"_energy_PDG2212").c_str(),
                                                   ("Particle energy when hitting "+trackerName+" (protons)").c_str(),
                                                   engNbins,0,beamEnergy);
        tracker_type_energy.back()[0]   = new TH1D((trackerName+"_energy_PDGother").c_str(),
                                                   ("Particle energy when hitting "+trackerName+" (other)").c_str(),
                                                   engNbins,0,beamEnergy);
        for (auto it : tracker_type_energy.back()) {
            it.second->GetXaxis()->SetTitle("Energy [MeV]");
        }

        tracker_type_cutoff_energy.push_back(std::map<G4int,TH1D*>());
        tracker_type_cutoff_energy.back()[11]  = new TH1D((trackerName+"_cutoff_energy_PDG11").c_str(),
                                                          ("Particle energy when hitting "+trackerName+" (electrons) (r < Rcut)").c_str(),
                                                          engNbins,0,beamEnergy);
        tracker_type_cutoff_energy.back()[-11] = new TH1D((trackerName+"_cutoff_energy_PDG-11").c_str(),
                                                          ("Particle energy when hitting "+trackerName+" (positrons) (r < Rcut)").c_str(),
                                                          engNbins,0,beamEnergy);
        tracker_type_cutoff_energy.back()[22]  = new TH1D((trackerName+"_cutoff_energy_PDG22").c_str(),
                                                          ("Particle energy when hitting "+trackerName+" (photons) (r < Rcut)").c_str(),
                                                          engNbins,0,beamEnergy);
        tracker_type_cutoff_energy.back()[2212]= new TH1D((trackerName+"_cutoff_energy_PDG2212").c_str(),
                                                          ("Particle energy when hitting "+trackerName+" (protons) (r < Rcut)").c_str(),
                                                          engNbins,0,beamEnergy);
        tracker_type_cutoff_energy.back()[0]   = new TH1D((trackerName+"_cutoff_energy_PDGother").c_str(),
                                                          ("Particle energy when hitting "+trackerName+" (other) (r < Rcut)").c_str(),
                                                          engNbins,0,beamEnergy);
        for (auto it : tracker_type_cutoff_energy.back()) {
            it.second->GetXaxis()->SetTitle("Energy [MeV]");
        }

        tracker_phasespaceX.push_back(
            new TH2D((trackerName+"_x").c_str(),
                     (trackerName+" phase space (x)").c_str(),
                     1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                     1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad) );
        tracker_phasespaceX.back()->GetXaxis()->SetTitle("X [mm]");
        tracker_phasespaceX.back()->GetYaxis()->SetTitle("X' [rad]");

        tracker_phasespaceY.push_back(
            new TH2D((trackerName+"_y").c_str(),
                     (trackerName+" phase space (y)").c_str(),
                    1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                    1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad) );
        tracker_phasespaceY.back()->GetXaxis()->SetTitle("Y [mm]");
        tracker_phasespaceY.back()->GetYaxis()->SetTitle("Y' [rad]");

        tracker_phasespaceXY.push_back(
            new TH2D((trackerName+"_xy").c_str(),
                     (trackerName+" phase space (x,y)").c_str(),
                     1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                     1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm) );
        tracker_phasespaceXY.back()->GetXaxis()->SetTitle("X [mm]");
        tracker_phasespaceXY.back()->GetYaxis()->SetTitle("Y [mm]");

        tracker_phasespaceX_cutoff.push_back(
            new TH2D((trackerName+"_cutoff_x").c_str(),
                     (trackerName+" phase space (x) (charged, energy > Ecut, r < Rcut)").c_str(),
                     1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                     1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad) );
        tracker_phasespaceX_cutoff.back()->GetXaxis()->SetTitle("X [mm]");
        tracker_phasespaceX_cutoff.back()->GetYaxis()->SetTitle("X' [rad]");
        
        tracker_phasespaceY_cutoff.push_back(
            new TH2D((trackerName+"_cutoff_y").c_str(),
                     (trackerName+" phase space (y) (charged, energy > Ecut, r < Rcut)").c_str(),
                     1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                     1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad) );
        tracker_phasespaceY_cutoff.back()->GetXaxis()->SetTitle("Y [mm]");
        tracker_phasespaceY_cutoff.back()->GetYaxis()->SetTitle("Y' [rad]");

        tracker_phasespaceXY_cutoff.push_back(
            new TH2D((trackerName+"_cutoff_xy").c_str(),
                     (trackerName+" phase space (x,y) (charged, energy > Ecut, r < Rcut)").c_str(),
                     1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                     1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm) );
        tracker_phasespaceXY_cutoff.back()->GetXaxis()->SetTitle("X [mm]");
        tracker_phasespaceXY_cutoff.back()->GetYaxis()->SetTitle("Y [mm]");
	
        tracker_phasespaceX_cutoff_PDG.push_back(std::map<G4int,TH2D*>());
        tracker_phasespaceY_cutoff_PDG.push_back(std::map<G4int,TH2D*>());
        tracker_phasespaceX_cutoff_PDG.back()[11]  = new TH2D((trackerName+"_cutoff_x_PDG11").c_str(),
                                                                  (trackerName+" phase space (x) (electrons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        tracker_phasespaceX_cutoff_PDG.back()[-11]  = new TH2D((trackerName+"_cutoff_x_PDG-11").c_str(),
                                                                  (trackerName+" phase space (x) (positrons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        tracker_phasespaceX_cutoff_PDG.back()[22]  = new TH2D((trackerName+"_cutoff_x_PDG22").c_str(),
                                                                  (trackerName+" phase space (x) (photons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        tracker_phasespaceX_cutoff_PDG.back()[2212]  = new TH2D((trackerName+"_cutoff_x_PDG2212").c_str(),
                                                                  (trackerName+" phase space (x) (protons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        tracker_phasespaceX_cutoff_PDG.back()[0]  = new TH2D((trackerName+"_cutoff_x_PDGother").c_str(),
                                                                  (trackerName+" phase space (x) (other, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        for (auto PDG : tracker_phasespaceX_cutoff_PDG.back()) {
            PDG.second->GetXaxis()->SetTitle("X [mm]");
            PDG.second->GetYaxis()->SetTitle("X' [rad]");
        }
	
        tracker_phasespaceY_cutoff_PDG.back()[11]  = new TH2D((trackerName+"_cutoff_y_PDG11").c_str(),
                                                                  (trackerName+" phase space (y) (electrons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        tracker_phasespaceY_cutoff_PDG.back()[-11]  = new TH2D((trackerName+"_cutoff_y_PDG-11").c_str(),
                                                                  (trackerName+" phase space (y) (positrons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        tracker_phasespaceY_cutoff_PDG.back()[22]  = new TH2D((trackerName+"_cutoff_y_PDG22").c_str(),
                                                                  (trackerName+" phase space (y) (photons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        tracker_phasespaceY_cutoff_PDG.back()[2212]  = new TH2D((trackerName+"_cutoff_y_PDG2212").c_str(),
                                                                  (trackerName+" phase space (y) (protons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        tracker_phasespaceY_cutoff_PDG.back()[0]  = new TH2D((trackerName+"_cutoff_y_PDGother").c_str(),
                                                                  (trackerName+" phase space (y) (other, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        for (auto PDG : tracker_phasespaceY_cutoff_PDG.back()) {
            PDG.second->GetXaxis()->SetTitle("Y [mm]");
            PDG.second->GetYaxis()->SetTitle("Y' [rad]");
        }

        tracker_phasespaceXY_cutoff_PDG.push_back(std::map<G4int,TH2D*>());
        tracker_phasespaceXY_cutoff_PDG.back()[11]  = new TH2D((trackerName+"_cutoff_xy_PDG11").c_str(),
                                                                  (trackerName+" phase space (x,y) (electrons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        tracker_phasespaceXY_cutoff_PDG.back()[-11]  = new TH2D((trackerName+"_cutoff_xy_PDG-11").c_str(),
                                                                  (trackerName+" phase space (x,y) (positrons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        tracker_phasespaceXY_cutoff_PDG.back()[22]  = new TH2D((trackerName+"_cutoff_xy_PDG22").c_str(),
                                                                  (trackerName+" phase space (x,y) (photons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        tracker_phasespaceXY_cutoff_PDG.back()[2212]  = new TH2D((trackerName+"_cutoff_xy_PDG2212").c_str(),
                                                                  (trackerName+" phase space (x,y) (protons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        tracker_phasespaceXY_cutoff_PDG.back()[0]  = new TH2D((trackerName+"_cutoff_xy_PDGother").c_str(),
                                                                  (trackerName+" phase space (x,y) (other, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        for (auto PDG : tracker_phasespaceXY_cutoff_PDG.back()) {
            PDG.second->GetXaxis()->SetTitle("X [mm]");
            PDG.second->GetYaxis()->SetTitle("Y [mm]");
        }

        // Tracker R position
        tracker_Rpos.push_back(std::map<G4int,TH1D*>());
        tracker_Rpos.back()[11]  = new TH1D((trackerName+"_rpos_PDG11").c_str(),
                                            (trackerName+" rpos (electrons)").c_str(),
                                            1000,0,minR);
        tracker_Rpos.back()[-11] = new TH1D((trackerName+"_rpos_PDG-11").c_str(),
                                            (trackerName+" rpos (positrons)").c_str(),
                                            1000,0,minR);
        tracker_Rpos.back()[22]  = new TH1D((trackerName+"_rpos_PDG22").c_str(),
                                            (trackerName+" rpos (photons)").c_str(),
                                            1000,0,minR);
        tracker_Rpos.back()[2212]= new TH1D((trackerName+"_rpos_PDG2212").c_str(),
                                            (trackerName+" rpos (protons)").c_str(),
                                            1000,0,minR);
        tracker_Rpos.back()[0]   = new TH1D((trackerName+"_rpos_PDGother").c_str(),
                                            (trackerName+" rpos (other)").c_str(),
                                            1000,0,minR);
        for (auto hist : tracker_Rpos.back()) {
            hist.second->GetXaxis()->SetTitle("R [mm]");
        }

        tracker_Rpos_cutoff.push_back(std::map<G4int,TH1D*>());
        tracker_Rpos_cutoff.back()[11]  = new TH1D((trackerName+"_rpos_cutoff_PDG11").c_str(),
                                                   (trackerName+" rpos (electrons, energy > E_cut)").c_str(),
                                                   1000,0,minR);
        tracker_Rpos_cutoff.back()[-11] = new TH1D((trackerName+"_rpos_cutoff_PDG-11").c_str(),
                                                   (trackerName+" rpos (positrons, energy > E_cut)").c_str(),
                                                   1000,0,minR);
        tracker_Rpos_cutoff.back()[22]  = new TH1D((trackerName+"_rpos_cutoff_PDG22").c_str(),
                                                   (trackerName+" rpos (photos, energy > E_cut)").c_str(),
                                                   1000,0,minR);
        tracker_Rpos_cutoff.back()[2212]= new TH1D((trackerName+"_rpos_cutoff_PDG2212").c_str(),
                                                   (trackerName+" rpos (protons, energy > E_cut)").c_str(),
                                                   1000,0,minR);
        tracker_Rpos_cutoff.back()[0]   = new TH1D((trackerName+"_rpos_cutoff_PDGother").c_str(),
                                                   (trackerName+" rpos (other, energy > E_cut)").c_str(),
                                                   1000,0,minR);
        for (auto hist : tracker_Rpos_cutoff.back()) {
            hist.second->GetXaxis()->SetTitle("R [mm]");
        }
    }

    init_phasespaceX   =
        new TH2D("init_x",
                 "Initial phase space (x)",
                 1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                 1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
    //init_phasespaceX->Sumw2();
    init_phasespaceX->GetXaxis()->SetTitle("X [mm]");
    init_phasespaceX->GetYaxis()->SetTitle("X' [rad]");
    init_phasespaceY   =
        new TH2D("init_y",
                 "Initial phase space (y)",
                 1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                 1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
    //init_phasespaceY->Sumw2();
    init_phasespaceY->GetXaxis()->SetTitle("Y [mm]");
    init_phasespaceY->GetYaxis()->SetTitle("Y' [rad]");
    init_phasespaceXY   =
        new TH2D("init_xy",
                 "Initial phase space (x,y)",
                 1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                 1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm);
    init_E =
        new TH1D("init_E",
                 "Initial particle energy",
                 1000, 0.0, max(beamEnergy*1.1,genAct->get_beam_energy_flatMax()));
    init_E->GetXaxis()->SetTitle("Energy [MeV]");

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

    //For counting the types of particles hitting the detectors (for magnets it is defined elsewhere)
    if (detCon->GetHasTarget()){
        typeCounter["target"]        = particleTypesCounter();
        typeCounter["target_cutoff"] = particleTypesCounter();
    }

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
                                        engNbins+1,0,beamEnergy*(1+1/engNbins)) );
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
        magnet_exit_phasespaceX.back()->GetYaxis()->SetTitle("X' [rad]");

        magnet_exit_phasespaceY.push_back
            ( new TH2D((magName+"_y").c_str(),
                       (magName+" phase space (y)").c_str(),
                       1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                       1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad)
              );
        magnet_exit_phasespaceY.back()->GetXaxis()->SetTitle("Y [mm]");
        magnet_exit_phasespaceY.back()->GetYaxis()->SetTitle("Y' [rad]");

        magnet_exit_phasespaceX_cutoff.push_back
            ( new TH2D((magName+"_cutoff_x").c_str(),
                       (magName+" phase space (x) (charged, energy > Ecut, r < Rcut)").c_str(),
                       1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                       1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad)
              );
        magnet_exit_phasespaceX_cutoff.back()->GetXaxis()->SetTitle("X [mm]");
        magnet_exit_phasespaceX_cutoff.back()->GetYaxis()->SetTitle("X' [rad]");

        magnet_exit_phasespaceY_cutoff.push_back
            ( new TH2D((magName+"_cutoff_y").c_str(),
                       (magName+" phase space (y) (charged, energy > Ecut, r < Rcut)").c_str(),
                       1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                       1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad)
              );
        magnet_exit_phasespaceY_cutoff.back()->GetXaxis()->SetTitle("Y [mm]");
        magnet_exit_phasespaceY_cutoff.back()->GetYaxis()->SetTitle("Y' [rad]");

        magnet_exit_phasespaceX_cutoff_PDG.push_back(std::map<G4int,TH2D*>());
        magnet_exit_phasespaceY_cutoff_PDG.push_back(std::map<G4int,TH2D*>());

        magnet_exit_phasespaceX_cutoff_PDG.back()[11]  = new TH2D((magName+"_cutoff_x_PDG11").c_str(),
                                                                  (magName+" phase space (x) (electrons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        magnet_exit_phasespaceX_cutoff_PDG.back()[-11]  = new TH2D((magName+"_cutoff_x_PDG-11").c_str(),
                                                                  (magName+" phase space (x) (positrons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        magnet_exit_phasespaceX_cutoff_PDG.back()[22]  = new TH2D((magName+"_cutoff_x_PDG22").c_str(),
                                                                  (magName+" phase space (x) (photons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        magnet_exit_phasespaceX_cutoff_PDG.back()[2212]  = new TH2D((magName+"_cutoff_x_PDG2212").c_str(),
                                                                  (magName+" phase space (x) (protons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        magnet_exit_phasespaceX_cutoff_PDG.back()[0]  = new TH2D((magName+"_cutoff_x_PDGother").c_str(),
                                                                  (magName+" phase space (x) (other, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        for (auto PDG : magnet_exit_phasespaceX_cutoff_PDG.back()) {
            PDG.second->GetXaxis()->SetTitle("X [mm]");
            PDG.second->GetYaxis()->SetTitle("X' [rad]");
        }


        magnet_exit_phasespaceY_cutoff_PDG.back()[11]  = new TH2D((magName+"_cutoff_y_PDG11").c_str(),
                                                                  (magName+" phase space (y) (electrons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        magnet_exit_phasespaceY_cutoff_PDG.back()[-11]  = new TH2D((magName+"_cutoff_y_PDG-11").c_str(),
                                                                  (magName+" phase space (y) (positrons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        magnet_exit_phasespaceY_cutoff_PDG.back()[22]  = new TH2D((magName+"_cutoff_y_PDG22").c_str(),
                                                                  (magName+" phase space (y) (photons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        magnet_exit_phasespaceY_cutoff_PDG.back()[2212]  = new TH2D((magName+"_cutoff_y_PDG2212").c_str(),
                                                                  (magName+" phase space (y) (protons, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        magnet_exit_phasespaceY_cutoff_PDG.back()[0]  = new TH2D((magName+"_cutoff_y_PDGother").c_str(),
                                                                  (magName+" phase space (y) (other, energy > Ecut, r < Rcut)").c_str(),
                                                                  1000, -phasespacehist_posLim/mm,phasespacehist_posLim/mm,
                                                                  1000, -phasespacehist_angLim/rad,phasespacehist_angLim/rad);
        for (auto PDG : magnet_exit_phasespaceY_cutoff_PDG.back()) {
            PDG.second->GetXaxis()->SetTitle("Y [mm]");
            PDG.second->GetYaxis()->SetTitle("Y' [rad]");
        }

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
                
                magnetExit[i] = new TTree((magName+"_ExitHits").c_str(), "MagnetExit tree");
                magnetExitBuffer[i] = trackerHitStruct();
                magnetExit[i]->Branch("magnetExitBranch",  &(magnetExitBuffer[i]), ("x/D:y:z:px:py:pz:E:PDG/I:charge:eventID"));
                
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
        G4int TargetEdep_CollID = SDman->GetCollectionID("target_edep");
        if (TargetEdep_CollID>=0){
            EdepHitsCollection* targetEdepHitsCollection = NULL;
            targetEdepHitsCollection = (EdepHitsCollection*) (HCE->GetHC(TargetEdep_CollID));
            if (targetEdepHitsCollection != NULL) {
                G4int nEntries = targetEdepHitsCollection->entries();
                G4double edep      = 0.0; // G4 units, normalized before Fill()
                G4double edep_NIEL = 0.0; // G4 units, normalized before Fill()
                G4double edep_IEL  = 0.0; // G4 units, normalized before Fill()
                for (G4int i = 0; i < nEntries; i++){
                    EdepHit* edepHit = (*targetEdepHitsCollection)[i];
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

                if (edep > 0.0) {
                    targetEdep->Fill(edep/MeV);
                    targetEdep_NIEL->Fill(edep_NIEL/keV);
                    targetEdep_IEL->Fill(edep_IEL/MeV);
                }
            }
            else {
                G4cout << "targetEdepHitsCollection was NULL!"<<G4endl;
            }
        }
        else {
            G4cout << "TargetEdep_CollID was " << TargetEdep_CollID << " < 0!"<<G4endl;
        }

        G4int TargetExitpos_CollID = SDman->GetCollectionID("target_exitpos");
        if (TargetExitpos_CollID>=0) {
            TrackerHitsCollection* targetExitposHitsCollection = NULL;
            targetExitposHitsCollection = (TrackerHitsCollection*) (HCE->GetHC(TargetExitpos_CollID));
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
                    target_exit_phasespaceXY->Fill(hitPos.x()/mm,hitPos.y()/mm);

                    if (charge != 0 and energy/MeV > beamEnergy*beamEnergy_cutoff and hitR/mm < position_cutoffR) {
                        target_exit_phasespaceX_cutoff->Fill(hitPos.x()/mm, momentum.x()/momentum.z());
                        target_exit_phasespaceY_cutoff->Fill(hitPos.y()/mm, momentum.y()/momentum.z());
                        target_exit_phasespaceXY_cutoff->Fill(hitPos.x()/mm,hitPos.y()/mm);
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
                        targetExitBuffer.y = hitPos.y()/mm;
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
            G4cout << "TargetExitpos_CollID was " << TargetExitpos_CollID << " < 0!"<<G4endl;
        }
    }

    //**Data from detectorTrackerSD**
    VirtualTrackerWorldConstruction* traCon = VirtualTrackerWorldConstruction::getInstance();
    for (int idx = 0; idx < traCon->getNumTrackers(); idx++) {

        std::string trackerName;
        if (traCon->getNumTrackers() == 1) { trackerName = "tracker"; }
        else                               { trackerName = std::string("tracker_") + std::to_string(idx+1); }
        
        G4int TrackerSD_CollID = SDman->GetCollectionID(std::string("tracker_") + std::to_string(idx+1)+"_exitpos"); //Internal name is always tracker_<idx>_exitpos

        if (TrackerSD_CollID>=0) {
            TrackerHitsCollection* trackerHitsCollection = NULL;
            trackerHitsCollection = (TrackerHitsCollection*) (HCE->GetHC(TrackerSD_CollID));
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
                    tracker_energy[idx]->Fill(energy/MeV);

                    if (tracker_type_energy[idx].find(PDG) != tracker_type_energy[idx].end()) {
                        tracker_type_energy[idx][PDG]->Fill(energy/MeV);
                    }
                    else {
                        tracker_type_energy[idx][0]->Fill(energy/MeV);
                    }

                    if (hitR/mm < position_cutoffR) {
                        if (tracker_type_cutoff_energy[idx].find(PDG) != tracker_type_cutoff_energy[idx].end()) {
                            tracker_type_cutoff_energy[idx][PDG]->Fill(energy/MeV);
                        }
                        else {
                            tracker_type_cutoff_energy[idx][0]->Fill(energy/MeV);
                        }
                    }

                    //Phase space
                    tracker_phasespaceX[idx]->Fill(hitPos.x()/mm, momentum.x()/momentum.z());
                    tracker_phasespaceY[idx]->Fill(hitPos.y()/mm, momentum.y()/momentum.z());
                    tracker_phasespaceXY[idx]->Fill(hitPos.x()/mm, hitPos.y()/mm);

                    if (energy/MeV > beamEnergy*beamEnergy_cutoff and hitR/mm < position_cutoffR) {
                        if (charge != 0) {
                            // All charged particles passing the cutoff
                            tracker_phasespaceX_cutoff[idx]->Fill(hitPos.x()/mm, momentum.x()/momentum.z());
                            tracker_phasespaceY_cutoff[idx]->Fill(hitPos.y()/mm, momentum.y()/momentum.z());
                            tracker_phasespaceXY_cutoff[idx]->Fill(hitPos.x()/mm, hitPos.y()/mm);
                        }

                        //Also separated by species
                        if(tracker_phasespaceX_cutoff_PDG[idx].find(PDG) != tracker_phasespaceX_cutoff_PDG[idx].end()) {
                            tracker_phasespaceX_cutoff_PDG[idx][PDG]->Fill(hitPos.x()/mm, momentum.x()/momentum.z());
                        }
                        else {
                            tracker_phasespaceX_cutoff_PDG[idx][0]->Fill(hitPos.x()/mm, momentum.x()/momentum.z());
                        }
			
                        if(tracker_phasespaceY_cutoff_PDG[idx].find(PDG) != tracker_phasespaceY_cutoff_PDG[idx].end()) {
                            tracker_phasespaceY_cutoff_PDG[idx][PDG]->Fill(hitPos.y()/mm, momentum.y()/momentum.z());
                        }
                        else {
                            tracker_phasespaceY_cutoff_PDG[idx][0]->Fill(hitPos.y()/mm, momentum.y()/momentum.z());
                        }

                        if(tracker_phasespaceXY_cutoff_PDG[idx].find(PDG) != tracker_phasespaceXY_cutoff_PDG[idx].end()) {
                            tracker_phasespaceXY_cutoff_PDG[idx][PDG]->Fill(hitPos.x()/mm, momentum.y()/mm);
                        }
                        else {
                            tracker_phasespaceXY_cutoff_PDG[idx][0]->Fill(hitPos.x()/mm, momentum.y()/mm);
                        }
                    }

                    //Particle type counting
                    FillParticleTypes(typeCounter[trackerName], PDG, type);
                    if (energy/MeV > beamEnergy*beamEnergy_cutoff and hitR/mm < position_cutoffR) {
                        FillParticleTypes(typeCounter[trackerName+"_cutoff"], PDG, type);
                    }

                    //R position
                    if (tracker_Rpos[idx].find(PDG) != tracker_Rpos[idx].end()) {
                        tracker_Rpos[idx][PDG]->Fill(hitR/mm);
                    }
                    else {
                        tracker_Rpos[idx][0]->Fill(hitR/mm);
                    }
                    if (energy/MeV > beamEnergy*beamEnergy_cutoff) {
                        if (tracker_Rpos_cutoff[idx].find(PDG) != tracker_Rpos_cutoff[idx].end()) {
                            tracker_Rpos_cutoff[idx][PDG]->Fill(hitR/mm);
                        }
                        else {
                            tracker_Rpos_cutoff[idx][0]->Fill(hitR/mm);
                        }
                    }

                    //Fill the TTree
                    if (not miniFile) {
                        trackerHitsBuffer.x = hitPos.x()/mm;
                        trackerHitsBuffer.y = hitPos.y()/mm;
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

                tracker_numParticles[idx]->Fill(nEntries);
            }
            else{
                G4cout << "trackerHitsCollection was NULL! for tracker '" + trackerName << "'" << G4endl;
            }
        }
        else{
            G4cout << "TrackerSD_CollID was " << TrackerSD_CollID << " < 0 for tracker '" << trackerName << "'" << G4endl;
        }
    }

    // Initial particle distribution
    init_phasespaceX->Fill(genAct->x/mm,genAct->xp/rad);
    init_phasespaceY->Fill(genAct->y/mm,genAct->yp/rad);
    init_phasespaceXY->Fill(genAct->x/mm,genAct->y/mm);
    init_E->Fill(genAct->E/MeV);

    //Fill the TTree
    if (not miniFile) {
        initPartsBuffer.x = genAct->x/mm;
        initPartsBuffer.y = genAct->y/mm;
        initPartsBuffer.z = genAct->z/mm; //Required to fill struct properly

        initPartsBuffer.px = genAct->xp;
        initPartsBuffer.py = genAct->yp;
        initPartsBuffer.pz = 0;

        initPartsBuffer.E = genAct->E / MeV;

        initPartsBuffer.PDG = genAct->PDG;
        initPartsBuffer.charge = genAct->PDG_Q;

        initPartsBuffer.eventID = eventCounter;

        initParts->Fill();
    }

    // *** Data from Magnets, which use a TargetSD ***
    size_t magIdx = -1;
    for (auto mag : detCon->magnets) {
        const G4String magName = mag->magnetName;
        magIdx++;

        //Edep data collection
        G4int MagnetEdep_CollID = SDman->GetCollectionID(magName+"_edep");
        if (MagnetEdep_CollID>=0){
            EdepHitsCollection* magnetEdepHitsCollection = NULL;
            magnetEdepHitsCollection = (EdepHitsCollection*) (HCE->GetHC(MagnetEdep_CollID));
            if (magnetEdepHitsCollection != NULL) {
                G4int nEntries = magnetEdepHitsCollection->entries();
                G4double edep      = 0.0;
                for (G4int i = 0; i < nEntries; i++){
                    EdepHit* edepHit = (*magnetEdepHitsCollection)[i];
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

                if (edep > 0.0) {
                    magnet_edep[magIdx]->Fill(edep/MeV);
                }

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
            G4cout << "MagnetEdep_CollID was " << MagnetEdep_CollID << " < 0 for '" << magName << "'!"<<G4endl;
        }


        // Exitpos data collection
        G4int MagnetExitpos_CollID = SDman->GetCollectionID(magName + "_exitpos");
        if (MagnetExitpos_CollID>=0) {
            TrackerHitsCollection* magnetExitposHitsCollection = NULL;
            magnetExitposHitsCollection = (TrackerHitsCollection*) (HCE->GetHC(MagnetExitpos_CollID));
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

                        if ( energy/MeV > beamEnergy*beamEnergy_cutoff and
                             hitR/mm < position_cutoffR
                             ) {
                            if(charge != 0) {
                                // All charged particles passing the cutoff
                                magnet_exit_phasespaceX_cutoff[magIdx]->
                                    Fill(hitPos.x()/mm, momentum.x()/momentum.z());
                                magnet_exit_phasespaceY_cutoff[magIdx]->
                                    Fill(hitPos.y()/mm, momentum.y()/momentum.z());
                            }

                            //Also separated by species
                            if(magnet_exit_phasespaceX_cutoff_PDG[magIdx].find(PDG) != magnet_exit_phasespaceX_cutoff_PDG[magIdx].end()) {
                                magnet_exit_phasespaceX_cutoff_PDG[magIdx][PDG]->Fill(hitPos.x()/mm, momentum.x()/momentum.z());
                            }
                            else {
                                magnet_exit_phasespaceX_cutoff_PDG[magIdx][0]->Fill(hitPos.x()/mm, momentum.x()/momentum.z());
                            }
                            if(magnet_exit_phasespaceY_cutoff_PDG[magIdx].find(PDG) != magnet_exit_phasespaceY_cutoff_PDG[magIdx].end()) {
                                magnet_exit_phasespaceY_cutoff_PDG[magIdx][PDG]->Fill(hitPos.y()/mm, momentum.y()/momentum.z());
                            }
                            else {
                                magnet_exit_phasespaceY_cutoff_PDG[magIdx][0]->Fill(hitPos.y()/mm, momentum.y()/momentum.z());
                            }
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

                    //Fill the TTree
                    if (not miniFile) {
                        magnetExitBuffer[magIdx].x = hitPos.x()/mm;
                        magnetExitBuffer[magIdx].y = hitPos.y()/mm;
                        magnetExitBuffer[magIdx].z = hitPos.z()/mm;

                        magnetExitBuffer[magIdx].px = momentum.x()/MeV;
                        magnetExitBuffer[magIdx].py = momentum.y()/MeV;
                        magnetExitBuffer[magIdx].pz = momentum.z()/MeV;

                        magnetExitBuffer[magIdx].E = energy / MeV;

                        magnetExitBuffer[magIdx].PDG = PDG;
                        magnetExitBuffer[magIdx].charge = charge;

                        magnetExitBuffer[magIdx].eventID = eventCounter;
                        
                        magnetExit[magIdx]->Fill();
                    }
                }
            }
            else {
                G4cout << "magnetExitposHitsCollection was NULL! for '" << magName << "'"<<G4endl;
            }
        }
        else {
            G4cout << "MagnetExitpos_CollID was " << MagnetExitpos_CollID << " < 0 for '" << magName << "'!"<<G4endl;
        }
    } // END loop over magnets
    if (not miniFile) {
        // Outside of the loop over magnets -- Fill to write all branches
        magnetEdeps->Fill();
    }
}
void RootFileWriter::finalizeRootFile() {

    //Needed for some of the processing
    G4RunManager*                      run  = G4RunManager::GetRunManager();
    DetectorConstruction*            detCon = (DetectorConstruction*)run->GetUserDetectorConstruction();
    VirtualTrackerWorldConstruction* traCon = VirtualTrackerWorldConstruction::getInstance();

    //Print out the particle types on all detector planes
    for (auto it : typeCounter) {
        PrintParticleTypes(it.second, it.first);
    }

    // ** Below cutoff **

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


    if (detCon->GetHasTarget()) {
        G4cout << G4endl
                << "Exit angle average (x) = " << exitangle_avg << " [deg]" << G4endl
                << "Exit angle RMS (x)     = " << exitangle_rms << " [deg]" << G4endl;
    }

    // ** Above cutoff **

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
           << beamEnergy*beamEnergy_cutoff <<" [MeV]) only:" << G4endl << G4endl;


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
        for (auto PDG : magnet_exit_phasespaceX_cutoff_PDG[magIdx]) {
            PrintTwissParameters(PDG.second);
        }
        for (auto PDG : magnet_exit_phasespaceY_cutoff_PDG[magIdx]) {
            PrintTwissParameters(PDG.second);
        }
    }
    for (int idx = 0; idx < traCon->getNumTrackers(); idx++) {
        PrintTwissParameters(tracker_phasespaceX[idx]);
        PrintTwissParameters(tracker_phasespaceY[idx]);
        PrintTwissParameters(tracker_phasespaceX_cutoff[idx]);
        PrintTwissParameters(tracker_phasespaceY_cutoff[idx]);
        for (auto PDG : tracker_phasespaceX_cutoff_PDG[idx]) {
            PrintTwissParameters(PDG.second);
        }
        for (auto PDG : tracker_phasespaceY_cutoff_PDG[idx]) {
            PrintTwissParameters(PDG.second);
        }
    }

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


    if (not quickmode) { //Write the 2D and 3D histograms to the ROOT file (slow)

        G4cout << "Writing 2D histograms..." << G4endl;

        init_phasespaceX->SetOption("COLZ");
        init_phasespaceX->Write();
        init_phasespaceY->SetOption("COLZ");
        init_phasespaceY->Write();
        init_phasespaceXY->SetOption("COLZ");
        init_phasespaceXY->Write();

        if (detCon->GetHasTarget()) {
            target_exit_phasespaceX->SetOption("COLZ");
            target_exit_phasespaceX->Write();
            target_exit_phasespaceY->SetOption("COLZ");
            target_exit_phasespaceY->Write();

            target_exit_phasespaceX_cutoff->SetOption("COLZ");
            target_exit_phasespaceX_cutoff->Write();
            target_exit_phasespaceY_cutoff->SetOption("COLZ");
            target_exit_phasespaceY_cutoff->Write();

            target_exit_phasespaceXY->SetOption("COLZ");
            target_exit_phasespaceXY->Write();
            target_exit_phasespaceXY_cutoff->SetOption("COLZ");
            target_exit_phasespaceXY_cutoff->Write();

            if (target_edep_rdens != NULL) {
                target_edep_rdens->Write();
            }
        }

        for (int idx = 0; idx < traCon->getNumTrackers(); idx++) {
            tracker_phasespaceX[idx]->SetOption("COLZ");
            tracker_phasespaceX[idx]->Write();
            tracker_phasespaceY[idx]->SetOption("COLZ");
            tracker_phasespaceY[idx]->Write();
            tracker_phasespaceXY[idx]->SetOption("COLZ");
	          tracker_phasespaceXY[idx]->Write();

            tracker_phasespaceX_cutoff[idx]->SetOption("COLZ");
            tracker_phasespaceX_cutoff[idx]->Write();
            tracker_phasespaceY_cutoff[idx]->SetOption("COLZ");
            tracker_phasespaceY_cutoff[idx]->Write();
            tracker_phasespaceXY_cutoff[idx]->SetOption("COLZ");
	          tracker_phasespaceXY_cutoff[idx]->Write();

            for (auto PDG : tracker_phasespaceX_cutoff_PDG[idx]) {
                PDG.second->SetOption("COLZ");
                PDG.second->Write();
            }
            for (auto PDG : tracker_phasespaceY_cutoff_PDG[idx]) {
                PDG.second->SetOption("COLZ");
                PDG.second->Write();
            }
            for (auto PDG : tracker_phasespaceXY_cutoff_PDG[idx]) {
                PDG.second->SetOption("COLZ");
                PDG.second->Write();
            }
        }

        for (auto it : magnet_edep_rdens) {
            if (it != NULL) {
                it->Write();
	    }
        }
        for (auto it : magnet_exit_phasespaceX) {
            it->SetOption("COLZ");
            it->Write();
        }
        for (auto it : magnet_exit_phasespaceY) {
            it->SetOption("COLZ");
            it->Write();
        }
        for (auto it : magnet_exit_phasespaceX_cutoff) {
            it->SetOption("COLZ");
            it->Write();
        }
        for (auto it : magnet_exit_phasespaceY_cutoff) {
            it->SetOption("COLZ");
            it->Write();
        }
        for (auto mag : magnet_exit_phasespaceX_cutoff_PDG) {
            for (auto PDG : mag) {
                PDG.second->SetOption("COLZ");
                PDG.second->Write();
            }
        }
        for (auto mag : magnet_exit_phasespaceY_cutoff_PDG) {
            for (auto PDG : mag) {
                PDG.second->SetOption("COLZ");
                PDG.second->Write();
            }
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
        initParts->Write(); //writing InitParts with other TTrees

        if (detCon->GetHasTarget()) {
            targetExit->Write();
        }
        trackerHits->Write();

        magnetEdeps->Write();
        delete magnetEdeps;
        magnetEdeps=NULL;

        for (auto mag : magnetExit) {
            mag.second->Write();
            delete mag.second;
        }
        magnetExit.clear();
        magnetExitBuffer.clear();

        if (magnetEdepsBuffer != NULL) {
            delete magnetEdepsBuffer;
            magnetEdepsBuffer = NULL;
        }
    }

    G4cout << "Writing 1D histograms..." << G4endl;

    init_E->Write();
    delete init_E; init_E = NULL;

    if (detCon->GetHasTarget()) {
        targetEdep->Write();
        delete targetEdep;      targetEdep      = NULL;
        targetEdep_NIEL->Write();
        delete targetEdep_NIEL; targetEdep_NIEL = NULL;
        targetEdep_IEL->Write();
        delete targetEdep_IEL;  targetEdep_IEL  = NULL;

        target_exitangle_hist->Write();
        delete target_exitangle_hist; target_exitangle_hist = NULL;

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

    for (int idx = 0; idx < traCon->getNumTrackers(); idx++) {
        tracker_numParticles[idx]->Write();
        delete tracker_numParticles[idx]; tracker_numParticles[idx] = NULL;

        tracker_energy[idx]->Write();
        delete tracker_energy[idx]; tracker_energy[idx] = NULL;

        for (auto it : tracker_type_energy[idx]) {
            it.second->Write();
            delete it.second;
        }
        tracker_type_energy[idx].clear();
        for (auto it : tracker_type_cutoff_energy[idx]) {
            it.second->Write();
            delete it.second;
        }
        tracker_type_cutoff_energy[idx].clear();

        for (auto it: tracker_Rpos[idx]) {
            it.second->Write();
            delete it.second;
        }
        tracker_Rpos[idx].clear();
        for (auto it: tracker_Rpos_cutoff[idx]) {
            it.second->Write();
            delete it.second;
        }   
        tracker_Rpos_cutoff[idx].clear();
    }
    tracker_numParticles.clear();
    tracker_energy.clear();
    tracker_type_energy.clear();
    tracker_type_cutoff_energy.clear();
    tracker_Rpos.clear();
    tracker_Rpos_cutoff.clear();

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

    for (auto mag : magnet_exit_phasespaceX_cutoff_PDG) {
        for (auto PDG : mag) {
            delete PDG.second;
        }
        mag.clear();
    }
    magnet_exit_phasespaceX_cutoff_PDG.clear();
    for (auto mag : magnet_exit_phasespaceY_cutoff_PDG) {
        for (auto PDG : mag) {
            delete PDG.second;
        }
        mag.clear();
    }
    magnet_exit_phasespaceY_cutoff_PDG.clear();

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

    //Now that we have plotted, delete 2D and 3D hists (even if we have disabled writing 2D and 3D histograms)

    delete init_phasespaceX; init_phasespaceX = NULL;
    delete init_phasespaceY; init_phasespaceY = NULL;
    delete init_phasespaceXY; init_phasespaceXY = NULL;

    if (detCon->GetHasTarget()) {
        delete target_exitangle_hist; target_exitangle_hist = NULL;

        delete target_exit_phasespaceX; target_exit_phasespaceX = NULL;
        delete target_exit_phasespaceY; target_exit_phasespaceY = NULL;

        delete target_exit_phasespaceX_cutoff; target_exit_phasespaceX_cutoff = NULL;
        delete target_exit_phasespaceY_cutoff; target_exit_phasespaceY_cutoff = NULL;

        delete target_exit_phasespaceXY; target_exit_phasespaceXY = NULL;
        delete target_exit_phasespaceXY_cutoff; target_exit_phasespaceXY_cutoff = NULL;
    }

    //2D tracker histos
    for (int idx = 0; idx < traCon->getNumTrackers(); idx++) {
        delete tracker_phasespaceX[idx]; tracker_phasespaceX[idx] = NULL;
        delete tracker_phasespaceY[idx]; tracker_phasespaceY[idx] = NULL;
        delete tracker_phasespaceXY[idx]; tracker_phasespaceXY[idx] = NULL;

        delete tracker_phasespaceX_cutoff[idx]; tracker_phasespaceX_cutoff[idx] = NULL;
        delete tracker_phasespaceY_cutoff[idx]; tracker_phasespaceY_cutoff[idx] = NULL;
        delete tracker_phasespaceXY_cutoff[idx]; tracker_phasespaceXY_cutoff[idx] = NULL;

        for (auto PDG : tracker_phasespaceX_cutoff_PDG[idx]) {
            delete PDG.second;
        }
        tracker_phasespaceX_cutoff_PDG[idx].clear();

        for (auto PDG : tracker_phasespaceY_cutoff_PDG[idx]) {
            delete PDG.second;
        }
        tracker_phasespaceY_cutoff_PDG[idx].clear();

        for (auto PDG : tracker_phasespaceXY_cutoff_PDG[idx]) {
            delete PDG.second;
        }
        tracker_phasespaceXY_cutoff_PDG[idx].clear();
    }
    tracker_phasespaceX.clear();
    tracker_phasespaceY.clear();
    tracker_phasespaceXY.clear();
    tracker_phasespaceX_cutoff.clear();
    tracker_phasespaceY_cutoff.clear();
    tracker_phasespaceXY_cutoff.clear();
    tracker_phasespaceX_cutoff_PDG.clear();
    tracker_phasespaceY_cutoff_PDG.clear();
    tracker_phasespaceXY_cutoff_PDG.clear();

    if (detCon->GetHasTarget()) {
        if (target_edep_dens != NULL) {
            delete target_edep_dens; target_edep_dens = NULL;
        }
        if (target_edep_rdens != NULL) {
            delete target_edep_rdens; target_edep_rdens = NULL;
        }

        delete target_exitangle_hist_cutoff; target_exitangle_hist_cutoff = NULL;
    }

    if(not miniFile) {
        delete trackerHits; trackerHits = NULL;
        if (detCon->GetHasTarget()) {
            delete targetExit; targetExit = NULL;
        }
        delete initParts; initParts = NULL;
    }

    histFile->Write();
    histFile->Close();
    delete histFile; histFile = NULL;
    G4cout << "Results written to ROOT file '" + rootFileName +"'." << G4endl;
    G4cout << G4endl;
}

void RootFileWriter::PrintTwissParameters(TH2D* phaseSpaceHist) {
    G4cout << "Stats for '" << phaseSpaceHist->GetTitle() << "':"  << G4endl;
    double stats[7];
    phaseSpaceHist->GetStats(stats);
    /* From docs: https://root.cern.ch/doc/master/classTH2.html#a7836d201c7ab24717a65a38779d2f243
     * Note: Here all weights = 1
     * stats[0] = sumw    (number of fills)
     * stats[1] = sumw2   (number of fills 1^1 = 1)
     * stats[2] = sumwx   (sum of X in fills)
     * stats[3] = sumwx2  (sum of X^2 in fills)
     * stats[4] = sumwy   (sum of Y^2 in fills)
     * stats[5] = sumwy2  (sum of Y^2 in fills)
     * stats[6] = sumwxy  (sum of X*Y in fills)
     */

    // Fill used [mm] and [rad].
    // These variance formulas are susceptible to catastrophic cancellation if mean far off-center compared to sigma,
    // but mean should be close to zero so we should be OK
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

    double det = posVar*angVar - coVar*coVar; // [mm^2 * rad^2]
    double epsG = sqrt(det);                  // [mm*rad]
    double epsN = epsG*beta_rel*gamma_rel;    // [mm*rad]
    double beta = posVar/epsG;                // [mm]
    double alpha = -coVar/epsG;               // [-]

    G4cout << "Geometrical emittance          = " << epsG*1e3 << " [um = mm*mrad]" << G4endl;
    G4cout << "Normalized emittance           = " << epsN*1e3 << " [um = mm*mrad]"
           << ", assuming beam kinetic energy = " << genAct->get_beam_energy() << " [MeV]"
           << ", and mass = " << genAct->get_beam_particlemass()/MeV << " [MeV/c^2]"
           << G4endl;
    G4cout << "Twiss beta  = " << beta*1e-3  << " [m]" << G4endl
           << "Twiss alpha = " << alpha      << " [-]"  << G4endl;

    G4cout << G4endl;

    // Write to root file
    TVectorD twissVector (8);
    twissVector[0] = epsN*1e3;  // [um = mm*mrad]
    twissVector[1] = beta*1e-3; // [m]
    twissVector[2] = alpha;     // [-]
    twissVector[3] = posAve;    // [mm]
    twissVector[4] = angAve;    // [rad]
    twissVector[5] = posVar;    // [mm^2]
    twissVector[6] = angVar;    // [rad^2]
    twissVector[7] = coVar;     // [mm*rad]
    twissVector.Write( (G4String(phaseSpaceHist->GetName())+"_TWISS").c_str() );
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
        G4String errormessage = "edepNbins must be > 0 (or 0 for auto)";
        G4Exception("RootFileWriter::setEngNbins()", "MSRootFile2000",FatalException,errormessage);
    }
}
