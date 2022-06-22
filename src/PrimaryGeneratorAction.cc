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

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4String.hh"

#include <iostream>
#include <cmath>
#include <string>

#include "TRandom1.h"

// -----------------------------------------------------------------------------------------

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* DC,
                                               G4double beam_energy_in,
                                               G4String beam_type_in,
                                               G4double beam_offset_in,
                                               G4double beam_zpos_in,
                                               G4bool   doBacktrack_in,
                                               G4double beam_angle_in,
                                               G4String covarianceString_in,
                                               G4double Rcut_in,
                                               G4int    rngSeed_in,
                                               G4double beam_energy_min_in,
                                               G4double beam_energy_max_in,
                                               G4String beam_loadFile_in) :
    Detector(DC),
    beam_energy(beam_energy_in),
    beam_type(beam_type_in),
    beam_offset(beam_offset_in*mm),
    beam_zpos(beam_zpos_in),
    doBacktrack(doBacktrack_in),
    beam_angle(beam_angle_in*deg),
    covarianceString(covarianceString_in),
    Rcut(Rcut_in*mm),
    rngSeed(rngSeed_in),
    beam_energy_min(beam_energy_min_in),
    beam_energy_max(beam_energy_max_in),
    beam_loadFile(beam_loadFile_in) {

    G4int n_particle = 1;
    particleGun  = new G4ParticleGun(n_particle);

    if (beam_zpos == 0.0) {
        beam_zpos = PrimaryGeneratorAction::GetDefaultZpos(Detector->getTargetThickness()/mm);
    }
    else {
        beam_zpos *= mm;

        //Sanity check
        if (beam_zpos >= - Detector->getTargetThickness() / 2.0) {
            G4cout << "Beam starting position = " << beam_zpos/mm
                   << " [mm] is not behind target back plane = "
                   << - Detector->getTargetThickness() / 2.0 / mm<< " [mm]"
                   << G4endl;
            exit(1);
        }
        if (beam_zpos <= - Detector->getWorldSizeZ()/2.0) {
            G4cout << "Beam starting position = " << beam_zpos/mm
                   << " [mm] is behind world back plane = "
                   << - Detector->getWorldSizeZ() / 2.0 / mm<< " [mm]"
                   << G4endl;
            exit(1);
        }
    }

    //Sanity check
    if (beam_loadFile != "") {
        if ( (beam_offset != 0.0) || (beam_angle != 0.0) || (covarianceString != "") || (Rcut != 0.0) || (beam_energy_min != -1) || (beam_energy_max != -1) ) {
            G4cerr << "Error, user specified a flag which is incompatible with --beamFile." << G4endl;
            exit(1);
        }
        G4String beam_loadFile_lower = beam_loadFile;
        beam_loadFile_lower.toLower();
        if (beam_loadFile.rfind(".csv") == (beam_loadFile.length()-4)) {
            beam_loadFile_csv = std::ifstream(beam_loadFile,std::ifstream::in);
            if (!beam_loadFile_csv.good()) {
                G4cerr << "Error when opening file '" + beam_loadFile + "', does the file exist?" << G4endl;
                exit(1);
            }
            G4cout << "Opened CSV file '" + beam_loadFile + "'" << G4endl;
        }
        else {
            G4cerr << "Error, the file type of '" + beam_loadFile + "' could not be determined." << G4endl;
            exit(1);
        }
    }

    //Sanity check
    if (abs(beam_offset) > Detector->getWorldSizeX()/2.0) {
        G4cout << "Beam offset = " << beam_offset/mm
                   << " [mm] is outside of world volume, max = world size / 2 = "
                   << Detector->getWorldSizeX() / 2.0 / mm<< " [mm]"
                   << G4endl;
            exit(1);
    }

    //Sanity check
    if(abs(beam_angle) != 0.0) {
        if (abs(beam_angle_in) >= 90.0) {
            G4cout << "Beam angles >= 90 degrees are not supported."
                   << G4endl;
            exit(1);
        }
        if (covarianceString != "") {
            G4cout << "Beam angle and covariance matrix cannot be used together."
                   << G4endl;
            exit(1);
        }
        if (Rcut != 0.0)  {
            G4cout << "Beam angle and Rcut cannot be used together."
                   << G4endl;
            exit(1);
        }
        if (beam_offset != 0.0) {
            G4cout << "Beam angle and beam offset cannot be used together."
                   << G4endl;
            exit(1);
        }
        if (doBacktrack) {
            G4cout << "Beam angle and backtracking cannot be used together."
                   << G4endl;
            exit(1);
        }
    }
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete particleGun;
}

G4double PrimaryGeneratorAction::GetDefaultZpos(G4double targetThickness_in) {
        G4double beam_zpos = targetThickness_in / 2.0;
        // Round up to nearest 10 mm
        G4double stepsize = 10.0*mm;
        beam_zpos = ceil(beam_zpos/stepsize)*stepsize;

        return -beam_zpos; //Negative!
    }

void PrimaryGeneratorAction::setupCovariance() {
    // Convert the covarianceString to a set of Twiss parameters,
    // then setup the covariance matrix and do the Cholesky Decomposition
    // Format: epsN[um]:beta[m]:alpha(::epsN_Y[um]:betaY[m]:alphaY)

    G4cout << G4endl;
    G4cout << "Initializing covariance matrices..." << G4endl;

    // Convert the string to relevant variables
    str_size startPos = 0;
    str_size endPos   = covarianceString.index(":",startPos);
    epsN_x            = convertColons(startPos,endPos, "epsN");

    startPos = endPos+1;
    endPos   = covarianceString.index(":",startPos);
    beta_x   = convertColons(startPos,endPos, "beta");

    if ( covarianceString.contains("::") ) {
        startPos = endPos+1;
        endPos   = covarianceString.index(":",startPos);
        alpha_x  = convertColons(startPos,endPos, "alpha");

        startPos = endPos+2;
        endPos   = covarianceString.index(":",startPos);
        epsN_y   = convertColons(startPos,endPos, "epsN_y");

        startPos = endPos+1;
        endPos   = covarianceString.index(":",startPos);
        beta_y   = convertColons(startPos,endPos, "beta_y");

        startPos = endPos+1;
        endPos   = covarianceString.length();
        alpha_y  = convertColons(startPos,endPos, "alpha_y");
    }
    else {
        startPos = endPos+1;
        endPos   = covarianceString.length();
        alpha_x  = convertColons(startPos,endPos, "alpha");

        epsN_y  = epsN_x;
        beta_y  = beta_x;
        alpha_y = alpha_x;
    }

    //Compute the geometrical emittance
    G4cout << "beam_energy       = " << beam_energy << " [MeV]" << G4endl;
    G4double beam_total_energy = beam_energy*MeV + particle->GetPDGMass();
    G4cout << "beam_total_energy = " << beam_total_energy/MeV << G4endl;
    G4double gamma_rel = beam_total_energy/particle->GetPDGMass();
    G4cout << "gamma_rel         = " << gamma_rel << G4endl;
    G4double beta_rel = sqrt(gamma_rel*gamma_rel - 1.0) / gamma_rel;
    G4cout << "beta_rel          = " << beta_rel << G4endl;

    epsG_x = epsN_x / (beta_rel*gamma_rel);
    epsG_y = epsN_y / (beta_rel*gamma_rel);

    G4cout << "Got Twiss parameters:" << G4endl;
    G4cout << "epsN_x  = " << epsN_x << "[um]" << G4endl;
    G4cout << "epsG_x  = " << epsG_x << "[um]" << G4endl;
    G4cout << "beta_x  = " << beta_x << "[m]"  << G4endl;
    G4cout << "alpha_x = " << alpha_x << G4endl;
    G4cout << "epsN_y  = " << epsN_y << "[um]" << G4endl;
    G4cout << "epsG_y  = " << epsG_y << "[um]" << G4endl;
    G4cout << "beta_y  = " << beta_y << "[m]"  << G4endl;
    G4cout << "alpha_y = " << alpha_y << G4endl;

    //Create covariance matrices
    covarX.ResizeTo(2,2);
    covarX[0][0] = beta_x;
    covarX[0][1] = -alpha_x;
    covarX[1][0] = -alpha_x;
    covarX[1][1] = (1+alpha_x*alpha_x)/beta_x;
    //covarX.Print(); // Raw matrix

    covarX *= (epsG_x*1e-6);

    G4cout << "Covariance matrix (X) [m^2, m * rad, rad^2]:" << G4endl;
    covarX.Print();

    covarY.ResizeTo(2,2);
    covarY[0][0] = beta_y;
    covarY[0][1] = -alpha_y;
    covarY[1][0] = -alpha_y;
    covarY[1][1] = (1+alpha_y*alpha_y)/beta_y;
    //covarY.Print(); // Raw matrix

    covarY *= (epsG_y*1e-6);

    G4cout << "Covariance matrix (Y) [m^2, m * rad, rad^2]:" << G4endl;
    covarY.Print();

    // Get the cholesky decomposition
    TDecompChol covarX_Utmp(covarX,1e-9);
    covarX_Utmp.Decompose();
    G4cout << "Decomposed upper-tridiagonal matrix for Y (to be transposed):" << G4endl;
    covarX_Utmp.Print();
    covarX_L.ResizeTo(covarX);
    covarX_L = covarX_Utmp.GetU();
    // Get the lower-tridiagonal that is needed for the generation
    covarX_L.Transpose(covarX_L);
    //covarX_L.Print();

    TDecompChol covarY_Utmp(covarY,1e-9);
    covarY_Utmp.Decompose();
    G4cout << "Decomposed upper-tridiagonal matrix for Y (to be transposed):" << G4endl;
    covarY_Utmp.Print();
    covarY_L.ResizeTo(covarY);
    covarY_L = covarY_Utmp.GetU();
    covarY_L.Transpose(covarY_L);

    G4cout << G4endl;
}
G4double PrimaryGeneratorAction::convertColons(str_size startPos, str_size endPos, G4String paramName) {
    if (endPos == std::string::npos) {
        G4cout << "PrimaryGeneratorAction::convertColons():" << G4endl
               << " Error while searching for " << paramName << " in '"
               << covarianceString << "'" << G4endl;
        exit(1);
    }
    G4String floatString = covarianceString(startPos,endPos-startPos);

    G4double floatData = 0.0;
    try {
        floatData = std::stod(std::string(floatString));
    }
    catch (const std::invalid_argument& ia) {
        G4cerr << "Invalid float in data '" << floatData << "'" << G4endl;
        exit(1);
    }

    /*
    G4cout << "got floatString = '" << floatString << "',"
           << " startPos = " << startPos << ", endPos = " << endPos << ","
           << " floatData = " << floatData << G4endl;
    */

    return floatData;
}

G4ParticleDefinition* PrimaryGeneratorAction::parseParticleName(G4String particleString) {
    const G4String ION = "ion";
    G4ParticleDefinition* particle_ret = NULL;
    if (particleString.compare(0, ION.length(), ION) == 0) {
        // Format: 'ion::Z,A'
        str_size ionZpos = particleString.index("::")+2;
        str_size ionApos = particleString.index(";")+1;
        if (ionZpos >= particleString.length() or ionApos >= particleString.length()) {
            G4cerr << "Error in parsing ion string; expected format: 'ion::Z,A'" << G4endl;
            exit(1);
        }
        G4int ionZ, ionA;
        try {
            ionZ = std::stoi(particleString(ionZpos,ionApos-ionZpos));
            ionA = std::stoi(particleString(ionApos,particleString.length()));
        }
        catch (const std::invalid_argument& ia) {
            G4cerr << "Error when extracting ionZ and ionA from string '" << particleString << "'" << G4endl;
            G4cerr << "ionZpos = " << ionZpos << G4endl;
            G4cerr << "ionApos = " << ionApos << G4endl;
            exit(1);
        }
        G4cout << "Initializing ion with Z = " << ionZ << ", A = " << ionA << G4endl;
        particle_ret = ionTable->GetIon(ionZ,ionA);
    }
    else {
        particle_ret = particleTable->FindParticle(particleString);
    }
    return particle_ret;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

    if (anEvent->GetEventID() == 0) {
        particleTable = G4ParticleTable::GetParticleTable();
        ionTable = G4IonTable::GetIonTable();
        
        particle = NULL;
        particle = parseParticleName(beam_type);
        if (particle == NULL) {
            G4cerr << "Error - particle named '" << beam_type << "' not found" << G4endl;
            //particleTable->DumpTable();
            exit(1);
        }
        PDG   = get_beam_particlePDG();
        PDG_Q = get_beam_particlecharge();

        particleGun->SetParticleDefinition(particle);

        G4cout << G4endl;
        G4cout << "Injecting beam at z0 = " << beam_zpos/mm << " [mm]" << G4endl;
        G4cout << "Distance to target   = " << (-beam_zpos - Detector->getTargetThickness()/2.0)/mm << "[mm]" << G4endl;
        G4cout << G4endl;

        //If needed, setup the covariance matrices for beam generation
        if (covarianceString != "") {
            hasCovariance = true;
            setupCovariance();
        }

        //If needed, setup the RNG
        if (covarianceString != "" or Rcut != 0.0 or (beam_energy_min >= 0.0 and beam_energy_max > 0.0)) {
            RNG = new TRandom1((UInt_t) rngSeed);
        }

	
    } //END first-event setup

    if (hasCovariance) {
        int loopCounter = 0;
        while(true) {
            G4double xn  = RNG->Gaus(0,1);
            G4double xpn = RNG->Gaus(0,1);;
            x  = (xn*covarX_L[0][0] + xpn*covarX_L[0][1])*m + beam_offset;
            xp = (xn*covarX_L[1][0] + xpn*covarX_L[1][1])*rad;

            G4double yn  = RNG->Gaus(0,1);
            G4double ypn = RNG->Gaus(0,1);
            y  = (yn*covarY_L[0][0] + ypn*covarY_L[0][1])*m;
            yp = (yn*covarY_L[1][0] + ypn*covarY_L[1][1])*rad;

            if (Rcut == 0.0) break;
            else {
                G4double init_r = sqrt(x*x + y*y)/mm;
                if (init_r <= Rcut) break;

                // No hit; do the anti-hang check.
                loopCounter++;
                if (loopCounter == 10000) {
                    G4cerr << "Warning in PrimaryGeneratorAction::GeneratePrimaries():" << G4endl
                           << " Number of iterations for the Rcut = " << Rcut/mm << " [mm]"
                           << " has reached 10k." << G4endl;
                }
                else if (loopCounter > 1000000) {
                    G4cerr << "Error in PrimaryGeneratorAction::GeneratePrimaries():" << G4endl
                           << " Number of iterations for the Rcut = " << Rcut/mm << " [mm]"
                           << " has reached 1M. Aborting." << G4endl;
                    exit(1);
                }
            }
        }
        z = beam_zpos;
    }
    else if (Rcut != 0.0) {
        G4double r = RNG->Uniform();
        G4double t = RNG->Uniform(2*M_PI);
        x = sqrt(r)*cos(t)*Rcut + beam_offset;
        y = sqrt(r)*sin(t)*Rcut;

        xp = 0.0;
        yp = 0.0;

        z = beam_zpos;
    }
    else if (beam_angle != 0.0) {
        x  = -beam_zpos*sin(beam_angle/rad);
        xp = tan(-beam_angle/rad);

        y  = 0.0;
        yp = 0.0;

        z = beam_zpos*cos(beam_angle/rad);
    }
    else if (beam_loadFile) {
        if (beam_loadFile_csv.is_open()) {
            if (beam_loadFile_csv.eof()) {
                G4cerr << "ERROR: Got EOF while reading file, event/line no. " << anEvent->GetEventID() << G4endl << G4endl;
                exit(1);
            }

            std::string line;
            std::getline(beam_loadFile_csv,line);

            std::stringstream lineStream(line);
            std::string word;
            try {
                G4cout << "Line='" << line << "' -> Words: '";

                std::getline(lineStream,word,',');
                G4String particle_name = word;
                particle = parseParticleName(particle_name);
                if (particle==NULL) {
                    G4cerr << "ERROR when parsing particle name '" + particle_name + "'." << G4endl;
                    exit(1);
                }
                PDG   = get_beam_particlePDG();
                PDG_Q = get_beam_particlecharge();
                particleGun->SetParticleDefinition(particle);
                

                std::getline(lineStream,word,',');
                x = std::stod(word)*mm;
                G4cout << word << "','";

                std::getline(lineStream,word,',');
                xp = std::stod(word);
                G4cout << word << "','";

                std::getline(lineStream,word,',');
                y = std::stod(word)*mm;
                G4cout << word << "','";

                std::getline(lineStream,word,',');
                yp = std::stod(word);
                G4cout << word << "','";

                std::getline(lineStream,word,',');
                z = std::stod(word)*mm;
                G4cout << word << "','";

                std::getline(lineStream,word,',');
                E = std::stod(word)*MeV;
                G4cout << word << "'" << G4endl;

                

            }
            catch (const std::invalid_argument& ia) {
                G4cerr << "Error when parsing line '" + line + "'" << G4endl;
                exit(1);
            }
        }
    }
    else {
        x  = beam_offset;
        xp = 0.0;

        y  = 0.0;
        yp = 0.0;

        z = beam_zpos;
    }

    if (doBacktrack) {
        //Bactrack from 0.0 to beam_zpos (<0.0)
        x -= (xp/rad)*(0.0 - beam_zpos);
        y -= (yp/rad)*(0.0 - beam_zpos);
    }

    particleGun->SetParticlePosition(G4ThreeVector(x,y,z));

    //Gives a non-unit vector, however it is normalized by
    //  SetParticleMomentumDirection, which fixes it.
    //  This works: With beam_angle, it hits (0,0,0) perfectly even at large angles.
    particleGun->SetParticleMomentumDirection(G4ThreeVector(xp,yp,1));

    if (not beam_loadFile) {
        if (beam_energy_min >= 0.0 and beam_energy_max > 0.0) {
            E = beam_energy_min+RNG->Uniform()*(beam_energy_max-beam_energy_min);
            E *= MeV;
        }
        else {
            E = beam_energy*MeV;
        }
    }
    particleGun->SetParticleEnergy(E); //Setting the kinetic energy (E>0 is valid)
    particleGun->GeneratePrimaryVertex(anEvent);

}

void PrimaryGeneratorAction::endOfRun() {
    if (beam_loadFile_csv) {
        beam_loadFile_csv.close();
    }
}