#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
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
                                               G4String covarianceString_in,
                                               G4double Rcut_in ) :
    Detector(DC),
    beam_energy(beam_energy_in),
    beam_type(beam_type_in),
    beam_offset(beam_offset_in),
    beam_zpos(beam_zpos_in),
    doBacktrack(doBacktrack_in),
    covarianceString(covarianceString_in),
    Rcut(Rcut_in) {

    G4int n_particle = 1;
    particleGun  = new G4ParticleGun(n_particle);

    if (beam_zpos == 0.0) {
        beam_zpos = - ( Detector->getTargetThickness() / 2.0 +
                        Detector->WorldSizeZ_buffer    / 2.0   );
    }
    else {
        beam_zpos *= mm;

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

}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete particleGun;
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
    G4double gamma_rel = beam_energy*MeV/particle->GetPDGMass();
    G4cout << "gamma_rel = " << gamma_rel << G4endl;
    G4double beta_rel = sqrt(gamma_rel*gamma_rel - 1.0) / gamma_rel;
    G4cout << "beta_rel = " << beta_rel << G4endl;

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

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

    if (anEvent->GetEventID() == 0) {
        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        G4IonTable* ionTable = G4IonTable::GetIonTable();
        G4String ION = "ion";
        if (beam_type.compare(0, ION.length(), ION) == 0) {
            // Format: 'ion::Z,A'
            str_size ionZpos = beam_type.index("::")+2;
            str_size ionApos = beam_type.index(",")+1;
            if (ionZpos >= beam_type.length() or ionApos >= beam_type.length()) {
                G4cerr << "Error in parsing ion string; expected format: 'ion::Z,A'" << G4endl;
                exit(1);
            }
            G4int ionZ = std::stoi(beam_type(ionZpos,ionApos-ionZpos));
            G4int ionA = std::stoi(beam_type(ionApos,beam_type.length()));
            G4cout << "Initializing ion with Z = " << ionZ << ", A = " << ionA << G4endl;
            particle = ionTable->GetIon(ionZ,ionA);
        }
        else {
            particle = particleTable->FindParticle(beam_type);
        }
        if (particle == NULL) {
            G4cerr << "Error - particle named '" << beam_type << "'not found" << G4endl;
            //particleTable->DumpTable();
            exit(1);
        }
        particleGun->SetParticleDefinition(particle);

        G4cout << G4endl;
        G4cout << "Injecting beam at z0 = " << beam_zpos/mm << " [mm]" << G4endl;
        G4cout << "Distance to target   = " << (-beam_zpos - Detector->getTargetThickness()/2.0)/mm << "[mm]" << G4endl;
        G4cout << G4endl;

        if (covarianceString != "") {
            hasCovariance = true;
            setupCovariance();
        }
        if (covarianceString != "" or Rcut != 0.0) {
            RNG = new TRandom1();
        }
    }

    if (hasCovariance) {
        int loopCounter = 0;
        while(true) {
            G4double xn  = RNG->Gaus(0,1);
            G4double xpn = RNG->Gaus(0,1);;
            x  = (xn*covarX_L[0][0] + xpn*covarX_L[0][1])*m + beam_offset*mm;
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
    }
    else if (Rcut != 0.0) {
        RNG->Circle(x,y,Rcut);
        x *= mm;
        y *= mm;

        xp = 0.0;
        yp = 0.0;
    }
    else {
        x  = beam_offset*mm;
        xp = 0.0;
        y  = 0.0;
        yp = 0.0;
    }

    if (doBacktrack) {
        //Bactrack from 0.0 to beam_zpos (<0.0)
        x -= (xp/rad)*(0.0 - beam_zpos);
        y -= (yp/rad)*(0.0 - beam_zpos);
    }

    particleGun->SetParticlePosition(G4ThreeVector(x,y,beam_zpos));

    //Technically not completely accurate but close enough for now
    particleGun->SetParticleMomentumDirection(G4ThreeVector(xp,yp,1));

    particleGun->SetParticleEnergy(beam_energy*MeV);
    particleGun->GeneratePrimaryVertex(anEvent);
}
