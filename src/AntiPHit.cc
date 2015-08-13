
//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id:$
//
// Jane Tinslay - adapted from A01 example
//
#include "AntiPHit.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4AttDef.hh"
#include "G4AttDefStore.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4ParticleGun.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "CLHEP/Units/PhysicalConstants.h"


G4Allocator<AntiPHit> AntiPHitAllocator;

AntiPHit::AntiPHit()
  :fCellID(-1),
   fDepositedEnergy(0),
   fHitPosition(),
   fParticleID(),
   fParticleName(),
   pLogicalVolume(0)
{}

AntiPHit::AntiPHit(G4int id)
  :fCellID(id),
   fDepositedEnergy(0),
   fHitPosition(),
   fParticleID(),
   fParticleName(),
   pLogicalVolume(0)
{}

AntiPHit::~AntiPHit() {}

void AntiPHit::Print()
{
}


