#include "MyEdepHit.hh"
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


G4Allocator<MyEdepHit> MyEdepHitAllocator;

MyEdepHit::MyEdepHit() :
   fDepositedEnergy(0.0), fDepositedEnergy_NIEL(0.0) {
}

MyEdepHit::~MyEdepHit() {
}

void MyEdepHit::Print() {
}


