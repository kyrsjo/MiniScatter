#include "MyEdepHit.hh"

G4Allocator<MyEdepHit> MyEdepHitAllocator;

MyEdepHit::MyEdepHit() :
   fDepositedEnergy(0.0), fDepositedEnergy_NIEL(0.0) {
}

MyEdepHit::~MyEdepHit() {
}

void MyEdepHit::Print() {
}


