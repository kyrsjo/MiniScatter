#include "MyTrackerHit.hh"

G4Allocator<MyTrackerHit> MyTrackerHitAllocator;

MyTrackerHit::MyTrackerHit() :
  trackEnergy(0.0), trackAngle(0.0), PDG(0) {
}

MyTrackerHit::MyTrackerHit(G4double energy, G4double angle, G4int id, G4ThreeVector momentum) :
  trackEnergy(energy), trackAngle(angle), PDG(id), trackMomentum(momentum) {
}


MyTrackerHit::~MyTrackerHit() {
}

void MyTrackerHit::Print() {
}


