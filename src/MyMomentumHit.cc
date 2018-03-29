#include "MyMomentumHit.hh"

G4Allocator<MyMomentumHit> MyMomentumHitAllocator;

MyMomentumHit::MyMomentumHit(G4ThreeVector momentum) :
  trackMomentum(momentum) {
}


MyMomentumHit::~MyMomentumHit() {
}

void MyMomentumHit::Print() {
}
