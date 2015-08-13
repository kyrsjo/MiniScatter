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
#ifndef AntiPHIT_HH
#define AntiPHIT_HH

#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4THitsCollection.hh"
#include "G4Transform3D.hh"
#include "G4VHit.hh"
#include "G4Track.hh"

class G4AttDef;
class G4AttValue;

class AntiPHit : public G4VHit {

public:
  
  // Constructors
  AntiPHit();
  AntiPHit(G4int id);

  // Destructor
  virtual ~AntiPHit();  
  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  // Methods
  virtual void Print();
  // Cell ID
  inline void  SetCellID(G4int id) { fCellID = id; }
  inline G4int GetCellID() { return fCellID; }
  // Deposited energy
  inline void AddDepositedEnergy(G4double energy) {fDepositedEnergy = energy;}
  inline G4double GetDepositedEnergy() const {return fDepositedEnergy;}
  // HitPosition
  inline void SetHitPosition(const G4ThreeVector& pos) { fHitPosition = pos; }
  inline G4ThreeVector& GetHitPosition() { return fHitPosition; }
  // Particle ID
  inline void SetParticleID(G4int partID){fParticleID = partID;}
  inline G4int GetParticleID() const {return fParticleID;}
  // Particle Name
  inline void SetParticleName(const G4String & partName){fParticleName = partName;}
  inline G4String GetParticleName() const {return fParticleName;}
  inline void SetHitTrack(G4Track* track){fTrack=track;}
  inline G4Track * GetHitTrack(){return fTrack;}

  // Logical volume
  inline void SetLogicalVolume(G4LogicalVolume* volume) {pLogicalVolume = volume;}
  inline const G4LogicalVolume* GetLogicalVolume() const {return pLogicalVolume;}
  
private:
  G4int fCellID;
  G4double fDepositedEnergy;   // Energy deposit
  G4ThreeVector fHitPosition;  // Hit Position
  G4int fParticleID;			// Particle ID
  G4String fParticleName;
  G4Track * fTrack;
  
  const G4LogicalVolume* pLogicalVolume;
  
};

typedef G4THitsCollection<AntiPHit> AntiPHitsCollection;

extern G4Allocator<AntiPHit> AntiPHitAllocator;

inline void* AntiPHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*)AntiPHitAllocator.MallocSingle();
  return aHit;
}

inline void AntiPHit::operator delete(void* aHit)
{
  AntiPHitAllocator.FreeSingle((AntiPHit*) aHit);
}

#endif


