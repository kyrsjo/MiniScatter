#include "TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "G4ParticleDefinition.hh"
#include <iostream>
#include <vector>
#include "G4SystemOfUnits.hh"

using namespace std;

TrackingAction::TrackingAction()
{;}

 void TrackingAction::PreUserTrackingAction(const G4Track* track){
  //Method only called once pr. track - no need to check if already added
  G4int trackID = track->GetTrackID();
  parentMap[trackID] = track->GetParentID();
  defMap[trackID] = track->GetDefinition();
  trackMap[trackID]=track;
  volumeMap[trackID]=track->GetVolume();
  energyMap[trackID]=track->GetKineticEnergy();
  if(childMap.find(GetParentID(trackID))==childMap.end()){
    vector <int > toAdd (1,trackID);    
    childMap[GetParentID(trackID)]=toAdd;
  }
  else {
    childMap[GetParentID(trackID)].push_back(trackID);
  }
}

 void TrackingAction::PostUserTrackingAction(const G4Track* track){
   G4cout<<"particle name "<<track->GetDefinition()->GetParticleName()<<"  track length "<<track->GetTrackLength()<<" name of volume "<<track->GetVolume()->GetName();
   G4int trackID = track->GetTrackID();
   lengthMap[trackID]=track->GetTrackLength();
  volumeMap[trackID]=track->GetVolume();
 }

void TrackingAction::ResetMaps(){
  std::map<G4int,G4int>::iterator it;
  parentMap.clear();
  defMap.clear(); 
  childMap.clear();
  lengthMap.clear();
  energyMap.clear();
  volumeMap.clear();
}
