#ifndef TrackingAction_h
#define TrackingAction_h 1
class G4Track;
#include "G4UserTrackingAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include <map>
#include <iostream>
#include <vector>
using namespace std;

class TrackingAction : public G4UserTrackingAction{
public:
  TrackingAction();

  virtual void PreUserTrackingAction(const G4Track*);
  virtual void PostUserTrackingAction(const G4Track*);

  //Clear the maps for a new event
  void ResetMaps();

  inline G4int GetParentID(G4int child) {
    return parentMap[child];
  }
  inline G4ParticleDefinition* GetDefByID(G4int ID) {
    return defMap[ID];
    }
  inline const G4Track* GetTrack(G4int ID){
    return trackMap[ID];
  }

  inline vector <int> GetKids(G4int ID){
    return childMap[ID];
  }
  
  inline G4VPhysicalVolume *GetVolume(G4int ID){
    return volumeMap[ID];
  }
  
  inline G4double GetLength(G4int ID){
    return lengthMap[ID];
  }

  inline G4double GetEnergy(G4int ID){
    return energyMap[ID];
  }


private:
  //What trackID <value> is parent to trackID <key> ?
  std::map<G4int, G4int> parentMap;
  // trackID <key>, particleDefinition<value>
  std::map<G4int, G4ParticleDefinition*> defMap;
  std::map<G4int, const G4Track*> trackMap;
  //std::map<G4int, std::vector<int> >childMap;
  std::map<G4int, std::vector<int> >childMap;
  std::map<G4int, G4VPhysicalVolume*> volumeMap;
  std::map<G4int, G4double> lengthMap;
  std::map<G4int, G4double> energyMap;
  
}; 


#endif

