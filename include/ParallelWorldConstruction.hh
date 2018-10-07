#ifndef ParallelWorldConstruction_hh
#define ParallelWolrdConstruction_hh 1

#include "G4VUserParallelWorld.hh"

#include "DetectorConstruction.hh"

// Parallel geometry in which the magnet sensitive detectors are defined;
// this gets around the problem of how to attach the SDs correctly onto complex magnet geometries.

class ParallelWorldConstruction : public G4VUserParallelWorld {
public:
    ParallelWorldConstruction(G4String worldName, DetectorConstruction* mainGeometryConstruction_in);
    virtual ~ParallelWorldConstruction();

    virtual void Construct();
    virtual void ConstructSD();

private:
    DetectorConstruction* mainGeometryConstruction = NULL;
    std::vector <G4VPhysicalVolume*> magnetDetectorPVs;
};

#endif
