/*
 * This file is part of MiniScatter.
 *
 *  MiniScatter is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  MiniScatter is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with MiniScatter.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef VirtualTrackerWorldConstruction_hh
#define VirtualTrackerWorldConstruction_hh 1

#include "G4VUserParallelWorld.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "DetectorConstruction.hh"

// Parallel geometry in which the tracker planes and their sensitive detectors are defined;
// this makes it possible to place the planes completely freely.

class VirtualTrackerWorldConstruction : public G4VUserParallelWorld {
public:
    VirtualTrackerWorldConstruction(G4String worldName, DetectorConstruction* mainGeometryConstruction_in,
        std::vector<G4double>* trackerDistances_in, G4double trackerAngle_in);
    virtual ~VirtualTrackerWorldConstruction();

    virtual void Construct();
    virtual void ConstructSD();

//    inline G4double getTrackerDistance() const {return TrackerDistance;};
    inline G4double getTrackerSizeX()    const {return TrackerSizeX;};
    inline G4double getTrackerSizeY()    const {return TrackerSizeY;};

    //Used to compute the length of the volume; all input in mm and deg
    static G4double ComputeMaxAbsZ(std::vector<G4double>* trackerDistances_in, G4double trackerAngle_in, G4double world_size) {
        if (trackerDistances_in->size() == 0) return 0.0; //If no trackers

        //Find largest offset from z=0
        G4double maxD = 0;
        for (auto d : *trackerDistances_in) {
            if (fabs(d) > maxD) maxD = fabs(d);
        }

        //From tilt
        if ( fabs(trackerAngle_in) > 89.0) {
            G4cerr << "abs(trackerAngle) too close to or above 90 degrees, current setting = " << trackerAngle_in << " [deg]" << G4endl;
            G4cerr << "This would result in a very or infinitely long volume or inverted left/right axis." << G4endl;
            exit(1);
        }
        maxD += (world_size/2.0) * fabs(tan(trackerAngle_in * M_PI/180.0));

        // In case we have negative offsets,
        // the volume will be between the tracking plane and the box wall
        maxD += trackerThickness;

        // Round up to nearest 10 mm
        return ceil(maxD/10.0)*10.0;
    }

private:
    DetectorConstruction* mainGeometryConstruction = NULL;

    std::vector<G4double> trackerDistances; //[G4units]
    G4double trackerAngle;                  //[G4units]

    static constexpr G4double trackerThickness = 1.0e-3; // Very thin [mm]


    G4double           TrackerSizeX;
    G4double           TrackerSizeY;
    // Distance from center of target to center of tracking 2D plane:
    G4double           TrackerThickness;
    //G4bool             TrackerRotated;
    
    std::vector <G4LogicalVolume*>   virtualTrackerLVs;
    std::vector <G4VPhysicalVolume*> virtualTrackerPVs;

public:
    int getNumTrackers() {return virtualTrackerLVs.size();};
    
    //Singleton
public:
    static VirtualTrackerWorldConstruction* getInstance() {
        assert (VirtualTrackerWorldConstruction::singleton != NULL);
        return VirtualTrackerWorldConstruction::singleton;
    }
private:
    static VirtualTrackerWorldConstruction* singleton;
};

#endif
