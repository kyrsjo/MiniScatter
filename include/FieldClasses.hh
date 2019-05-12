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

#ifndef FieldClasses_h$
#define FieldClasses_h$

#include "G4SystemOfUnits.hh"

#include "G4MagneticField.hh"
#include "G4Navigator.hh"

class FieldBase : public G4MagneticField {
    // Class that has the navigator and transform
public:
    FieldBase(G4ThreeVector centerPoint_in, G4LogicalVolume* fieldLV_in) :
        centerPoint(centerPoint_in), fieldLV(fieldLV_in) {};
    virtual void PostInitialize() {
        SetupTransform();
    }
private:
    G4ThreeVector centerPoint;
    G4LogicalVolume* fieldLV;
    static G4Navigator* fNavigator;
protected:
    void SetupTransform();
    G4AffineTransform fGlobalToLocal;

    G4double gradient; // [T/m]
};

class FieldPLASMA1 : public FieldBase {
public:
    FieldPLASMA1(G4double current_in, G4double radius_in,
                 G4ThreeVector centerPoint_in, G4LogicalVolume* fieldLV_in);
    virtual void GetFieldValue(const G4double point[4], G4double field[6]) const;

private:
    G4double plasmaTotalCurrent; // [A]
    G4double capRadius;  // [G4 units]

};

#endif
