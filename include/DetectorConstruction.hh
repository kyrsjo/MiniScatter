// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.hh,v 1.1 2010/10/18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4UniformMagField.hh"
//#include "DetectorMessenger.hh"
//#include "globals.hh"
// class G4Box;
// class G4LogicalVolume;
// class G4VPhysicalVolume;
// class G4Material;
// class G4UniformMagField;
// class DetectorMessenger;
//--------------------------------------------------------------------------------

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
  void SetSiliconMaterial (G4String); 
  void SetMagField(G4double);
  G4VPhysicalVolume* Construct();
public:
  G4double GetWorldSizeX()           {G4cout<<"returning "<<WorldSizeX<<G4endl; return WorldSizeX;}; 
  G4double GetWorldSizeYZ()          {return WorldSizeYZ;};
  
  G4double GetSiliconSizeX()          {return SiliconSizeX;};
  G4double GetSiliconSizeY()          {return SiliconSizeY;};
  G4int GetNXpixels(){return xpixels;};
  G4int GetNYpixels(){return ypixels;};

  G4int GetNbOfLayers()              {return NbOfLayers;}; 
  
  G4Material* GetSiliconMaterial()  {return SiliconMaterial;};  
  const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
  const G4VPhysicalVolume* GetSilicon()   {return physiSilicon;};
  
private:
  bool               background;
  G4Material*        SiliconMaterial;
  G4Material*        ThermalOMaterial;
  G4Material*        AlMaterial;
  G4Material*        CMaterial;
  G4Material*        CuMaterial;
  G4Material*        PbMaterial;
  G4Material*        TiMaterial;
  G4Material*        StainlessSteel; 
  G4Material*        ceramic;


  G4int              NbOfLayers;
  G4double           LayerThickness;
  
  G4double           SiliconSizeX;
  G4double           SiliconSizeY;
  G4int xpixels;
  G4int ypixels;
  G4double xpixel_pitch;
  G4double ypixel_pitch;

  G4double           SiliconThickness;

  G4Material*        defaultMaterial;
  G4double           WorldSizeYZ;
  G4double           WorldSizeX;
  
  G4Box*             solidWorld;    //pointer to the solid World 
  G4LogicalVolume*   logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* physiWorld;    //pointer to the physical World
  
  G4Box*             solidSilicon; //pointer to the solid Silicon
  G4LogicalVolume*   logicSilicon; //pointer to the logical Silicon
  G4VPhysicalVolume* physiSilicon; //pointer to the physical Silicon
  
  
  G4UniformMagField* magField;      //pointer to the magnetic field
  
   
private:
  void DefineMaterials();
  //  void ComputeWorldParameters();
  G4VPhysicalVolume* ConstructSilicon();     
};

//--------------------------------------------------------------------------------

//--------------------------------------------------------------------------------

#endif

