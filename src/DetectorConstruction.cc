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
// $Id: DetectorConstruction.cc,v 1.1 2010/10/18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "AntiPSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4UniformMagField.hh"
#include "G4Element.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"



//------------------------------------------------------------------------------

DetectorConstruction::DetectorConstruction() :
  AlMaterial(0), TargetMaterial(0),
  solidWorld(0),logicWorld(0),physiWorld(0),
  solidTarget(0),logicTarget(0),physiTarget(0),
  magField(0) {

  WorldSizeX  = 200*cm;
  WorldSizeYZ = 200*cm;

  TargetSizeX     = 1.408*cm;
  TargetSizeY	  = 1.408*cm;
  TargetThickness = 0.023*cm;
  
  DetectorSizeX     = 150*cm;
  DetectorSizeY     = 150*cm;
  DetectorThickness = 1*um;
  
  DetectorDistance = 5*cm;

  // materials
  DefineMaterials();
  SetTargetMaterial("G4_Cu");
  DetectorMaterial = vacuumMaterial;
}

//------------------------------------------------------------------------------

G4VPhysicalVolume* DetectorConstruction::Construct() {
  // Clean old geometry, if any
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  // World
  solidWorld = new G4Box("World",				//its name
			 (WorldSizeX/2)*mm,(WorldSizeYZ/2)*mm,(WorldSizeYZ/2)*mm);	//its size, divide by two to get the right size
  
  logicWorld = new G4LogicalVolume(solidWorld,     //its solid
				   vacuumMaterial, //its material
				   "World");       //its name
  
  physiWorld = new G4PVPlacement(0,			//no rotation
				 G4ThreeVector(),	//at (0,0,0)
				 logicWorld,		//its logical volume
				 "World",		//its name
				 0,			//its mother  volume
				 false,			//no boolean operation
				 0);			//copy number
  
  //constructing the target
  
  solidTarget = new G4Box("Target",
			  TargetThickness/2, TargetSizeY/2,TargetSizeX/2);
  logicTarget = new G4LogicalVolume(solidTarget,    //its solid
				    TargetMaterial, //its material
				    TargetMaterial->GetName()); //name
  
  physiTarget = new G4PVPlacement(0,		   //no rotation
				  G4ThreeVector(0.0,0.0,0.0),  //its position
				  logicTarget,     //its logical volume		    
				  TargetMaterial->GetName(), //its name
				  logicWorld,        //its mother
				  false,             //no boulean operat
				  0);                //copy number

  solidDetector = new G4Box("Detector",
			  DetectorThickness/2, DetectorSizeY/2,DetectorSizeX/2);
  logicDetector = new G4LogicalVolume(solidDetector,    //its solid
				      DetectorMaterial, //its material
				      DetectorMaterial->GetName()); //name
  
  physiDetector = new G4PVPlacement(0,		   //no rotation
				  G4ThreeVector(0.0,0.0,DetectorDistance),  //its position
				  logicDetector,     //its logical volume		    
				  DetectorMaterial->GetName(), //its name
				  logicWorld,        //its mother
				  false,             //no boulean operat
				  0);                //copy number

  
  G4VSensitiveDetector* detector = new AntiPSD("/mydet/Target");
  // Get pointer to detector manager                                                     
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  // Register detector with manager                                                      
  SDman->AddNewDetector(detector);
  // Attach detector to volume defining calorimeter cells                                
  logicTarget->SetSensitiveDetector(detector);
                                         
  //Visualization attributes
  
  //logicWorld->SetVisAttributes (G4VisAttributes::Invisible);
  //cellLogical->SetVisAttributes (G4VisAttributes::Invisible);
  // logicTarget->SetVisAttributes(G4VisAttributes::Invisible);
  //  logicThermalO->SetVisAttributes(G4VisAttributes::Invisible);
  //  logicDopSi->SetVisAttributes(G4VisAttributes::Invisible);
  //  logicAl->SetVisAttributes(G4VisAttributes::Invisible);
  
  return physiWorld;
}

//------------------------------------------------------------------------------

void DetectorConstruction::DefineMaterials() {
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);

  G4Material* Al = man->FindOrBuildMaterial("G4_Al");
  G4Material* C = man->FindOrBuildMaterial("G4_C");
  G4Material* Cu = man->FindOrBuildMaterial("G4_Cu");
  G4Material* Pb = man->FindOrBuildMaterial("G4_Pb");
  G4Material* Ti =man->FindOrBuildMaterial("G4_Ti");  
  
  G4Material* Vacuum = man->FindOrBuildMaterial("G4_Galactic");
  
  //default materials 
  vacuumMaterial   = Vacuum;
  AlMaterial       = Al;
  CMaterial        = C;
  CuMaterial       = Cu;
  PbMaterial       = Pb;
  TiMaterial       = Ti;
}

//------------------------------------------------------------------------------

void DetectorConstruction::SetTargetMaterial(G4String materialChoice) {
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) TargetMaterial = pttoMaterial;
}
/*
void DetectorConstruction::SetDetectorMaterial(G4String materialChoice) {
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) DetectorMaterial = pttoMaterial;
}
*/
//------------------------------------------------------------------------------
