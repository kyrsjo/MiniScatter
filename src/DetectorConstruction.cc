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

DetectorConstruction::DetectorConstruction()
  :SiliconMaterial(0),ThermalOMaterial(0),AlMaterial(0),
   solidWorld(0),logicWorld(0),physiWorld(0),
   solidSilicon(0),logicSilicon(0),physiSilicon(0),
   solidThermalO(0),logicThermalO(0),physiThermalO(0),
  solidDopSi(0),logicDopSi(0),physiDopSi(0),
  solidAl(0),logicAl(0),physiAl(0),
  magField(0)
{
  // default parameter values of the Silicon
  SiliconThickness = 0.023*cm;
  //  SiliconThickness = 2.0*cm;

  ThermalOThickness = 0.0*mm;
  DopSiThickness = 0.0*mm;
  AlThickness = 0.0*mm;
  ThermalOThickness = 0.0;
  DopSiThickness = 0.0;
  AlThickness = 0.0;
  SiliconSizeX       = 1.408*cm;
  SiliconSizeY	= 1.408*cm;
  xpixels = 256;
  ypixels=256;
  xpixel_pitch=SiliconSizeX/(xpixels*1.0);
  ypixel_pitch=SiliconSizeY/(ypixels*1.0);
  // materials
  DefineMaterials();
  SetSiliconMaterial("Silicon");
}

//------------------------------------------------------------------------------

DetectorConstruction::~DetectorConstruction()
{}

//------------------------------------------------------------------------------

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructSilicon();
}

//------------------------------------------------------------------------------

void DetectorConstruction::DefineMaterials()
{ 
  
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);
  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;   
  // define Elements  
  new G4Material("Silicon", z=14., a= 28.09*g/mole, density=2.3290*g/cm3);
  G4Material* Al = man->FindOrBuildMaterial("G4_Al");
  G4Material* C = man->FindOrBuildMaterial("G4_C");
  G4Material* Cu = man->FindOrBuildMaterial("G4_Cu");
  G4Material* Pb = man->FindOrBuildMaterial("G4_Pb");
  G4Material* Ti =man->FindOrBuildMaterial("G4_Ti");  
  G4Material* SiO2 = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4Material* Vacuum = new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
				      kStateGas, 2.73*kelvin, 3.e-18*pascal);  
  //default materials 
  defaultMaterial  = Vacuum;
  ThermalOMaterial = SiO2;
  AlMaterial = Al;
  CMaterial = C;
  CuMaterial = Cu;
  PbMaterial = Pb;
  TiMaterial =Ti;
}

//------------------------------------------------------------------------------

G4VPhysicalVolume* DetectorConstruction::ConstructSilicon()
{
  
  // Clean old geometry, if any
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  // complete the World parameters definition
  ComputeWorldParameters();  
  //     
  // World
  //
  solidWorld = new G4Box("World",				//its name
			 (WorldSizeX/2)*mm,(WorldSizeYZ/2)*mm,(WorldSizeYZ/2)*mm);	//its size, divide by two to get the right size
  
  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
				   defaultMaterial,	//its material
				   "World");		//its name
  
  physiWorld = new G4PVPlacement(0,			//no rotation
				 G4ThreeVector(),	//at (0,0,0)
				 logicWorld,		//its logical volume				 
				 "World",		//its name
				 0,			//its mother  volume
				 false,			//no boolean operation
				 0);			//copy number
  
  //constructing the detector itself
  solidSilicon=0; logicSilicon=0; physiSilicon=0;  
  if (SiliconThickness > 0.) 
    { 
      solidSilicon = new G4Box("Silicon",		//its name
			       SiliconThickness/2, SiliconSizeY/2,SiliconSizeX/2); 
      logicSilicon = new G4LogicalVolume(solidSilicon,    //its solid
					 SiliconMaterial, //its material
					 SiliconMaterial->GetName()); //name
      
      physiSilicon = new G4PVPlacement(0,		   //no rotation
				       G4ThreeVector(0.0,0.0,0.0*cm),  //its position
				       logicSilicon,     //its logical volume		    
				       SiliconMaterial->GetName(), //its name
				       logicWorld,        //its mother
				       false,             //no boulean operat
				       0);                //copy number
    }
  
  if (ThermalOThickness > 0.) 
    {
      solidThermalO = new G4Box("ThermalOLayer",		//its name
				ThermalOThickness/2, SiliconSizeY/2,SiliconSizeX/2); 
      
      logicThermalO = new G4LogicalVolume(solidThermalO,    //its solid
					  ThermalOMaterial, //its material
					  "ThermalOLayer"); //name
      
      //tenker oss forand som negavie x retning. Da legger vi oksidasjonslaget forand silicon laget
      physiThermalO = new G4PVPlacement(0,		   //no rotation
					G4ThreeVector(-((SiliconThickness/2)+(ThermalOThickness/2)),0.,0.0*cm),  //its position
					logicThermalO,     //its logical volume		    
					"ThermalOLayer", //its name
					logicWorld,        //its mother
					false,             //no boulean operat
					0);                //copy number  
    }
  
  if (DopSiThickness > 0.) 
    {
      solidDopSi = new G4Box("DopedSiliconLayer",		//its name
			     DopSiThickness/2, SiliconSizeY/2,SiliconSizeX/2); 
      
      logicDopSi = new G4LogicalVolume(solidDopSi,    //its solid
				       SiliconMaterial, //its material
				       "DopedSilicon"); //name
      
      physiDopSi = new G4PVPlacement(0,		   //no rotation
				     G4ThreeVector(-((SiliconThickness/2)+(ThermalOThickness)+(DopSiThickness/2))  ,0.,0.0*cm),  //its position
				     logicDopSi,     //its logical volume		    
				     "DopedSilicon", //its name
				     logicWorld,        //its mother
				     false,             //no boulean operat
				     0);                //copy number
    }
  
  if (AlThickness > 0.) 
    {
      solidAl = new G4Box("AlLayer",		//its name
			  AlThickness/2, SiliconSizeY/2,SiliconSizeX/2); 
      
      logicAl = new G4LogicalVolume(solidAl,    //its solid
				    AlMaterial, //its material
				    "Al"); //name
      
      physiAl = new G4PVPlacement(0,		   //no rotation
				  G4ThreeVector(-((SiliconThickness/2)+(ThermalOThickness)+(DopSiThickness)+(AlThickness/2))  ,0.,0.0*cm),  //its position
				  logicAl,     //its logical volume		    
				  "Al", //its name
				  logicWorld,        //its mother
				  false,             //no boulean operat
				  0);                //copy number      
    }

  
   
   G4VSolid* cellSolid = new G4Box("Cell_Solid", // Name                                                             
  				  (SiliconThickness/2),         // x half length                                                    
  				  (ypixel_pitch/2.0),         // y half length                                                    
  				  (xpixel_pitch/2.00));      // z half length                                                    
 

  G4LogicalVolume* cellLogical   
    = new G4LogicalVolume(cellSolid,       // Solid                                                                 
  			  SiliconMaterial,             // Material                                                              
                          "Cell_Logical"); // Name                                                                  


  
   G4VSolid* stripSolid = new G4Box("strip_Solid", // Name                                                             
  				  (SiliconThickness/2),         // x half length                                                    
  				  (SiliconSizeX/2.0),         // y half length                                                    
  				  (xpixel_pitch/2.00));      // z half length                                                    
 

  G4LogicalVolume* stripLogi
    = new G4LogicalVolume(stripSolid,       // Solid                                                                 
  			  SiliconMaterial,             // Material                                                              
                          "Strip_Logical"); // Name                                                                  
 
  // //G4VPVParameterisation* cellParam = new AntiPTestCellParameterisation();

  G4VPhysicalVolume* stripPhysi=new G4PVReplica("Cell_Physical",    // Name                                                                 
  			stripLogi,                                                          
  			physiSilicon,                                                        
  			kZAxis,                                                               
  			xpixels,                                                                    
  			1.408*cm/xpixels);         
  
  // //her legges alle cellene inn i silicon volumet
  G4VPhysicalVolume* physiStrip=new G4PVReplica("Strip_Physical",    // Name                                                                 
  			cellLogical,                                                          
  			stripPhysi,                                                        
  			kYAxis,                                                               
  			xpixels,                                                                    
  			1.408*cm/xpixels);         
  

  
  G4VSensitiveDetector* detector = new AntiPSD("/mydet/Silicon");
  // Get pointer to detector manager                                                     
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  // Register detector with manager                                                      
  SDman->AddNewDetector(detector);
  // Attach detector to volume defining calorimeter cells                                
  cellLogical->SetSensitiveDetector(detector);
                                         
  //Visualization attributes
  
  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);
  //cellLogical->SetVisAttributes (G4VisAttributes::Invisible);
  // logicSilicon->SetVisAttributes(G4VisAttributes::Invisible);
  //  logicThermalO->SetVisAttributes(G4VisAttributes::Invisible);
  //  logicDopSi->SetVisAttributes(G4VisAttributes::Invisible);
  //  logicAl->SetVisAttributes(G4VisAttributes::Invisible);
  
  return physiWorld;
}

//------------------------------------------------------------------------------

void DetectorConstruction::SetSiliconMaterial(G4String materialChoice)
{
	// search the material by its name   
	G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
	if (pttoMaterial) SiliconMaterial = pttoMaterial;
}

//------------------------------------------------------------------------------

void DetectorConstruction::SetSiliconThickness(G4double val)
{
	// change Silicon thickness and recompute the calorimeter parameters
	SiliconThickness = val;
}

//------------------------------------------------------------------------------
void DetectorConstruction::ComputeWorldParameters()
{
  WorldSizeX = 200*cm; WorldSizeYZ = 200*cm;
}


//------------------------------------------------------------------------------
