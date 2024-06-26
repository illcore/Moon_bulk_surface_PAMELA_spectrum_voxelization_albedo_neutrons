//
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
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
     
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

//     
// World
//  
G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4Box* solidWorld =    
    new G4Box("WorldBox",                       //its name
      3*m, 3*m,3*m);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "LogicalWorld");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                      
//Moon surface
G4Material* Oxygen = nist->FindOrBuildMaterial("G4_O");
G4Material* Sodium = nist->FindOrBuildMaterial("G4_Na");
G4Material* Magnesium = nist->FindOrBuildMaterial("G4_Mg");
G4Material* Aluminum = nist->FindOrBuildMaterial("G4_Al");
G4Material* Silicon = nist->FindOrBuildMaterial("G4_Si");       
G4Material* Calcium = nist->FindOrBuildMaterial("G4_Ca");       
G4Material* Titanium = nist->FindOrBuildMaterial("G4_Ti");       
G4Material* Manganese = nist->FindOrBuildMaterial("G4_Mn");       
G4Material* Iron = nist->FindOrBuildMaterial("G4_Fe"); 
G4String name;
G4int ncomponents;
G4double fracMass;
G4Material* bulk_mat = new G4Material(name = "bulk_mat", 1.93*g/cm3, ncomponents = 9);
bulk_mat->AddMaterial(Oxygen, fracMass = 43.7 * perCent);
bulk_mat->AddMaterial(Sodium, fracMass = 0.3 * perCent);
bulk_mat->AddMaterial(Magnesium, fracMass = 5.6 * perCent);         
bulk_mat->AddMaterial(Aluminum, fracMass = 9 * perCent);         
bulk_mat->AddMaterial(Silicon, fracMass = 21.1 * perCent);         
bulk_mat->AddMaterial(Calcium, fracMass = 8.5 * perCent);         
bulk_mat->AddMaterial(Titanium, fracMass = 1.5 * perCent);
bulk_mat->AddMaterial(Manganese, fracMass = 0.1 * perCent);
bulk_mat->AddMaterial(Iron, fracMass = 10.2 * perCent);
//200-cm bulk, density=1.93
G4Box* Layer_box =    
    new G4Box("Layer_box",                       //its name
      1.5*m, 1*cm,1.5*m);     //its size

G4LogicalVolume* Bulk =                         
    new G4LogicalVolume(Layer_box,          //its solid
                        bulk_mat,           //its material
                        "Bulk_logic");            //its name
                                   
for (int i = 1; i < 101; i++){
                         new G4PVPlacement(0, G4ThreeVector(0*m,(-100+(2*i))*cm, 0*m), Bulk, "Bulk", logicWorld, false, i, checkOverlaps);                      
  }
//Pseudo detector
G4ThreeVector Pos_surface = G4ThreeVector(0*m, 1.011*m, 0*m);
G4Box* Surface_detector_box =    
    new G4Box("Surface_detector_box",                      
      1.5*m, 1*mm,1.5*m);     
      
G4LogicalVolume* Surface_detector =                         
    new G4LogicalVolume(Surface_detector_box,         
                        world_mat,          
                        "Surface_detector_logic");         
                                   
    new G4PVPlacement(0,                    
                      Pos_surface,      
                      Surface_detector,           
                      "Surface_detector",               
                      logicWorld,                     
                      false,                 
                      0,                     
                      checkOverlaps);  
G4VisAttributes* visAttributesbulk = new G4VisAttributes(G4Colour(0.65, 0.16, 0.16));
visAttributesbulk->SetVisibility(true);
Bulk->SetVisAttributes(visAttributesbulk);                            
  //always return the physical World
  //  
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
