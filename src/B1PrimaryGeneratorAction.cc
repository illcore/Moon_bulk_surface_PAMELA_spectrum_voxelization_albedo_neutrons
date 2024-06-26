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
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <iostream>
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="proton");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,-1,0));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
//69-bin PAMELA proton spectrum between 0.4 and 115 GeV by AIT
double A[]= {0, 0.0567745244255567,
0.055625387941295,
0.054835356608365,
0.0532912044576383,
0.0518547838523111,
0.0503465422167176,
0.0484432849146591,
0.046432296067201,
0.0442776651592103,
0.0420512132209532,
0.0399684033432287,
0.037813772435238,
0.0357668730726467,
0.0335045106192564,
0.0312780586809993,
0.0291952488032749,
0.02722017047095,
0.0252450921386252,
0.0234352021759129,
0.0216109480071474,
0.0199482911564812,
0.0183143627179215,
0.0169497631428607,
0.0155959367223398,
0.0142385192503057,
0.0130139706842642,
0.011904335766649,
0.01088806818838,
0.00991130217675753,
0.0089848108863215,
0.00815168693523174,
0.00739038401440833,
0.00672244843293119,
0.00577081978190194,
0.00467554907033997,
0.0037813772435238,
0.00304521168329362,
0.0024419150290562,
0.00195640486445562,
0.00155600262072067,
0.00123532172058137,
0.000979638852833136,
0.00077638533717934,
0.000610478757264051,
0.000479405377027946,
0.000377419514049716,
0.000297339065302726,
0.000232341032911671,
0.000179911680817229,
0.0001400510090194,
0.000109527071156197,
8.50001893202358E-05,
6.56444216634521E-05,
5.09929314891148E-05,
3.97529402524297E-05,
3.04880273480694E-05,
2.34495663819662E-05,
1.82784522027883E-05,
1.42564745078722E-05,
1.0809065055087E-05,
8.36715002603082E-06,
6.39207169370594E-06,
5.0633826337783E-06,
3.96092981918969E-06,
2.98416380756721E-06,
2.27313560793026E-06,
1.64829264461294E-06,
1.1204080721552E-06}; 
//Discrete energies
double B[]= {0.0, 440,
480,
520,
560,
600,
650,
700,
750,
800,
860,
920,
990,
1060,
1130,
1210,
1290,
1380,
1470,
1570,
1670,
1780,
1890,
2010,
2130,
2270,
2400,
2550,
2700,
2860,
3030,
3210,
3390,
3590,
3900,
4350,
4840,
5380,
5980,
6640,
7360,
8150,
9010,
9970,
11010,
12160,
13410,
14790,
16310,
17960,
19780,
21780,
23970,
26370,
29000,
31880,
35050,
38520,
42320,
46490,
51070,
56080,
61580,
67610,
74220,
81470,
89420,
100410,
115420};
double Grid[69];
double sum=0;
  for(int x=0; x < sizeof(Grid)/sizeof(Grid[0]); x++){
  sum=sum+A[x];
  Grid[x]=sum;
  std::ofstream GridFile;
  GridFile.open("Probability_grid.txt", std::ios::app);
  GridFile <<  Grid[x] << G4endl;
  GridFile.close();
  }  
  for (int n_particle = 1; n_particle < 100000; n_particle++){
  G4double y0 = 3*m;
  G4double z0 = 10*cm;
  G4double x0 = 10*cm; 
  x0 = -x0+2*x0*G4UniformRand();
  z0 = -z0+2*z0*G4UniformRand();
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  G4double Energy; //Just for initialization
  G4double pseudo=G4UniformRand();
  for (int i=0; i < sizeof(Grid)/sizeof(Grid[0]); i++){
  if(pseudo > Grid[i] && pseudo <= Grid[i+1]){
  Energy=B[i+1];
  std::ofstream EnergyFile;
  EnergyFile.open("Energy.txt", std::ios::app);
  EnergyFile <<  Energy << G4endl;
  EnergyFile.close();
  }
  }   
  fParticleGun->SetParticleEnergy(Energy);
  fParticleGun->GeneratePrimaryVertex(anEvent);
 }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

