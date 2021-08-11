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
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4IntersectionSolid.hh"
#include "HexParameterisation.hh"
#include "G4PVParameterised.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreadLocal
G4GlobalMagFieldMessenger* B1DetectorConstruction::fMagFieldMessenger = 0;



B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0),
  fCheckOverlaps(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{
    G4String name, symbol;
    G4double z, a, density,fractionmass;
    G4int ncomponents, natoms;
    G4NistManager* manager = G4NistManager::Instance();
    //vaccum
    G4Material* Vaccum = manager->FindOrBuildMaterial("G4_Galactic");
    G4Element* elSi = manager->FindOrBuildElement(14, 28); //silicon
    G4Element* elNa = manager->FindOrBuildElement(11, 23); //Na
    G4Element* elCs = manager->FindOrBuildElement(55, 133); // Cs
    G4Element* elBa = manager->FindOrBuildElement(56, 138); //Ba
    G4Element* elGd152 = manager->FindOrBuildElement(64, 152); //Gd
    G4Element* elGd154 = manager->FindOrBuildElement(64, 154);
    G4Element* elGd155 = manager->FindOrBuildElement(64, 155);
    G4Element* elGd156 = manager->FindOrBuildElement(64, 156);
    G4Element* elGd157 = manager->FindOrBuildElement(64, 157);
    G4Element* elGd158 = manager->FindOrBuildElement(64, 158);
    G4Element* elGd160 = manager->FindOrBuildElement(64, 160);
    G4Element* elPb = manager->FindOrBuildElement(82, 207); //Pb
    G4Element* elBi = manager->FindOrBuildElement(83, 209); //Bi
    G4Element* elTi = manager->FindOrBuildElement(22, 48); //Ti
    G4Element* elO = manager->FindOrBuildElement(8, 16); //O
    density = 4.29*g/cm3;
    G4Material* LeadGlass = new G4Material(name="LeadGlass",density,ncomponents = 15);
    LeadGlass->AddElement(elSi, fractionmass=16.55*perCent);
    LeadGlass->AddElement(elNa, fractionmass=0.96*perCent);
    LeadGlass->AddElement(elCs, fractionmass=2.45*perCent);
    LeadGlass->AddElement(elBa, fractionmass=5.02*perCent);
    LeadGlass->AddElement(elGd152, fractionmass=0.01718*perCent);
    LeadGlass->AddElement(elGd154, fractionmass=0.18039*perCent);
    LeadGlass->AddElement(elGd155, fractionmass=1.27132*perCent);
    LeadGlass->AddElement(elGd156, fractionmass=1.76954*perCent);
    LeadGlass->AddElement(elGd157, fractionmass=1.34863*perCent);
    LeadGlass->AddElement(elGd158, fractionmass=2.13032*perCent);
    LeadGlass->AddElement(elGd160, fractionmass=1.87262*perCent);
    LeadGlass->AddElement(elPb, fractionmass=37.6*perCent);
    LeadGlass->AddElement(elBi, fractionmass=0.45*perCent);
    LeadGlass->AddElement(elTi, fractionmass=2.64*perCent);
    LeadGlass->AddElement(elO, fractionmass=25.75*perCent);

    //defining world environment
    G4double worldSizeXYZ = 1.*mm;
    G4Box* solidWorld =
            new G4Box("World", 0.5*worldSizeXYZ, 0.5*worldSizeXYZ, 0.5*worldSizeXYZ);
    G4LogicalVolume* worldLog = new G4LogicalVolume(solidWorld, Vaccum, "World");

    /*G4VisAttributes *worldAttributes = new G4VisAttributes;
    worldAttributes->SetVisibility(false);
    worldLog->SetVisAttributes(worldAttributes);*/

    G4VPhysicalVolume* physWorld =
            new G4PVPlacement(0, G4ThreeVector(), worldLog, "World", 0, false, 0, fCheckOverlaps);

    //defining first mcp
    G4double innerRadius = 0.*um;
    G4double outerRadius = 1000.*um;
    G4double mcp_h = 540. * um;             //mcp thickness unit:mm
    G4double startAngle = 0.*deg;
    G4double spanningAngle = 360.*deg;

    G4Tubs* mcpTube = new G4Tubs("mcp", innerRadius, outerRadius, mcp_h/2, startAngle, spanningAngle);
    G4LogicalVolume* mcpLog = new G4LogicalVolume(mcpTube, LeadGlass, "mcp");

    /*G4VisAttributes *mcpAttributes = new G4VisAttributes;
    mcpAttributes->SetVisibility(true);
    mcpLog->SetVisAttributes(mcpAttributes);*/

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0., 0., 0.),
                      mcpLog,
                      "mcp",
                      worldLog,
                      false,
                      0,
                      fCheckOverlaps);

    //defining pore
    G4double poreRadius = 5.*um;
    G4double poreL = 620.*um;
    G4double biasAngle = 8.*deg;
    G4Tubs* poreTube = new G4Tubs("poreTube", innerRadius, poreRadius, poreL/2, startAngle, spanningAngle);
    G4RotationMatrix* rm = new G4RotationMatrix();
    rm->rotateX(biasAngle);
    G4VSolid* pore = new G4IntersectionSolid("pore",
                                             mcpTube, poreTube, rm, G4ThreeVector(0.,0.,0.));
    G4LogicalVolume* poreLog = new G4LogicalVolume(pore, Vaccum, "pore");

    G4double pore_pitch = 11.8427*um;
    G4int ncolumns = 100;
    G4int nrow = 100;
    G4int ncells =  ncolumns * nrow;
    G4double first_position_X = - pore_pitch * (ncolumns - 1) / 2;
    G4double first_position_Y = - pore_pitch * sqrt(3) / 2 * (nrow - 1) / 2;
    G4double fisrt_position_Z = 0;
    G4VPVParameterisation* poreHexParameterisation =
            new HexParameterisation (ncells,
                                     ncolumns,
                                     nrow,
                                     first_position_X,
                                     first_position_Y,
                                     fisrt_position_Z,
                                     pore_pitch);

    new G4PVParameterised("poreArray",
                          poreLog,
                          mcpLog,
                          kUndefined,
                          ncells,
                          poreHexParameterisation,
                          0);




  fScoringVolume = mcpLog;

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::ConstructSDandField()
{
	// Sensitive detectors

	// Create global magnetic field messenger.
	// Uniform magnetic field is then created automatically if
	// the field value is not zero.
	G4ThreeVector fieldValue = G4ThreeVector(0.,0.,-1*tesla);
	fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
	fMagFieldMessenger->SetVerboseLevel(1);
	// Register the field messenger for deleting
	G4AutoDelete::Register(fMagFieldMessenger);
}