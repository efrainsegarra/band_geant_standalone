// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBDetectorConstruction.cc
/// \brief Implementation of the NeutronHallBDetectorConstruction class

#include "NeutronHallBDetectorConstruction.h"
#include "NeutronHallBDetectorMessenger.h"
#include "NeutronHallBTrackerSD.h"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4GenericMessenger.hh"

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"

#include "G4Sphere.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Trd.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"

#include "G4SystemOfUnits.hh"



G4ThreadLocal G4GlobalMagFieldMessenger* NeutronHallBDetectorConstruction::fMagFieldMessenger = 0;

NeutronHallBDetectorConstruction::NeutronHallBDetectorConstruction(void * t)
 :G4VUserDetectorConstruction(),
  fNbOfChambers(0),
  fLogicTarget(NULL), fLogicChamber(NULL),
  fTargetMaterial(NULL), fChamberMaterial(NULL),
  fStepLimit(NULL),
  fVisAttributes(),
  fCheckOverlaps(true){
  fMessenger = new NeutronHallBDetectorMessenger(this);
  
  fNbOfChambers = 5;
  fLogicChamber = new G4LogicalVolume*[fNbOfChambers];
  
  treePtr = t;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
NeutronHallBDetectorConstruction::~NeutronHallBDetectorConstruction(){
  delete [] fLogicChamber;
  delete fStepLimit;
  delete fMessenger;
  for (auto visAttributes: fVisAttributes) {
    delete visAttributes;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* NeutronHallBDetectorConstruction::Construct(){
  // Define materials
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHallBDetectorConstruction::DefineMaterials(){
    // Material definition
    
    G4NistManager* man = G4NistManager::Instance();
    
    // Define StainlessSteel not in NIST, following http://hypernews.slac.stanford.edu/HyperNews/geant4/get/geometry/915.html?inline=1
    G4Element* C  = man->FindOrBuildElement("C");
    G4Element* Si = man->FindOrBuildElement("Si");
    G4Element* Cr = man->FindOrBuildElement("Cr");
    G4Element* Mn = man->FindOrBuildElement("Mn");
    G4Element* Fe = man->FindOrBuildElement("Fe");
    G4Element* Ni = man->FindOrBuildElement("Ni");
    
    G4Material* StainlessSteel = new G4Material( "StainlessSteel", 8.06*g/cm3, 6 );
    StainlessSteel->AddElement(C, 0.001);
    StainlessSteel->AddElement(Si, 0.007);
    StainlessSteel->AddElement(Cr, 0.18);
    StainlessSteel->AddElement(Mn, 0.01);
    StainlessSteel->AddElement(Fe, 0.712);
    StainlessSteel->AddElement(Ni, 0.09);
    
    G4Material* Air = man->FindOrBuildMaterial("G4_AIR");
    G4Material* vacuum = new G4Material("Vacuum", 1.e-25*g/cm3, 1, kStateGas, 2.e-2*bar, CLHEP::STP_Temperature);
    vacuum -> AddMaterial(Air, 1.);
    
    
    
    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* NeutronHallBDetectorConstruction::DefineVolumes(){
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* vacuum          = nist->FindOrBuildMaterial("Vacuum");
    G4Material* default_mat     = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* SS_mat          = nist->FindOrBuildMaterial("StainlessSteel");
    G4Material* LeadWall_mat    = nist->FindOrBuildMaterial("G4_Pb");
    G4Material* BAND_mat        = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    
    G4double length;
    G4double shift;

    // ------------------ Defining World ---------------- //
    G4double world_sizeX = 5000.*cm;
    G4double world_sizeY = 5000.*cm;
    G4double world_sizeZ = 5000.*cm;

  
    G4Box* solidWorld =    
        new G4Box("World",                                     // its name
                  0.5*world_sizeX, 
                  0.5*world_sizeY, 
                  0.5*world_sizeZ);                            // its size
      
    G4LogicalVolume* logicWorld =                         
        new G4LogicalVolume(solidWorld,                        // its solid
                            vacuum,                            // its world material
                            "World");                          // its name
                                   
    G4VPhysicalVolume* physWorld = 
        new G4PVPlacement(0,                                   // no rotation
                          G4ThreeVector(),                     // at (0,0,0)
                          logicWorld,                          // its logical volume
                          "World",                             // its name
                          0,                                   // its mother  volume
                          false,                               // no boolean operation
                          0,                                   // copy number
                          fCheckOverlaps);                     // overlaps checking
  // ------------------ End Defining World ---------------- //

  
  // ------------------ Defining BAND Detector - in World ---------------- //

    const int BAND_wallReplicas =1;
    const int BAND_wallGroups = 1;
    
    G4double BAND_barCrossSection = 7.4*cm;
    G4double BAND_groupE_length = 2018.35*mm; // temp
    
    G4double BAND_Z_coord_center = -2805.*mm;
    G4double BAND_Z_thickness = BAND_barCrossSection * BAND_wallReplicas;
    G4double BAND_Z_start = BAND_Z_coord_center + (BAND_Z_thickness/2.)*mm;


    G4int numBars;
    G4double wall_Z_offset;
    // create 5 box walls that will hold the BAND array
    for (int wall = 0; wall < BAND_wallReplicas; wall++){
      wall_Z_offset = (BAND_Z_start - BAND_barCrossSection*(1./2+wall));

      G4Box* solidDetectorBarsE =
                new G4Box("DetectorBars_E",
                          0.5*BAND_groupE_length,
                          0.5*BAND_barCrossSection,
                          0.5*BAND_barCrossSection);
      logicDetectorBarsE[wall] = 
          new G4LogicalVolume(solidDetectorBarsE,
                              BAND_mat,
                              "DetectorBarsLV_E"); 
             

      // now for each wall, need to parameterize it into the 5 groups.
      for (int group = 0; group < BAND_wallGroups; group++){
        if    (group == 0){
          numBars = 1;
          for(int barNo = 0; barNo < numBars; barNo++){
                new G4PVPlacement(0, G4ThreeVector(0,0, wall_Z_offset), logicDetectorBarsE[wall],                 
                            "BarsA", logicWorld, false, barNo+24*wall, fCheckOverlaps);
          }
        }
      }
    }
  // ------------------ END Defining BAND Detector  ---------------- //
    
  
    // Always return the physical world
    return physWorld;
}

void NeutronHallBDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors
  auto sdManager = G4SDManager::GetSDMpointer();
  G4String SDname;

  NeutronHallBTrackerSD* bandSD = new NeutronHallBTrackerSD(treePtr,SDname="band");
  sdManager->AddNewDetector(bandSD);

  for (int wall = 0; wall < 1; wall++) {
      logicDetectorBarsE[wall]->SetSensitiveDetector(bandSD);        
  }

  // Electro-magnetic fields (currently zero)  
  G4ThreeVector fieldValue = G4ThreeVector();
  G4AutoDelete::Register(fMagFieldMessenger);
}

void NeutronHallBDetectorConstruction::SetTargetMaterial(G4String materialName)
{
    G4NistManager* nistManager = G4NistManager::Instance();
    
    G4Material* pttoMaterial =
    nistManager->FindOrBuildMaterial(materialName);
    
    if (fTargetMaterial != pttoMaterial) {
        if ( pttoMaterial ) {
            fTargetMaterial = pttoMaterial;
            if (fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
            //        G4cout << "\n----> The target is made of " << materialName << G4endl;
        } else {
            //        G4cout << "\n-->  WARNING from SetTargetMaterial : "
            //               << materialName << " not found" << G4endl;
        }
    }
}

void NeutronHallBDetectorConstruction::SetChamberMaterial(G4String materialName)
{
    G4NistManager* nistManager = G4NistManager::Instance();
    
    G4Material* pttoMaterial =
    nistManager->FindOrBuildMaterial(materialName);
    
    if (fChamberMaterial != pttoMaterial) {
        if ( pttoMaterial ) {
            fChamberMaterial = pttoMaterial;
            for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++) {
                if (fLogicChamber[copyNo]) fLogicChamber[copyNo]->
                    SetMaterial(fChamberMaterial);
            }
            //        G4cout << "\n----> The chambers are made of " << materialName << G4endl;
        } else {
            //        G4cout << "\n-->  WARNING from SetChamberMaterial : "
            //               << materialName << " not found" << G4endl;
        }
    }
}


void NeutronHallBDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

void NeutronHallBDetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}  
