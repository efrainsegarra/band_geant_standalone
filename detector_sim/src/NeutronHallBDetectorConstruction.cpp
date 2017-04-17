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
    G4Material* BAND_mat        = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    
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

  

  // ------------------ Defining Dubna NeuLand Detector - in World ---------------- //

    G4double DUBNA_offset = 22.5*cm;
    G4double DUBNA_length = 250.*cm;
    G4double DUBNA_width = 250.*cm;
    G4double DUBNA_barCrossSection = 5.*cm;
    G4double DUBNA_thickness = 40*cm;


    bool doSegment = true;

    if(doSegment == true){

      const int DUBNA_wallReplicas = 8;
      G4int DUBNA_barReplicas = 50;
      
      G4Box                 *solidDetector[DUBNA_wallReplicas];                    // solid 
      G4VPhysicalVolume     *physDetector[DUBNA_wallReplicas];                     // the physical volume
      G4ThreeVector         posDetector[DUBNA_wallReplicas];                       // the vector for the position of the boxes
      
      G4Colour Detec_Color( 0.75 , 0.6 , 0.75 );
      G4VisAttributes* DetecVisAttributes = new G4VisAttributes( true , Detec_Color );
      fVisAttributes.push_back(DetecVisAttributes);

      for (int wall = 0; wall < DUBNA_wallReplicas; wall++) {

          posDetector[wall]     = G4ThreeVector(0,0,(-DUBNA_offset - DUBNA_barCrossSection*wall));

          solidDetector[wall] =                                      // solid volume for each box
              new G4Box("Detector",                                  // name
                        0.5*DUBNA_length, 
                        0.5*DUBNA_width, 
                        0.5*DUBNA_barCrossSection);                                  // ending phi 


          logicDetectorWall[wall] = 
              new G4LogicalVolume(solidDetector[wall],
                              BAND_mat,
                              "DetectorWallLV");
        
          physDetector[wall]=                                        // physical volume for TOF
              new G4PVPlacement(0,                               // rotation
                                G4ThreeVector(0,0,(-DUBNA_offset - DUBNA_barCrossSection*wall) ),                 // at 0,0,0
                                logicDetectorWall[wall],                 // logical volume of each TOF
                                "DetectorWallPhys",                           // name of each TOF
                                logicWorld,                      // mother volume for TOF
                                false,                           // no boolean operation
                                0,                               // copy number
                                fCheckOverlaps);                 // checking overlap
      
            if(wall%2==0){ 
              // dubna long-y extent bars (50 bars in x direction, defining xy look)
              G4Box* solidDetectorBars =
                  new G4Box("DetectorBars",
                            0.5*DUBNA_barCrossSection,
                            0.5*DUBNA_width,
                            0.5*DUBNA_barCrossSection);

              logicDetectorBars[wall] = 
                  new G4LogicalVolume(solidDetectorBars,
                                      BAND_mat,
                                      "DetectorBarsLV");

                new G4PVReplica("DetectorBarsPhys",logicDetectorBars[wall],logicDetectorWall[wall],kXAxis,DUBNA_barReplicas,DUBNA_barCrossSection); 
            }
            else{
              // dubna long-y extent bars (50 bars in x direction, defining xy look)
              G4Box* solidDetectorBars =
                  new G4Box("DetectorBars",
                            0.5*DUBNA_length,
                            0.5*DUBNA_barCrossSection,
                            0.5*DUBNA_barCrossSection);

              logicDetectorBars[wall] = 
                  new G4LogicalVolume(solidDetectorBars,
                                      BAND_mat,
                                      "DetectorBarsLV");
                  
                new G4PVReplica("DetectorBarsPhys",logicDetectorBars[wall],logicDetectorWall[wall],kYAxis,DUBNA_barReplicas,DUBNA_barCrossSection); 
            }
          
          logicDetectorBars[wall] -> SetVisAttributes(DetecVisAttributes);
      }
    }
    else{
      G4Box* solidDetector =                                      // solid volume for each box
              new G4Box("Detector",                                  // name
                        0.5*DUBNA_length, 
                        0.5*DUBNA_width, 
                        0.5*DUBNA_thickness);                                  // ending phi 


      logicDetector = 
          new G4LogicalVolume(solidDetector,
                          BAND_mat,
                          "DetectorLV");


                                            // physical volume for TOF
              new G4PVPlacement(0,                               // rotation
                                G4ThreeVector(0,0,-DUBNA_thickness ),                 // at 0,0,0
                                logicDetector,                 // logical volume of each TOF
                                "DetectorPhys",                           // name of each TOF
                                logicWorld,                      // mother volume for TOF
                                false,                           // no boolean operation
                                0,                               // copy number
                                fCheckOverlaps);                 // checking overlap
    }
    
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

  bool doSegment = true;

  if(doSegment == true){

    for (int wall = 0; wall < 8; wall++) {
      logicDetectorBars[wall]->SetSensitiveDetector(bandSD);
    }
  }
  else{
    logicDetector->SetSensitiveDetector(bandSD);

  }


  //G4String nameSD = "BAND";
  //NeutronHallBTrackerSD* bandSD = new NeutronHallBTrackerSD(treePtr,nameSD);
  //SetSensitiveDetector("DetectorLV" , bandSD , true);
  
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
