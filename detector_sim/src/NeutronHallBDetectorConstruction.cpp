// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBDetectorConstruction.cc
/// \brief Implementation of the NeutronHallBDetectorConstruction class

#include "NeutronHallBDetectorConstruction.h"
#include "NeutronHallBDetectorMessenger.h"
#include "NeutronHallBTrackerSD.h"

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"

#include "G4Sphere.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"

#include "G4SystemOfUnits.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"

G4ThreadLocal G4GlobalMagFieldMessenger* NeutronHallBDetectorConstruction::fMagFieldMessenger = 0;

NeutronHallBDetectorConstruction::NeutronHallBDetectorConstruction(void * t)
 :G4VUserDetectorConstruction(),
  fNbOfChambers(0),
  fLogicTarget(NULL), fLogicChamber(NULL),
  fTargetMaterial(NULL), fChamberMaterial(NULL),
  fStepLimit(NULL),
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
    
    bool DoTarget           = false;
    bool DoSSTube           = true;
    bool DoElectronicBoxes  = true;
    bool DoLeadWall         = true;
    bool DoBandAsFullCircle = true;
    bool DoBandAsProposed   = !DoBandAsFullCircle;
    
    
    //--------- World --------------//
    G4double worldXYZ =  8000.*cm;
    G4Box* solidWorld = new G4Box("World", 0.5*worldXYZ, 0.5*worldXYZ, 0.5*worldXYZ );
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, vacuum, "World");
    G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld,  "World", 0,  false,  0, fCheckOverlaps);
    
    
    //--------- Target --------------//
    if (DoTarget){
        G4Sphere* solid_target = new G4Sphere ("target", 0 , 1.*cm , 0., CLHEP::twopi, 0., CLHEP::twopi/2.);
        G4LogicalVolume* logic_target = new G4LogicalVolume(solid_target, default_mat , "target");
        new G4PVPlacement( 0, G4ThreeVector(),  logic_target, "target", logicWorld, false,  0,  fCheckOverlaps);
        logic_target -> SetVisAttributes (G4Colour( 0.65 , 0.82 , 0.77 ));
    }
    
    //--------- Stainless-Steel --------------//
    if (DoSSTube){
        G4double length;
        G4double shift;
        //Numbers come from counting pixels in Hall-B sketches
        
        length = 658.*mm;
        G4Tubs* solid_SSTube_P1 = new G4Tubs("SSTube_P1", 210.*mm , 213.*mm, length , 0., CLHEP::twopi);
        G4LogicalVolume* logic_SSTube_P1 = new G4LogicalVolume(solid_SSTube_P1, SS_mat , "SSTube_P1");
        new G4PVPlacement( 0, G4ThreeVector(),  logic_SSTube_P1, "SSTube_P1", logicWorld, false,  0,  fCheckOverlaps);
        logic_SSTube_P1 -> SetVisAttributes (G4Colour( 0.99 , 0.88 , 0.66 ));
        
        if (DoElectronicBoxes){
        // 1 cm plastic coating around the pipe to represent the cables going to the electronic box
        G4Tubs* solid_SSTube_P1_coating = new G4Tubs("SSTube_P1_coating", 214.*mm , 224.*mm, length , 0., CLHEP::twopi);
        G4LogicalVolume* logic_SSTube_P1_coating = new G4LogicalVolume(solid_SSTube_P1_coating, BAND_mat , "SSTube_P1_coating");
        new G4PVPlacement( 0, G4ThreeVector(),  logic_SSTube_P1_coating, "SSTube_P1_coating", logicWorld, false,  0,  fCheckOverlaps);
        logic_SSTube_P1_coating -> SetVisAttributes (G4Colour( 0.75 , 0.6 , 0.75 ));
        }
        
        
        
        length = 24.*mm;
        shift = 658.*mm;
        G4Tubs* solid_SSTube_Connector1 = new G4Tubs("SSTube_Connector1", 213.*mm , 240.*mm, length , 0., CLHEP::twopi);
        G4LogicalVolume* logic_SSTube_Connector1 = new G4LogicalVolume(solid_SSTube_Connector1, SS_mat , "SSTube_Connector1");
        new G4PVPlacement( 0, G4ThreeVector(0,0,length+shift),  logic_SSTube_Connector1, "SSTube_Connector1", logicWorld, false,  0,  fCheckOverlaps);
        logic_SSTube_Connector1 -> SetVisAttributes (G4Colour( 0.99 , 0.88 , 0.66 ));

        if (DoElectronicBoxes){
        // 1 cm plastic coating around the pipe to represent the cables going to the electronic box
        G4Tubs* solid_SSTube_Connector1_coating = new G4Tubs("SSTube_Connector1_coating", 241.*mm , 251.*mm, length , 0., CLHEP::twopi);
        G4LogicalVolume* logic_SSTube_Connector1_coating = new G4LogicalVolume(solid_SSTube_Connector1_coating, BAND_mat , "SSTube_Connector1_coating");
        new G4PVPlacement( 0, G4ThreeVector(0,0,length+shift),  logic_SSTube_Connector1_coating, "SSTube_Connector1_coating", logicWorld, false,  0,  fCheckOverlaps);
        logic_SSTube_Connector1_coating -> SetVisAttributes (G4Colour( 0.75 , 0.6 , 0.75 ));
        }
        
        
        
        
        length = 34.*mm;
        shift = 658.*mm;
        G4Tubs* solid_SSTube_Connector2 = new G4Tubs("SSTube_Connector2", 240.*mm , 258*mm, length , 0., CLHEP::twopi);
        G4LogicalVolume* logic_SSTube_Connector2 = new G4LogicalVolume(solid_SSTube_Connector2, SS_mat , "SSTube_Connector2");
        new G4PVPlacement( 0, G4ThreeVector(0,0,length+shift),  logic_SSTube_Connector2, "SSTube_Connector2", logicWorld, false,  0,  fCheckOverlaps);
        logic_SSTube_Connector2 -> SetVisAttributes (G4Colour( 0.99 , 0.88 , 0.66 ));
        
        if (DoElectronicBoxes){
            // 1 cm plastic coating around the pipe to represent the cables going to the electronic box
            G4Tubs* solid_SSTube_Connector2_coating = new G4Tubs("SSTube_Connector2_coating", 259.*mm , 269*mm, length , 0., CLHEP::twopi);
            G4LogicalVolume* logic_SSTube_Connector2_coating = new G4LogicalVolume(solid_SSTube_Connector2_coating, BAND_mat , "SSTube_Connector2_coating");
            new G4PVPlacement( 0, G4ThreeVector(0,0,length+shift),  logic_SSTube_Connector2_coating, "SSTube_Connector2_coating", logicWorld, false,  0,  fCheckOverlaps);
            logic_SSTube_Connector2_coating -> SetVisAttributes (G4Colour( 0.75 , 0.6 , 0.75 ));
        }
        
        
        length = 1000.*mm;
        shift = 682.*mm;
        G4Tubs* solid_SSTube_P2 = new G4Tubs("SSTube_P2", 237.*mm , 240*mm, length , 0., CLHEP::twopi);
        G4LogicalVolume* logic_SSTube_P2 = new G4LogicalVolume(solid_SSTube_P2, SS_mat , "SSTube_P2");
        new G4PVPlacement( 0, G4ThreeVector(0,0,length+shift),  logic_SSTube_P2, "SSTube_P2", logicWorld, false,  0,  fCheckOverlaps);
        logic_SSTube_P2 -> SetVisAttributes (G4Colour( 0.99 , 0.88 , 0.66 ));
        
        if (DoElectronicBoxes){
            // 1 cm plastic coating around the pipe to represent the cables going to the electronic box
            
            G4Tubs* solid_SSTube_P2_coating = new G4Tubs("SSTube_P2_coating", 241.*mm , 251*mm, length , 0., CLHEP::twopi);
            G4LogicalVolume* logic_SSTube_P2_coating = new G4LogicalVolume(solid_SSTube_P2_coating, BAND_mat , "SSTube_P2_coating");
            new G4PVPlacement( 0, G4ThreeVector(0,0,length+shift),  logic_SSTube_P2_coating, "SSTube_P2", logicWorld, false,  0,  fCheckOverlaps);
            logic_SSTube_P2_coating -> SetVisAttributes (G4Colour( 0.75 , 0.6 , 0.75 ));
        }
    
    }
    
        //--------- electronics boxes --------------//
    if (DoElectronicBoxes){
    
        G4double BoxLength = 20*cm , BoxWidth = 25*cm , BoxThickness = 10*cm;
        G4double BoxZ = 150*cm , BoxDisFromAxis = 240.*mm + BoxWidth/2. + 5*mm , BoxX , BoxY  , BoxAngle;
        const int NBoxes = 6;
        G4Box             * SolidElBox[NBoxes];
        G4LogicalVolume   * LogicElBox[NBoxes];
        G4VPhysicalVolume * PhysElBox[NBoxes];
        G4ThreeVector       posElBox[NBoxes];
        G4VisAttributes   * ElBoxVisAtt = new G4VisAttributes( true , G4Colour( 0.75 , 0.6 , 0.75 ));
        G4RotationMatrix    rot_mat;
        G4Transform3D       transform;
        
        for (int i_box = 0; i_box < NBoxes; i_box++) {
            
            BoxAngle            = i_box * 60 * deg;
            BoxX                = BoxDisFromAxis * std::cos(BoxAngle);
            BoxY                = BoxDisFromAxis * std::sin(BoxAngle);
            rot_mat = G4RotationMatrix();
            rot_mat.rotateZ( BoxAngle );
            posElBox[i_box]     = G4ThreeVector( BoxX , BoxY , BoxZ);
            transform           = G4Transform3D( rot_mat , posElBox[i_box] );

            SolidElBox[i_box]   = new G4Box("ElBox" , BoxLength/2. , BoxWidth/2. , BoxThickness/2.);
            LogicElBox[i_box]   = new G4LogicalVolume(SolidElBox[i_box], BAND_mat, "ElBoxLV" );
            LogicElBox[i_box]   -> SetVisAttributes (ElBoxVisAtt);
            PhysElBox[i_box]    = new G4PVPlacement(transform, LogicElBox[i_box],  "ElBox" , logicWorld,  false,  0, fCheckOverlaps);
        }
        
    
    }
    
   
    G4double RinDetector,  RoutDetector;
    G4double BANDThickness = 24.*cm;
    if (DoBandAsFullCircle){
        //--------- Full circle Detector for scattering angle --------------//
        RinDetector = 0.*cm;  RoutDetector = 2.*m;
    }
    if (DoBandAsProposed){
        //--------- BAND proposed Detector ---------------------------------//
        RinDetector = 612.*mm;  RoutDetector = 1292.;
    }
    G4Tubs* solidDetector = new G4Tubs("Detector", RinDetector , RoutDetector , BANDThickness/2. , 0., CLHEP::twopi);
    G4LogicalVolume* logicDetector = new G4LogicalVolume(solidDetector,BAND_mat,"DetectorLV");
    new G4PVPlacement(0,G4ThreeVector(0,0,3304*mm+BANDThickness/2.),logicDetector, "Detector",logicWorld, false,  0, fCheckOverlaps);
    G4Colour Detec_Color( 0.75 , 0.6 , 0.75 );
    G4VisAttributes* DetecVisAttributes = new G4VisAttributes( true , Detec_Color );
    logicDetector -> SetVisAttributes(DetecVisAttributes);
    
    
    if (DoLeadWall){
        //--------- Full circle Detector for scattering angle--------------//
        G4Tubs* solidLeadWall = new G4Tubs("LeadWall", RinDetector , RoutDetector , 3.*cm/2. , 0., CLHEP::twopi);
        G4LogicalVolume* logicLeadWall = new G4LogicalVolume(solidLeadWall,LeadWall_mat,"LeadWall");
        new G4PVPlacement(0,G4ThreeVector(0,0,3264*mm),logicLeadWall, "LeadWall",logicWorld, false,  0, fCheckOverlaps);
        G4Colour LeadWall_Color( 0.25 , 0.38 , 0.25 );
        G4VisAttributes* LeadWallVisAttributes = new G4VisAttributes( true , LeadWall_Color );
        logicLeadWall -> SetVisAttributes(LeadWallVisAttributes);
        
    }
    
    // Always return the physical world
    return physWorld;
}

void NeutronHallBDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors
  G4String nameSD = "BAND";
  NeutronHallBTrackerSD* bandSD = new NeutronHallBTrackerSD(treePtr,nameSD);
  SetSensitiveDetector("DetectorLV" , bandSD , true);
  
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
