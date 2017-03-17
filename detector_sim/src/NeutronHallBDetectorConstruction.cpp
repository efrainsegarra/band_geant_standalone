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
    
    G4double length;
    G4double shift;

    // ------------------ Defining World ---------------- //
    G4double world_sizeX = 5000.*cm;
    G4double world_sizeY = 5000.*cm;
    G4double world_sizeZ = 5000.*cm;

  
    G4Box* solidWorld =    
        new G4Box("World",              // its name
                  0.5*world_sizeX, 
                  0.5*world_sizeY, 
                  0.5*world_sizeZ);     // its size
      
    G4LogicalVolume* logicWorld =                         
        new G4LogicalVolume(solidWorld,          // its solid
                            vacuum,           // its world material
                            "World");            // its name
                                   
    G4VPhysicalVolume* physWorld = 
        new G4PVPlacement(0,                     // no rotation
                          G4ThreeVector(),       // at (0,0,0)
                          logicWorld,            // its logical volume
                          "World",               // its name
                          0,                     // its mother  volume
                          false,                 // no boolean operation
                          0,                     // copy number
                          fCheckOverlaps);        // overlaps checking
  // ------------------ End Defining World ---------------- //

  // ------------------ Defining Target - in World ---------------- //
    G4double targetRadius = 10.*mm;

    G4Sphere* solid_target = 
        new G4Sphere("target",                    // name
                      0.,                         // inner radius 
                      targetRadius,               // outer radius
                      0.,                         // starting phi angle
                      CLHEP::twopi,               // ending phi angle
                      0.,                         // starting theta angle
                      CLHEP::twopi/2.);           // ending theta angle

    G4LogicalVolume* logic_target = 
        new G4LogicalVolume(solid_target,         // its solid
                            default_mat,           // its material
                            "target");            // its name

        new G4PVPlacement(0,                      // no rotation
                          G4ThreeVector(),          // at 0,0,0
                          logic_target,             // its logical volume
                          "target",                 // its name
                          logicWorld,               // its mother volume
                          false,                    // no boolean operation
                          0,                        // copy number
                          fCheckOverlaps);           // overlap checkign

    logic_target -> SetVisAttributes (G4Colour::G4Colour( 0.65 , 0.82 , 0.77 ));
  // ------------------ End Defining Target ---------------- //

  // ------------------ Defining BAND Detector - in World ---------------- //
    G4double RinDetector = 612.*mm;
    G4double RoutDetector = 1292.*mm;
    G4double BANDThickness = 24.*cm;
    G4double BAND_Z_coord = -(3304*mm+BANDThickness/2.);
    G4Tubs* solidDetector = 
        new G4Tubs("Detector",                    // name
                    RinDetector,                  // inner radius
                    RoutDetector,                 // outer radius
                    BANDThickness/2.,             // z extent
                    0.,                           // starting phi angle
                    CLHEP::twopi/2.);             // ending phi angle

    G4LogicalVolume* logicDetector =            
        new G4LogicalVolume(solidDetector,        // its solid
                            BAND_mat,             // material
                            "DetectorLV");        // name

        new G4PVPlacement(0,                                   // no rotation
                          G4ThreeVector(0,0,BAND_Z_coord),     // (0,0,z)
                          logicDetector,                       // logical volume
                          "Detector",                          // name
                          logicWorld,                          // mother volume
                          false,                               // no boolean operation
                          0,                                   // copy number
                          fCheckOverlaps);                      // overlap checking

    G4Colour Detec_Color( 0.75 , 0.6 , 0.75 );
    G4VisAttributes* DetecVisAttributes = new G4VisAttributes( true , Detec_Color );
    logicDetector -> SetVisAttributes(DetecVisAttributes);
  // ------------------ END Defining BAND Detector  ---------------- //

// ------------------ Defining Lead Wall  - in World ---------------- //
    G4double lead_Z_coord = -3264*mm;

    G4Tubs* solidLeadWall = 
      new G4Tubs("LeadWall",                                // name
                  RinDetector,                              // inner radius
                  RoutDetector,                             // outer raidus
                  3.*cm/2.,                                 // z extent
                  0.,                                       // starting phi angle
                  CLHEP::twopi/2.);                            // ending phi angle

    G4LogicalVolume* logicLeadWall = 
      new G4LogicalVolume(solidLeadWall,                    // its solid
                          LeadWall_mat,                     // material
                          "LeadWall");                      // its name

      new G4PVPlacement(0,                                  // no rotation
                        G4ThreeVector(0, 0,lead_Z_coord),    // its (0,0,z)
                        logicLeadWall,                      // logical volume
                        "LeadWall",                         // name
                        logicWorld,                         // mother volume
                        false,                              // no boolean operation
                        0,                                  // copy number
                        fCheckOverlaps);                     // overlap checking

    G4Colour LeadWall_Color( 0.25 , 0.38 , 0.25 );
    G4VisAttributes* LeadWallVisAttributes = new G4VisAttributes( true , LeadWall_Color );
    logicLeadWall -> SetVisAttributes(LeadWallVisAttributes);
    // ------------------ END Defining Lead Wall  - in World ---------------- //

    // ------------------ Defining Solenoid ---------------- //
    length = 480.*mm; // length of the first tube pipe

    G4Tubs* solid_solenoid =     
        new G4Tubs("Solenoid",           // name
                    476.*mm ,             // inner radius
                    1200.*mm,              // outer radius
                    length ,              // z-extent
                    0.,                   // starting phi
                    CLHEP::twopi);        // ending phi

    G4LogicalVolume* logic_solenoid =      
        new G4LogicalVolume(solid_solenoid,  // solid logic volume
                            BAND_mat,           // material
                            "Solenoid");     // name
  
        new G4PVPlacement(0,                  // rotation
                          G4ThreeVector(),    // at (0,0,0)
                          logic_solenoid,    // logical volume
                          "Solenoid",        // name
                          logicWorld,         // mother volume
                          false,              // no boolean operation
                          0,                  // copy number
                          fCheckOverlaps);     // checking overlap

    logic_solenoid -> SetVisAttributes (G4Colour::G4Colour( 0.99 , 0.88 , 0.66 ));


    const int NTOFs = 4;                          // Creating 6 boxes that will go around the beam line
    G4Tubs                *SolidTOF[NTOFs];    // solid box
    G4LogicalVolume       *LogicTOF[NTOFs];    // the logical volumes
    G4VPhysicalVolume     *PhysTOF[NTOFs];     // the physical volume
    G4ThreeVector         posTOF[NTOFs];       // the vector for the position of the boxes
  
    for (int i_tof = 0; i_tof < NTOFs; i_tof++) {
        posTOF[i_tof]     = G4ThreeVector();             // defining position vector for box

        SolidTOF[i_tof] =                                     // solid volume for each box
            new G4Tubs("TOF",                                    // name
                (252.+50.8*i_tof)*mm,             // inner radius
                (252.+50.8*(i_tof+1))*mm,              // outer radius
                length ,              // z-extent
                0.,                   // starting phi
                CLHEP::twopi);        // ending phi                         // thickness of box

        LogicTOF[i_tof] =                                     // logical volume for each box
            new G4LogicalVolume(SolidTOF[i_tof],                // box is solid
                                BAND_mat,                         // material for each box
                                "TOF" );                      // name for each logical volume box

        LogicTOF[i_tof]->SetVisAttributes (G4Colour::G4Colour( 0.75+(i_tof/10.) , 0.6-(i_tof/10.) , 0.75 ));      // setting color
      
        PhysTOF[i_tof]=                                       // physical volume for box
            new G4PVPlacement(0,   
                              G4ThreeVector(),                       // rotation for each box
                              LogicTOF[i_tof],                  // logical volume of each box
                              "TOF",                            // name of each box
                              logicWorld,                         // mother volume for box
                              false,                              // no boolean operation
                              0,                                  // copy number
                              fCheckOverlaps);                     // checking overlap
    }

    length = 1200/2.*mm;
    shift = 480*mm;
      // 30 degrees for 3 CND and 25 degrees for 1 TOF
    const int slanted_TOFs = 4;                          // Creating 6 boxes that will go around the beam line

    G4Cons                *Solid_slantedTOF_front[slanted_TOFs];    // solid box
    G4LogicalVolume       *Logic_slantedTOF_front[slanted_TOFs];    // the logical volumes
    G4VPhysicalVolume     *Phys_slantedTOF_front[slanted_TOFs];     // the physical volume
    G4ThreeVector         pos_slantedTOF_front[slanted_TOFs];       // the vector for the position of the boxes

    for (int i_tof = 0; i_tof < slanted_TOFs; i_tof++) {
        pos_slantedTOF_front[i_tof]     = G4ThreeVector();             // defining position vector for box

        if(i_tof==0){
            Solid_slantedTOF_front[i_tof] =                                     // solid volume for each box
                new G4Cons("slanted_TOF_front",                                    // name
                           (560+252.)*mm,
                           (560+252.+50.8)*mm,
                           (252.)*mm,
                           (252.+50.8)*mm,
                           length ,              // z-extent
                           0.,                   // starting phi
                           CLHEP::twopi);        // ending phi                         // thickness of box
        }
        else{
            Solid_slantedTOF_front[i_tof] =                                     // solid volume for each box
                new G4Cons("slanted_TOF_front",                                    // name
                           (695+252.+50.8*(i_tof))*mm,
                           (695+252.+50.8*(i_tof+1))*mm,
                           (252.+50.8*(i_tof))*mm,
                           (252.+50.8*(i_tof+1))*mm,
                           length ,              // z-extent
                           0.,                   // starting phi
                           CLHEP::twopi);        // ending phi                         // thickness of box
        }
        Logic_slantedTOF_front[i_tof] =                                     // logical volume for each box
            new G4LogicalVolume(Solid_slantedTOF_front[i_tof],                // box is solid
                                BAND_mat,                         // material for each box
                                "slanted_TOF_front" );                      // name for each logical volume box

        Logic_slantedTOF_front[i_tof]->SetVisAttributes(G4Colour::G4Colour( 0.75 , 0.6+(i_tof/10.) , 0.75-(i_tof/10.) ));      // setting color
          
        Phys_slantedTOF_front[i_tof]=                                       // physical volume for box
            new G4PVPlacement(0,   
                              G4ThreeVector(0,0,-1.*(shift+length)),                       // rotation for each box
                              Logic_slantedTOF_front[i_tof],                  // logical volume of each box
                              "slanted_TOF_front",                            // name of each box
                              logicWorld,                         // mother volume for box
                              false,                              // no boolean operation
                              0,                                  // copy number
                              fCheckOverlaps);                     // checking overlap
    }
    // ------------------ END Defining Solenoid  - in World ---------------- //

    // ------------------ Defining SS Tubing around beam pipe ---------------- //
    length = 658.*mm; // length of the first tube pipe

    G4Tubs* solid_SSTube_P1 =     
        new G4Tubs("SSTube_P1",           // name
                    210.*mm ,             // inner radius
                    213.*mm,              // outer radius
                    length ,              // z-extent
                    0.,                   // starting phi
                    CLHEP::twopi);        // ending phi

    G4LogicalVolume* logic_SSTube_P1 =      
        new G4LogicalVolume(solid_SSTube_P1,  // solid logic volume
                            SS_mat,           // material
                            "SSTube_P1");     // name
      
        new G4PVPlacement(0,                  // rotation
                          G4ThreeVector(),    // at (0,0,0)
                          logic_SSTube_P1,    // logical volume
                          "SSTube_P1",        // name
                          logicWorld,         // mother volume
                          false,              // no boolean operation
                          0,                  // copy number
                          fCheckOverlaps);     // checking overlap

    logic_SSTube_P1 -> SetVisAttributes (G4Colour::G4Colour( 0.99 , 0.88 , 0.66 ));

        //  ----------1 cm plastic coating around tub to represent cables ---------- //
    G4Tubs* solid_SSTube_P1_coating = 
        new G4Tubs("SSTube_P1_coating",                   // name
                    214.*mm ,                             // inner radius -- starting at outer radius of 1st pipe
                    224.*mm,                              // outer radius
                    length ,                              // z extent
                    0.,                                   // starting phi
                    CLHEP::twopi);                        // ending phi

    G4LogicalVolume* logic_SSTube_P1_coating =      
        new G4LogicalVolume(solid_SSTube_P1_coating,      // solid logic volume
                            BAND_mat ,                    // material
                            "SSTube_P1_coating");         // name

        new G4PVPlacement(0,                              // no rotation
                          G4ThreeVector(),                // at (0,0,0)
                          logic_SSTube_P1_coating,        // logical volume
                          "SSTube_P1_coating",            // name
                          logicWorld,                     // mother volume
                          false,                          // no boolean operation
                          0,                              // copy number
                          fCheckOverlaps);                 // checking overlaps

    logic_SSTube_P1_coating -> SetVisAttributes (G4Colour::G4Colour( 0.75 , 0.6 , 0.75 ));
  
        //  ---------- SS connector from target tube to beam tube ---------- //
    length = 24.*mm;
    shift = 658.*mm;

    G4Tubs* solid_SSTube_Connector1 =         
        new G4Tubs("SSTube_Connector1",     // name
                    213.*mm ,               // inner radius
                    240.*mm,                // outer radius
                    length ,                // z extent
                    0.,                     // starting phi
                    CLHEP::twopi);          // ending phi

    G4LogicalVolume* logic_SSTube_Connector1 = 
        new G4LogicalVolume(solid_SSTube_Connector1,    // solid object
                            SS_mat ,                    // material
                            "SSTube_Connector1");       // name

        new G4PVPlacement(0,                                    // no rotation
                          G4ThreeVector(0,0,-(length+shift)),      // at (x,y,z)
                          logic_SSTube_Connector1,              // logical volume
                          "SSTube_Connector1",                  // name
                          logicWorld,                           // mother volume
                          false,                                // boolean operation  
                          0,                                    // copy number
                          fCheckOverlaps);                       // checking overlaps
    
    logic_SSTube_Connector1 -> SetVisAttributes (G4Colour::G4Colour( 0.99 , 0.88 , 0.66 ));
  
      //  ---------- SS connector coating to represent cables ---------- //
    G4Tubs* solid_SSTube_Connector1_coating = 
        new G4Tubs("SSTube_Connector1_coating", 
                   241.*mm, 
                   251.*mm, 
                   length, 
                   0.,
                   CLHEP::twopi);

    G4LogicalVolume* logic_SSTube_Connector1_coating = 
        new G4LogicalVolume(solid_SSTube_Connector1_coating, 
                            BAND_mat , 
                            "SSTube_Connector1_coating");

        new G4PVPlacement(0, 
                          G4ThreeVector(0,0,-(length+shift)),
                          logic_SSTube_Connector1_coating,
                          "SSTube_Connector1_coating",
                          logicWorld,
                          false,
                          0,
                          fCheckOverlaps);

    logic_SSTube_Connector1_coating -> SetVisAttributes (G4Colour::G4Colour( 0.75 , 0.6 , 0.75 ));
     
  
      //  ---------- SS connector #2 from target tube to beam tube ---------- //
    length = 34.*mm;
    shift = 658.*mm;

    G4Tubs* solid_SSTube_Connector2 = 
        new G4Tubs("SSTube_Connector2",
                   240.*mm,
                   258*mm,
                   length,
                   0.,
                   CLHEP::twopi);

    G4LogicalVolume* logic_SSTube_Connector2 = 
        new G4LogicalVolume(solid_SSTube_Connector2, 
                            SS_mat, 
                            "SSTube_Connector2");

        new G4PVPlacement(0, 
                          G4ThreeVector(0,0,-(length+shift)),
                          logic_SSTube_Connector2, 
                          "SSTube_Connector2",
                          logicWorld, 
                          false,  
                          0,  
                          fCheckOverlaps);

    logic_SSTube_Connector2 -> SetVisAttributes (G4Colour::G4Colour( 0.99 , 0.88 , 0.66 ));
     
      //  ---------- SS connector #2 coating to represent cables ---------- //
    G4Tubs* solid_SSTube_Connector2_coating = 
        new G4Tubs("SSTube_Connector2_coating", 
                   259.*mm, 
                   269*mm, 
                   length, 
                   0., 
                   CLHEP::twopi);

    G4LogicalVolume* logic_SSTube_Connector2_coating = 
        new G4LogicalVolume(solid_SSTube_Connector2_coating, 
                            BAND_mat, 
                            "SSTube_Connector2_coating");

        new G4PVPlacement(0, 
                          G4ThreeVector(0,0,-(length+shift)),
                          logic_SSTube_Connector2_coating, 
                          "SSTube_Connector2_coating", 
                          logicWorld, 
                          false,  
                          0,  
                          fCheckOverlaps);

    logic_SSTube_Connector2_coating -> SetVisAttributes (G4Colour::G4Colour( 0.75 , 0.6 , 0.75 ));
              
      //  ---------- SS tube from connector #2 to BAND detector ---------- //
    length = 3304*mm+BANDThickness/2.*mm-1750*mm;
    shift = 682.*mm;

    G4Tubs* solid_SSTube_P2 = 
        new G4Tubs("SSTube_P2", 
                   237.*mm, 
                   240*mm, 
                   length, 
                   0., 
                   CLHEP::twopi);

    G4LogicalVolume* logic_SSTube_P2 = 
        new G4LogicalVolume(solid_SSTube_P2, 
                            SS_mat, 
                            "SSTube_P2");

        new G4PVPlacement(0, 
                          G4ThreeVector(0,0,-(length+shift)),
                          logic_SSTube_P2, 
                          "SSTube_P2", 
                          logicWorld, 
                          false,  
                          0,  
                          fCheckOverlaps);

    logic_SSTube_P2 -> SetVisAttributes (G4Colour::G4Colour( 0.99 , 0.88 , 0.66 ));
  
      //  ---------- SS tube coatings from connector #2 to BAND detector ---------- //
      
    G4Tubs* solid_SSTube_P2_coating = 
        new G4Tubs("SSTube_P2_coating", 
                   241.*mm, 
                   251*mm, 
                   length, 
                   0., 
                   CLHEP::twopi);

    G4LogicalVolume* logic_SSTube_P2_coating = 
        new G4LogicalVolume(solid_SSTube_P2_coating, 
                            BAND_mat, 
                            "SSTube_P2_coating");

        new G4PVPlacement(0, 
                          G4ThreeVector(0,0,-(length+shift)),  
                          logic_SSTube_P2_coating, 
                          "SSTube_P2", 
                          logicWorld, 
                          false,  
                          0,  
                          fCheckOverlaps);

    logic_SSTube_P2_coating -> SetVisAttributes (G4Colour::G4Colour( 0.75 , 0.6 , 0.75 ));
  // ------------------ END SS Tubing around beam pipe ---------------- //


// ------------------ Defining Electronic Boxes  - in World ---------------- //
    G4double BoxLength = 20*cm;
    G4double BoxWidth = 25*cm;
    G4double BoxThickness = 10*cm;
    G4double BoxZ = -150*cm;
    G4double BoxDisFromAxis = 240.*mm + BoxWidth/2. + 5*mm;

    G4double BoxX;
    G4double BoxY;
    G4double BoxAngle;

    const int NBoxes = 6;                          // Creating 6 boxes that will go around the beam line
    G4Box                 *SolidElBox[NBoxes];    // solid box
    G4LogicalVolume       *LogicElBox[NBoxes];    // the logical volumes
    G4VPhysicalVolume     *PhysElBox[NBoxes];     // the physical volume
    G4ThreeVector         posElBox[NBoxes];       // the vector for the position of the boxes
    G4VisAttributes       *ElBoxVisAtt = new G4VisAttributes( true , G4Colour::G4Colour( 0.75 , 0.6 , 0.75 ));   // setting color
    G4RotationMatrix      rot_mat;                // rotation matrix to rotate each box
    G4Transform3D         transform;              // the rotation that will be applied to the placement
    
    for (int i_box = 0; i_box < NBoxes; i_box++) {
        BoxAngle            = i_box * 60 * deg;                               // this box's angle around beam line
        BoxX                = BoxDisFromAxis * std::cos(BoxAngle);            // the x-placement based on angle
        BoxY                = BoxDisFromAxis * std::sin(BoxAngle);            // the y-placement based on angle
        rot_mat             = G4RotationMatrix();                             // initializing rotation matrix
        rot_mat.rotateZ(BoxAngle);                                            // defining rotation matrix around z-beamLine for box
        posElBox[i_box]     = G4ThreeVector( BoxX , BoxY , BoxZ);             // defining position vector for box
        transform           = G4Transform3D( rot_mat , posElBox[i_box] );     // creating rotation of each box 

        SolidElBox[i_box] =                                     // solid volume for each box
            new G4Box("ElBox",                                    // name
                     BoxLength/2.,                              // length of box
                     BoxWidth/2.,                               // width of box
                     BoxThickness/2.);                          // thickness of box

        LogicElBox[i_box] =                                     // logical volume for each box
            new G4LogicalVolume(SolidElBox[i_box],                // box is solid
                              BAND_mat,                         // material for each box
                              "ElBoxLV" );                      // name for each logical volume box

        LogicElBox[i_box]->SetVisAttributes (ElBoxVisAtt);      // setting color
        
        PhysElBox[i_box]=                                       // physical volume for box
            new G4PVPlacement(transform,                          // rotation for each box
                            LogicElBox[i_box],                  // logical volume of each box
                            "ElBox",                            // name of each box
                            logicWorld,                         // mother volume for box
                            false,                              // no boolean operation
                            0,                                  // copy number
                            fCheckOverlaps);                     // checking overlap
    }
  // ------------------ END Defining Electronic Boxes  ---------------- //
    
    
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