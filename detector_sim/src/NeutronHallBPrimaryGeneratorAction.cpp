// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBPrimaryGeneratorAction.cc
/// \brief Implementation of the NeutronHallBPrimaryGeneratorAction class

#include "NeutronHallBPrimaryGeneratorAction.h"
#include "NeutronHallBTrackerSD.h"
#include "NeutronHallBAnalysis.h"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include "G4UnitsTable.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHallBPrimaryGeneratorAction::NeutronHallBPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction()
{
    G4int nofParticles = 1;
    fParticleGun = new G4ParticleGun(nofParticles);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHallBPrimaryGeneratorAction::~NeutronHallBPrimaryGeneratorAction()
{
    delete fParticleGun;
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHallBPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // This function is called at the begining of event
    
    // In order to avoid dependence of PrimaryGeneratorAction
    // on DetectorConstruction class we get world volume
    // from G4LogicalVolumeStore.
    
    G4double worldZHalfLength = 0;
    G4LogicalVolume* worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
    G4Box* worldBox = NULL;
    if ( worldLV ) worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
    if ( worldBox ) worldZHalfLength = worldBox->GetZHalfLength();
    else  { G4cerr << "World volume of box not found, The gun will be place in the center." << G4endl; }
    
    
    
    
    
    //-----------Particle--------------//
    G4String ParticleType = "neutron";
    //---------------------------------//
    G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(ParticleType);
    fParticleGun -> SetParticleDefinition(particleDefinition);
    
    
    
    
    //----------- Beam --------------//
    fParticleGun->SetParticlePosition(G4ThreeVector(0,0,0));
    G4double theta , phi;
    bool beam_shape = false;
    while (!beam_shape){
        theta = 180.*G4UniformRand();
        phi = 360.*G4UniformRand();
        if ( ((theta > 10.5) && (theta < 20.5)) ){
//        if ( ((theta > 0) && (theta < 30)) ){
            beam_shape = true;
        }
    }
    theta = theta*(CLHEP::twopi/360.);
    phi = phi*(CLHEP::twopi/360.);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)));
    
    if (ParticleType=="neutron")
        fParticleGun -> SetParticleMomentum( 500.*MeV * G4UniformRand() + 0.*MeV ) ;
    else if (ParticleType=="gamma")
        fParticleGun -> SetParticleEnergy( 97.*MeV * G4UniformRand() + 3.*MeV ) ;
    
    fParticleGun -> GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
