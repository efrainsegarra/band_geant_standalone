// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBTrackerSD.cc
/// \brief Implementation of the NeutronHallBTrackerSD class

#include "NeutronHallBTrackerSD.h"
#include "NeutronHallBAnalysis.h"
#include "NeutronHallBPrimaryGeneratorAction.h"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHallBTrackerSD::NeutronHallBTrackerSD(const G4String& name, const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHallBTrackerSD::~NeutronHallBTrackerSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHallBTrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection = new NeutronHallBTrackerHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool NeutronHallBTrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{  
    NeutronHallBTrackerHit* newHit = new NeutronHallBTrackerHit();
    G4int Track_ID          = aStep ->  GetTrack() ->  GetTrackID();
    G4String particle_Name  = aStep ->  GetTrack() ->  GetDefinition() ->  GetParticleName();
   
    
    G4int particle_ID       = aStep ->  GetTrack() ->  GetDefinition() ->  GetPDGEncoding();
    G4double E_kin          = aStep ->  GetPreStepPoint()  ->  GetKineticEnergy();//GetTrack() -> GetKineticEnergy()
    G4double E_tot          = aStep ->  GetTrack() -> GetTotalEnergy();//GetPostStepPoint()  ->  GetTotalEnergy();
    G4double p              = (aStep -> GetTrack() -> GetMomentum()).mag();
    G4double Edep           = aStep ->  GetTotalEnergyDeposit() ;
    G4ThreeVector Position  = aStep ->  GetPreStepPoint()  ->  GetPosition();
    G4double Theta          = (360./CLHEP::twopi)*Position.getTheta();
    G4double HitTime        = aStep -> GetPreStepPoint() -> GetGlobalTime();
    
    
    newHit -> SetParticleCode       (particle_ID);
    newHit -> SetParticleName       (particle_Name);
    newHit -> SetKineticEnergy      (E_kin);
    newHit -> SetMomentum           (p);
    newHit -> SetTotalEnergy        (E_tot);
    newHit -> SetTrackID            (Track_ID);
    newHit -> SetEdep               (Edep);
    newHit -> SetPos                (Position);
    newHit -> SetTheta              (Theta);
    newHit -> SetHitTime            (HitTime);
    fHitsCollection -> insert( newHit );
    
    
    
//    G4double dose = aStep -> GetTotalEnergyDeposit()/MassTarget;
//    Run->AddDose(dose);
    
    
    if (Track_ID != 1)
        return false;
    else
        return true;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHallBTrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     G4int nofHits = fHitsCollection->entries();
     G4cout << "\n-------->Hits Collection: in this event they are " << nofHits
            << " hits in the tracker chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
