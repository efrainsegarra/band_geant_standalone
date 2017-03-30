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

#include "TTree.h"

#include "prop_tree.h"

using namespace CLHEP;

NeutronHallBTrackerSD::NeutronHallBTrackerSD(void * treePtr, const G4String& name) : G4VSensitiveDetector(name)
{
  hitList = new BAND_Event;
  outTree = (TTree*) treePtr;
  outTree->Branch(name.data(),&hitList);
}

NeutronHallBTrackerSD::~NeutronHallBTrackerSD() 
{}

void NeutronHallBTrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Clear the last event
  hitList->hits.clear();
}

G4bool NeutronHallBTrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{ 
  // only enters this loop 709 times but saves 133k events in the root file...?




  // create a new band hit and add to the vector!
  if (aStep->GetTotalEnergyDeposit() < 1*eV){
    return false;
  }
  
  // asks only if the particle in BAND is a neutron
  if (aStep->GetTrack()->GetDefinition()->GetPDGEncoding() != 2112){ // particleID
    return false;
  }
  BAND_Hit newHit;

  // Hit position
        // MUST USE PRESTEPPOINT FOR GEOMETRY
  G4ThreeVector worldPos = aStep->GetPreStepPoint()->GetPosition();
  newHit.pos = TVector3(worldPos.x() / mm,worldPos.y() / mm,worldPos.z() / mm);

  // Time
  newHit.time = aStep->GetPreStepPoint()->GetGlobalTime() / ns;

  // Energy deposition
  newHit.E_dep = aStep->GetTotalEnergyDeposit()/MeV;
  
  // Track ID
  newHit.track = aStep->GetTrack()->GetTrackID();
  
  G4cout << aStep->GetTrack()->GetCurrentStepNumber() << G4endl;
  // push the hit into the vector
  hitList->hits.push_back(newHit);

  //G4cout << aStep->GetTrack()->GetTrackID() << G4endl;
  //G4cout << newHit.track << G4endl;
  //G4cout << newHit << G4endl;


  /*
  // Stuff that Erez was doing that is not useful.
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
  */

  return true;
}

void NeutronHallBTrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 )
    { 
      G4int nofHits = hitList->hits.size();
      G4cout << "\n-------->Hits Collection: in this event they were " << nofHits << G4endl;
    }
}
