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


using namespace CLHEP;

NeutronHallBTrackerSD::NeutronHallBTrackerSD(void * treePtr, const G4String& name)
 : G4VSensitiveDetector(name)
   //fEventAction(eventAction)
{
  //fEventAction(eventAction)
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
  barHits.clear();
  vecOfbarHits.clear();

  vecOfbarHits.push_back(barHits);
  vecOfbarHits.push_back(barHits);
  vecOfbarHits.push_back(barHits);
  vecOfbarHits.push_back(barHits);
  vecOfbarHits.push_back(barHits);

  eff_energy = 0.;
  eff_time = 0.;
  eff_pos *= 0;
}

G4bool NeutronHallBTrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{ 
  G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable()); 
  hitBarNo = touchable->GetReplicaNumber(0);

  hitE = aStep->GetTotalEnergyDeposit()/MeV;
  if (hitE == 0.) return true;

  hitTime = aStep->GetPreStepPoint()->GetGlobalTime()/ns;
  hitPos = aStep->GetPreStepPoint()->GetPosition();

  if( abs( hitPos.z() / mm + 2620. + (1.) * 74/2. ) < 74/2. ) ind = 0;
  else if( abs( hitPos.z() / mm + 2620. + (3.) * 74/2. ) < 74/2. ) ind = 1;
  else if( abs( hitPos.z() / mm + 2620. + (5.) * 74/2. ) < 74/2. ) ind = 2;
  else if( abs( hitPos.z() / mm + 2620. + (7.) * 74/2. ) < 74/2. ) ind = 3;
  else if( abs( hitPos.z() / mm + 2620. + (9.) * 74/2. ) < 74/2. ) ind = 4;
  else return true;

  if(vecOfbarHits[ind].count(hitBarNo) == 0){
    BAND_Hit newHit;

    newHit.E_dep = hitE;
    newHit.time = hitTime * hitE;
    newHit.pos = TVector3(hitPos.x() / mm,hitPos.y() / mm,hitPos.z() / mm) * hitE;
    newHit.barNo = hitBarNo;

    vecOfbarHits[ind][hitBarNo] = newHit;
  }
  else{
    vecOfbarHits[ind][hitBarNo].E_dep += hitE;
    vecOfbarHits[ind][hitBarNo].time += hitTime * hitE;
    vecOfbarHits[ind][hitBarNo].pos += TVector3(hitPos.x() / mm,hitPos.y() / mm,hitPos.z() / mm) * hitE;
  }

  
  return true;
}

void NeutronHallBTrackerSD::EndOfEvent(G4HCofThisEvent*)
{

  for(int i=0; i<5;i++)
  {
    for(std::map<int,BAND_Hit>::iterator it=vecOfbarHits[i].begin() ; it != vecOfbarHits[i].end() ; it++){

      eventE = it->second.E_dep;

      it->second.time *= (1./eventE);
      it->second.pos *= (1./eventE);
      hitList->hits.push_back(it->second);

    }
  }


}
