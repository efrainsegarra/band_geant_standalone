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
{
  hitList = new BAND_Event;
  outTree = (TTree*) treePtr;
  outTree->Branch(name.data(),&hitList);
}

NeutronHallBTrackerSD::~NeutronHallBTrackerSD() 
{}

void NeutronHallBTrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // clear last event
  hitList->hits.clear();
  barHits.clear();
}

G4bool NeutronHallBTrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{ 

  G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable()); 
  hitBarNo = touchable->GetReplicaNumber(0);


  hitE = aStep->GetTotalEnergyDeposit()/MeV;
  if (hitE == 0.) return true;

  hitTime = aStep->GetPreStepPoint()->GetGlobalTime()/ns;
  hitPos = aStep->GetPreStepPoint()->GetPosition();

  if(barHits.count(hitBarNo) == 0){
    BAND_Hit newHit;

    newHit.E_dep = hitE;
    newHit.time = hitTime * hitE;
    newHit.pos = TVector3(hitPos.x() / mm,hitPos.y() / mm,hitPos.z() / mm) * hitE;
    newHit.barNo = hitBarNo;

    barHits[hitBarNo] = newHit;
  }
  else{
    barHits[hitBarNo].E_dep += hitE;
    barHits[hitBarNo].time += hitTime * hitE;
    barHits[hitBarNo].pos += TVector3(hitPos.x() / mm,hitPos.y() / mm,hitPos.z() / mm) * hitE;
  }
  

  return true;
}

void NeutronHallBTrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  for(std::map<int,BAND_Hit>::iterator it=barHits.begin() ; it != barHits.end() ; it++){

    eventE = it->second.E_dep;

    it->second.time *= (1./eventE);
    it->second.pos *= (1./eventE);
    hitList->hits.push_back(it->second);

  }

  
}
