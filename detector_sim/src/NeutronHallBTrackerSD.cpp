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

}

G4bool NeutronHallBTrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{ 
  G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable()); 
  hitBarNo = touchable->GetReplicaNumber(0);

  hitE = aStep->GetTotalEnergyDeposit()/MeV;

  bool testGeometryEffect = true; // looks only at neutrons hitting the face
  // of the detector to see the effect geometry had on getting them there
  if (testGeometryEffect){
    G4String hitParticleName = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
    G4String wantedParticle = "neutron";
    if (hitParticleName != wantedParticle) return true;
    if( std::abs(aStep->GetPreStepPoint()->GetPosition().z()/mm +2620.) > 0.001) return true;
    if (hitE != 0.) return true;

    // I'm not sure how to get the information he wants with this chain -- need to talk 
    // to Axel. We want to track for a given event, what is the information of that event
    // right when the neutron gets to (if if does...) the FACE of the detector, and we want
    // to skim those events based on if that event then caused A SINGLE hit in A SINGLE
    // bar above threshold. I'm not sure what the right way to track this is because I can't
    // use parentID number, because neutrons are destroyed/recreated at interactions, and we
    // want to see the effect of the geometry that lies BEFORE the detector...

  }
  else{
    if (hitE == 0.) return true;
  }
  
  hitTime = aStep->GetPreStepPoint()->GetGlobalTime()/ns;
  hitPos = aStep->GetPreStepPoint()->GetPosition();

  if(barHits.count(hitBarNo) == 0){
    BAND_Hit newHit;

    newHit.E_dep = hitE;
    newHit.time = hitTime; //* hitE;
    newHit.pos = TVector3(hitPos.x() / mm,hitPos.y() / mm,hitPos.z() / mm); //* hitE;
    newHit.barNo = hitBarNo;

    barHits[hitBarNo] = newHit;
  }
  else{
    barHits[hitBarNo].E_dep += hitE;
    if (hitTime < barHits[hitBarNo].time){
      barHits[hitBarNo].time = hitTime;
      barHits[hitBarNo].pos = TVector3(hitPos.x() / mm,hitPos.y() / mm,hitPos.z() / mm); //* hitE;
    }
  }

  return true;
}

void NeutronHallBTrackerSD::EndOfEvent(G4HCofThisEvent*)
{

  for(std::map<int,BAND_Hit>::iterator it=barHits.begin() ; it != barHits.end() ; it++){
    //eventE = it->second.E_dep;
    //it->second.time *= (1./eventE);
    //it->second.pos *= (1./eventE);
    hitList->hits.push_back(it->second);

  }


}
