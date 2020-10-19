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
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
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

	G4String hitParticleName = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
	G4double hitParticleCharge = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetPDGCharge();
	if( hitE > 0 ){
		G4StepPoint * preStep = (G4StepPoint*) aStep->GetPreStepPoint();
		G4StepPoint * postStep = (G4StepPoint*) aStep->GetPostStepPoint();
		double Tin = preStep->GetKineticEnergy()/MeV;
		double Tout = postStep->GetKineticEnergy()/MeV;
		//	note that hitE is the same as Tin-Tout but we do the conversion on Tin and Tout separately

		double Tin_e, Tout_e, dT_e;
		Tin_e 	= MeVtoMeVee( Tin , hitParticleName , hitParticleCharge );
		Tout_e 	= MeVtoMeVee( Tout , hitParticleName , hitParticleCharge );
		dT_e = Tin_e - Tout_e;
		if( dT_e < 0 ){ return true; }


		hitE = dT_e;
		//G4Track * track = (G4Track*) aStep->GetTrack();
		//const G4VProcess * creatorProcess = track->GetCreatorProcess();
		//creatorProcess->DumpInfo();
	}	
	else{
		return true;
	}

	hitTime = aStep->GetPreStepPoint()->GetGlobalTime()/ns;
	hitPos = aStep->GetPreStepPoint()->GetPosition();

	// For each event with multiple hits, if there are multiple hits
	// IN THE SAME BAR, take only the earliest time / position in the bar.
	// But, for each event, with multiple bars fired, save all bars fired.
	if(barHits.count(hitBarNo) == 0){
		BAND_Hit newHit;

		newHit.E_dep = hitE;
		newHit.time = hitTime; //* hitE;
		newHit.pos = TVector3(hitPos.x() / mm,hitPos.y() / mm,hitPos.z() / mm); //* hitE;
		newHit.barNo = hitBarNo;

		barHits[hitBarNo] = newHit;
	}
	else{
		//return true;
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
		hitList->hits.push_back(it->second);

	}

}

double NeutronHallBTrackerSD::MeVtoMeVee( double E_MeV, G4String name , G4double Z){
	double E_MeVee = 0;
	double a1, a2, a3, a4;
	if(name == "e-" || name == "e+" || name == "mu-" || name == "mu+" || name == "gamma" ){
		a1 = 1;
		a2 = 0;
		a3 = 0;
		a4 = 0;
	}
	else if(name == "proton"){
		a1 = 0.902713 ;
		a2 = 7.55009  ;
		a3 = 0.0990013;
		a4 = 0.736281 ;
	}
	else if(name == "deuteron"){
		a1 = 0.891575 ;
		a2 = 12.2122  ;
		a3 = 0.0702262;
		a4 = 0.782977 ;
	}
	else if(name == "triton"){
		a1 = 0.881489 ;
		a2 = 15.9064  ;
		a3 = 0.0564987;
		a4 = 0.811916 ;
	}
	else if(name == "He3"){
		a1 = 0.803919 ;
		a2 = 34.4153  ;
		a3 = 0.0254322;
		a4 = 0.894859 ;
	}
	else if(name == "alpha" || Z == 2){
		a1 = 0.781501 ;
		a2 = 39.3133  ;
		a3 = 0.0217115;
		a4 = 0.910333 ;
	}
	else if( Z == 3 ){
		a1 = 0.613491 ;
		a2 = 57.1372  ;
		a3 = 0.0115948;
		a4 = 0.951875 ;
	}
	else if( Z == 4 ){
		a1 = 0.435772 ;
		a2 = 45.538   ;
		a3 = 0.0104221;
		a4 = 0.916373 ;
	}
	else if( Z == 6 ){
		a1 = 0.298394 ;
		a2 = 25.5679  ;
		a3 = 0.0130345;
		a4 = 0.908512 ;
	}
	else if( Z == 5 ){
		a1 = 0.350273 ;
		a2 = 34.4664  ;
		a3 = 0.0112395;
		a4 = 0.912711 ;
	}
	else{ // Treat as C12 and call it good
		a1 = 0.298394 ;
		a2 = 25.5679  ;
		a3 = 0.0130345;
		a4 = 0.908512 ;
	}
	E_MeVee = a1 * E_MeV  -  a2 * (1-exp(-a3*pow(E_MeV, a4)));

	return E_MeVee;
}
