// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBTrackerHit.cc
/// \brief Implementation of the NeutronHallBTrackerHit class

#include "NeutronHallBTrackerHit.h"
#include "NeutronHallBAnalysis.h"

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<NeutronHallBTrackerHit>* NeutronHallBTrackerHitAllocator=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHallBTrackerHit::NeutronHallBTrackerHit()
: G4VHit(),
fTrackID(-10000),
fEdep(-10000.),
fPos(G4ThreeVector(-10000.,-10000.,-10000.)),
fTotalEnergy(-10000.),
fKineticEnergy(-10000.),
fMomentum(-10000.),
fParticleCode(-10000),
fTheta(-10000.),
fX(-10000.),
fY(-10000.),
fZ(-10000.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHallBTrackerHit::~NeutronHallBTrackerHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHallBTrackerHit::NeutronHallBTrackerHit(const NeutronHallBTrackerHit& right)
: G4VHit()
{
    fTrackID        = right.fTrackID;
    fEdep           = right.fEdep;
    fPos            = right.fPos;
    fTotalEnergy    = right.fTotalEnergy;
    fKineticEnergy  = right.fKineticEnergy;
    fMomentum       = right.fMomentum;
    fParticleCode   = right.fParticleCode;
    fParticleName   = right.fParticleName;
    fTheta          = right.fTheta;
    fX              = right.fPos.x();
    fY              = right.fPos.y();
    fZ              = right.fPos.z();
    fHitTime        = right.fHitTime;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const NeutronHallBTrackerHit& NeutronHallBTrackerHit::operator=(const NeutronHallBTrackerHit& right)
{
    fTrackID        = right.fTrackID;
    fEdep           = right.fEdep;
    fPos            = right.fPos;
    fTotalEnergy    = right.fTotalEnergy;
    fKineticEnergy  = right.fKineticEnergy;
    fMomentum       = right.fMomentum;
    fParticleCode   = right.fParticleCode;
    fParticleName   = right.fParticleName;
    fTheta          = right.fTheta;
    fX              = right.fPos.x();
    fY              = right.fPos.y();
    fZ              = right.fPos.z();
    fHitTime        = right.fHitTime;
    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int NeutronHallBTrackerHit::operator==(const NeutronHallBTrackerHit& right) const
{
    return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHallBTrackerHit::Draw()
{
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager)
    {
        G4Circle circle(fPos);
        circle.SetScreenSize(4.);
        circle.SetFillStyle(G4Circle::filled);
        G4Colour colour(1.,0.,0.);
        G4VisAttributes attribs(colour);
        circle.SetVisAttributes(attribs);
        pVVisManager->Draw(circle);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHallBTrackerHit::Print()
{
    std::cout << "\n---------------- Hit - start --------------------"
    << "\n Particle Code: " << fParticleCode
    << "\n Particle Name: " << fParticleName
    << "\n Momentum: "  << G4BestUnit(fMomentum , "Energy" )
    << "\n Kinetic Energy: " << std::setw(7) << G4BestUnit( fKineticEnergy , "Energy" )
    << "\n Theta =  " << fTheta
    << "\n X =  " << fX
    << "\n Y =  " << fY
    << "\n Z =  " << fZ
    << "\n hit time =  " << fHitTime
    << "\n-------------------- Hit - End ---------------------\n"
    << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
