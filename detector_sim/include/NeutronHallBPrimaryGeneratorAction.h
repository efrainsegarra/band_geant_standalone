// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBPrimaryGeneratorAction.hh
/// \brief Definition of the NeutronHallBPrimaryGeneratorAction class

#ifndef NeutronHallBPrimaryGeneratorAction_h
#define NeutronHallBPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the Tracker 
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class 
/// (see the macros provided with this example).

class NeutronHallBPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    NeutronHallBPrimaryGeneratorAction();    
    virtual ~NeutronHallBPrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event* );

    G4ParticleGun* GetParticleGun() {return fParticleGun;}
  
    // Set methods
    void SetRandomFlag(G4bool );
    
    
    G4double E_kin_g;
    G4int particle_ID_g;
    G4double p_g;
    G4int Event_Number_g;

  private:
    G4ParticleGun*          fParticleGun; // G4 particle gun
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
