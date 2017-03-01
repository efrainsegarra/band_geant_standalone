// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBActionInitialization.cc
/// \brief Implementation of the NeutronHallBActionInitialization class

#include "NeutronHallBActionInitialization.h"
#include "NeutronHallBPrimaryGeneratorAction.h"
#include "NeutronHallBRunAction.h"
#include "NeutronHallBEventAction.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHallBActionInitialization::NeutronHallBActionInitialization()
 : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NeutronHallBActionInitialization::~NeutronHallBActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHallBActionInitialization::BuildForMaster() const
{
  SetUserAction(new NeutronHallBRunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NeutronHallBActionInitialization::Build() const
{
  SetUserAction(new NeutronHallBPrimaryGeneratorAction);
  SetUserAction(new NeutronHallBRunAction);
  SetUserAction(new NeutronHallBEventAction);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
