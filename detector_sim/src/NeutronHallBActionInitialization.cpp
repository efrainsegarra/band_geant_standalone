// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBActionInitialization.cc
/// \brief Implementation of the NeutronHallBActionInitialization class

#include "NeutronHallBActionInitialization.h"
#include "NeutronHallBPrimaryGeneratorAction.h"
#include "NeutronHallBRunAction.h"
#include "NeutronHallBEventAction.h"


NeutronHallBActionInitialization::NeutronHallBActionInitialization(void * tp)
 : G4VUserActionInitialization()
{
  treePtr = tp;
}

NeutronHallBActionInitialization::~NeutronHallBActionInitialization()
{}

void NeutronHallBActionInitialization::BuildForMaster() const
{
  SetUserAction(new NeutronHallBRunAction);
}

void NeutronHallBActionInitialization::Build() const
{
  SetUserAction(new NeutronHallBPrimaryGeneratorAction(treePtr));
  SetUserAction(new NeutronHallBRunAction);
  SetUserAction(new NeutronHallBEventAction);
}  
