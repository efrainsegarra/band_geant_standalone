// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBEventAction.hh
/// \brief Definition of the NeutronHallBEventAction class

#ifndef NeutronHallBEventAction_h
#define NeutronHallBEventAction_h 1

#include "G4UserEventAction.hh"
#include "TRandom2.h"

#include "globals.hh"

/// Event action class

class NeutronHallBEventAction : public G4UserEventAction
{
  public:
    TRandom2 * rand2;
    NeutronHallBEventAction();
    virtual ~NeutronHallBEventAction();

    virtual void  BeginOfEventAction(const G4Event* );
    virtual void    EndOfEventAction(const G4Event* );
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
