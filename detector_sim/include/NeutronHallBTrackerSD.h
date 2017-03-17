// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBTrackerSD.hh
/// \brief Definition of the NeutronHallBTrackerSD class

#ifndef NeutronHallBTrackerSD_h
#define NeutronHallBTrackerSD_h 1

#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;
class BAND_Event;
class TTree;

/// NeutronHallBTracker sensitive detector class
///
/// The hits are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step. A hit is created for each step and stored
/// in the output branch.

class NeutronHallBTrackerSD : public G4VSensitiveDetector
{
 public:
  NeutronHallBTrackerSD(void * treePtr, const G4String& name);
  virtual ~NeutronHallBTrackerSD();
  
  // methods from base class
  virtual void   Initialize(G4HCofThisEvent* hitCollection);
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
  virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

 private:
  BAND_Event * hitList;
  TTree * outTree;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
