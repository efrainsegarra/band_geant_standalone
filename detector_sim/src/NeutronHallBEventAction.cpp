// ********************************************************************
// NeutronHallB simulation
// ********************************************************************
/// \file NeutronHallBEventAction.cc
/// \brief Implementation of the NeutronHallBEventAction class

#include "NeutronHallBEventAction.h"
#include "NeutronHallBRunAction.h"
#include "NeutronHallBAnalysis.h"
#include "NeutronHallBTrackerHit.h"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"


#include <G4VHit.hh>
#include <G4VHitsCollection.hh>
#include "Randomize.hh"

#include "G4UnitsTable.hh"
#include <iomanip>
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"

#include "TTree.h"

NeutronHallBEventAction::NeutronHallBEventAction(void * treePtr)
  : G4UserEventAction()
{
  rand2  = new TRandom2();
  outTree = (TTree*) treePtr;
}

NeutronHallBEventAction::~NeutronHallBEventAction()
{}

void NeutronHallBEventAction::BeginOfEventAction(const G4Event*){
//    G4cout << "starting event...." << G4endl;
}

void NeutronHallBEventAction::EndOfEventAction(const G4Event* event)
{
  outTree->Fill();
  /*
  // get number of stored trajectories
  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
    
  // periodic printing
  G4int eventID = event -> GetEventID();
  //bool PrintData = (eventID%100 == 0)? true : false;
  PrintData = false;
  
  //generated event
  G4double rInnerSection      = 3304*tan(17.3*(CLHEP::twopi/360.));   //10.5 + 6.8 = 17.3 is the angle of the inner section

  G4PrimaryParticle * primary_particle = event -> GetPrimaryVertex() -> GetPrimary();
  G4int       particle_ID_g       = primary_particle -> GetPDGcode();
  G4double    mass                = primary_particle -> GetMass ();
  G4double    p_g                 = (primary_particle -> GetMomentum()).mag();
  G4double    E_kin_g             = mass*( sqrt(p_g*p_g + mass*mass)/mass - 1 );
  G4double    x_g                 = event -> GetPrimaryVertex() -> GetX0();
  G4double    y_g                 = event -> GetPrimaryVertex() -> GetY0();
  G4double    theta_g             = ((primary_particle -> GetMomentum()).theta())*(360./CLHEP::twopi);
  G4int       Sector_g            = ( (theta_g > 10.5) && (theta_g < 17.3) ) ? 1 : 2 ;  //10.5+6.8=17.3 is the angle of the inner section
  
    
    //event after the wall
    G4VHitsCollection* hc           = event->GetHCofThisEvent()->GetHC(0);
    size_t      HitsSize            = hc -> GetSize();
    G4int       NeutronReachedDet   = 0;
    G4double    particle_ID         = -10000;
    G4double    E_kin               = -10000;
    G4double    p                   = -10000;
    G4double    Theta               = -10000;
    G4double    x                   = -10000;
    G4double    y                   = -10000;
    G4double    z                   = -10000;
    G4int       Sector              = -10000;
    G4int       TrackID             = -10000;
    G4String    PrtclName           = "";
    G4double    MomentumLoss        = -10000;
    G4double    TotEnergyDep        = 0;
    G4double    HitTime             = 0;
    G4double    l_reconstructed     = 0;    // reconstructed path lenght
    G4double    v_reconstructed     = 0;    // reconstructed velocity from rec. path length + TOF
    G4double    gamma_reconstrc     = 0;
    G4double    p_reconstructed     = 0;    // corresponding momentum
    G4double    delta_p             = 0;    // (p(rec)-p(gen) / p(gen)) [%]
    bool        ParticleHitDetector = false;
    

    
    
    if (HitsSize > 0){
        
        for (size_t j = 0 ; j < HitsSize ; j++ ){
            
            NeutronHallBTrackerHit * track =  (NeutronHallBTrackerHit*)(event->GetHCofThisEvent()->GetHC(0)->GetHit(j)) ;
            
            TotEnergyDep += track -> GetEdep();
            if (PrintData){
                G4cout << track -> GetParticleName() << "(" << track -> GetParticleCode() << ") Hit the Detector at Hit " << j
                << ",  TotEnergyDep = "  << std::setw(7) << G4BestUnit( TotEnergyDep , "Energy" )
                << ",  KineticEnergy : " << std::setw(7) << G4BestUnit( track -> GetKineticEnergy() , "Energy" )
                << G4endl;
            }

            
            if ( (track -> GetParticleCode()) == particle_ID_g){
                if (!ParticleHitDetector){
                    // take the hit iformation from the first primary particle hit
                    x               = track -> GetX();
                    y               = track -> GetY();
                    Sector          = (sqrt(x*x + y*y) < rInnerSection) ? 1 : 2 ;
                    z               = track -> GetZ();
                    Theta           = track -> GetTheta();
                    TrackID         = track -> GetTrackID();
                    PrtclName       = track -> GetParticleName();
                    particle_ID     = track -> GetParticleCode();
                    HitTime         = rand2 -> Gaus( track -> GetHitTime() , 0.200*CLHEP::ns ) ; //  ‘smear’ each TOF by the resolution: a gaussian of width of 200ps
                    l_reconstructed = sqrt(x*x + y*y + z*z);
                    v_reconstructed = l_reconstructed/HitTime;
                    gamma_reconstrc = 1./sqrt( 1 - pow( (v_reconstructed / CLHEP::c_light ),2) );
                    p_reconstructed = (CLHEP::neutron_mass_c2)*(v_reconstructed/CLHEP::c_light)*gamma_reconstrc; // = mv𝞬
                    delta_p         = 100*(p_reconstructed - p_g)/p_g;
                    NeutronReachedDet = 1;
                    ParticleHitDetector = true;
                }
                
                // take the energy, momentum, and energy deposition information from the traveling within the detector
                E_kin               = track -> GetKineticEnergy();
                p                   = track -> GetMomentum();
            }
        }
        MomentumLoss    = p_g - p;
    }
    if ( PrintData ){
        std::cout << "\n\n**************************** Event: " << eventID << " ********************************************\n\n"
        << "\n Generated Particle Code: "       << particle_ID_g
        << "\n Generated Momentum: "            << G4BestUnit( p_g , "Energy" )
        << "\n Generated Kinetic Energy: "      << std::setw(7) << G4BestUnit( E_kin_g , "Energy" )
        << "\n Generated x: "                   << std::setw(7) << G4BestUnit( x_g , "Length" )
        << "\n Generated y: "                   << std::setw(7) << G4BestUnit( y_g , "Length" )
        << "\n theta_g: "                       << ((primary_particle -> GetMomentum()).theta())*(360./CLHEP::twopi)
        << "\n Was suppose to hit sector: "     << Sector_g
        << "\n Has the particle hit the Detector?: " << NeutronReachedDet
        << "\n Hit sector: "                    << Sector
        << "\n Recieved Particle Code: "        << particle_ID
        << "\n Particle Name: "                 << PrtclName
        << "\n Momentum leaving detector: "     << G4BestUnit( p , "Energy" )
        << "\n Kinetic Energy leaving det: "    << std::setw(7) << G4BestUnit( E_kin , "Energy" )
        << "\n Detected x: "                    << std::setw(7) << G4BestUnit( x , "Length" )
        << "\n Detected y: "                    << std::setw(7) << G4BestUnit( y , "Length" )
        << "\n Detected z: "                    << std::setw(7) << G4BestUnit( z , "Length" )
        << "\n Momentum Loss: "                 << std::setw(7) << G4BestUnit( MomentumLoss , "Energy" )
        << "\n Energy Deposition: "             << std::setw(7) << G4BestUnit( TotEnergyDep , "Energy" )
        << "\n Hit Time: "                      << std::setw(7) << G4BestUnit( HitTime , "Time" )
        << "\n l_reconstructed: "               << std::setw(7) << G4BestUnit( l_reconstructed , "Length" )
//        << "\n v_reconstructed: "               << std::setw(7) << G4BestUnit( v_reconstructed , "Speed" )
        << "\n beta: "                          << std::setw(7) << v_reconstructed / CLHEP::c_light
        << "\n gamma_reconstrc: "               << std::setw(7) << gamma_reconstrc
        << "\n reconstructed momentum: "        << std::setw(7) << G4BestUnit( p_reconstructed , "Energy" )
        << "\n (p(rec) - p(gen)) / p(gen): "    << std::setw(7) << (p_reconstructed - p_g)/p_g
        << "\n-------------------- Fill  - End ---------------------\n"
        << std::endl;
    }
    
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man -> FillNtupleIColumn( 0 , particle_ID_g );
    man -> FillNtupleDColumn( 1 , E_kin_g/CLHEP::MeV  );
    man -> FillNtupleDColumn( 2 , p_g/CLHEP::MeV  );
    man -> FillNtupleDColumn( 3 , x_g/CLHEP::mm  );
    man -> FillNtupleDColumn( 4 , y_g/CLHEP::mm  );
    man -> FillNtupleDColumn( 5 , theta_g  );
    man -> FillNtupleIColumn( 6 , Sector_g  );
    man -> FillNtupleIColumn( 7 , NeutronReachedDet  );
    man -> FillNtupleIColumn( 8 , particle_ID );
    man -> FillNtupleDColumn( 9 , E_kin/CLHEP::MeV  );
    man -> FillNtupleDColumn( 10 , p/CLHEP::MeV  );
    man -> FillNtupleDColumn( 11 , Theta  );
    man -> FillNtupleDColumn( 12 , x/CLHEP::mm  );
    man -> FillNtupleDColumn( 13 , y/CLHEP::mm  );
    man -> FillNtupleDColumn( 14 , z/CLHEP::mm  );
    man -> FillNtupleIColumn( 15 , Sector  );
    man -> FillNtupleDColumn( 16 , MomentumLoss/CLHEP::MeV  );
    man -> FillNtupleDColumn( 17 , TotEnergyDep/CLHEP::MeV  );
    man -> FillNtupleDColumn( 18 , HitTime/CLHEP::ns  );
    man -> FillNtupleDColumn( 19 , p_reconstructed/CLHEP::MeV  );
    man -> FillNtupleDColumn( 20 , delta_p  );
    man -> AddNtupleRow();
  */
}

