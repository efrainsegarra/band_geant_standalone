#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorT.h"

#include "constants.h"

#include "gen_tree.h"

using namespace std;

// Foam ranges
//const double min_theta_e = 5.*M_PI/180.;
//const double max_theta_e = 35.*M_PI/180.;
//const double min_p_e = 2.;
//const double max_p_e = 8.;

const double min_theta_r =M_PI;
const double max_theta_r =M_PI;
const double min_p_r =0.25;
const double max_p_r =0.25;
const double min_phi_r=0.;
const double max_phi_r=2*M_PI;

//const double min_phi_er=0.;
//const double max_phi_er=2.*M_PI;

int main(int argc, char *argv[])
{
  // Read in arguments
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead use: [nEvents] [/path/to/output/file]\n";
      exit(-1);
    }

  const int nEvents=atoi(argv[1]);
  TFile * outputFile = new TFile(argv[2],"RECREATE");

  // Create a tree
  TTree * outputTree = new TTree("MCout","Generator Output");
  
  // Initialize the branches
  Gen_Event * thisEvent = new Gen_Event;
  outputTree->Branch("event",&thisEvent);

  // Random number generator
  TRandom3 * rand = new TRandom3(0);


  for (int i=0 ; i<nEvents ; i++)
    {
      if (nEvents >= 10000)
  if (i%10000==0)
    cerr << "Working on event " << i << "\n";


      // Extract useful quantities
      double theta_r = min_theta_r + (max_theta_r - min_theta_r) * rand->Rndm();
      double p_r = min_p_r + (max_p_r - min_p_r) * rand->Rndm();
      double phi_r = min_phi_r + (max_phi_r - min_phi_r) * rand->Rndm();


      // Write to tree
      thisEvent->particles.clear();
      //Gen_Particle electron;
      //electron.type="e-";
      //electron.momentum.SetMagThetaPhi(p_e,theta_e,phi_e);
      //thisEvent->particles.push_back(electron);
      Gen_Particle neutron;
      neutron.type="neutron";
      neutron.momentum.SetMagThetaPhi(p_r,theta_r,phi_r);
      thisEvent->particles.push_back(neutron);
      outputTree->Fill();
    }
   
  outputTree->Write();


  outputFile->Close();

  
  delete rand;

  return 0;
}

