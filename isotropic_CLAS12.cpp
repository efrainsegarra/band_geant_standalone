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
const double min_theta_r = 170.*M_PI/180.;//154.3*M_PI/180.;
const double max_theta_r = 170.*M_PI/180.;//175.8*M_PI/180.;
const double min_p_r =0.55;
const double max_p_r =0.55;
const double min_phi_r=M_PI/2.;
const double max_phi_r=M_PI/2.;

int checkBAND(double theta_r,double phi_r){
  // ONLY SAVE THE EVENTS THAT WILL MAKE IT INTO BAND
  double numbars,prevbars;
  double max_x, min_x,max_x_2,min_x_2;
  double max_y = 7.4*13;
  double min_y = -7.4*5;

  double r = -262. / cos(theta_r) ;
  double x_traj = r * cos(phi_r) * sin(theta_r) ;
  double y_traj = r * sin(phi_r) * sin(theta_r) ;
  if ((y_traj > max_y) || (y_traj < min_y)) return -1;

  // first group
  prevbars = 0;
  numbars = 2;
  if ((y_traj > min_y + prevbars*7.4) && (y_traj < min_y + (numbars+prevbars)*7.4)){
    max_x = 2018.35 / 10 / 2;
    min_x = -max_x;
    if ((x_traj > max_x) || (x_traj < min_x)) return -1;
  }
  // second group
  prevbars += numbars;
  numbars = 6;
  if ((y_traj > min_y + prevbars*7.4) && (y_traj < min_y + (numbars+prevbars)*7.4)){
    max_x = (509.17/10 + 500/10);
    min_x = (509.17/10 );
    max_x_2 = -(509.17/10 );
    min_x_2 = -(509.17/10 + 500/10);
    if ((x_traj > max_x) || (x_traj < min_x_2)) return -1;
    if (( x_traj < min_x ) && (x_traj > max_x_2) ) return -1;
  }
  // third group
  prevbars += numbars;
  numbars = 4;
  if ((y_traj > min_y + prevbars*7.4) && (y_traj < min_y + (numbars+prevbars)*7.4)){
    max_x = 2018.35 / 10 / 2;
    min_x = -max_x;
    if ((x_traj > max_x) || (x_traj < min_x)) return -1;
  }
  // fourth group
  prevbars += numbars;
  numbars = 3;
  if ((y_traj > min_y + prevbars*7.4) && (y_traj < min_y + (numbars+prevbars)*7.4)){
    max_x = 1946.53 / 10 / 2;
    min_x = -max_x;
    if ((x_traj > max_x) || (x_traj < min_x)) return -1;
  }
  // five group
  prevbars += numbars;
  numbars = 3;
  if ((y_traj > min_y + prevbars*7.4) && (y_traj < min_y + (numbars+prevbars)*7.4)){
    max_x = 1634.58 / 10 / 2;
    min_x = -max_x;
    if ((x_traj > max_x) || (x_traj < min_x)) return -1;
  }
  cout << x_traj << " " << y_traj << "\n";
  return 0;
}

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

      int pass = checkBAND(theta_r,phi_r);
      if (pass < 0) continue;
      
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

