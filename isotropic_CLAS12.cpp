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

const double max_y = 7.4*13;
const double min_y = -7.4*5;
const double max_z = -262.;
const double min_z = -262.-5*7.4;
const double glob_max_x = 2018.35 / 10 / 2;
const double glob_min_x = -2018.35 / 10 / 2;
const double min_p_r =0.25;
const double max_p_r =0.25;

int checkBAND(double x_traj,double y_traj){

  // ONLY SAVE THE EVENTS THAT WILL MAKE IT INTO BAND
  double numbars,prevbars;
  double max_x, min_x,max_x_2,min_x_2;
 
  // first group
  prevbars = 0;
  numbars = 2;
  if ((y_traj >= min_y + prevbars*7.4) && (y_traj < min_y + (numbars+prevbars)*7.4)){
    max_x = 2018.35 / 10 / 2;
    min_x = -max_x;

    if ((x_traj > max_x) || (x_traj < min_x)) return 0;

  }
  // second group
  prevbars += numbars;
  numbars = 6;
  if ((y_traj >= min_y + prevbars*7.4) && (y_traj < min_y + (numbars+prevbars)*7.4)){
    max_x = (509.17/10 + 500/10);
    min_x = (509.17/10 );
    min_x_2 = -(509.17/10 + 500/10);
    max_x_2 = -(509.17/10 );

    if ((x_traj > max_x) || (x_traj < min_x_2)) return 0;
    if ((x_traj < min_x) && (x_traj > max_x_2)) return 0;

  }

  // third group
  prevbars += numbars;
  numbars = 4;
  if ((y_traj >= min_y + prevbars*7.4) && (y_traj < min_y + (numbars+prevbars)*7.4)){
    max_x = 2018.35 / 10 / 2;
    min_x = -max_x;

    if ((x_traj > max_x) || (x_traj < min_x)) return 0;

  }
  // fourth group
  prevbars += numbars;
  numbars = 3;
  if ((y_traj >= min_y + prevbars*7.4) && (y_traj < min_y + (numbars+prevbars)*7.4)){
    max_x = 1946.53 / 10 / 2;
    min_x = -max_x;

    if ((x_traj > max_x) || (x_traj < min_x)) return 0;

  }
  // five group
  prevbars += numbars;
  numbars = 3;
  if ((y_traj >= min_y + prevbars*7.4) && (y_traj < min_y + (numbars+prevbars)*7.4)){
    max_x = 1634.58 / 10 / 2;
    min_x = -max_x;

    if ((x_traj > max_x) || (x_traj < min_x)) return 0;

  }
  
  return 1;
}


int main(int argc, char *argv[])
{
  // Read in arguments
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead use: [nEvents] [/path/to/output/file]\n";
      exit(-1);
    }

  double willHitBand = 0.;

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
      if (nEvents >= 10000) if (i%10000==0) cerr << "Working on event " << i << "\n";

      double z_r = min_z + (max_z - min_z) * rand->Rndm();
      double y_r = min_y + (max_y - min_y) * rand->Rndm();
      double x_r = glob_min_x + (glob_max_x - glob_min_x) * rand->Rndm();
      int pass = checkBAND(x_r,y_r);
      if (pass == 0) continue;
      willHitBand++;

      double Costheta_r = z_r / sqrt( pow(x_r,2) + pow(y_r,2) + pow(z_r,2) );
      double theta_r = acos(Costheta_r);
      double phi_r = atan2(y_r , x_r);

      if (phi_r < 0) phi_r += 2.*M_PI;
      if (phi_r > 2.*M_PI) phi_r -= 2.*M_PI;


      //cout << x_r << " " << y_r << " " << z_r << " " << theta_r << " " << phi_r << "\n";

      // Extract useful quantities
      double p_r = min_p_r + (max_p_r - min_p_r) * rand->Rndm();

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
  cout << "Percent of Generated Hitting BAND: " << willHitBand / float(nEvents) << endl;
  
  delete rand;

  return 0;
}

