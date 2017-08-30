#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TFoam.h"
#include "TFoamIntegrand.h"
#include "TRandom3.h"
#include "TVectorT.h"

#include "constants.h"
#include "deuteronwf.h"
#include "crossdis.h"

#include "gen_tree.h"

using namespace std;

int checkBAND(double theta_r,double phi_r){
  // ONLY SAVE THE EVENTS THAT WILL MAKE IT INTO BAND
  double numbars,prevbars;
  double max_x, min_x,max_x_2,min_x_2;
  double max_y = 7.4*13;
  double min_y = -7.4*5;

  double r = -262. / cos(theta_r) ;
  double x_traj = r * cos(phi_r) * sin(theta_r) ;
  double y_traj = r * sin(phi_r) * sin(theta_r) ;

  double r_2 = (-262.-7.4*5) / cos(theta_r) ;
  double x_traj_2 = r_2 * cos(phi_r) * sin(theta_r) ;
  double y_traj_2 = r_2 * sin(phi_r) * sin(theta_r) ;

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
    min_x_2 = -(509.17/10 + 500/10);
    if ((x_traj > max_x) || (x_traj < min_x_2)) return -1;
  }
  if ((y_traj_2 > min_y + prevbars*7.4) && (y_traj_2 < min_y + (numbars+prevbars)*7.4)){
    min_x = (509.17/10 );
    max_x_2 = -(509.17/10 );
    if (( x_traj_2 < min_x ) && (x_traj_2 > max_x_2) ) return -1;
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
  return 1;
}


// Settings for cross section: 
double sigmainput=40.;          //sigma parameter in rescattering amplitude [mb], default value
double betainput=8.;            //beta parameter in rescattering amplitude [GeV^-2], default value
double epsinput=-0.5;           //epsilon parameter in rescattering amplitude [], default value
double lambdainput=1.2;
double betaoffinput=8.;
int offshellset=3;              // 0=mass diff suppression, 1=dipole suppression, 2=beta parameterization, 3=no offshell, 4=full off-shell
int symm=0;                     //symmetric FSI in inclusive reaction
int phiavg=0;                   //average cross section over phi
int F_param=0;                  //structure functions parametrization 0=SLAC, 1=Christy&Bosted, 2=Alekhin et al leading twist

extern "C"{
    extern struct{
        char homedir[128];
        char a09file1[128];
        char a09file2[128];
    } dir_;
}//shared fortran variables

// Foam ranges
const double min_theta_e = 5.*M_PI/180.;
const double max_theta_e = 35.*M_PI/180.;
const double min_p_e = 2.;
const double max_p_e = 8.;

const double min_theta_r =150*M_PI/180.;
const double max_theta_r =178*M_PI/180.;
const double min_p_r =0.200;
const double max_p_r =0.600;

const double min_phi_er=0.;
const double max_phi_er=2.*M_PI;

// Foam Integrand
class disCS : public TFoamIntegrand
{
public:
  double Density(int nDim, double * args);
};


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

  // Set external variables
  if(offshellset==0) epsinput=-0.5;
  if(offshellset==1) lambdainput=0;
  if(offshellset==2) betaoffinput=16.8;
  if(offshellset==3) epsinput=-0.5;
  if(offshellset==4) epsinput=-0.5;
 
  strcpy(dir_.homedir,getenv("HOME"));
  strcpy(dir_.a09file1,getenv("HOME"));
  strcpy(dir_.a09file2,getenv("HOME"));
  strcat(dir_.a09file1,"/.deuteron_dis/grids/a09.sfs_lNNC");
  strcat(dir_.a09file2,"/.deuteron_dis/grids/a09.dsfs_lNNC");

  // Create a tree
  TTree * outputTree = new TTree("MCout","Generator Output");
  
  // Initialize the branches
  Gen_Event * thisEvent = new Gen_Event;
  outputTree->Branch("event",&thisEvent);

  // Random number generator
  TRandom3 * rand = new TRandom3(0);

  // Initialize the foam
  TFoam * csFoam = new TFoam("csFoam");
  disCS * csTotal = new disCS();
  csFoam->SetkDim(5);
  csFoam->SetRho(csTotal);
  csFoam->SetPseRan(rand);
  // optional
  csFoam->SetnCells(10000);
  csFoam->SetnSampl(1000);
  // initialize
  csFoam->Initialize();

  // Create memory for each event
  double * eventData = new double[5];
  int willHitBAND = 0;

  for (int i=0 ; i<nEvents ; i++)
    {
      if (nEvents >= 10000)
	if (i%10000==0)
	  cerr << "Working on event " << i << "\n";

      csFoam->MakeEvent();
      csFoam->GetMCvect(eventData);

      // Extract useful quantities
      double theta_e = min_theta_e + eventData[0]*(max_theta_e - min_theta_e);
      double p_e = min_p_e + eventData[1]*(max_p_e - min_p_e);

      double theta_r = min_theta_r + eventData[2]*(max_theta_r - min_theta_r);
      double p_r = min_p_r + eventData[3]*(max_p_r - min_p_r);

      double phi_er = min_phi_er + eventData[4]*(max_phi_er - min_phi_er);

      // Handle phi (** DANGER ** )
      double phi_e = 2.*M_PI * rand->Rndm();
      double phi_r = phi_e + phi_er;
      if (phi_r > 2.*M_PI) phi_r -= 2.*M_PI;

      int pass = checkBAND(theta_r,phi_r);
      if (pass < 0) continue;
      willHitBAND++;

      // Write to tree
      thisEvent->particles.clear();

      Gen_Particle neutron;
      neutron.type="neutron";
      neutron.momentum.SetMagThetaPhi(p_r,theta_r,phi_r);
      thisEvent->particles.push_back(neutron);
      
      Gen_Particle electron;
      electron.type="e-";
      electron.momentum.SetMagThetaPhi(p_e,theta_e,phi_e);
      thisEvent->particles.push_back(electron);

      outputTree->Fill();
    }
  
  cout << "Of generated, the percent that have initial path to hit BAND: "<< float(willHitBAND)/float(nEvents)*100. << endl;

  outputTree->Write();

  // Store the total cross section (and error) in the output file in a TVectorT
  TVectorT<double> totalCS(2);
  csFoam->GetIntegMC(totalCS[0],totalCS[1]);
  totalCS.Write("totalCS");

  outputFile->Close();

  // Clean up
  delete csFoam;
  delete csTotal;
  delete rand;

  return 0;
}

inline double sq(double x) {return x*x;};

double disCS::Density(int nDim, double *args)
{
  // Parameters
  int proton = 1; // (0 for DIS on neutron with spectator proton, 1 for DIS on proton with spectator neutron).
  int which_wave =1;     // Non relativistic deuteron wf: 0 - Paris wf, 1 - AV18, 2 - CD Bonn, 3 - AV18*
  int decay=0;
  int num_res=1;

  // Variables (there should be 5 of them)
  double theta_e = min_theta_e + args[0]*(max_theta_e - min_theta_e);
  double p_e = min_p_e + args[1]*(max_p_e - min_p_e);
  double theta_r = min_theta_r + args[2]*(max_theta_r - min_theta_r);
  double p_r = min_p_r + args[3]*(max_p_r - min_p_r);
  double phi_er = min_phi_er + args[4]*(max_phi_er - min_phi_er);

  //cout << theta_e*180/M_PI << " " << p_e << " " << theta_r*180./M_PI << " " << p_r << " " << phi_er*180./M_PI << "\n";

  // Develop derived quantities
  double E_r = sqrt(sq(mP) + sq(p_r));
  double Q2 = 4.*E1*p_e * sq(sin(0.5*theta_e));
  double nu = E1 - p_e;
  double x = Q2 / (2.*mN*nu);
  double q = sqrt(Q2 + sq(nu));
  double theta_q = acos((E1 - p_e*cos(theta_e))/q);
  double theta_rq = acos( cos(theta_r)*cos(theta_q) - sin(theta_r)*sin(theta_q)*cos(phi_er));

  //cout << E_r << " " << Q2 << " " << nu << " " << x << " " << q << " " << theta_q*180./M_PI << " " << theta_rq*180./M_PI << "\n";

  // Sigma-input parameter
  double W_primeSq = sq(mD) - Q2 + sq(mN) + 2.*mD*(nu-E_r) -2.* nu * E_r + 2.*q*p_r*cos(theta_rq);
  if (W_primeSq < 0.)
    return 0.;

  double W_prime = sqrt(W_primeSq);
  sigmainput = (25.3+53*(W_prime-mN))/(Q2);

  double crosstotal1 = calc_cross(E1, Q2, x, p_r, theta_rq, phi_er, proton, which_wave, decay, num_res, 0);
  double jacobian = x*E1*p_e*p_r/(M_PI*nu);
  double differential_e = (max_theta_e - min_theta_e)*(2.*M_PI)*(max_p_e - min_p_e)*sin(theta_e);
  double differential_p = (max_theta_r - min_theta_r)*(max_phi_er - min_phi_er)*(max_p_r - min_p_r)*sin(theta_r)*p_r/E_r;
  
  double result = crosstotal1 * jacobian * differential_e * differential_p;

  if (result < 0.)
    return 0.;
  return result;
}

