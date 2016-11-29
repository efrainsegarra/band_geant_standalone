#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <complex>
#include <fstream>
#include <cstdarg>
#include <time.h>

#include "TFile.h"
#include "TTree.h"
#include "TFoam.h"
#include "TRandom3.h"

using namespace std;

#include "constants.h"
#include "deuteronwf.h"
#include "crossdis.h"

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
double csTotal(int nDim, double *args);
const double min_theta_e = 5.*M_PI/180.;
const double max_theta_e = 35.*M_PI/180.;
const double min_p_e = 2.;
const double max_p_e = 8.;
const double min_theta_r =160.*M_PI/180.;
const double max_theta_r =170.*M_PI/180.;
const double min_p_r =0.275;
const double max_p_r =0.600;
const double min_phi_er=0.;
const double max_phi_er=2.*M_PI;

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
  double mom_e[3], mom_r[3];
  outputTree->Branch("x_e",&(mom_e[0]),"x_e/D");
  outputTree->Branch("y_e",&(mom_e[1]),"y_e/D");
  outputTree->Branch("z_e",&(mom_e[2]),"z_e/D");
  outputTree->Branch("x_r",&(mom_r[0]),"x_r/D");
  outputTree->Branch("y_r",&(mom_r[1]),"y_r/D");
  outputTree->Branch("z_r",&(mom_r[2]),"z_r/D");

  // Random number generator
  TRandom3 * rand = new TRandom3(0);

  // Initialize the foam
  TFoam * csFoam = new TFoam("csFoam");
  csFoam->SetkDim(5);
  csFoam->SetRhoInt(csTotal);
  csFoam->SetPseRan(rand);
  csFoam->Initialize();

  // Create memory for each event
  double * eventData = new double[5];

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
      if (phi_r > 2.*M_PI)
	phi_r -= 2.*M_PI;

      // Write to tree
      mom_e[0] = p_e*sin(theta_e)*cos(phi_e);
      mom_e[1] = p_e*sin(theta_e)*sin(phi_e);
      mom_e[2] = p_e*cos(theta_e);
      mom_r[0] = p_r*sin(theta_r)*cos(phi_r);
      mom_r[1] = p_r*sin(theta_r)*sin(phi_r);
      mom_r[2] = p_r*cos(theta_r);
      outputTree->Fill();
    }
   
  outputTree->Write();
  outputFile->Close();

  // Clean up
  delete csFoam;
  delete rand;

  return 0;
}

inline double sq(double x) {return x*x;};

double csTotal(int nDim, double *args)
{
  // Parameters
  double Ein = 10.9; // In GeV
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

  // Develop derived quantities
  double Q2 = 4.*Ein*p_e * sq(sin(0.5*theta_e));
  double nu = Ein - p_e;
  double x = Q2 / (2.*MASSN*nu);
  double q = sqrt(Q2 + sq(nu));
  double theta_q = acos((Ein - p_e*cos(theta_e))/q);
  double theta_rq = acos( cos(theta_r)*cos(theta_q) - sin(theta_r)*sin(theta_q)*cos(phi_er));
  //printf("Derived: %g %g %g %g %g %g\n",Q2,nu,x,q,theta_q,theta_rq);

  // Sigma-input parameter
  double W_prime=sqrt( MASSD*MASSD-
		       Q2+MASSN*MASSN+2.*MASSD*((Q2/(2*MASSN*x))-sqrt(p_r*p_r+MASSP*MASSP))-
		       2.*(Q2/(2*MASSN*x))*sqrt(p_r*p_r+MASSP*MASSP)+
		       2.*sqrt(Q2+(Q2/(2*MASSN*x))*(Q2/(2*MASSN*x)))*p_r*cos(theta_rq));

  double result=0.; // This will be the final returned cross section
  if (W_prime > 2.)
    {      
      sigmainput = (25.3+53*(W_prime-MASSN))/(Q2);
      
      double crosstotal1 = calc_cross(Ein, Q2, x, p_r, theta_rq, phi_er, proton, which_wave, decay, num_res, 0);
      double jacobian = x*Ein*p_e*p_r/(M_PI*nu);

      result = crosstotal1 * jacobian;
    }

  //printf("%g %g %g %g %g %g\n",theta_e*180./M_PI, p_e, theta_r*180./M_PI, p_r, phi_er*180./M_PI,result);

  return result;
}

