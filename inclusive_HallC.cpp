#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "TFile.h"
#include "TTree.h"
#include "TFoam.h"
#include "TFoamIntegrand.h"
#include "TRandom3.h"
#include "TVectorT.h"

using namespace std;

#include "constants.h"
#include "deuteronwf.h"
#include "crossdis.h"
#include "crossincl.h"
#include "gen_tree.h"

// Settings for cross section: 
double sigmainput=40.;          //sigma parameter in rescattering amplitude [mb], default value
double betainput=8.;            //beta parameter in rescattering amplitude [GeV^-2], default value
double epsinput=-0.5;           //epsilon parameter in rescattering amplitude [], default value
double lambdainput=1.2;
double betaoffinput=8.;
int offshellset=3;              // 0=mass diff suppression, 1=dipole suppression, 2=beta parameterization, 3=no offshell, 4=full off-shell
int symm=0;                     //symmetric FSI in inclusive reaction
int phiavg=1;                   //average cross section over phi
int F_param=0;                  //structure functions parametrization 0=SLAC, 1=Christy&Bosted, 2=Alekhin et al leading twist

// Constant settings
const int calc=0;
const int which_wave=1;
const int this_decay=0;
const int num_res=1;

// Memory
double *c, *d, *m;
int c_length=0;

extern "C"{
    extern struct{
        char homedir[128];
        char a09file1[128];
        char a09file2[128];
    } dir_;
}//shared fortran variables

// Foam ranges
const double min_theta_e=hallc_min_theta;
const double max_theta_e=hallc_max_theta;
double max_phi_e;
double min_phi_e;
double min_p_e;
double max_p_e;

// Foam Integrand
class InclusiveCS : public TFoamIntegrand
{
public:
  double Density(int nDim, double * args);
};

int main(int argc, char *argv[])
{
  // Read in arguments
  if (argc != 4)
    {
      cerr << "Wrong number of arguments. Instead use: [nEvents] [/path/to/output/file] [0 for HMS, 1 for SHMS] \n";
      exit(-1);
    }
  const int nEvents=atoi(argv[1]);
  TFile * outputFile = new TFile(argv[2],"RECREATE");

  // Set ranges
  if (atoi(argv[3])==0)
    {
      min_p_e=0.4;
      max_p_e=9.0;
      min_phi_e=M_PI-hallc_phi_range;
      max_phi_e=M_PI+hallc_phi_range;
    }
  else if (atoi(argv[3])==1)
    {
      min_p_e=1.8;
      max_p_e=10.;
      min_phi_e=-hallc_phi_range;
      max_phi_e= hallc_phi_range;
    }

  // Run sanity check on parameters
  sanitycheck(calc,0,F_param,which_wave,offshellset);

  // Get some parameters
  get_wf_param(&c, &d, &m, c_length, which_wave);

  // Read in data files
  strcpy(dir_.homedir,getenv("HOME"));
  strcpy(dir_.a09file1,getenv("HOME"));
  strcpy(dir_.a09file2,getenv("HOME"));
  strcat(dir_.a09file1,"/.deuteron_dis/grids/a09.sfs_lNNC");
  strcat(dir_.a09file2,"/.deuteron_dis/grids/a09.dsfs_lNNC");

  // Create a tree
  TTree * outputTree = new TTree("MCout","Generator Output");
  
  // Initialize the branches
  Gen_Event * thisEvent = new Gen_Event;
  double ze;
  outputTree->Branch("event",&thisEvent);
  outputTree->Branch("ze",&ze,"ze/D");

  // Random number generator
  TRandom3 * rand = new TRandom3(0);

  // Initialize the foam
  TFoam * csFoam = new TFoam("csFoam");
  InclusiveCS * csTotal = new InclusiveCS();
  csFoam->SetkDim(2);
  csFoam->SetRho(csTotal);
  csFoam->SetPseRan(rand);
  // optional
  csFoam->SetnCells(100);
  csFoam->SetnSampl(100);
  // initialize
  csFoam->Initialize();

  // Create memory for each event
  double * eventData = new double[2];

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

      // Handle phi
      double phi_e = min_phi_e + (max_phi_e-min_phi_e)*rand->Rndm();

      // Vertex
      ze = lad_target_z*(rand->Rndm()-0.5);

      // Write to tree
      thisEvent->particles.clear();
      Gen_Particle electron;
      electron.type="e-";
      electron.momentum.SetMagThetaPhi(p_e,theta_e,phi_e);
      electron.t0=0.;
      thisEvent->particles.push_back(electron);
      outputTree->Fill();
    }
   
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

double InclusiveCS::Density(int nDim, double *args)
{
  // Variables (there should be 2 of them)
  double theta_e = min_theta_e + args[0]*(max_theta_e - min_theta_e);
  double p_e = min_p_e + args[1]*(max_p_e - min_p_e);

  // Develop derived quantities
  double Q2 = 4.*E1*p_e * sq(sin(0.5*theta_e));
  double nu = E1 - p_e;
  double x = Q2 / (2.*mN*nu);

  // Some memory for the cross section calculation results
  double QEpw, DISpw, DISfsi, DISfsi2;
  
  // Calculate the cross section
  calc_inclusive2(QEpw, DISpw, DISfsi, DISfsi2, E1, Q2,x, c, d, m, c_length, which_wave, offshellset, this_decay, num_res, calc);
  double dsigma=DISpw+QEpw-DISfsi;
  double jacobian = x*E1*p_e/(M_PI*nu);
  double differential = (max_theta_e - min_theta_e)*(max_phi_e - min_phi_e)*(max_p_e - min_p_e)*sin(theta_e);

  return dsigma*jacobian*differential;
}

