#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorT.h"
#include "TH1.h"

#include "constants.h"
#include "gen_tree.h"

using namespace std;

const double minMomR = 0.2;
const double maxMomR = 1.0;
const double maxThetaR = lad_max_theta_deg*M_PI/180.;
const double minThetaR = lad_min_theta_deg*M_PI/180.;
const double max_lad_phi=20.*M_PI/180.;
const double csP = 500.; // nb
const double timeWindow = 20.; // ns

inline double sq(double x){ return x*x; }
inline double beta(double p){ return 1./sqrt(1. + sq(mP/p)); }

int main(int argc, char ** argv)
{
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead use:\n\t/path/to/inclusive/file /path/to/outputfile\n";
      exit(-1);
    }

  // Prep
  TRandom3 * myRand = new TRandom3(0);
  const double minCosThetaR = cos(maxThetaR);
  const double maxCosThetaR = cos(minThetaR);
  
  // Fill the proton momentum distribution histogram
  TH1D * protonMomDist = new TH1D("protonMomDist","Proton Momentum Distribution;Momentum [GeV];Prob",725,0.275,1.);
  for (int i=1 ; i<=275 ; i++)
    {
      double momMeV = 1000.*protonMomDist->GetBinCenter(i);
      double probPerMeV = 1400. - 4.18*momMeV + 0.0041*momMeV*momMeV - 1.3E-6*momMeV*momMeV*momMeV;
      if (probPerMeV < 0.) probPerMeV = 0.;
      protonMomDist->SetBinContent(i,probPerMeV*1000.);
    }

  // Open the files
  TFile *infile = new TFile(argv[1]);
  TFile *outfile = new TFile(argv[2],"RECREATE");

  // Get the cross section information
  TVectorT<double> * csVec = (TVectorT<double>*)infile->Get("totalCS");

  // Get the tree and assign branches
  TTree * inTree = (TTree*) infile->Get("MCout");
  Gen_Event * inEvent = NULL;
  double ze, zr;
  inTree->SetBranchAddress("event",&inEvent);
  inTree->SetBranchAddress("ze",&ze);

  // Produce an output tree
  TTree * outTree = new TTree("MCout","Random Coincidence Output");
  Gen_Event * outEvent = new Gen_Event;
  outTree->Branch("event",&outEvent);
  outTree->Branch("ze",&ze,"ze/D");
  outTree->Branch("zr",&zr,"zr/D");
  double t0;
  outTree->Branch("t0",&t0,"t0/D");

  // Loop over the events
  const int nEvents = inTree->GetEntries();
  for (int i=0 ; i<nEvents ; i++)
    {
      inTree->GetEntry(i);

      // Require that there be one particle in the event
      if (inEvent->particles.size() != 1)
	{
	  cerr << "Event " << i << " is not an inclusive event! It has " << inEvent->particles.size() << "particles!!!\n";
	  continue;
	}

      // Generate random neutron
      t0 = (myRand->Rndm()-0.5)*timeWindow;
      double momR = protonMomDist->GetRandom();
      double cosThetaR = minCosThetaR + myRand->Rndm()*(maxCosThetaR-minCosThetaR);
      double thetaR = acos(cosThetaR);
      double phiR = 2.*max_lad_phi * (myRand->Rndm()-0.5);

      // Vertex
      zr = lad_target_z*(myRand->Rndm()-0.5);

      // Write tree
      outEvent->particles.clear();
      outEvent->particles.push_back(inEvent->particles[0]);
      Gen_Particle proton;
      proton.type="proton";
      proton.momentum.SetMagThetaPhi(momR,thetaR,phiR);
      outEvent->particles.push_back(proton);
      outTree->Fill();
    }
 
  // Update cross section info
  double protonCS = csP * (2.*max_lad_phi)*(maxCosThetaR - minCosThetaR);
  TVectorT<double> csSqVec(3);
  csSqVec[0]=protonCS*timeWindow*(*csVec)[0]; // Units of nb^2 * s
  csSqVec[1]=protonCS*timeWindow*(*csVec)[1];
  csSqVec[2]=timeWindow;
  csSqVec.Write("totalCSSq");

  // Clean-up
  gFile = outfile;
  outTree->Write();
  delete protonMomDist;

  return 0;
}
