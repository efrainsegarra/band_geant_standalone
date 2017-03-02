#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorT.h"

#include "constants.h"
#include "gen_tree.h"

using namespace std;

const double minMomR = 0.275;
const double maxMomR = 0.700;
const double maxThetaR = 170.*M_PI/180.;
const double minThetaR = 160.*M_PI/180.;
const double bandZ = -250.; // cm
const double csN = 4.E5 / (0.1 * (6.E36 * 1.E-33)); // nb/sr

inline double sq(double x){ return x*x; };
inline double beta(double p){ return 1./sqrt(1. + sq(mN/p)); };

int main(int argc, char ** argv)
{
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead use:\n\t/path/to/inclusive/file /path/to/outputfile\n";
      exit(-1);
    }

  // Prep
  TRandom3 * myRand = new TRandom3(0);
  const double maxBetaR = beta(maxMomR);
  const double minBetaR = beta(minMomR);
  const double timeWindow = (fabs(bandZ)/(cAir))*((1./minBetaR) - (1./maxBetaR));
  const double minCosThetaR = cos(maxThetaR);
  const double maxCosThetaR = cos(minThetaR);
  
  // Open the files
  TFile *infile = new TFile(argv[1]);
  TFile *outfile = new TFile(argv[2],"RECREATE");

  // Get the cross section information
  TVectorT<double> * csVec = (TVectorT<double>*)infile->Get("totalCS");

  // Get the tree and assign branches
  TTree * inTree = (TTree*) infile->Get("MCout");
  Gen_Event * inEvent = NULL;
  inTree->SetBranchAddress("event",&inEvent);

  // Produce an output tree
  TTree * outTree = new TTree("MCout","Random Coincidence Output");
  Gen_Event * outEvent = new Gen_Event;
  outTree->Branch("event",&outEvent);

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
      double tR = myRand->Rndm()*timeWindow + fabs(bandZ)/(cAir*maxBetaR);
      double betaR = fabs(bandZ)/(cAir*tR);
      double momR = mN /sqrt(1./sq(betaR) - 1.);
      double cosThetaR = minCosThetaR + myRand->Rndm()*(maxCosThetaR-minCosThetaR);
      double thetaR = acos(cosThetaR);
      double phiR = 2.*M_PI * myRand->Rndm();

      // Write tree
      outEvent->particles.clear();
      outEvent->particles.push_back(inEvent->particles[0]);
      Gen_Particle neutron;
      neutron.type="neutron";
      neutron.momentum.SetMagThetaPhi(momR,thetaR,phiR);
      outEvent->particles.push_back(neutron);
      outTree->Fill();
    }
 
  // Update cross section info
  double neutronCS = csN * (2.*M_PI)*(maxCosThetaR - minCosThetaR);
  TVectorT<double> csSqVec(3);
  csSqVec[0]=neutronCS*timeWindow*(*csVec)[0]; // Units of nb^2 * s
  csSqVec[1]=neutronCS*timeWindow*(*csVec)[1];
  csSqVec[2]=timeWindow;
  csSqVec.Write("totalCSSq");

  // Clean-up
  gFile = outfile;
  outTree->Write();
}
