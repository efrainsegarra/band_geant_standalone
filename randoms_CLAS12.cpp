#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorT.h"

#include "constants.h"

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
  double mom_e[3];
  TTree * inTree = (TTree*) infile->Get("MCout");
  inTree->SetBranchAddress("x_e",&(mom_e[0]));
  inTree->SetBranchAddress("y_e",&(mom_e[1]));
  inTree->SetBranchAddress("z_e",&(mom_e[2]));

  // Produce an output tree
  TTree * outTree = new TTree("MCout","Random Coincidence Output");
  double mom_r[3];
  outTree->Branch("x_e",&(mom_e[0]),"x_e/D");
  outTree->Branch("y_e",&(mom_e[1]),"y_e/D");
  outTree->Branch("z_e",&(mom_e[2]),"z_e/D");
  outTree->Branch("x_r",&(mom_r[0]),"x_r/D");
  outTree->Branch("y_r",&(mom_r[1]),"y_r/D");
  outTree->Branch("z_r",&(mom_r[2]),"z_r/D");

  // Loop over the events
  const int nEvents = inTree->GetEntries();
  for (int i=0 ; i<nEvents ; i++)
    {
      inTree->GetEntry(i);

      // Generate random neutron
      double tR = myRand->Rndm()*timeWindow + fabs(bandZ)/(cAir*maxBetaR);
      double betaR = fabs(bandZ)/(cAir*tR);
      double momR = mN /sqrt(1./sq(betaR) - 1.);
      double cosThetaR = minCosThetaR + myRand->Rndm()*(maxCosThetaR-minCosThetaR);
      double thetaR = acos(cosThetaR);
      double phiR = 2.*M_PI * myRand->Rndm();

      // Fill momentum vector
      mom_r[0] = momR*sin(thetaR)*cos(phiR);
      mom_r[1] = momR*sin(thetaR)*sin(phiR);
      mom_r[2] = momR*cosThetaR;

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
