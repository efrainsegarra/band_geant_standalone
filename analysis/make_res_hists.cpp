#include <iostream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TVector3.h"

using namespace std;

double sq(double x){ return x*x;};
const double mN=0.93956541;

int main(int argc, char ** argv)
{
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead use\n"
	   << "\t/path/to/smeared/file /path/to/output/file\n";
      exit(-1);
    }

  TFile * f = new TFile(argv[1]);
  TTree * t = (TTree*)f->Get("ReconTree");
  
  // Set memory addresses for tree
  double trueXp, trueWp, reconXp, reconWp;
  double mom_e[3], mom_r[3];
  t->SetBranchAddress("reconPe_x",&(mom_e[0]));
  t->SetBranchAddress("reconPe_y",&(mom_e[1]));
  t->SetBranchAddress("reconPe_z",&(mom_e[2]));
  t->SetBranchAddress("reconPr_x",&(mom_r[0]));
  t->SetBranchAddress("reconPr_y",&(mom_r[1]));
  t->SetBranchAddress("reconPr_z",&(mom_r[2]));
  t->SetBranchAddress("trueXp",&trueXp);
  t->SetBranchAddress("trueWp",&trueWp);
  t->SetBranchAddress("reconXp",&reconXp);
  t->SetBranchAddress("reconWp",&reconWp);

  // Prep histograms
  int nXpBins=20;
  int nWpBins=26;
  TH2D * xPrimeRes[5];
  TH2D * WPrimeRes[5];
  double minAlpha[5];
  double maxAlpha[5];

  // Create 2D histograms based on alpha bins
  for (int i=0 ; i<5 ; i++)
    {
      minAlpha[i] = 1.30 + 0.05*i;
      maxAlpha[i] = 1.35 + 0.05*i;
      char name[100];
      char title[100];

      sprintf(name,"xPrimeRes_%d",i);
      sprintf(title,"Recon Q^2 > 2., W_prime > 1.8, %g < alpha_s < %g;x_prime;delta x_prime;Counts",minAlpha[i],maxAlpha[i]);
      xPrimeRes[i]=new TH2D(name,title,nXpBins,0.2,0.6,80,-0.2,.2);

      sprintf(name,"WPrimeRes_%d",i);
      sprintf(title,"Recon Q^2 > 2., W_prime > 1.8, %g < alpha_s < %g;W_prime;delta W_prime;Counts",minAlpha[i],maxAlpha[i]);
      WPrimeRes[i] = new TH2D(name,title,nWpBins,1.8,3.1,80,-0.4,0.4);
    }

  // Loop over events
  int nCounts[5]={0,0,0,0,0};
  const int nEvents = t->GetEntries();
  for (int i=0 ; i < nEvents ; i++)
    //for (int i=0 ; i < 1000000 ; i++)
    {
      if (i%10000==0)
	cerr << "Working on event " << i << "\n";

      t->GetEvent(i);

      double p_e = sqrt(sq(mom_e[0]) +sq(mom_e[1]) + sq(mom_e[2]));
      double cosTheta = mom_e[2]/p_e;
      double QSq = 2.*(10.9)*p_e*(1.-cosTheta);
      // Make kinematic cuts
      if ((QSq > 2.) && (reconWp > 1.8))
	{
	  // Reconstruct alpha_s
	  double nu = 10.9 - p_e;
	  TVector3 qVec(-mom_e[0],-mom_e[1],10.9-mom_e[2]);
	  TVector3 rVec(mom_r[0],mom_r[1],mom_r[2]);
	  double E_r = sqrt(sq(mN) + rVec.Mag2());
	  double alpha_s = (E_r - rVec.Dot(qVec)/qVec.Mag())/mN; // ERROR IN THIS LINE

	  // Figure out what bin
	  int aBin = floor((alpha_s - minAlpha[0])/0.05);

	  //cout << "alpha_s is " << alpha_s << " and falls in bin " << aBin <<"\n";

	  if ((aBin<0)||(aBin>4))
	    continue;
	  
	  nCounts[aBin]++;
	  xPrimeRes[aBin]->Fill(trueXp,reconXp-trueXp);
	  WPrimeRes[aBin]->Fill(trueWp,reconWp-trueWp);
	}
    }

  TFile *outRoot = new TFile(argv[2],"RECREATE");
  
  for (int aBin=0 ; aBin <5 ; aBin++)
    {
      xPrimeRes[aBin]->Write();
      WPrimeRes[aBin]->Write();
      cout << "There were " << nCounts[aBin] << " in bin " << aBin << ".\n";
    }

  outRoot->Close();  
  f->Close();
}
