#include <iostream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TVector3.h"

#include "constants.h"

using namespace std;

double sq(double x){ return x*x;};

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
  double mom_e[3], true_mom_e[3], mom_r[3], true_mom_r[3];
  t->SetBranchAddress("reconPe_x",&(mom_e[0]));
  t->SetBranchAddress("reconPe_y",&(mom_e[1]));
  t->SetBranchAddress("reconPe_z",&(mom_e[2]));
  t->SetBranchAddress("reconPr_x",&(mom_r[0]));
  t->SetBranchAddress("reconPr_y",&(mom_r[1]));
  t->SetBranchAddress("reconPr_z",&(mom_r[2]));
  t->SetBranchAddress("truePe_x",&(true_mom_e[0]));
  t->SetBranchAddress("truePe_y",&(true_mom_e[1]));
  t->SetBranchAddress("truePe_z",&(true_mom_e[2]));
  t->SetBranchAddress("truePr_x",&(true_mom_r[0]));
  t->SetBranchAddress("truePr_y",&(true_mom_r[1]));
  t->SetBranchAddress("truePr_z",&(true_mom_r[2]));
  t->SetBranchAddress("trueXp",&trueXp);
  t->SetBranchAddress("trueWp",&trueWp);
  t->SetBranchAddress("reconXp",&reconXp);
  t->SetBranchAddress("reconWp",&reconWp);

  // Prep histograms
  int nXpBins=20;
  int nWpBins=26;
  TH2D * trueXpAsFine =new TH2D("trueXpAsFine","Recon Q^2 > 2., W_prime > 1.8;x';a_s;Counts",100,0.,1.,100,1.,2.);
  TH2D * prRes = new TH2D("prRes","Recon Q^2 > 2., W_prime > 1.8;p_r;delta p_r;Counts",20,.28,.68,80,-0.1,0.1);
  TH2D * asRes = new TH2D("asRes","Recon Q^2 > 2., W_prime > 1.8;a_s;delta a_s;Counts",20,1.3,1.55,80,-0.2,0.2);
  TH2D * reconXpAsSame = new TH2D("reconXpAsSame","Stayed in same bin;Recon x_p;Recon a_s;Counts",5,0.15,0.65,5,1.3,1.55);
  TH2D * reconXpAsDiff = new TH2D("reconXpAsDiff","Originated in other bin;Recon x_p;Recon a_s;Counts",5,0.15,0.65,5,1.3,1.55);
  TH2D * trueXpAsDiff = new TH2D("trueXpAsDiff","Left the bin;True x_p;True a_s;Counts",5,0.15,0.65,5,1.3,1.55);
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
      double QSq = 2.*E1*p_e*(1.-cosTheta);
      // Make kinematic cuts
      if ((QSq > 2.) && (reconWp > 1.8))
	{
	  // Recoil momentum
	  TVector3 rVec_true(true_mom_r[0],true_mom_r[1],true_mom_r[2]);
	  TVector3 rVec(mom_r[0],mom_r[1],mom_r[2]);
	  prRes->Fill(rVec_true.Mag(),rVec.Mag() - rVec_true.Mag());

	  // Reconstruct alpha_s
	  double nu = E1 - p_e;
	  TVector3 qVec(-mom_e[0],-mom_e[1],E1-mom_e[2]);
	  double E_r = sqrt(sq(mN) + rVec.Mag2());
	  double alpha_s = (E_r - rVec.Dot(qVec)/qVec.Mag())/mN; 

	  double true_pe = sqrt(sq(true_mom_e[0]) +sq(true_mom_e[1]) + sq(true_mom_e[2]));
	  double true_nu = E1 - true_pe;
	  TVector3 true_qVec(-true_mom_e[0],-true_mom_e[1],E1-true_mom_e[2]);
	  double true_E_r = sqrt(sq(mN) + rVec_true.Mag2());
	  double true_alpha_s = (true_E_r - rVec_true.Dot(true_qVec)/true_qVec.Mag())/mN;

	  trueXpAsFine->Fill(trueXp,alpha_s);
	  asRes->Fill(true_alpha_s,alpha_s-true_alpha_s);

	  // Determine if there is bin migration
	  if (reconXpAsSame->FindBin(trueXp,true_alpha_s) == reconXpAsSame->FindBin(reconXp,alpha_s))
	    {
	      reconXpAsSame->Fill(reconXp,alpha_s);
	    }
	  else
	    {
	      reconXpAsDiff->Fill(reconXp,alpha_s);
	      trueXpAsDiff->Fill(trueXp,true_alpha_s);
	    }

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

  prRes->Write();
  asRes->Write();
  trueXpAsFine->Write();
  reconXpAsSame->Write();
  reconXpAsDiff->Write();
  trueXpAsDiff->Write();
  for (int aBin=0 ; aBin <5 ; aBin++)
    {
      xPrimeRes[aBin]->Write();
      WPrimeRes[aBin]->Write();
      cout << "There were " << nCounts[aBin] << " in bin " << aBin << ".\n";
    }

  outRoot->Close();  
  f->Close();
}
