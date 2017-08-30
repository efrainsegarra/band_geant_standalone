#include <iostream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TVector3.h"
#include "TF1.h"

#include "constants.h"

using namespace std;

double sq(double x){ return x*x;};
void doProj( TH2D * hist,int bin, int itBins,double &sig, bool write);
double gauss(double *x, double *p);
double gaussAndBkg(double *x, double *p);
void fitMeanSig(TH1D * hist, double &mean, double &sig);


int main(int argc, char ** argv)
{
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead use\n"
	   << "\t/path/to/kinematic/file /path/to/hists/file\n";
      exit(-1);
    }

  TFile * inFile = new TFile(argv[1]);
  TTree * inTree = (TTree*)inFile->Get("ResTree");

  // Memory for the output tree
  double trueWp, trueXp, trueAs;
  double reconWp, reconXp, reconAs;
  double truePe, truePn, reconPe, reconPn;
  double truereconMom_e[3],trueMom_n[3];
  double reconreconMom_e[3],reconMom_n[3];

  double trueQSq, trueNu, trueX, true_q, truePhi_q, trueTheta_q;
  double reconQSq, reconNu, reconX, recon_q, reconPhi_q, reconTheta_q;
  double truePhi_e, truePhi_n, truePhi_en;
  double reconPhi_e, reconPhi_n, reconPhi_en;
  double trueTheta_n, reconTheta_n;
  double trueTheta_e, reconTheta_e;
  double trueEn, reconEn;
  double trueCosTheta_qn, reconCosTheta_qn;

  inTree->SetBranchAddress("trueWp",&trueWp);
  inTree->SetBranchAddress("trueXp",&trueXp);
  inTree->SetBranchAddress("trueAs",&trueAs);

  inTree->SetBranchAddress("reconWp",&reconWp);
  inTree->SetBranchAddress("reconXp",&reconXp);
  inTree->SetBranchAddress("reconAs",&reconAs);

  inTree->SetBranchAddress("truePe",&truePe);
  inTree->SetBranchAddress("truePn",&truePn);
  inTree->SetBranchAddress("reconPe",&reconPe);
  inTree->SetBranchAddress("reconPn",&reconPn);

  inTree->SetBranchAddress("truePe_x",&(truereconMom_e[0]));
  inTree->SetBranchAddress("truePe_y",&(truereconMom_e[1]));
  inTree->SetBranchAddress("truePe_z",&(truereconMom_e[2]));

  inTree->SetBranchAddress("truePn_x",&(truereconMom_e[0]));
  inTree->SetBranchAddress("truePn_y",&(truereconMom_e[1]));
  inTree->SetBranchAddress("truePn_z",&(truereconMom_e[2]));

  inTree->SetBranchAddress("reconPe_x",&(reconreconMom_e[0]));
  inTree->SetBranchAddress("reconPe_y",&(reconreconMom_e[1]));
  inTree->SetBranchAddress("reconPe_z",&(reconreconMom_e[2]));

  inTree->SetBranchAddress("reconPn_x",&(reconMom_n[0]));
  inTree->SetBranchAddress("reconPn_y",&(reconMom_n[1]));
  inTree->SetBranchAddress("reconPn_z",&(reconMom_n[2]));

  inTree->SetBranchAddress("trueQSq",&(trueQSq));
  inTree->SetBranchAddress("trueNu",&(trueNu));
  inTree->SetBranchAddress("trueX",&(trueX));
  inTree->SetBranchAddress("true_q",&(true_q));
  inTree->SetBranchAddress("truePhi_q",&(truePhi_q));
  inTree->SetBranchAddress("trueTheta_q",&(trueTheta_q));

  inTree->SetBranchAddress("reconQSq",&(reconQSq));
  inTree->SetBranchAddress("reconNu",&(reconNu));
  inTree->SetBranchAddress("reconX",&(reconX));
  inTree->SetBranchAddress("recon_q",&(recon_q));
  inTree->SetBranchAddress("reconPhi_q",&(reconPhi_q));
  inTree->SetBranchAddress("reconTheta_q",&(reconTheta_q));

  inTree->SetBranchAddress("truePhi_e",&(truePhi_e));
  inTree->SetBranchAddress("truePhi_n",&(truePhi_n));
  inTree->SetBranchAddress("truePhi_en",&(truePhi_en));

  inTree->SetBranchAddress("reconPhi_e",&(reconPhi_e));
  inTree->SetBranchAddress("reconPhi_n",&(reconPhi_n));
  inTree->SetBranchAddress("reconPhi_en",&(reconPhi_en));

  inTree->SetBranchAddress("trueTheta_n",&(trueTheta_n));
  inTree->SetBranchAddress("reconTheta_n",&(reconTheta_n));
  inTree->SetBranchAddress("trueTheta_e",&(trueTheta_e));
  inTree->SetBranchAddress("reconTheta_e",&(reconTheta_e));

  inTree->SetBranchAddress("trueEn",&(trueEn));
  inTree->SetBranchAddress("reconEn",&(reconEn));

  inTree->SetBranchAddress("trueCosTheta_qn",&(trueCosTheta_qn));
  inTree->SetBranchAddress("reconCosTheta_qn",&(reconCosTheta_qn));

  TFile *outFile = new TFile(argv[2],"RECREATE");
  //cout << "\tLoaded input root file..." << "\n";

  // Prep resolution 2D histograms
  int Pn_nBin = 150;
  int As_nBin = 200;
  int Wp_nBin = 200;
  int Xp_nBin = 250;
  TH2D * PnRes = new TH2D("PnRes","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;Pn;delta Pn / Pn;Counts",Pn_nBin,0.25,0.55,100,-0.2,0.2);
  TH2D * AsRes = new TH2D("AsRes","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;As;delta As / As;Counts",As_nBin,1.25,1.65,100,-0.1,0.1);
  TH2D * WpRes = new TH2D("WpRes","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;Wp;delta Wp / Wp;Counts",Wp_nBin,1.5,3.5,100,-0.2,0.2);
  TH2D * XpRes = new TH2D("XpRes","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;Xp;delta Xp / Xp;Counts",Xp_nBin,0.2,0.7,100,-0.2,0.2);

  // Loop over events to fill histograms
  const int nEvents = inTree->GetEntries();
  //cout << "Filling 2D histograms..." << endl;
  for (int i=0 ; i < nEvents ; i++){
    //if (i % 10000 == 0) cout << "Event " << i << "\n";
    inTree->GetEvent(i);

    // Fill histograms
	  PnRes->Fill(truePn,(reconPn - truePn)/truePn);
	  AsRes->Fill(trueAs,(reconAs-trueAs)/trueAs);
    WpRes->Fill(trueWp,(reconWp-trueWp)/trueWp);
    XpRes->Fill(trueXp,(reconXp-trueXp)/trueXp);
  }

  // Now loop over the histograms and get projections to fit resolution
  //cout << "Looping over a second time to fit 2D histgrams..." << "\n";

  cout << "\tPn Resolution: " << "\n";
  int numBins = 10;
  int PnStartBin = PnRes->GetXaxis()->FindBin(0.25);
  int PnEndBin = PnRes->GetXaxis()->FindBin(0.5);
  int itBins = (PnEndBin-PnStartBin)/numBins;
  for (int bin=PnStartBin ; bin<= PnEndBin; bin += itBins){
    bool write = false;
    double PnSigGuess = 0.02;
    doProj( PnRes, bin,  itBins,PnSigGuess,write);
  }

  cout << "\tAs Resolution: " << "\n";
  numBins = 10;
  int AsStartBin = AsRes->GetXaxis()->FindBin(1.3);
  int AsEndBin = AsRes->GetXaxis()->FindBin(1.5);
  itBins = (AsEndBin-AsStartBin)/numBins;
  for (int bin=AsStartBin ; bin<= AsEndBin; bin += itBins){
    bool write = false;
    double AsSigGuess = 0.01;
    doProj( AsRes, bin,  itBins,AsSigGuess,write);
  }

  cout << "\tWp Resolution: " << "\n";
  numBins = 10;
  int WpStartBin = WpRes->GetXaxis()->FindBin(1.8);
  int WpEndBin = WpRes->GetXaxis()->FindBin(3.);
  itBins = (WpEndBin-WpStartBin)/numBins;
  for (int bin=WpStartBin ; bin<= WpEndBin; bin += itBins){
    bool write = false;
    double WpSigGuess = 0.03;
    doProj( WpRes, bin,  itBins,WpSigGuess,write);
  }

  cout << "\tXp Resolution: " << "\n";
  numBins = 10;
  int XpStartBin = XpRes->GetXaxis()->FindBin(0.2);
  int XpEndBin = XpRes->GetXaxis()->FindBin(0.6);
  itBins = (XpEndBin-XpStartBin)/numBins;
  for (int bin=XpStartBin ; bin<= XpEndBin; bin += itBins){
    bool write = false;
    double XpSigGuess = 0.02;
    doProj( XpRes, bin,  itBins,XpSigGuess,write);
  }
  

  PnRes->Write();
  AsRes->Write();
  WpRes->Write();
  XpRes->Write();

  outFile->Close();  
  inFile->Close();
}










double gauss(double *x, double *p)
{
  double var = *x;

  double amp = fabs(p[0]);
  double mean = p[1];
  double sig = p[2];

  double temp = (var - mean)/sig;
  return amp*exp(-0.5*temp*temp);
}

double gaussAndBkg(double *x, double *p)
{
  double bkg = fabs(p[3]);
  return bkg + gauss(x,p);
}

void fitMeanSig(TH1D * hist, double &mean, double &sig)
{
  // fill the mean based on peak of the TDC spectrum
  mean = hist->GetBinCenter(hist->GetMaximumBin());
  
  // Create model of gaussian + background signal
  TF1 * model = new TF1("fitRes",gauss,-1.,1.,3); // (name,function,min,max,num of parameters)
  
  model->SetParameter(1,mean); // set the mean
  model->SetParameter(2,sig);  // set the initial guess width
  model->SetParameter(0,hist->GetMaximum()); // set amplitude of gaussian ontop of background

  hist->Fit("fitRes","QES","",mean-3*sig,mean+3*sig);

  mean = model->GetParameter(1);
  sig  = fabs(model->GetParameter(2)); // update sigma that is fitted and output

  //delete model;
}

void doProj( TH2D * hist,int bin, int itBins, double &sig, bool write){
  
  char temp[100];
  sprintf(temp,"slice_%d",bin);

  TH1D * pj = hist->ProjectionY(temp,bin,bin+itBins-1);
  double mean;
  fitMeanSig(pj,mean,sig);
  double err =  pj->GetFunction("fitRes")->GetParError(2);
  double binCenter = 0.5*(hist->GetXaxis()->GetBinCenter(bin) + hist->GetXaxis()->GetBinCenter(bin+itBins-1));
  if(write) pj->Write();
  cout << binCenter << " " << sig*100 << " " << err*100 << "\n";

}
