#include <iostream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TVector3.h"
#include "TVectorT.h"

#include "recon_tree.h"
#include "constants.h"

using namespace std;


double gauss(double *x, double *p);
double gaussAndBkg(double *x, double *p);
void fitMeanSig(TH1D * hist, double &mean, double &sig);

int main(int argc, char ** argv)
{
  if (argc != 4){
      cerr << "Wrong number of arguments. Instead use\n"
	   << "\tmakeHists /path/to/reconstructions_noGeo/file \n"
     << "\t /path/to/reconstructions_fullGeo/file /path/to/output/file\n";
      exit(-1);
  }

  TFile * inFile_nogeo = new TFile(argv[1]);
  TFile * inFile_geo = new TFile(argv[2]);
  TFile * outFile = new TFile(argv[3],"RECREATE");

  TTree * inTree_nogeo = (TTree*)inFile_nogeo->Get("ReconTree");
  TTree * inTree_geo = (TTree*)inFile_geo->Get("ReconTree");

  Recon_Event * reconEvent_nogeo = NULL;
  Recon_Event * reconEvent_geo = NULL;

  inTree_nogeo-> SetBranchAddress("event",&reconEvent_nogeo);
  inTree_geo-> SetBranchAddress("event",&reconEvent_geo);

  int firstBinVal = 100;
  int lastBinVal  = 600;
  int numBins     = 500;

  // Prep histograms
  TH1D * pr_nogeo = new TH1D("pr_nogeo","Pr Magnitude reconstruction\n from first hit above threshold, No Geo",numBins,firstBinVal,lastBinVal);
  TH1D * pr_geo = new TH1D("pr_fullgeo","Pr Magnitude reconstruction\n from first hit above threshold, Full Geo",numBins,firstBinVal,lastBinVal);

  double pr;
  int numEvents;

  numEvents = inTree_nogeo->GetEntries();
  for(int i =0; i<numEvents; i++){
    if (i % 10000 == 0) cout << "Event " << i << endl;
    inTree_nogeo->GetEntry(i);
    for(int j=0; j<reconEvent_nogeo->particles.size(); j++){
      pr = reconEvent_nogeo->particles[j].momRecon.Mag();
      pr_nogeo->Fill(pr);
    }
  }

  numEvents = inTree_geo->GetEntries();
  for(int i =0; i<numEvents; i++){
    if (i % 10000 == 0) cout << "Event " << i << endl;
    inTree_geo->GetEntry(i);
    for(int j=0; j<reconEvent_geo->particles.size(); j++){
      pr = reconEvent_geo->particles[j].momRecon.Mag();
      pr_geo->Fill(pr);
    }
  }

  // Now fit peak of nogeo to gaussian to get momentum interval
  double nogeo_sigma = 10.;
  double nogeo_mean;
  fitMeanSig(pr_nogeo,nogeo_mean,nogeo_sigma);

  double minBin, maxBin;

  double pr_peak_low = nogeo_mean - 2*nogeo_sigma;
  double pr_peak_high = nogeo_mean + 2*nogeo_sigma;

  // peak integration
  maxBin = pr_nogeo->GetXaxis()->FindBin(pr_peak_high);
  minBin = pr_nogeo->GetXaxis()->FindBin(pr_peak_low);
  double peak_nogeo = pr_nogeo->Integral(minBin,maxBin,"width");
  double peak_geo = pr_geo->Integral(minBin,maxBin,"width");

  // tail integration
  maxBin = pr_nogeo->GetXaxis()->FindBin(pr_peak_low)-1;
  minBin = pr_nogeo->GetXaxis()->FindBin(200);
  double tail_nogeo = pr_nogeo->Integral(minBin,maxBin,"width");
  double tail_geo = pr_geo->Integral(minBin,maxBin,"width");

  cout << peak_nogeo/tail_nogeo << endl;
  cout << peak_geo/tail_geo << endl;

  pr_nogeo->Write();
  outFile->Close();  
  inFile_geo->Close();
  inFile_nogeo->Close();
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

void fitMeanSig(TH1D * hist, double &mean, double &sig)
{ 
  double SigmaFromMean = 5;
  // fill the mean based on peak of the TDC spectrum
  mean = hist->GetBinCenter(hist->GetMaximumBin());
  
  // Create model of gaussian + background signal
  TF1 * model = new TF1("fitGausPeak",gauss,200,600,3); // (name,function,min,max,num of parameters)
  
  model->SetParameter(1,mean); // set the mean
  model->SetParameter(2,sig);  // set the initial guess width
  model->SetParameter(0,hist->GetMaximum()); // set amplitude of gaussian ontop of background

  hist->Fit("fitGausPeak","QES","",mean-SigmaFromMean*sig,mean+SigmaFromMean*sig);

  mean = model->GetParameter(1);
  sig  = fabs(model->GetParameter(2)); // update sigma that is fitted and output

  //delete model;
}
