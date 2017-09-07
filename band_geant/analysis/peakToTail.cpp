#include <iostream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TVector3.h"
#include "TVectorT.h"

#include "recon_tree.h"
#include "constants.h"

using namespace std;

double sq(double x){ return x*x;};



int main(int argc, char ** argv)
{
  if (argc != 4){
      cerr << "Wrong number of arguments. Instead use\n"
	   << "\tmakeHists /path/to/reconstructions_fullGeo/file \n"
     << "\t /path/to/reconstructions_noGeo/file /path/to/output/file\n";
      exit(-1);
  }

  TFile * inFile_Geo = new TFile(argv[1]);
  TFile * inFile_NoGeo = new TFile(argv[2]);
  TFile * outFile = new TFile(argv[3],"RECREATE");

  TTree * inTree_Geo = (TTree*)inFile_Geo->Get("ReconTree");
  TTree * inTree_NoGeo = (TTree*)inFile_NoGeo->Get("ReconTree");

  Recon_Event * reconEvent_noGeo = NULL;
  Recon_Event * reconEvent_Geo = NULL;

  inTree_NoGeo-> SetBranchAddress("event",&reconEvent_noGeo);
  inTree_Geo-> SetBranchAddress("event",&reconEvent_Geo);

  int firstBinVal = 100;
  int lastBinVal  = 600;
  int numBins     = 500;

  // Prep histograms
  TH1D * pr_fullGeo = new TH1D("pr_fullGeo","Pr Magnitude reconstruction\n from first hit above threshold",numBins,firstBinVal,lastBinVal);
  TH1D * pr_noGeo = new TH1D("pr_noGeo","Pr Magnitude reconstruction\n from first hit above threshold",numBins,firstBinVal,lastBinVal);
  double pr,pr_peak;
  int numEvents;

  numEvents = inTree_NoGeo->GetEntries();
  for(int i =0; i<numEvents; i++){
    
    //if(i%100000==0)cerr << "Working on event " << i << "\n";
    inTree_NoGeo->GetEntry(i);

    //cout << "Event " << i << endl;

    for(int j=0; j<reconEvent_noGeo->particles.size(); j++){
      
      pr = reconEvent_noGeo->particles[j].momRecon.Mag();
      pr_noGeo->Fill(pr);
    
    }
  }
  // This is the peak momentum value -- i.e. the momentum that we shot at the detector
  //pr_peak = firstBinVal + (pr_fullGeo->GetMaximumBin())*(lastBinVal-firstBinVal)/numBins;

  //double peak, tail;
  //int maxBin, minBin;

  // Defines the peak as +/- 1.5% of the delta peak value
  //maxBin = pr_fullGeo->GetXaxis()->FindBin(pr_peak*1.015);
  //minBin = pr_fullGeo->GetXaxis()->FindBin(pr_peak*0.985);
  //peak = pr_fullGeo->Integral(minBin,maxBin,"width");

  // Defines the tail as all events from 100 MeV/c to the peak
  //maxBin = pr_fullGeo->GetXaxis()->FindBin(pr_peak*0.985)-1;
  //minBin = pr_fullGeo->GetXaxis()->FindBin(100);
  //tail = pr_fullGeo->Integral(minBin,maxBin,"width");
  
  //cout << "Neutron Momentum [MeV/c]: " << pr_peak << "\nPeak-to-tail ratio: " << peak/tail << endl;

  pr_noGeo->Write();

  outFile->Close();  
  inFile_Geo->Close();
  inFile_NoGeo->Close();
}
