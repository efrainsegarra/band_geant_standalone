#include <iostream>
#include <cmath>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TVector3.h"
#include "TVectorT.h"

#include "prop_tree.h"
#include "constants.h"

using namespace std;

double sq(double x){ return x*x;};



int main(int argc, char ** argv)
{
  if (argc != 3){
      cerr << "Wrong number of arguments. Instead use\n"
	   << "\tmakeHists /path/to/reconstructions_fullGeo/file \n"
     << "\t /path/to/output/file\n";
      exit(-1);
  }

  TFile * inFile = new TFile(argv[1]);
  TFile * outFile = new TFile(argv[2],"RECREATE");
  TTree * inTree = (TTree*)inFile->Get("ReconTree");

  BAND_Event * reconEvent = NULL;
  inTree-> SetBranchAddress("band",&reconEvent);

  int firstBinVal = 100;
  int lastBinVal  = 600;
  int numBins     = 500;

  // Prep histograms
  TH1D * pr_fullGeo = new TH1D("pr_fullGeo","Pr Magnitude at face of BAND, full geometry",numBins,firstBinVal,lastBinVal);
  double pr,pr_peak;
  int numEvents;

  numEvents = inTree->GetEntries();
  for(int i =0; i<numEvents; i++){
    
    //if(i%100000==0)cerr << "Working on event " << i << "\n";
    inTree->GetEntry(i);

    //cout << "Event " << i << endl;

    for(int j=0; j<reconEvent->hits.size(); j++){
      
      pr = reconEvent->hits[j].mom.Mag();
      pr_fullGeo->Fill(pr);
    
    }
  }
  // This is the peak momentum value -- i.e. the momentum that we shot at the detector
  pr_peak = firstBinVal + (pr_fullGeo->GetMaximumBin())*(lastBinVal-firstBinVal)/numBins;

  double peak, tail;
  int maxBin, minBin;

  // Defines the peak as +/- 1.5% of the delta peak value
  maxBin = pr_fullGeo->GetXaxis()->FindBin(pr_peak*1.015);
  minBin = pr_fullGeo->GetXaxis()->FindBin(pr_peak*0.985);
  peak = pr_fullGeo->Integral(minBin,maxBin,"width");

  // Defines the tail as all events from 100 MeV/c to the peak
  maxBin = pr_fullGeo->GetXaxis()->FindBin(pr_peak*0.985)-1;
  minBin = pr_fullGeo->GetXaxis()->FindBin(100);
  tail = pr_fullGeo->Integral(minBin,maxBin,"width");
  
  cout << "Neutron Momentum [MeV/c]: " << pr_peak << "\nPeak-to-tail ratio: " << peak/tail << endl;

  pr_fullGeo->Write();

  outFile->Close();  
  inFile->Close();
}
