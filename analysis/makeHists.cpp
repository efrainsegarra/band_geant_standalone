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
  if (argc != 7){
      cerr << "Wrong number of arguments. Instead use\n"
	   << "\tmakeHists /path/to/reconstructions_noGeo/file \n"
     << "\t /path/to/reconstructions_fullGeo/file \n"
     << "\t /path/to/output/file firstBinVal secondBinVal numBins\n";
      exit(-1);
  }

  TFile * inFile1 = new TFile(argv[1]);
  TFile * inFile2 = new TFile(argv[2]);
  TFile * outFile = new TFile(argv[3],"RECREATE");
  TTree * inTree1 = (TTree*)inFile1->Get("ReconTree");
  TTree * inTree2 = (TTree*)inFile2->Get("ReconTree");

  BAND_Event * reconEvent_noGeo = NULL;
  BAND_Event * reconEvent_fullGeo = NULL;
  inTree1-> SetBranchAddress("band",&reconEvent_noGeo);
  inTree2-> SetBranchAddress("band",&reconEvent_fullGeo);

  double firstBinVal = atof(argv[4]);
  double lastBinVal  = atof(argv[5]);
  double numBins     = atof(argv[6]);

  // Prep histograms
  TH1D * pr_noGeo = new TH1D("pr_noGeo","Pr Magnitude at face of BAND, no geometry",numBins,firstBinVal,lastBinVal); 
  TH1D * pr_fullGeo = new TH1D("pr_fullGeo","Pr Magnitude at face of BAND, full geometry",numBins,firstBinVal-numBins/2,lastBinVal-numBins/2);
  double pr,pr_peak;
  int numEvents;

  // Loop over events
  numEvents = inTree1->GetEntries();
  for(int i =0; i<numEvents; i++){
    
    //if(i%100000==0)cerr << "Working on event " << i << "\n";
    inTree1->GetEntry(i);

    //cout << "Event " << i << endl;

    for(int j=0; j<reconEvent_noGeo->hits.size(); j++){
      
      pr = reconEvent_noGeo->hits[j].mom.Mag();
      pr_noGeo->Fill(pr);
    
    }
  }
  pr_peak = firstBinVal+(pr_noGeo->GetMaximumBin())*(lastBinVal-firstBinVal)/numBins;

  numEvents = inTree2->GetEntries();
  for(int i =0; i<numEvents; i++){
    
    //if(i%100000==0)cerr << "Working on event " << i << "\n";
    inTree2->GetEntry(i);

    //cout << "Event " << i << endl;

    for(int j=0; j<reconEvent_fullGeo->hits.size(); j++){
      
      pr = reconEvent_fullGeo->hits[j].mom.Mag();

      if (abs(pr - pr_peak) < 1.) continue;

      pr_fullGeo->Fill(pr);
    
    }
  }

  double peak, tail;
  peak = pr_noGeo->GetMaximum();
  tail = pr_fullGeo->GetEntries();
  cout << "Neutron Momentum [MeV/c]: " << pr_peak << "\nPeak-to-tail ratio: " << peak/tail << endl;

  pr_fullGeo->Write();
  pr_noGeo->Write();

  outFile->Close();  
  inFile1->Close();
  inFile2->Close();
}
