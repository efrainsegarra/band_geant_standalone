#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"
#include "TF1.h"
#include "TRandom.h"
#include "TH2.h"
#include "prop_tree.h"

int main(int argc, char** argv){
	if (argc != 3)
    {
      cerr << "Wrong number of arguments! Instead use:\n"
	   << "\tstudy_geantThreshold /path/to/inputGeant/root/file threshold[MeV-ee]\n\n";
      exit(-1);
    }

  TFile *inFile   = new TFile(argv[1]);
  TTree *inTree   = (TTree*)inFile -> Get("PropTree");

  TFile * outFile   = new TFile("tomography_250.root","RECREATE");
  TH2D * tomography = new TH2D("tomography_250","",240,-120,120,160,110,-50);
  tomography->GetXaxis()->SetTitle("x position [cm]");
  tomography->GetYaxis()->SetTitle("y position [cm]");
  tomography->SetStats(0);

	double threshold = atof(argv[2]);  
	BAND_Event * event = NULL;
	inTree-> SetBranchAddress("band",&event);

	int events_generated = inTree->GetEntries(); // total number of particles shot in given phase space

  for( int i = 0; i < events_generated ; i++){
    if (i%10000 == 0) cout << "Working on event: "<< i <<"\n";
    inTree->GetEntry(i);

    //for ( int j = 0; j < event->hits.size() ; j++){
      
      //double en_MeVee = 0.83 * en - 2.82 * ( 1 - exp( - 0.25 * ( pow(en,0.93)) ) );

      double x = event->hits[0].pos.X()/10.;
      double y = event->hits[0].pos.Y()/10.;
      double z = event->hits[0].pos.Z()/10.;
      if ((abs(z + 262)) < 0.2) tomography->Fill(x,y);
      
    //}
  }
  tomography->Write();
  outFile->Close();
  inFile->Close();
  

	return 0;
}
