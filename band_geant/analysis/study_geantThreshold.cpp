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

  TFile * outFile   = new TFile("MeVee_spectrum.root","RECREATE");
  TH1D * spectrum = new TH1D("MeVee_spectrum","Energy Deposited in BAND for all Neutron Hits in Single Event",2000,0,20);
  spectrum->GetXaxis()->SetTitle("Energy Deposited [MeVee]");
  spectrum->GetYaxis()->SetTitle("Entries");
  spectrum->SetStats(0);

	double threshold = atof(argv[2]);  
	BAND_Event * event = NULL;
	inTree-> SetBranchAddress("band",&event);

	int events_generated = inTree->GetEntries(); // total number of particles shot in given phase space
  int tot_particles = 0; // number of particles that pass through BAND or deposit energy in BAND
  int int_particles = 0; // number of interacting particles in BAND
  int it = 0;
  bool counted;

  for( int i = 0; i < events_generated ; i++){
    //if (i%10000 == 0) cout << "Working on event: "<< i <<"\n";
    inTree->GetEntry(i);

    if ( event->hits.size() > 0){
      counted = false;
      // Then we had a particle either pass through BAND or deposit energy in BAND
      tot_particles++;


      // A single event may have many hits in different bars -- look for first hit
      // that is above threshold, and call that an interacting particle
      for( int j = 0; j < event->hits.size() ; j++){

        double en = event->hits[j].E_dep; // in MeV
        double en_MeVee = 0.83 * en - 2.82 * ( 1 - exp( - 0.25 * ( pow(en,0.93)) ) );

        if (en_MeVee > 0){
          spectrum->Fill(en_MeVee);
          it++;
        } 

        if ((en_MeVee > threshold) && (counted==false)){
          int_particles++;
          counted = true;
          break;
        }

      }
      
    }
  }

  cout << events_generated << " " << tot_particles << " " << int_particles << " " << threshold << "\n";
  spectrum->Scale(1./it);
  spectrum->Write();
  outFile->Close();

	return 0;
}
