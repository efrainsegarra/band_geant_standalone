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

// assumes that the detector is located in negative z-direction, upstream
// the target with the following parameters:


//
int main(int argc, char** argv){
	if (argc != 3)
    {
      cerr << "Wrong number of arguments! Instead use:\n"
	   << "\tavg_barFires /path/to/inputGeant/root/file threshold[MeV-ee]\n\n";
      exit(-1);
    }
	double threshold = atof(argv[2]);  
	double num_barFired = 0;
	double num_events_barFired = 0;

	TFile *inFile 	= new TFile(argv[1]);
	TTree *inTree   = (TTree*)inFile -> Get("PropTree");

	BAND_Event * event = NULL;
	inTree-> SetBranchAddress("band",&event);

	int numEvents = inTree->GetEntries();
	
	for(int i =0; i<numEvents; i++){
		
		double barFires = 0;
		inTree->GetEntry(i);

		for(int j=0; j<event->hits.size(); j++){
			// get energy deposit, covert to MeVee energy equivalent
			// deposit and then ask if that MeVee > some threshold set
			// in MeVee
			// this conversion from http://shop-pdp.net/efhtml/NIM_151_1978_445-450_Madey.pdf
			double hitE = event->hits[j].E_dep;
			double hitE_MeVee = 0.83 * hitE - 2.82 * ( 1 - exp( 0.25 * ( pow(hitE,0.93)) ) );
			if(abs(hitE_MeVee)<threshold) continue;
			
			barFires++;


		}
		// For that event, get the total number of bars fired in the event
		// and save that to a running total, and also check if that is > 0
		// to count that event
		
		num_barFired += barFires;
		if (barFires > 0) num_events_barFired++;

	}
	
	cout << threshold  << " " << (num_events_barFired)/numEvents << " " << num_barFired/num_events_barFired << endl;
	
	
	return 0;
}
