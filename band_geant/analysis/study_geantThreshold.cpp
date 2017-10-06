#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <math.h>
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

int checkBAND(double theta_r,double phi_r){
  // ONLY SAVE THE EVENTS THAT WILL MAKE IT INTO BAND
  double numbars,prevbars;
  double max_x, min_x,max_x_2,min_x_2;
  double max_y = 7.4*13;
  double min_y = -7.4*5;

  double r = -262. / cos(theta_r) ;
  double x_traj = r * cos(phi_r) * sin(theta_r) ;
  double y_traj = r * sin(phi_r) * sin(theta_r) ;

  double r_2 = (-262.-7.4*5) / cos(theta_r) ;
  double x_traj_2 = r_2 * cos(phi_r) * sin(theta_r) ;
  double y_traj_2 = r_2 * sin(phi_r) * sin(theta_r) ;

  if ((y_traj > max_y) || (y_traj < min_y)) return -1;

  // first group
  prevbars = 0;
  numbars = 2;
  if ((y_traj > min_y + prevbars*7.4) && (y_traj < min_y + (numbars+prevbars)*7.4)){
    max_x = 2018.35 / 10 / 2;
    min_x = -max_x;
    if ((x_traj > max_x) || (x_traj < min_x)) return -1;
  }
  // second group
  prevbars += numbars;
  numbars = 6;
  if ((y_traj > min_y + prevbars*7.4) && (y_traj < min_y + (numbars+prevbars)*7.4)){
    max_x = (509.17/10 + 500/10);
    min_x_2 = -(509.17/10 + 500/10);
    if ((x_traj > max_x) || (x_traj < min_x_2)) return -1;
  }
  if ((y_traj_2 > min_y + prevbars*7.4) && (y_traj_2 < min_y + (numbars+prevbars)*7.4)){
    min_x = (509.17/10 );
    max_x_2 = -(509.17/10 );
    if (( x_traj_2 < min_x ) && (x_traj_2 > max_x_2) ) return -1;
  }

  // third group
  prevbars += numbars;
  numbars = 4;
  if ((y_traj > min_y + prevbars*7.4) && (y_traj < min_y + (numbars+prevbars)*7.4)){
    max_x = 2018.35 / 10 / 2;
    min_x = -max_x;
    if ((x_traj > max_x) || (x_traj < min_x)) return -1;
  }
  // fourth group
  prevbars += numbars;
  numbars = 3;
  if ((y_traj > min_y + prevbars*7.4) && (y_traj < min_y + (numbars+prevbars)*7.4)){
    max_x = 1946.53 / 10 / 2;
    min_x = -max_x;
    if ((x_traj > max_x) || (x_traj < min_x)) return -1;
  }
  // five group
  prevbars += numbars;
  numbars = 3;
  if ((y_traj > min_y + prevbars*7.4) && (y_traj < min_y + (numbars+prevbars)*7.4)){
    max_x = 1634.58 / 10 / 2;
    min_x = -max_x;
    if ((x_traj > max_x) || (x_traj < min_x)) return -1;
  }
  return 1;
}



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
			double hitE_MeVee = 0.83 * hitE - 2.82 * ( 1 - exp( - 0.25 * ( pow(hitE,0.93)) ) );
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