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
#include "gen_tree.h"
#include "TLatex.h"

// Units of GeV
const double mN = 0.93956536;
// Units of cm / ns
const double cAir = 29.9792458;


int main(int argc, char** argv){
	if (argc != 4)
    {
      cerr << "Wrong number of arguments! Instead use:\n"
	   << "\tstudy_geometry /path/to/inputGeant/root/file /path/to/generator/root/file threshold[MeV-ee]\n\n";
      exit(-1);
    }

  TFile *inFileGeant   = new TFile(argv[1]);
  TTree *inTreeGeant   = (TTree*)inFileGeant -> Get("PropTree");

  TFile *inFileGen   = new TFile(argv[2]);
  TTree *inTreeGen   = (TTree*)inFileGen -> Get("MCout");

  TFile * outFile   = new TFile("geoEffects.root","RECREATE");
  TH2D * thetaN_pN = new TH2D("thetaN_pN","Phase Space of Neutrons that Hit BAND",1800,0,180,1800*2,-180,180);
  thetaN_pN->GetXaxis()->SetTitle("#theta_{n}");
  thetaN_pN->GetYaxis()->SetTitle("#phi_{n}");

  TH2D * thetaN_pN_gen = new TH2D("thetaN_pN_gen","Phase Space of Generated Neutrons",1800,0,180,1800*2,-180,180);
  thetaN_pN_gen->GetXaxis()->SetTitle("#theta_{n}");
  thetaN_pN_gen->GetYaxis()->SetTitle("#phi_{n}");


	double threshold = atof(argv[3]);  
	BAND_Event * geantevent = NULL;
	inTreeGeant-> SetBranchAddress("band",&geantevent);

  Gen_Event * generatedEvent = NULL;
  inTreeGen-> SetBranchAddress("event",&generatedEvent);



  for( int i = 0; i < inTreeGeant->GetEntries() ; i++){
    if ( i % 100000 == 0) cout << "Working on entry " << i << "\n";

    inTreeGen->GetEntry(i);
    double gen_phi = generatedEvent->particles[1].momentum.Phi() * 180./M_PI; // in MeV/c
    double gen_theta = generatedEvent->particles[1].momentum.Theta() * 180/M_PI; // in degrees
    thetaN_pN_gen->Fill(gen_theta,gen_phi);


    inTreeGeant->GetEntry(i);

    if( geantevent->hits.size() <= 0) continue;

    for( int j = 0; j<geantevent->hits.size(); j++){
      double en = geantevent->hits[j].E_dep;
      double MeVee = 0.83 * en - 2.82 * ( 1 - exp( - 0.25 * ( pow(en,0.93)) ) );

      //if(abs(MeVee)<threshold) continue;

      double time = geantevent->hits[j].time;
      TVector3 pos( geantevent->hits[j].pos.x()/10.,geantevent->hits[j].pos.y()/10.,geantevent->hits[j].pos.z()/10. ); // in cm
      double x = pos[0];
      double y = pos[1];
      double z = pos[2];

      double path = sqrt( pow(x,2) + pow(y,2) + pow(z,2) );
      double beta = path/(cAir*time);
      double momMag = mN / sqrt( pow(1./beta,2) - 1.) * 1.E3; // in MeV/c

      double CosTheta = z/path;
      double theta = acos(CosTheta) * 180./M_PI;
      double phi = atan2(y,x) * 180./M_PI;

      // Only fill for one hit -- right now the first hit, although option to take
      // first hit above threshold
      thetaN_pN->Fill(theta,phi);
      //cout << gen_theta << " " << gen_momMag << " " << theta << " " << momMag << "\n";
      break;
    }
    
  }


  
  thetaN_pN->Write();
  thetaN_pN_gen->Write();
  outFile->Close();

	return 0;
}
