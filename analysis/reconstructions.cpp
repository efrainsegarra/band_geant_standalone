#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorT.h"
#include "TVector3.h"

#include "prop_tree.h"
#include "constants.h"

using namespace std;

const double cScint = 15.; // in cm/ns
const double barCrossSection = 7.4;// in cm
const double BAND_offset = -262; // in cm
const int BAND_numLayers = 5;

int main(int argc, char** argv){

	if (argc != 5)
    {
      cerr << "\nWrong number of arguments! Instead use:\n"
	   << "\treconstructions /path/to/inputGeant/root/file /path/to/output/root/file\n"
	   << "\t\t[threshold (MeV-ee)] [PMTtimeRes(ps)]\n\n";
      exit(-1);
    }
    const double threshold = atof(argv[3]); // in MeV-ee
    const double tResPMT = 1.E-3 * atof(argv[4]); // in ns

    // Set up input and output root files
    TFile * inFile = new TFile(argv[1]);
  	TFile * outFile = new TFile(argv[2],"RECREATE");
  	TTree * inTree = (TTree*)inFile->Get("PropTree");
  	TTree * outTree = new TTree("ReconTree","Reconstructed Values");
	 
	TRandom3 * myRand = new TRandom3(0);

	// Initialize the branches
	BAND_Event * trueEvent = NULL;
	BAND_Event * reconEvent = new BAND_Event;

	inTree-> SetBranchAddress("band",&trueEvent);
  	outTree->Branch("band",&reconEvent);

	// Loop over number of events and all the hits in the event
	int numEvents = inTree->GetEntries();
	
	for(int i =0; i<numEvents; i++){
		
		if(i%100000==0)cerr << "Working on event " << i << "\n";
		inTree->GetEntry(i);

		cout << "Event " << i << endl;
		cout << trueEvent->hits.size() << endl;

		for(int j=0; j<trueEvent->hits.size(); j++){
			cout << "   bar: " << trueEvent->hits[j].barNo << endl;
			//if(abs(trueEvent->hits[j].pos.z()/10. - BAND_offset + barCrossSection/2. ) > barCrossSection/2.) continue;

			// get energy deposit, covert to MeVee energy equivalent
			// deposit and then ask if that MeVee > some threshold set
			// in MeVee
			// this conversion from http://shop-pdp.net/efhtml/NIM_151_1978_445-450_Madey.pdf
			
			double trueE = trueEvent->hits[j].E_dep; // in MeV
			double trueE_MeVee = 0.83 * trueE - 2.82 * ( 1 - exp( -0.25 * ( pow(trueE,0.93)) ) );
			if(abs(trueE_MeVee)<threshold) continue;

			// A single neutron can hit many bars, but for each bar, there is
			// one 'hit'. If it did 'hit' more than once in a bar, we already
			// energy-averaged those hits in geant. So we now have multiple 
			// TOF & position measurements for one neutron if it his many bars.

			// get all the true values from the hit; energy is above
			double trueT = trueEvent->hits[j].time; // in ns
			int trueBarNo = trueEvent->hits[j].barNo;
			TVector3 truePos( trueEvent->hits[j].pos.x()/10.,trueEvent->hits[j].pos.y()/10.,trueEvent->hits[j].pos.z()/10. ); // in cm

			// reconstruct the true values based on resolutions

			// time smearing based on PMT resolution
			double reconT = myRand->Gaus(trueT,tResPMT/sqrt(2.));

			// x smearing based on PMT resolution for time differences
			double reconX = myRand->Gaus(truePos[0],cScint * tResPMT/sqrt(2.));

			// y smearing just assumed to be in middle of bar that registered hit
			int recBarNo = floor(truePos[1] / barCrossSection);
			double reconY = ( float(recBarNo) + 0.5)* barCrossSection;

			// z smearing just assumed to be in middle of bar that registered hit
			int recZLayer = floor( (BAND_offset - truePos[2]) / barCrossSection );
			double reconZ = BAND_offset - ( float(recZLayer) + 0.5)* barCrossSection;

			TVector3 reconPos( reconX ,reconY ,reconZ ); // in cm

			// now reconstruct momentum & energy from TOF and path length
			double reconPath = sqrt( pow(reconX,2) + pow(reconY,2) + pow(reconZ,2) );
			double reconBeta = reconPath/(cAir*reconT);
			double reconMomMag = mN / sqrt( pow(1./reconBeta,2) - 1.) * 1.E3; // in MeV/c
			double reconE = sqrt( pow(reconMomMag,2) + pow(mN * 1.E3 ,2) ) ; // in MeV

			// reconstruct momentum vector
	      	double reconCosTheta = reconZ/reconPath;
	      	double reconTheta = acos(reconCosTheta);
	      	double reconPhi = atan2(reconY,reconX);
	      	double reconMomX = reconMomMag * sin(reconTheta) * cos(reconPhi);
	      	double reconMomY = reconMomMag * sin(reconTheta) * sin(reconPhi);
	      	double reconMomZ = reconMomMag * reconCosTheta;  

	      	TVector3 reconMom( reconMomX ,reconMomY ,reconMomZ ); // in MeV/c

	      	// Now fill the outTree
	      	reconEvent->hits.clear();

	      	BAND_Hit reconHit;
	      	reconHit.E_dep = reconE; // in MeV
	      	reconHit.time = reconT;  // in ns
	      	reconHit.pos = reconPos; // in cm
	      	reconHit.mom = reconMom; // in MeV/c

	      	reconEvent->hits.push_back(reconHit);

	      	outTree->Fill();

		}
	}
		
	outTree->Write();
  	outFile->Close();
  	inFile->Close();
  	delete myRand;

	return 0;
}