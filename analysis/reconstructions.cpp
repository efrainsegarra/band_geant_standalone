// TOOK OUT THE SMEARING IN THE FILE FOR RIGHT NOW TO DO SOME BENCHMARK
// CHECKS

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

// Parameters
const double cScint = 15.; // in cm/ns
const double barCrossSection = 7.4;// in cm
const double BAND_Zoffset = -262; // in cm
const int BAND_numLayers = 5;
const int numBANDbars = 24; // bars in a single layer of the BAND array
const double BAND_Yoffset = barCrossSection * 5;

// Reconstruction options
const bool DoFirstHit = false;
const bool smearingOn = false;

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

		for(int j=0; j<trueEvent->hits.size(); j++){
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


			// get true time, location, and bar number from the event
			double trueT = trueEvent->hits[j].time; // in ns
			TVector3 truePos( trueEvent->hits[j].pos.x()/10.,trueEvent->hits[j].pos.y()/10.,trueEvent->hits[j].pos.z()/10. ); // in cm
			int trueBarNo = trueEvent->hits[j].barNo;


			// Looking at only first bars hit or front of detector hit
			if (DoFirstHit == true){
				if (abs(truePos[2]+262)>0.5) continue;
				//if (trueBarNo > 23) continue;
			}
			// DO RECONSTRUCTIONS:
			double reconT, reconX, reconY, reconZ;
			if (smearingOn == true){
				// time reconstruction based on PMT resolution gaussian
				reconT = myRand->Gaus(trueT,tResPMT/sqrt(2.));

				// x reconstruction based on PMT resolution for PMT time differences
				reconX = myRand->Gaus(truePos[0],cScint * tResPMT/sqrt(2.));

				// y reconstruction based on the bar that registered a hit, taken to be in
				// middle of the bar
					// first get the z-layer that the bar fired in and do y-reconstruction
				int recZLayer = floor(trueBarNo/numBANDbars);
					// take out the z-layer effect in the bar numbering
				int recYBarNo = trueBarNo - numBANDbars * recZLayer; 
					// do the Y reconstruction based on middle of the bar, but note that
					// bars 8-13 have the same y position as 2-7, so move the rest to level
				if (recYBarNo > 7) recYBarNo -= 6;
				reconY = (float(recYBarNo) + 0.5)* barCrossSection - BAND_Yoffset;

				// z smearing just assumed to be in middle of bar that registered hit
				reconZ = BAND_Zoffset - ( float(recZLayer) + 0.5)* barCrossSection;
			}
			else{
				//reconT = trueT;
				reconT = myRand->Gaus(trueT,tResPMT/sqrt(2.));
				reconX = truePos[0];
				reconY = truePos[1];
				reconZ = truePos[2];
			}
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
	      	reconHit.barNo = trueBarNo;

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