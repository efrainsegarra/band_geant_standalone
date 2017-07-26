#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorT.h"
#include "TVector3.h"

#include "prop_tree.h"
#include "gen_tree.h"
#include "constants.h"

using namespace std;

double getCLAS12_PRes(double theta, double p);
const double clasFac=1.;
// Parameters
const double cScint = 15.; // in cm/ns
const double barCrossSection = 7.4;// in cm
const double BAND_Zoffset = -262; // in cm
const int BAND_numLayers = 5;
const int numBANDbars = 24; // bars in a single layer of the BAND array
const double BAND_Yoffset = barCrossSection * 5;

// ********************************************************************************
// Reconstruction options
const bool smearingOn = true;
// ********************************************************************************

int main(int argc, char** argv){

	if (argc != 6)
    {
      cerr << "\nWrong number of arguments! Instead use:\n"
	   << "\treconstructions /path/to/inputGeant/root/file /path/to/inputGenerator/root/file\n"\
	   << "\t\t /path/to/output/root/file [threshold (MeV-ee)] [PMTtimeRes(ps)]\n";
      exit(-1);
    }
    const double threshold = atof(argv[4]); // in MeV-ee
    const double tResPMT = 1.E-3 * atof(argv[5]); // in ns

    // Set up input and output root files
    TFile * inFile_neutrons = new TFile(argv[1]);
    TTree * inTree_neutrons = (TTree*)inFile_neutrons->Get("PropTree");

    TFile * inFile_electrons = new TFile(argv[2]);
    TTree * inTree_electrons = (TTree*)inFile_electrons->Get("MCout");

  	TFile * outFile = new TFile(argv[3],"RECREATE");
  	TTree * outTree = new TTree("ReconTree","Reconstructed Values");
	 
	TRandom3 * myRand = new TRandom3(0);

	// Initialize the branches
	BAND_Event * trueEvent_neutrons = NULL;
	Gen_Event * trueEvent_electrons = NULL;
	BAND_Event * reconEvent = new BAND_Event;

	inTree_neutrons-> SetBranchAddress("band",&trueEvent_neutrons);
	inTree_electrons-> SetBranchAddress("event",&trueEvent_electrons);
  	outTree->Branch("band",&reconEvent);

  	// DIGITIZATION
	int numEvents_neutrons = inTree_neutrons->GetEntries();
	for(int i =0; i<numEvents_neutrons; i++){
		
		//if(i%100000==0)cerr << "Working on event " << i << "\n";
		inTree_neutrons->GetEntry(i);

		cout << "Event " << i << endl;
		
		// first find the earliest hit in the event that is above threshold
		double minTime = -1;
		int indexOfMinTime =-1;
		for(int j=0; j<trueEvent_neutrons->hits.size(); j++){

			// First check if the hit is above threshold
			// get energy deposit, covert to MeVee energy equivalent
			// deposit and then ask if that MeVee > some threshold set
			// in MeVee
			// this conversion from http://shop-pdp.net/efhtml/NIM_151_1978_445-450_Madey.pdf
			double trueE = trueEvent_neutrons->hits[j].E_dep; // in MeV
			double trueE_MeVee = 0.83 * trueE - 2.82 * ( 1 - exp( -0.25 * ( pow(trueE,0.93)) ) );
			if(abs(trueE_MeVee)<threshold) continue;

			// If hit above threshold, save the hit time to
			// find the earliest hit time
			double trueT = trueEvent_neutrons->hits[j].time; // in ns

			if ((trueT < minTime) || (indexOfMinTime < 0)){
				minTime = trueT;
				indexOfMinTime = j;
			}
		}
		if (indexOfMinTime >= 0) {
			int j = indexOfMinTime;
			double trueE = trueEvent_neutrons->hits[j].E_dep;
			double trueE_MeVee = 0.83 * trueE - 2.82 * ( 1 - exp( -0.25 * ( pow(trueE,0.93)) ) );

			// get true time, location, and bar number from the event
			double trueT = trueEvent_neutrons->hits[j].time; // in ns
			TVector3 truePos( trueEvent_neutrons->hits[j].pos.x()/10.,trueEvent_neutrons->hits[j].pos.y()/10.,trueEvent_neutrons->hits[j].pos.z()/10. ); // in cm
			int trueBarNo = trueEvent_neutrons->hits[j].barNo;
			

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
				reconT = trueT;
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

	      	std::string type = "neutron";
	      	// Now fill the outTree
	      	reconEvent->hits.clear();

	      	BAND_Hit reconHit;
	      	reconHit.E_dep = reconE; // in MeV
	      	reconHit.time = reconT;  // in ns
	      	reconHit.pos = reconPos; // in cm
	      	reconHit.mom = reconMom; // in MeV/c
	      	reconHit.barNo = trueBarNo;
	      	reconHit.type = type;
	      	reconEvent->hits.push_back(reconHit);


	      	// Now read the matching electron information from the generator tree:
	      	inTree_electrons->GetEntry(i);
	      		// for both the particle types in the generator tree:
	      	for (unsigned int p=0 ; p<trueEvent_electrons->particles.size() ; p++){
	      		std::string particleType = trueEvent_electrons->particles[p].type;
	      		if (particleType != "e-") continue;

	      		// Get the true momentum and true angles:
	      		// recall that the momentum in the generator is in GeV and in Geant it is MeV!!!
	      		double trueMomMag = trueEvent_electrons->particles[p].momentum.Mag() * 1000;
	      		double trueMomTheta = trueEvent_electrons->particles[p].momentum.Theta();
	      		double trueMomPhi = trueEvent_electrons->particles[p].momentum.Phi();

	      		double reconMomMag = myRand->Gaus(trueMomMag,clasFac*getCLAS12_PRes(trueMomTheta,trueMomMag));
	      		double reconMomTheta = myRand->Gaus(trueMomTheta,clasFac*0.001);
	      		double reconMomPhi = myRand->Gaus(trueMomPhi,clasFac*0.001/sin(trueMomTheta));

	      		if (reconMomPhi > M_PI) reconMomPhi -=2.*M_PI;
		    	if (reconMomPhi < -M_PI) reconMomPhi +=2.*M_PI;

				reconMomX = reconMomMag * sin(reconMomTheta) * cos(reconMomPhi);
				reconMomY = reconMomMag * sin(reconMomTheta) * sin(reconMomPhi);
				reconMomZ = reconMomMag * cos(reconMomTheta);
				TVector3 reconMom( reconMomX ,reconMomY ,reconMomZ ); // in MeV/c

	      		BAND_Hit reconHit;
	      		reconHit.type = particleType;
	      		reconHit.mom = reconMom;
	      		reconEvent->hits.push_back(reconHit);

	      	}
		
	      	outTree->Fill();
	    }
	  			
	}
	// ******************************************************************************************

		
	outTree->Write();
  	outFile->Close();
  	inFile_neutrons->Close();
  	inFile_electrons->Close();
  	delete myRand;

	return 0;
}

double getCLAS12_PRes(double theta, double p)
{
  const double a=0.0549141;
  const double b=-0.00332126;
  const double c=-0.052163;

  return p*(a+b*sin(theta)+c*cos(theta));  
}