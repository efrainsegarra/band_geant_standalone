#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorT.h"
#include "TVector3.h"

#include "recon_tree.h"
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
const int BAND_numLayers = 1;
const int numBANDbars = 1; // bars in a single layer of the BAND array
const double BAND_Yoffset = barCrossSection * BAND_numLayers;

// ********************************************************************************
// Reconstruction options
const bool smearingOn = false;
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
    TFile * inFile_geant = new TFile(argv[1]);
    TTree * inTree_geant = (TTree*)inFile_geant->Get("PropTree");

    TFile * inFile_generator = new TFile(argv[2]);
    TTree * inTree_generator = (TTree*)inFile_generator->Get("MCout");

  	TFile * outFile = new TFile(argv[3],"RECREATE");
  	TTree * outTree = new TTree("ReconTree","Reconstructed Values");
	 
	TRandom3 * myRand = new TRandom3(0);

	// Initialize the branches
	BAND_Event * geantEvent = NULL;
	Gen_Event * generatedEvent = NULL;
	Recon_Event * reconEvent = new Recon_Event;

	inTree_geant-> SetBranchAddress("band",&geantEvent);
	inTree_generator-> SetBranchAddress("event",&generatedEvent);
  	outTree->Branch("event",&reconEvent);

  	// DIGITIZATION
	int numEvents_geant = inTree_geant->GetEntries();
	for(int i =0; i<numEvents_geant; i++){
		
		//if(i%100000==0)cerr << "Working on event " << i << "\n";
		inTree_geant->GetEntry(i);
		if (i % 10000 == 0) cout << "Event " << i << endl;
		
		// first find the earliest hit in the event that is above threshold
		double minTime = -1;
		int indexOfMinTime =-1;
		for(int j=0; j<geantEvent->hits.size(); j++){

			// First check if the hit is above threshold
			// get energy deposit, covert to MeVee energy equivalent
			// deposit and then ask if that MeVee > some threshold set
			// in MeVee
			// this conversion from http://shop-pdp.net/efhtml/NIM_151_1978_445-450_Madey.pdf
			double trueE = geantEvent->hits[j].E_dep; // in MeV
			double trueE_MeVee = 0.95 * trueE - 8.0 * ( 1 - exp( 0.10 * ( pow(trueE,0.90)) ) );

			if(abs(trueE_MeVee)<threshold) continue;

			// If hit above threshold, save the hit time to
			// find the earliest hit time
			double trueT = geantEvent->hits[j].time; // in ns

			if ((trueT < minTime) || (indexOfMinTime < 0)){
				minTime = trueT;
				indexOfMinTime = j;
			}
		}
		if (indexOfMinTime >= 0) {
			int j = indexOfMinTime;
			double trueE = geantEvent->hits[j].E_dep;
			double trueE_MeVee = 0.95 * trueE - 8.0 * ( 1 - exp( 0.10 * ( pow(trueE,0.90)) ) );

			// get true time, location, and bar number from the event
			double trueT = geantEvent->hits[j].time; // in ns
			TVector3 truePos( geantEvent->hits[j].pos.x()/10.,geantEvent->hits[j].pos.y()/10.,geantEvent->hits[j].pos.z()/10. ); // in cm
			int trueBarNo = geantEvent->hits[j].barNo;
			

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
	      	reconEvent->particles.clear();

	      	Recon_Particle reconHit_n;
	      	reconHit_n.E_dep = reconE; // in MeV
	      	reconHit_n.time = reconT;  // in ns
	      	reconHit_n.pos = reconPos; // in cm
	      	reconHit_n.momRecon = reconMom; // in MeV/c
	      	reconHit_n.barNo = trueBarNo;
	      	reconHit_n.type = type;

	      	// Now read the matching electron information and the true neutron momentum!
	      	inTree_generator->GetEntry(i);
	      		// for both the particle types in the generator tree:
	      	for (unsigned int p=0 ; p<generatedEvent->particles.size() ; p++){
	      		double trueMomMag, trueMomTheta, trueMomPhi;
	      		double trueMomX, trueMomY, trueMomZ;

	      		std::string particleType = generatedEvent->particles[p].type;
	      		if (particleType == "neutron"){ // Finish writing the true momentum & push back event
	      			trueMomMag = generatedEvent->particles[p].momentum.Mag() * 1000;
		      		trueMomTheta = generatedEvent->particles[p].momentum.Theta();
		      		trueMomPhi = generatedEvent->particles[p].momentum.Phi();

		      		trueMomX = trueMomMag * sin(trueMomTheta) * cos(trueMomPhi);
					trueMomY = trueMomMag * sin(trueMomTheta) * sin(trueMomPhi);
					trueMomZ = trueMomMag * cos(trueMomTheta);
					TVector3 trueMom( trueMomX ,trueMomY ,trueMomZ ); // in MeV/c

		      		reconHit_n.momTrue = trueMom;
		      		reconEvent->particles.push_back(reconHit_n);
		      	}
	      		else{
	      			double reconMomMag, reconMomTheta, reconMomPhi;
	      			double reconMomX, reconMomY, reconMomZ;
		      		// Get the true momentum and true angles:
		      		// recall that the momentum in the generator is in GeV and in Geant it is MeV!!!
		      		trueMomMag = generatedEvent->particles[p].momentum.Mag() * 1000;
		      		trueMomTheta = generatedEvent->particles[p].momentum.Theta();
		      		trueMomPhi = generatedEvent->particles[p].momentum.Phi();

		      		trueMomX = trueMomMag * sin(trueMomTheta) * cos(trueMomPhi);
					trueMomY = trueMomMag * sin(trueMomTheta) * sin(trueMomPhi);
					trueMomZ = trueMomMag * cos(trueMomTheta);
					TVector3 trueMom( trueMomX ,trueMomY ,trueMomZ ); // in MeV/c

					if (smearingOn == true){
			      		reconMomMag = myRand->Gaus(trueMomMag,clasFac*getCLAS12_PRes(trueMomTheta,trueMomMag));
			      		reconMomTheta = myRand->Gaus(trueMomTheta,clasFac*0.001);
			      		reconMomPhi = myRand->Gaus(trueMomPhi,clasFac*0.001/sin(trueMomTheta));
			      	}
			      	else{
			      		reconMomMag = trueMomMag;
			      		reconMomTheta = trueMomTheta;
			      		reconMomPhi = trueMomPhi;
			      	}

		      		if (reconMomPhi > M_PI) reconMomPhi -=2.*M_PI;
			    	if (reconMomPhi < -M_PI) reconMomPhi +=2.*M_PI;

					reconMomX = reconMomMag * sin(reconMomTheta) * cos(reconMomPhi);
					reconMomY = reconMomMag * sin(reconMomTheta) * sin(reconMomPhi);
					reconMomZ = reconMomMag * cos(reconMomTheta);
					TVector3 reconMom( reconMomX ,reconMomY ,reconMomZ ); // in MeV/c

		      		Recon_Particle reconHit_e;
		      		reconHit_e.type = particleType;
		      		reconHit_e.momRecon = reconMom;
		      		reconHit_e.momTrue = trueMom;
		      			// now i don't want floating 0's in the tree, so make them negative
		      		TVector3 reconPos(-2147483648,-2147483648,-2147483648); // in mm
		      		reconHit_e.pos = reconPos;
					reconHit_e.time = -2147483648; // in ns
					reconHit_e.E_dep = -2147483648; // in MeV
					reconHit_e.barNo = -2147483648; // indicating what bar
		      		reconEvent->particles.push_back(reconHit_e);
		      	}

	      	}
	      	outTree->Fill();

	    }
	  			
	}
	// ******************************************************************************************

		
	outTree->Write();
  	outFile->Close();
  	inFile_geant->Close();
  	inFile_generator->Close();
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