#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorT.h"
#include "TVector3.h"

#include "recon_tree.h"
#include "constants.h"
using namespace std;
inline double sq(double x){ return x*x;};

int main(int argc, char** argv){

	if (argc != 3)
    {
      cerr << "\nWrong number of arguments! Instead use:\n"
	   << "\treconstructions /path/to/digitized/file /path/to/output/file\n";
      exit(-1);
    }

    // Set up input and output root files
    TFile * inFile_digitized = new TFile(argv[1]);
    TTree * inTree_digitized = (TTree*)inFile_digitized->Get("ReconTree");

  	TFile * outFile = new TFile(argv[2],"RECREATE");
  	TTree * outTree = new TTree("ResTree","Resolutions");
	 
	// Initialize the branches
	Recon_Event * reconEvent = NULL;
	inTree_digitized-> SetBranchAddress("event",&reconEvent);

	// Memory for the output tree
	double trueWp, trueXp, trueAs;
	double reconWp, reconXp, reconAs;
	double truePe, truePn, reconPe, reconPn;
	double trueMom_e[3],trueMom_n[3];
	double reconMom_e[3],reconMom_n[3];

	double trueQSq, trueNu, trueX, true_q, truePhi_q, trueTheta_q;
	double reconQSq, reconNu, reconX, recon_q, reconPhi_q, reconTheta_q;
	double truePhi_e, truePhi_n, truePhi_en;
	double reconPhi_e, reconPhi_n, reconPhi_en;
	double trueTheta_n, reconTheta_n;
	double trueTheta_e, reconTheta_e;
	double trueEn, reconEn;
	double trueCosTheta_qn, reconCosTheta_qn;

	outTree->Branch("trueWp",&trueWp,"trueWp/D");
	outTree->Branch("trueXp",&trueXp,"trueXp/D");
	outTree->Branch("trueAs",&trueAs,"trueAs/D");

	outTree->Branch("reconWp",&reconWp,"reconWp/D");
	outTree->Branch("reconXp",&reconXp,"reconXp/D");
	outTree->Branch("reconAs",&reconAs,"reconAs/D");

	outTree->Branch("truePe",&truePe,"truePe/D");
	outTree->Branch("truePn",&truePn,"truePn/D");
	outTree->Branch("reconPe",&reconPe,"reconPe/D");
	outTree->Branch("reconPn",&reconPn,"reconPn/D");

	outTree->Branch("truePe_x",&(trueMom_e[0]),"truePe_x/D");
	outTree->Branch("truePe_y",&(trueMom_e[1]),"truePe_y/D");
	outTree->Branch("truePe_z",&(trueMom_e[2]),"truePe_z/D");

	outTree->Branch("truePn_x",&(trueMom_e[0]),"truePn_x/D");
	outTree->Branch("truePn_y",&(trueMom_e[1]),"truePn_y/D");
	outTree->Branch("truePn_z",&(trueMom_e[2]),"truePn_z/D");

	outTree->Branch("reconPe_x",&(reconMom_e[0]),"reconPe_x/D");
	outTree->Branch("reconPe_y",&(reconMom_e[1]),"reconPe_y/D");
	outTree->Branch("reconPe_z",&(reconMom_e[2]),"reconPe_z/D");

	outTree->Branch("reconPn_x",&(reconMom_n[0]),"reconPn_x/D");
	outTree->Branch("reconPn_y",&(reconMom_n[1]),"reconPn_y/D");
	outTree->Branch("reconPn_z",&(reconMom_n[2]),"reconPn_z/D");

	outTree->Branch("trueQSq",&(trueQSq),"trueQSq/D");
	outTree->Branch("trueNu",&(trueNu),"trueNu/D");
	outTree->Branch("trueX",&(trueX),"trueX/D");
	outTree->Branch("true_q",&(true_q),"true_q/D");
	outTree->Branch("truePhi_q",&(truePhi_q),"truePhi_q/D");
	outTree->Branch("trueTheta_q",&(trueTheta_q),"trueTheta_q/D");

	outTree->Branch("reconQSq",&(reconQSq),"reconQSq/D");
	outTree->Branch("reconNu",&(reconNu),"reconNu/D");
	outTree->Branch("reconX",&(reconX),"reconX/D");
	outTree->Branch("recon_q",&(recon_q),"recon_q/D");
	outTree->Branch("reconPhi_q",&(reconPhi_q),"reconPhi_q/D");
	outTree->Branch("reconTheta_q",&(reconTheta_q),"reconTheta_q/D");

	outTree->Branch("truePhi_e",&(truePhi_e),"truePhi_e/D");
	outTree->Branch("truePhi_n",&(truePhi_n),"truePhi_n/D");
	outTree->Branch("truePhi_en",&(truePhi_en),"truePhi_en/D");

	outTree->Branch("reconPhi_e",&(reconPhi_e),"reconPhi_e/D");
	outTree->Branch("reconPhi_n",&(reconPhi_n),"reconPhi_n/D");
	outTree->Branch("reconPhi_en",&(reconPhi_en),"reconPhi_en/D");

	outTree->Branch("trueTheta_n",&(trueTheta_n),"trueTheta_n/D");
	outTree->Branch("reconTheta_n",&(reconTheta_n),"reconTheta_n/D");
	outTree->Branch("trueTheta_e",&(trueTheta_e),"trueTheta_e/D");
	outTree->Branch("reconTheta_e",&(reconTheta_e),"reconTheta_e/D");

	outTree->Branch("trueEn",&(trueEn),"trueEn/D");
	outTree->Branch("reconEn",&(reconEn),"reconEn/D");

	outTree->Branch("trueCosTheta_qn",&(trueCosTheta_qn),"trueCosTheta_qn/D");
	outTree->Branch("reconCosTheta_qn",&(reconCosTheta_qn),"reconCosTheta_qn/D");


	int numEvents = inTree_digitized->GetEntries();
	for(int i =0; i<numEvents; i++){
		inTree_digitized->GetEntry(i);
		if (i % 10000 == 0) cout << "Event " << i << endl;


		
		for (unsigned int p=0 ; p<reconEvent->particles.size() ; p++){

			double momTrueMag = reconEvent->particles[p].momTrue.Mag()/1000.;
			double momTrueTheta = reconEvent->particles[p].momTrue.Theta();
			double momTruePhi = reconEvent->particles[p].momTrue.Phi();

			double momReconMag = reconEvent->particles[p].momRecon.Mag()/1000.;
			double momReconTheta = reconEvent->particles[p].momRecon.Theta();
			double momReconPhi = reconEvent->particles[p].momRecon.Phi();

			
			if (reconEvent->particles[p].type == "e-"){
				truePe = momTrueMag;
				reconPe = momReconMag;

				trueMom_e[0] = momTrueMag * sin(momTrueTheta) * cos(momTruePhi);
				trueMom_e[1] = momTrueMag * sin(momTrueTheta) * sin(momTruePhi);
				trueMom_e[2] = momTrueMag * cos(momTrueTheta);

				reconMom_e[0] = momReconMag * sin(momReconTheta) * cos(momReconPhi);
				reconMom_e[1] = momReconMag * sin(momReconTheta) * sin(momReconPhi);
				reconMom_e[2] = momReconMag * cos(momReconTheta);

				// Kinematical Variables -- true
				trueQSq = 2.*E1*momTrueMag*(1.-cos(momTrueTheta));
				trueNu = E1 - momTrueMag;
				trueX = trueQSq / (2.*mP*trueNu);
				true_q = sqrt(trueQSq + sq(trueNu));
				truePhi_q = momTruePhi + M_PI;
				if (truePhi_q > M_PI) truePhi_q -= 2.*M_PI;
				trueTheta_q = acos((E1 - momTrueMag*cos(momTrueTheta))/true_q);
				// Kinematical Variables -- reconstructed
				reconQSq = 2.*E1*momReconMag*(1.-cos(momReconTheta));
				reconNu = E1 - momReconMag;
				reconX = reconQSq / (2.*mP*reconNu);
				recon_q = sqrt(reconQSq + sq(reconNu));
				reconPhi_q = momReconPhi + M_PI;
				if (reconPhi_q > M_PI) reconPhi_q -= 2.*M_PI;
				reconTheta_q = acos((E1 - momReconMag*cos(momReconTheta))/recon_q);

				truePhi_e = momTruePhi;
				reconPhi_e = momReconPhi;

				trueTheta_e = momTrueTheta;
				reconTheta_e = momReconTheta;

			}
			else{
				truePn = momTrueMag;
				reconPn = momReconMag;
				trueTheta_n = momTrueTheta;
				reconTheta_n = momReconTheta;
				truePhi_n = momTruePhi;
				reconPhi_n = momReconPhi;

				trueEn = sqrt(sq(momTrueMag) + sq(mN));
				reconEn = sqrt(sq(momReconMag) + sq(mN));

				trueMom_n[0] = momTrueMag * sin(momTrueTheta) * cos(momTruePhi);
				trueMom_n[1] = momTrueMag * sin(momTrueTheta) * sin(momTruePhi);
				trueMom_n[2] = momTrueMag * cos(momTrueTheta);

				reconMom_n[0] = momReconMag * sin(momReconTheta) * cos(momReconPhi);
				reconMom_n[1] = momReconMag * sin(momReconTheta) * sin(momReconPhi);
				reconMom_n[2] = momReconMag * cos(momReconTheta);

			}
			
		}
		// Kinematical Variables that require both e / n information
		truePhi_en = truePhi_n - truePhi_e;
		if (truePhi_en < -M_PI) truePhi_en += 2.*M_PI;
		if (truePhi_en > M_PI) truePhi_en -= 2.*M_PI;
		reconPhi_en = reconPhi_n - reconPhi_e;
		if (reconPhi_en < -M_PI) reconPhi_en += 2.*M_PI;
		if (reconPhi_en > M_PI) reconPhi_en -= 2.*M_PI;

		trueCosTheta_qn = cos(trueTheta_n)*cos(trueTheta_q) - sin(trueTheta_n)*sin(trueTheta_q)*cos(truePhi_en);
		reconCosTheta_qn = cos(reconTheta_n)*cos(reconTheta_q) - sin(reconTheta_n)*sin(reconTheta_q)*cos(reconPhi_en);

		trueXp = trueQSq/(2.*( trueNu*(mD-trueEn) + truePn*true_q*trueCosTheta_qn));
		reconXp = reconQSq/(2.*( reconNu*(mD-reconEn) + reconPn*recon_q*reconCosTheta_qn));

		// W prime
		trueWp = sqrt(sq(mD) - trueQSq + sq(mP) + 2.*mD*(trueNu-trueEn) -2.* trueNu * trueEn + 2.*true_q*truePn*trueCosTheta_qn);
		reconWp = sqrt(sq(mD) - reconQSq + sq(mP) + 2.*mD*(reconNu-reconEn) -2.* reconNu * reconEn + 2.*recon_q*reconPn*reconCosTheta_qn);

		// alpha_s
		trueAs = (trueEn - truePn*trueCosTheta_qn)/mN;
		reconAs = (reconEn - reconPn*reconCosTheta_qn)/mN;

		//if (acos(reconCosTheta_qn)*180./M_PI < 110.) continue;
		//if (reconQSq < 2) continue;
		//if (reconWp < 1.8) continue;

		outTree->Fill();
	}



	gFile=outFile;
	outTree->Write();
	// Clean up
	inFile_digitized->Close();
	outFile->Close();

}