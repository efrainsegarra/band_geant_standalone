#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorT.h"
#include "TVector3.h"

using namespace std;

int main(int argc, char** argv){

	if (argc < 3)
    {
      cerr << "\nWrong number of arguments! Instead use:\n"
	   << "\tcombineFiles /path/to/combined/root/ /path/to/allInput/root/ \n";
      exit(-1);
    }

    // Setting up the output file
    TFile * outFile = new TFile(argv[1],"RECREATE");
  	TTree * outTree = new TTree("ResTree","CombinedHistograms");

  	double total_trueWp, total_trueXp, total_trueAs;
	double total_reconWp, total_reconXp, total_reconAs;
	double total_truePe, total_truePn, total_reconPe, total_reconPn;
	double total_trueMom_e[3],total_trueMom_n[3];
	double total_reconMom_e[3],total_reconMom_n[3];

	double total_trueQSq, total_trueNu, total_trueX, total_true_q, total_truePhi_q, total_trueTheta_q;
	double total_reconQSq, total_reconNu, total_reconX, total_recon_q, total_reconPhi_q, total_reconTheta_q;
	double total_truePhi_e, total_truePhi_n, total_truePhi_en;
	double total_reconPhi_e, total_reconPhi_n, total_reconPhi_en;
	double total_trueTheta_n, total_reconTheta_n;
	double total_trueTheta_e, total_reconTheta_e;
	double total_trueEn, total_reconEn;
	double total_trueCosTheta_qn, total_reconCosTheta_qn;

	outTree->Branch("trueWp",&total_trueWp,"trueWp/D");
	outTree->Branch("trueXp",&total_trueXp,"trueXp/D");
	outTree->Branch("trueAs",&total_trueAs,"trueAs/D");

	outTree->Branch("reconWp",&total_reconWp,"reconWp/D");
	outTree->Branch("reconXp",&total_reconXp,"reconXp/D");
	outTree->Branch("reconAs",&total_reconAs,"reconAs/D");

	outTree->Branch("truePe",&total_truePe,"truePe/D");
	outTree->Branch("truePn",&total_truePn,"truePn/D");
	outTree->Branch("reconPe",&total_reconPe,"reconPe/D");
	outTree->Branch("reconPn",&total_reconPn,"reconPn/D");

	outTree->Branch("truePe_x",&(total_trueMom_e[0]),"truePe_x/D");
	outTree->Branch("truePe_y",&(total_trueMom_e[1]),"truePe_y/D");
	outTree->Branch("truePe_z",&(total_trueMom_e[2]),"truePe_z/D");

	outTree->Branch("truePn_x",&(total_trueMom_e[0]),"truePn_x/D");
	outTree->Branch("truePn_y",&(total_trueMom_e[1]),"truePn_y/D");
	outTree->Branch("truePn_z",&(total_trueMom_e[2]),"truePn_z/D");

	outTree->Branch("reconPe_x",&(total_reconMom_e[0]),"reconPe_x/D");
	outTree->Branch("reconPe_y",&(total_reconMom_e[1]),"reconPe_y/D");
	outTree->Branch("reconPe_z",&(total_reconMom_e[2]),"reconPe_z/D");

	outTree->Branch("reconPn_x",&(total_reconMom_n[0]),"reconPn_x/D");
	outTree->Branch("reconPn_y",&(total_reconMom_n[1]),"reconPn_y/D");
	outTree->Branch("reconPn_z",&(total_reconMom_n[2]),"reconPn_z/D");

	outTree->Branch("trueQSq",&(total_trueQSq),"trueQSq/D");
	outTree->Branch("trueNu",&(total_trueNu),"trueNu/D");
	outTree->Branch("trueX",&(total_trueX),"trueX/D");
	outTree->Branch("true_q",&(total_true_q),"true_q/D");
	outTree->Branch("truePhi_q",&(total_truePhi_q),"truePhi_q/D");
	outTree->Branch("trueTheta_q",&(total_trueTheta_q),"trueTheta_q/D");

	outTree->Branch("reconQSq",&(total_reconQSq),"reconQSq/D");
	outTree->Branch("reconNu",&(total_reconNu),"reconNu/D");
	outTree->Branch("reconX",&(total_reconX),"reconX/D");
	outTree->Branch("recon_q",&(total_recon_q),"recon_q/D");
	outTree->Branch("reconPhi_q",&(total_reconPhi_q),"reconPhi_q/D");
	outTree->Branch("reconTheta_q",&(total_reconTheta_q),"reconTheta_q/D");

	outTree->Branch("truePhi_e",&(total_truePhi_e),"truePhi_e/D");
	outTree->Branch("truePhi_n",&(total_truePhi_n),"truePhi_n/D");
	outTree->Branch("truePhi_en",&(total_truePhi_en),"truePhi_en/D");

	outTree->Branch("reconPhi_e",&(total_reconPhi_e),"reconPhi_e/D");
	outTree->Branch("reconPhi_n",&(total_reconPhi_n),"reconPhi_n/D");
	outTree->Branch("reconPhi_en",&(total_reconPhi_en),"reconPhi_en/D");

	outTree->Branch("trueTheta_n",&(total_trueTheta_n),"trueTheta_n/D");
	outTree->Branch("reconTheta_n",&(total_reconTheta_n),"reconTheta_n/D");
	outTree->Branch("trueTheta_e",&(total_trueTheta_e),"trueTheta_e/D");
	outTree->Branch("reconTheta_e",&(total_reconTheta_e),"reconTheta_e/D");

	outTree->Branch("trueEn",&(total_trueEn),"trueEn/D");
	outTree->Branch("reconEn",&(total_reconEn),"reconEn/D");

	outTree->Branch("trueCosTheta_qn",&(total_trueCosTheta_qn),"trueCosTheta_qn/D");
	outTree->Branch("reconCosTheta_qn",&(total_reconCosTheta_qn),"reconCosTheta_qn/D");


    int numFiles = argc - 2;
    cout << "Number of Files to Combine: " << numFiles << "\n";
    for( int i = 0 ; i < numFiles ; ++i){
    	cout << "\tWorking on file: " << argv[i+2] << "\n";
    	
    	TFile * inFile = new TFile(argv[i+2]);
   		TTree * inTree = (TTree*)inFile->Get("ResTree");

   		// Setting address for which branches to read in the input files
		double reconXp;
		inTree->SetBranchAddress("reconXp",&reconXp);

		// Looping over events in the input file
		const int nEvents = inTree->GetEntries();
		//cout << "Filling 2D histograms..." << endl;
		for (int j = 0 ; j < nEvents ; ++j){
			if (j % 100 == 0) cout << "\t\tEvent " << j << "\n";
			inTree->GetEvent(j);

			total_reconXp = reconXp;

			outTree->Fill();
		}

		inFile->Close();
	}

	outFile->cd();
  	outTree->Write();
  	outFile->Close();
  	   
}