#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorT.h"
#include "TVector3.h"

using namespace std;

const double luminosity = 1e35*1e-33; 
const double runtime = 40*24*60*60; 
const double azim_CLAS12 = 0.5;

const int numFiles = 200;
const int numGenEvents = 50000;

int main(int argc, char** argv){

	if (argc < 4)
    {
      cerr << "\nWrong number of arguments! Instead use:\n"
	   << "\tcombineFiles /path/to/combined/output /path/to/generator/allInput/ /path/to/kinVar/allInput/ \n";
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

	outTree->Branch("truePn_x",&(total_trueMom_n[0]),"truePn_x/D");
	outTree->Branch("truePn_y",&(total_trueMom_n[1]),"truePn_y/D");
	outTree->Branch("truePn_z",&(total_trueMom_n[2]),"truePn_z/D");

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


    int kinVar_startInd = 2+numFiles;
    int gen_startInd = 2;
    cout << "Number of Files to Combine: " << numFiles << "\n";
    for( int i = 0 ; i < numFiles ; ++i){
    	cout << "\tWorking on kin file   : " << argv[i+kinVar_startInd] << "\n";
    	cout << "\t Corresponding gen file: " << argv[i+gen_startInd] << "\n";
    	cout << "\t i: " << i << " " << argc << " " << i+kinVar_startInd << " " <<  i+gen_startInd << "\n";


		// Create the weight unique to this generation
		TFile * inGen = new TFile(argv[i+gen_startInd]);
		TTree * inTreeGen = (TTree*)inGen->Get("MCout");
		double numEvents_sim = inTreeGen->GetEntries(); 
		if (numGenEvents != numEvents_sim ){
			cerr << "*** GENERATED BATCH DOESN'T MATCH WITH NUM EVENTS ***\n"; exit(-1);
		}

		double total_num_events;
		TVectorT<double> *CSVec_dis = NULL;
		TVectorT<double> *CSVec_rand = NULL;
		CSVec_dis  = (TVectorT<double>*)inGen->Get("totalCS");
		CSVec_rand = (TVectorT<double>*)inGen->Get("totalCSSq");

		if (CSVec_dis){
			double CS = (*CSVec_dis)[0];
		    total_num_events = CS * luminosity * runtime;
		}
		else{
			double CSSqDt = (*CSVec_rand)[0] * 1e-9;   // CSSqDt is product of the cross-sections and coincidence time window
			total_num_events = CSSqDt * pow(luminosity, 2) * runtime;
		}

		// Calculate weighting ratio
		double acceptance = azim_CLAS12;
		double num_events_detected = total_num_events * acceptance;
		double weighting = (num_events_detected / numEvents_sim) / (numFiles);

		// Now import the kinematic root tree and combine into one
    	
    	TFile * inFile = new TFile(argv[i+kinVar_startInd]);
   		TTree * inTree = (TTree*)inFile->Get("ResTree");

   		// Setting address for which branches to read in the input files
		double trueWp, trueXp, trueAs;
		double reconWp, reconXp, reconAs;
		double truePe, truePn, reconPe, reconPn;
		double truereconMom_e[3],trueMom_n[3];
		double reconreconMom_e[3],reconMom_n[3];

		double trueQSq, trueNu, trueX, true_q, truePhi_q, trueTheta_q;
		double reconQSq, reconNu, reconX, recon_q, reconPhi_q, reconTheta_q;
		double truePhi_e, truePhi_n, truePhi_en;
		double reconPhi_e, reconPhi_n, reconPhi_en;
		double trueTheta_n, reconTheta_n;
		double trueTheta_e, reconTheta_e;
		double trueEn, reconEn;
		double trueCosTheta_qn, reconCosTheta_qn;

		inTree->SetBranchAddress("trueWp",&trueWp);
		inTree->SetBranchAddress("trueXp",&trueXp);
		inTree->SetBranchAddress("trueAs",&trueAs);

		inTree->SetBranchAddress("reconWp",&reconWp);
		inTree->SetBranchAddress("reconXp",&reconXp);
		inTree->SetBranchAddress("reconAs",&reconAs);

		inTree->SetBranchAddress("truePe",&truePe);
		inTree->SetBranchAddress("truePn",&truePn);
		inTree->SetBranchAddress("reconPe",&reconPe);
		inTree->SetBranchAddress("reconPn",&reconPn);

		inTree->SetBranchAddress("truePe_x",&(truereconMom_e[0]));
		inTree->SetBranchAddress("truePe_y",&(truereconMom_e[1]));
		inTree->SetBranchAddress("truePe_z",&(truereconMom_e[2]));

		inTree->SetBranchAddress("truePn_x",&(trueMom_n[0]));
		inTree->SetBranchAddress("truePn_y",&(trueMom_n[1]));
		inTree->SetBranchAddress("truePn_z",&(trueMom_n[2]));

		inTree->SetBranchAddress("reconPe_x",&(reconreconMom_e[0]));
		inTree->SetBranchAddress("reconPe_y",&(reconreconMom_e[1]));
		inTree->SetBranchAddress("reconPe_z",&(reconreconMom_e[2]));

		inTree->SetBranchAddress("reconPn_x",&(reconMom_n[0]));
		inTree->SetBranchAddress("reconPn_y",&(reconMom_n[1]));
		inTree->SetBranchAddress("reconPn_z",&(reconMom_n[2]));

		inTree->SetBranchAddress("trueQSq",&(trueQSq));
		inTree->SetBranchAddress("trueNu",&(trueNu));
		inTree->SetBranchAddress("trueX",&(trueX));
		inTree->SetBranchAddress("true_q",&(true_q));
		inTree->SetBranchAddress("truePhi_q",&(truePhi_q));
		inTree->SetBranchAddress("trueTheta_q",&(trueTheta_q));

		inTree->SetBranchAddress("reconQSq",&(reconQSq));
		inTree->SetBranchAddress("reconNu",&(reconNu));
		inTree->SetBranchAddress("reconX",&(reconX));
		inTree->SetBranchAddress("recon_q",&(recon_q));
		inTree->SetBranchAddress("reconPhi_q",&(reconPhi_q));
		inTree->SetBranchAddress("reconTheta_q",&(reconTheta_q));

		inTree->SetBranchAddress("truePhi_e",&(truePhi_e));
		inTree->SetBranchAddress("truePhi_n",&(truePhi_n));
		inTree->SetBranchAddress("truePhi_en",&(truePhi_en));

		inTree->SetBranchAddress("reconPhi_e",&(reconPhi_e));
		inTree->SetBranchAddress("reconPhi_n",&(reconPhi_n));
		inTree->SetBranchAddress("reconPhi_en",&(reconPhi_en));

		inTree->SetBranchAddress("trueTheta_n",&(trueTheta_n));
		inTree->SetBranchAddress("reconTheta_n",&(reconTheta_n));
		inTree->SetBranchAddress("trueTheta_e",&(trueTheta_e));
		inTree->SetBranchAddress("reconTheta_e",&(reconTheta_e));

		inTree->SetBranchAddress("trueEn",&(trueEn));
		inTree->SetBranchAddress("reconEn",&(reconEn));

		inTree->SetBranchAddress("trueCosTheta_qn",&(trueCosTheta_qn));
		inTree->SetBranchAddress("reconCosTheta_qn",&(reconCosTheta_qn));

		// Looping over events in the input file
		const int nEvents = inTree->GetEntries();
		//cout << "Filling 2D histograms..." << endl;
		for (int j = 0 ; j < nEvents ; ++j){
			if (j % 100 == 0) cout << "\t\tEvent " << j << "\n";
			inTree->GetEvent(j);

			// Now adding each branch to the full total output tree
			total_trueWp = trueWp;
			total_trueXp = trueXp;
			total_trueAs = trueAs;
			total_reconWp = reconWp;
			total_reconXp = reconXp;
			total_reconAs = reconAs;
			total_truePe = truePe;
			total_truePn = truePn;
			total_reconPe = reconPe;
			total_reconPn = reconPn;
			memcpy(total_trueMom_e, truereconMom_e, sizeof(truereconMom_e));
			memcpy(total_trueMom_n, trueMom_n, sizeof(trueMom_n));
			memcpy(total_reconMom_e, reconreconMom_e, sizeof(reconreconMom_e));
			memcpy(total_reconMom_n, reconMom_n, sizeof(reconMom_n));

			total_trueQSq = trueQSq;
			total_trueNu = trueNu;
			total_trueX = trueX;
			total_true_q = true_q;
			total_truePhi_q = truePhi_q;
			total_trueTheta_q = trueTheta_q;
			total_reconQSq = reconQSq;
			total_reconNu = reconNu;
			total_reconX = reconX;
			total_recon_q = recon_q;
			total_reconPhi_q = reconPhi_q;
			total_reconTheta_q = reconTheta_q;
			total_truePhi_e = truePhi_e;
			total_truePhi_n = truePhi_n;
			total_truePhi_en = truePhi_en;
			total_reconPhi_e = reconPhi_e;
			total_reconPhi_n = reconPhi_n;
			total_reconPhi_en = reconPhi_en;
			total_trueTheta_n = trueTheta_n;
			total_reconTheta_n = reconTheta_n;
			total_trueTheta_e = trueTheta_e;
			total_reconTheta_e = reconTheta_e;
			total_trueEn = trueEn;
			total_reconEn = reconEn;
			total_trueCosTheta_qn = trueCosTheta_qn;
			total_reconCosTheta_qn = reconCosTheta_qn;

			outTree->SetWeight(weighting);
			outTree->Fill();


		}

		inFile->Close();
	}


	
	

	


	outFile->cd();
  	outTree->Write();
  	outFile->Close();
  	   
}