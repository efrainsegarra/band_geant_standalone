#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TVector3.h"
#include "TF1.h"
#include "TVectorT.h"

using namespace std;

int main(int argc, char** argv){

	if (argc < 3){
	    cerr << "Wrong number of arguments. Instead use\n"
	    << "\t./combineHists /path/to/output/file /path/to/allInput/file\n";
	    exit(-1);
	}

	TFile * outFile = new TFile(argv[1],"RECREATE");

	TH1D * total_Xp = new TH1D("Xp","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;xp;xp;Counts",14,0.1,0.8);

	TH1D * total_QSq_highX = new TH1D("QSq_highX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;QSq;QSq;Counts",42,2,8);
	TH1D * total_Wp_highX = new TH1D("Wp_highX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;Wp;Wp;Counts",13,1.8,3.1);
	TH1D * total_As_highX = new TH1D("As_highX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;As;As;Counts",5,1.3,1.55);
	TH1D * total_Pn_highX = new TH1D("Pn_highX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;Pn;Pn;Counts",6,0.28,0.42);
	TH1D * total_Theta_qn_highX = new TH1D("Theta_qn_highX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;Theta_qn;Theta_qn;Counts",40,140,180);

	TH1D * total_QSq_lowX = new TH1D("QSq_lowX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;QSq;QSq;Counts",42,2,8);
	TH1D * total_Wp_lowX = new TH1D("Wp_lowX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;Wp;Wp;Counts",13,1.8,3.1);
	TH1D * total_As_lowX = new TH1D("As_lowX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;As;As;Counts",5,1.3,1.55);
	TH1D * total_Pn_lowX = new TH1D("Pn_lowX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;Pn;Pn;Counts",6,0.28,0.42);
	TH1D * total_Theta_qn_lowX = new TH1D("Theta_qn_lowX","Recon Q^2 > 2., W_prime > 1.8, Recon Pn > 0.25;Theta_qn;Theta_qn;Counts",40,140,180);

  	
  	const int numFiles = argc - 2;
  	for(int i = 0; i < numFiles; ++i){

  		TFile * inFile = new TFile(argv[i+2]);

  		TH1D  * inXp = (TH1D*)inFile->Get("Xp");

  		TH1D  * in_QSq_highX = (TH1D*)inFile->Get("QSq_highX");
  		TH1D  * in_Wp_highX = (TH1D*)inFile->Get("Wp_highX");
  		TH1D  * in_As_highX = (TH1D*)inFile->Get("As_highX");
  		TH1D  * in_Pn_highX = (TH1D*)inFile->Get("Pn_highX");
  		TH1D  * in_Theta_qn_highX = (TH1D*)inFile->Get("Theta_qn_highX");

  		TH1D  * in_QSq_lowX = (TH1D*)inFile->Get("QSq_lowX");
  		TH1D  * in_Wp_lowX = (TH1D*)inFile->Get("Wp_lowX");
  		TH1D  * in_As_lowX = (TH1D*)inFile->Get("As_lowX");
  		TH1D  * in_Pn_lowX = (TH1D*)inFile->Get("Pn_lowX");
  		TH1D  * in_Theta_qn_lowX = (TH1D*)inFile->Get("Theta_qn_lowX");

  		total_Xp->Add(inXp);

  		total_QSq_highX->Add(in_QSq_highX);
  		total_Wp_highX->Add(in_Wp_highX);
  		total_As_highX->Add(in_As_highX);
  		total_Pn_highX->Add(in_Pn_highX);
  		total_Theta_qn_highX->Add(in_Theta_qn_highX);

  		total_QSq_lowX->Add(in_QSq_lowX);
  		total_Wp_lowX->Add(in_Wp_lowX);
  		total_As_lowX->Add(in_As_lowX);
  		total_Pn_lowX->Add(in_Pn_lowX);
  		total_Theta_qn_lowX->Add(in_Theta_qn_lowX);

  		inFile->Close();

  	}
	total_Xp->Scale(1./numFiles);

	total_QSq_highX->Scale(1./numFiles);
	total_Wp_highX->Scale(1./numFiles);
	total_As_highX->Scale(1./numFiles);
	total_Pn_highX->Scale(1./numFiles);
	total_Theta_qn_highX->Scale(1./numFiles);

	total_QSq_lowX->Scale(1./numFiles);
	total_Wp_lowX->Scale(1./numFiles);
	total_As_lowX->Scale(1./numFiles);
	total_Pn_lowX->Scale(1./numFiles);
	total_Theta_qn_lowX->Scale(1./numFiles);

	
  	outFile->cd();
  	total_QSq_highX->Write();
  	outFile->Close();


}
