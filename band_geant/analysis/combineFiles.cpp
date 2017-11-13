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
  	double total_reconXp;
  	outTree->Branch("reconXp",&total_reconXp,"reconXp/D");


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
  	outTree->Write();
  	outFile->Close();
  	   
}