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

	// need to read in signal and background file
	// and get the Xp histogram and collect
	// events in the high X and low X region

	if (argc < 4){
	    cerr << "Wrong number of arguments. Instead use\n"
	    << "\t./makeF2 /path/to/output/file /path/to/signal/ /path/to/background/\n";
	    exit(-1);
	}


	TFile * inFileSig = new TFile(argv[2]);
	TTree * inTreeSig = (TTree*)inFileSig->Get("ResTree");
	
	TFile * inFileBac = new TFile(argv[3]);
	TTree * inTreeBac = (TTree*)inFileBac->Get("ResTree");

	double reconXp_sig,reconXp_bac;
	double reconAs_sig,reconAs_bac;
	inTreeSig->SetBranchAddress("reconXp",&reconXp_sig);
	inTreeSig->SetBranchAddress("reconAs",&reconAs_sig);
	inTreeBac->SetBranchAddress("reconXp",&reconXp_bac);
	inTreeBac->SetBranchAddress("reconAs",&reconAs_bac);


	// We have 5 bins we want in alphaS
	double signal_hi[5], background_hi[5];
	double signal_lo[5], background_lo[5];

	const int nEvents_Sig = inTreeSig->GetEntries();
	for (int i = 0 ; i < nEvents_Sig ; ++i){
		inTreeSig->GetEvent(i);

		if ((reconXp_sig > 0.5) && (reconXp_sig < 1.5)){
			if (abs(reconAs_sig - 1.325)<0.025) signal_hi[0]++;
			if (abs(reconAs_sig - 1.375)<0.025) signal_hi[1]++;
			if (abs(reconAs_sig - 1.425)<0.025) signal_hi[2]++;
			if (abs(reconAs_sig - 1.475)<0.025) signal_hi[3]++;
			if (abs(reconAs_sig - 1.525)<0.025) signal_hi[4]++;
		}
		if (abs(reconXp_sig - 0.3) < 0.05){
			if (abs(reconAs_sig - 1.325)<0.025) signal_lo[0]++;
			if (abs(reconAs_sig - 1.375)<0.025) signal_lo[1]++;
			if (abs(reconAs_sig - 1.425)<0.025) signal_lo[2]++;
			if (abs(reconAs_sig - 1.475)<0.025) signal_lo[3]++;
			if (abs(reconAs_sig - 1.525)<0.025) signal_lo[4]++;
		}

	}

	const int nEvents_Bac = inTreeBac->GetEntries();
	for (int i = 0 ; i < nEvents_Bac ; ++i){
		inTreeBac->GetEvent(i);

		if ((reconXp_bac > 0.5) && (reconXp_bac < 1.5)){
			if (abs(reconAs_bac - 1.325)<0.025) background_hi[0]++;
			if (abs(reconAs_bac - 1.375)<0.025) background_hi[1]++;
			if (abs(reconAs_bac - 1.425)<0.025) background_hi[2]++;
			if (abs(reconAs_bac - 1.475)<0.025) background_hi[3]++;
			if (abs(reconAs_bac - 1.525)<0.025) background_hi[4]++;
		}
		if (abs(reconXp_bac - 0.3) < 0.05){
			if (abs(reconAs_bac - 1.325)<0.025) background_lo[0]++;
			if (abs(reconAs_bac - 1.375)<0.025) background_lo[1]++;
			if (abs(reconAs_bac - 1.425)<0.025) background_lo[2]++;
			if (abs(reconAs_bac - 1.475)<0.025) background_lo[3]++;
			if (abs(reconAs_bac - 1.525)<0.025) background_lo[4]++;
		}

	}

	for(int i=0; i<5; ++i){
		double delta = pow((signal_hi[i] + background_hi[i]) / (pow(signal_hi[i],2)) + (signal_lo[i] + background_lo[i]) / (pow(signal_lo[i],2)),0.5);

		cout << delta << " " << signal_hi[i] << " " << background_hi[i] << " " << signal_lo[i] << " " << background_lo[i] << "\n";
	}

}