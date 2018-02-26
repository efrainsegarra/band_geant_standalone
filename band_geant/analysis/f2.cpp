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

	if (argc < 3){
	    cerr << "Wrong number of arguments. Instead use\n"
	    << "\t./makeF2 /path/to/signal/ /path/to/background/\n";
	    exit(-1);
	}

	TFile * inFileSig = new TFile(argv[1]);
	TH1F * As_highX_sig = new TH1F("As_highX_sig","Recon Q^2 > 2., W_prime > 1.8, 0.6 > Recon Pn > 0.25, Recon Xp < 1.;As;As;Counts",5,1.3,1.55);
	As_highX_sig = (TH1F*)inFileSig->Get("As_highX");

	TH1F * As_lowX_sig = new TH1F("As_lowX_sig","Recon Q^2 > 2., W_prime > 1.8, 0.6 > Recon Pn > 0.25, Recon Xp < 1.;As;As;Counts",5,1.3,1.55);
	As_lowX_sig = (TH1F*)inFileSig->Get("As_lowX");


	TFile * inFileBac = new TFile(argv[2]);
	TH1F * As_highX_bac = new TH1F("As_highX_bac","Recon Q^2 > 2., W_prime > 1.8, 0.6 > Recon Pn > 0.25, Recon Xp < 1.;As;As;Counts",5,1.3,1.55);
	As_highX_bac = (TH1F*)inFileBac->Get("As_highX");

	TH1F * As_lowX_bac = new TH1F("As_lowX_bac","Recon Q^2 > 2., W_prime > 1.8, 0.6 > Recon Pn > 0.25, Recon Xp < 1.;As;As;Counts",5,1.3,1.55);
	As_lowX_bac = (TH1F*)inFileBac->Get("As_lowX");

	// We have 5 bins we want in alphaS
	double signal_hi[5], background_hi[5];
	double signal_lo[5], background_lo[5];
	double bin_centers[5];

	for(int i=0; i<5;++i){
		bin_centers[i] = 1.325 + 0.05*i;
	}


	
	
	for (int i=1; i<6;i++){
		signal_lo[i-1] = As_lowX_sig->GetBinContent(i);
		signal_hi[i-1] = As_highX_sig->GetBinContent(i);
		background_lo[i-1] = As_lowX_bac->GetBinContent(i);
		background_hi[i-1] = As_highX_bac->GetBinContent(i);

	}

	cout << "# Uncertainty\tAlphaBinCenter\tSignal_hiX\tBackground_hiX\tSignal_lowX\tBackground_lowX\n";
	for(int i=0; i<5; ++i){
		double delta = pow((signal_hi[i] + background_hi[i]) / (pow(signal_hi[i],2)) + (signal_lo[i] + background_lo[i]) / (pow(signal_lo[i],2)),0.5);
		cout << delta << " " << bin_centers[i] << " " << signal_hi[i] << " " << background_hi[i] << " " << signal_lo[i] << " " << background_lo[i] << "\n";
	}


	inFileSig->Close();
	inFileBac->Close();


}