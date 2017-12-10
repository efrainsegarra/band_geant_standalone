#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TVectorT.h"
#include "TVector3.h"
#include "constants.h"  // contains constants E1
using namespace std;

const double luminosity = 1e35*1e-33;  // convert luminosity from cm^-2*s-1 -> nb-1*s-1
const double runtime = 40*24*60*60; // convert 40 day runtime to s
const double azim_CLAS12 = 0.5;  // 50% azimuthal electron coverage for CLAS12
const double azim_BAND = 0.6;    // 60% neutron coverage for BAND

int main(int argc, char **argv)
{
    // CLI arguments handling
    if (argc != 4)
    {
        cerr << "Wrong number of arguments, try:\n\t"
             << "ajalysis  path/to/inputfile  /path/to/outputfile  "
             << "Threshold_eeEDep" << endl;
        return -1;
    }
        
    // Unpack user input for input file, output file and band detection efficiency
    TFile *infile = new TFile(argv[1]);
    TFile *outfile = new TFile(argv[2], "RECREATE");
    double thresh_eeEDep = atof(argv[3]);
    
    // Get values from ReconTree
    TTree *inTree = (TTree*)infile -> Get("ReconTree");

    int num_events_sim = inTree->GetEntries();  // number of events generated in simulation
    double truePe[3], truePn[3], reconPe[3], reconPn[3];
    double reconWp, reconAs, reconQSq, reconXp, recon_eeEDep;
    inTree->SetBranchAddress("truePe_x", &(truePe[0]));
    inTree->SetBranchAddress("truePe_y", &(truePe[1]));
    inTree->SetBranchAddress("truePe_z", &(truePe[2]));
    inTree->SetBranchAddress("truePr_x", &(truePn[0]));
    inTree->SetBranchAddress("truePr_y", &(truePn[1]));
    inTree->SetBranchAddress("truePr_z", &(truePn[2]));
    inTree->SetBranchAddress("reconPe_x", &(reconPe[0]));
    inTree->SetBranchAddress("reconPe_y", &(reconPe[1]));
    inTree->SetBranchAddress("reconPe_z", &(reconPe[2]));
    inTree->SetBranchAddress("reconPr_x", &(reconPn[0]));
    inTree->SetBranchAddress("reconPr_y", &(reconPn[1]));
    inTree->SetBranchAddress("reconPr_z", &(reconPn[2]));
    inTree->SetBranchAddress("reconQSq", &(reconQSq));
    inTree->SetBranchAddress("reconXp", &(reconXp)); 
    inTree->SetBranchAddress("reconWp", &(reconWp));
    inTree->SetBranchAddress("reconAs", &(reconAs));
    inTree->SetBranchAddress("recon_eeEDep", &(recon_eeEDep));

    // Attempt to pull CS data from input root file
    TVectorT<double> *CSVec = NULL;
    CSVec = (TVectorT<double>*)infile->Get("totalCS");
    
    // Test if input root file from signal-generated run
    double tot_num_events;
    if (CSVec)
    {
        double CS = (*CSVec)[0];    // cross section in units of nb  
        tot_num_events = CS * luminosity * runtime;
    }
    // or from background run
    else
    {
         CSVec = (TVectorT<double>*)infile->Get("totalCSSq");
         
         // Error handling if no CS found in file
         if (!CSVec)
         {
            cerr << "No cross section was found in the file";
            return -2;
         }

         double CSSqDt = (*CSVec)[0] * 1e-9;   // CSSqDt is product of the cross-sections and coincidence time window
         tot_num_events = CSSqDt * pow(luminosity, 2) * runtime;
    }

    // Calculate weighting ratio
    double acceptance = azim_CLAS12 * azim_BAND;
    double num_events_detected = tot_num_events * acceptance;
    double weighting = num_events_detected / num_events_sim;

    // Hist output 
    TH1D *lx = new TH1D("lx", "Low x' Events; #alpha_{s}; Counts", 5, 1.30, 1.55);
    TH1D *hx = new TH1D("hx", "High x' Events; #alpha_{s}; Counts", 5, 1.30, 1.55);
    TH1D *Xp_dist = new TH1D("Xp_dist",
    "x' Distribution; x'; Counts", 14, 0.1, 0.8);
    TH1D *As_dist = new TH1D("As_dist",
    "#alpha_{s} Distribution; #alpha_{s}; Counts", 10, 1, 2);
    TH2D *Xp_As_dist = new TH2D("Xp_As",
	"#alpha_{s} vs x' Distribution; x'; #alpha_{s}; Counts", 14, 0.1, 0.8, 20, 1, 2);
    TH2D *the_pe_dist = new TH2D("the_pe_dist", 
    "#theta_{e} vs. #vec{p}_{e}; #theta_{e} [deg.]; #vec{p}_{e} #left[#frac{GeV}{c}#right]", 
    100, 5, 35, 100, 2, 7);
    TH2D *thrq_pr_dist = new TH2D("thrq_pr_dist", 
    "#theta_{rq} vs. #vec{p}_{r} Distribution; #theta_{rq} [deg.]; #vec{p}_{r} #left[#frac{GeV}{c}#right]",
    100, 140, 180, 100, 0.28, 0.52);
    TH2D *Xp_QSq_dist = new TH2D("Xp_QSq_dist", 
    "x' vs. Q^{2}; x'; Q^{2} #left[#frac{GeV}{c^{2}}#right]", 
    100, 0.1, 0.8, 100, 2, 8);

    // Counts for BAND efficiency
    int num_src = 0;
    int num_src_hits = 0;

    for(int iii = 0; iii < num_events_sim; ++iii)
    {
        // if (iii % 50000 == 0) cout << "Working on event " << iii << endl;
        inTree -> GetEvent(iii);

        TVector3 rePn(reconPn[0], reconPn[1], reconPn[2]);
        TVector3 rePe(reconPe[0], reconPe[1], reconPe[2]);
        TVector3 q = TVector3(0, 0, E1) - rePe;
        double theta_nq = rePn.Angle(q) * 180/TMath::Pi();
        double pn_mag = rePn.Mag();

        // DIS Selection Cuts
        if (reconQSq < 2) continue;
        if (reconWp < 1.8) continue;        
        if (theta_nq < 110) continue;
        if (pn_mag > 0.6 ) continue;
        if (reconXp > 1) continue;      // Get rid of pesky high x' background events
        ++num_src;

        // Band Selection Cut
        if (recon_eeEDep < thresh_eeEDep) continue;
        ++num_src_hits;
        
        // Fill histograms
        Xp_dist->Fill(reconXp, weighting);
        As_dist->Fill(reconAs, weighting);
        Xp_As_dist->Fill(reconXp, reconAs, weighting);

        if (reconXp > 0.25 && reconXp < 0.35) {
            lx -> Fill(reconAs, weighting);
        }
        else if (reconXp > 0.5) {
            hx -> Fill(reconAs, weighting);
        }
        
        double theta_e = rePe.Theta() * 180/TMath::Pi();
        double rePe_mag = rePe.Mag();
        the_pe_dist -> Fill(theta_e, rePe_mag, weighting);

        thrq_pr_dist -> Fill(theta_nq, pn_mag, weighting);

        Xp_QSq_dist -> Fill(reconXp, reconQSq, weighting);
    }
    
    // Calculate efficiency of BAND
    double band_eff = static_cast<double>(num_src_hits) / num_src;

    // Root file I/O
    infile -> Close();
    CSVec -> Write("CSdata");
    Xp_dist -> Write();
    As_dist -> Write();
    Xp_As_dist -> Write();
    lx -> Write();
    hx -> Write();
    the_pe_dist -> Write();
    thrq_pr_dist -> Write();
    Xp_QSq_dist -> Write();
    
    // Data (text) filenames string format
    ofstream txtfile;
    string root_filename(argv[2]);
    size_t raw_index = root_filename.find_last_of(string("."));
    string filename = root_filename.substr(0, raw_index);
    string dat_filename = filename + "_dat.txt";
    string cts_filename = filename + "_cts.txt";

    // Write efficiency data
    txtfile.open(dat_filename.c_str());
    txtfile << thresh_eeEDep << "\t" << band_eff << "\t" 
            << num_src << endl;
    txtfile.close();

    // Write signal and background counts for the As and x' bins
    txtfile.open(cts_filename.c_str());
    txtfile << thresh_eeEDep << "\t";

    TH1D* histlist[2] = {lx,hx};
    for (int h=0; h<2; ++h)
    {
        for (int bin=1; bin<=5; ++bin)
        {
            txtfile << histlist[h] -> GetBinContent(bin) << "\t";
        }
    }

    txtfile << "\n";
    txtfile.close();
    outfile -> Close();

    return 0;
}
