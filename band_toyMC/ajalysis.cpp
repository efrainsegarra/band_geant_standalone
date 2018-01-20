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
             << "Threshold(in MeVee)" << endl;
        return -1;
    }
        
    // Unpack user input for input file, output file and band detection efficiency
    TFile *infile = new TFile(argv[1]);
    TFile *outfile = new TFile(argv[2], "RECREATE");
    double thresh_eeEDep = atof(argv[3]);
    
    // Get values from ReconTree
    TTree *inTree = (TTree*)infile -> Get("ReconTree");

    int num_events_sim = inTree->GetEntries();  // number of events generated in simulation
    double truePe[3], truePr[3], reconPe[3], reconPr[3];
    double reconWp, reconAs, reconQSq, reconXp, ErDep;
    inTree->SetBranchAddress("truePe_x", &(truePe[0]));
    inTree->SetBranchAddress("truePe_y", &(truePe[1]));
    inTree->SetBranchAddress("truePe_z", &(truePe[2]));
    inTree->SetBranchAddress("truePr_x", &(truePr[0]));
    inTree->SetBranchAddress("truePr_y", &(truePr[1]));
    inTree->SetBranchAddress("truePr_z", &(truePr[2]));
    inTree->SetBranchAddress("reconPe_x", &(reconPe[0]));
    inTree->SetBranchAddress("reconPe_y", &(reconPe[1]));
    inTree->SetBranchAddress("reconPe_z", &(reconPe[2]));
    inTree->SetBranchAddress("reconPr_x", &(reconPr[0]));
    inTree->SetBranchAddress("reconPr_y", &(reconPr[1]));
    inTree->SetBranchAddress("reconPr_z", &(reconPr[2]));
    inTree->SetBranchAddress("reconQSq", &(reconQSq));
    inTree->SetBranchAddress("reconXp", &(reconXp)); 
    inTree->SetBranchAddress("reconWp", &(reconWp));
    inTree->SetBranchAddress("reconAs", &(reconAs));
    inTree->SetBranchAddress("ErDep", &(ErDep));

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
    // Reach histograms
    TH1D *Xp_dist = new TH1D("Xp_dist",
    "x' Distribution; x'; Counts", 14, 0.1, 0.8);
    TH1D *As_dist = new TH1D("As_dist",
    "#alpha_{s} Distribution; #alpha_{s}; Counts", 10, 1, 2);
    TH1D *lx_As = new TH1D("lx_As", "#alpha_{s} Distribution for Low x' Events; #alpha_{s}; Counts", 5, 1.30, 1.55);
    TH1D *hx_As = new TH1D("hx_As", "#alpha_{s} Distribution for High x' Events; #alpha_{s}; Counts", 5, 1.30, 1.55);
    TH1D *lx_qsq = new TH1D("lx_qsq", "Q^{2} Distribution for Low x' Events; Q^{2} [GeV}]^{2}; Counts", 60, 2, 8);
    TH1D *hx_qsq = new TH1D("hx_qsq", "Q^{2} Distribution for High x' Events; Q^{2} [GeV]^{2}; Counts", 60, 2, 8);
    TH1D *lx_theqr = new TH1D("lx_theqr", "#theta_{qr} Distribution for Low x' Events; #theta_{qr} [deg.]; Counts", 80, 140, 180);
    TH1D *hx_theqr = new TH1D("hx_theqr", "#theta_{qr} Distribution for High x' Events; #theta_{qr} [deg.]; Counts", 80, 140, 180);
    TH1D *lx_wp = new TH1D("lx_wp", "W' Distribution for Low x' Events; W' [GeV]; Counts", 14, 1.8, 3.2);
    TH1D *hx_wp = new TH1D("hx_wp", "W' Distribution for High x' Events; W' [GeV]; Counts", 14, 1.8, 3.2);
    TH1D *lx_pr = new TH1D("lx_pr", "#vec{p}_{r} Distribution for Low x' Events; |#vec{p}_{r}|; Counts", 14, 0.28, 0.42);
    TH1D *hx_pr = new TH1D("hx_pr", "#vec{p}_{r} Distribution for High x' Events; |#vec{p}_{r}|; Counts", 14, 0.28, 0.42);
    
    // Kinematic distributions histograms
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
    int num_dis = 0;
    int num_dis_hits = 0;

    for(int iii = 0; iii < num_events_sim; ++iii)
    {
        // if (iii % 50000 == 0) cout << "Working on event " << iii << endl;
        inTree -> GetEvent(iii);

        // Calculate values from reconstructed momentum vectors
        TVector3 rePn(reconPr[0], reconPr[1], reconPr[2]);
        double rePn_mag = rePn.Mag();

        TVector3 rePe(reconPe[0], reconPe[1], reconPe[2]);
        double rePe_mag = rePe.Mag();
        TVector3 q = TVector3(0, 0, E1) - rePe;
        double theta_e = rePe.Theta() * 180/TMath::Pi();

        double theta_qn = rePn.Angle(q) * 180/TMath::Pi();

        // Fill neutron histograms for low and high x'
        if (reconXp > 0.25 && reconXp < 0.35)
        {
            lx_qsq -> Fill(reconQSq, weighting);
            lx_theqr -> Fill(theta_qn, weighting);
            lx_wp -> Fill(reconWp, weighting);
            lx_pr -> Fill(rePn_mag, weighting);
        }

        if (reconXp > 0.5 && reconXp < 1)
        {
            hx_qsq -> Fill(reconQSq, weighting);
            hx_theqr -> Fill(theta_qn, weighting);
            hx_wp -> Fill(reconWp, weighting);
            hx_pr -> Fill(rePn_mag, weighting);
        }

        // Perform DIS Selection Cuts
        if (reconQSq < 2) continue;
        if (reconWp < 1.8) continue;        
        if (theta_qn < 110) continue;
        ++num_dis;

        // Perform BAND Selection Cut
        if (ErDep < thresh_eeEDep) continue;
        ++num_dis_hits;
        
        // Fill reach histograms
        Xp_dist->Fill(reconXp, weighting);
        As_dist->Fill(reconAs, weighting);
        Xp_As_dist->Fill(reconXp, reconAs, weighting);

        if (reconXp > 0.25 && reconXp < 0.35) lx_As -> Fill(reconAs, weighting);
        else if (reconXp > 0.5 && reconXp < 1) hx_As -> Fill(reconAs, weighting);
        
        // Fill kinematic distribution of DIS events hists        
        the_pe_dist -> Fill(theta_e, rePe_mag, weighting);
        thrq_pr_dist -> Fill(theta_qn, rePn_mag, weighting);
        Xp_QSq_dist -> Fill(reconXp, reconQSq, weighting);

    }
    
    // Calculate efficiency of BAND
    double band_eff = static_cast<double>(num_dis_hits) / num_dis;

    // Root file I/O
    infile -> Close();
    CSVec -> Write("CSdata");
    Xp_dist -> Write();
    As_dist -> Write();
    Xp_As_dist -> Write();
    the_pe_dist -> Write();
    thrq_pr_dist -> Write();
    Xp_QSq_dist -> Write();
    lx_As -> Write();
    hx_As -> Write();
    lx_qsq -> Write();
    hx_qsq -> Write();
    lx_theqr -> Write();
    hx_theqr -> Write();
    lx_wp -> Write();
    hx_wp -> Write();
    lx_pr -> Write();
    hx_pr -> Write();
    
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
            << num_dis << endl;
    txtfile.close();

    // Write signal and background counts for the As and x' bins
    txtfile.open(cts_filename.c_str());
    txtfile << thresh_eeEDep << "\t";

    TH1D* histlist[2] = {lx_As,hx_As};
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
