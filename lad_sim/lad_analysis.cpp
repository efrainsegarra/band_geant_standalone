#include <iostream>
#include <cmath>
#include <cstdlib>

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVector3.h"
#include "TVectorT.h"

#include "lad_helpers.h"
#include "gen_tree.h"
#include "constants.h"

using namespace std;

int main(int argc, char ** argv)
{
  if (argc != 5)
    {
      cerr << "Wrong number of arguments. Instead use\n\tlad_analysis /path/to/gen/file /path/to/prop/file /path/to/digi/file /path/to/output/file\n";
      exit(-1);
    }

  TFile * genFile = new TFile(argv[1]);
  TFile * proFile = new TFile(argv[2]);
  TFile * digFile = new TFile(argv[3]);
  TFile * outFile = new TFile(argv[4],"RECREATE");

  // Copy the cross section
  TVectorT<double> * csVec;
  const double lumi=2000.; // 2E36 in /nb/s
  const double duration=3600.; // one hour of run time
  double weight=0.;
  if (genFile->GetListOfKeys()->Contains("totalCS"))
    {
      csVec = (TVectorT<double>*)genFile->Get("totalCS");
      gFile=outFile;
      csVec->Write("totalCS");
      weight = (*csVec)[0]*lumi*duration;
    }
  else if (genFile->GetListOfKeys()->Contains("totalCSSq"))
    {
      csVec = (TVectorT<double>*)genFile->Get("totalCSSq");
      gFile=outFile;
      csVec->Write("totalCSSq");
      weight = (*csVec)[0]*lumi*lumi*duration;
    }
  else
    {
      cerr << "No cross section found in the file. Exiting...\n";
      exit(-3);
    }

  // Prep some histograms
  TH1D * histZRecon = new TH1D("z_res","Vertex resolution;z_recon - z_gen [cm];Counts",100,-20.,20.);
  TH2D * histMomRes = new TH2D("mom_res","Recoil protons;true mom [GeV];Recon - true [GeV];Counts",100,0.,1.,100,-0.05,0.05);

  TH3D * hmsHist = new TH3D("hms","hms_acceptance;x_prime;angle;momentum;Counts",20,0.,1.,20,4.,24.,200,0.,10.0);
  TH3D * shmsHist = new TH3D("shms","shms_acceptance;x_prime;angle;momentum;Counts",20,0.,1.,200,4.,24.,200,0.,10.0);

  // Get the trees
  TTree * genTree = (TTree*) genFile->Get("MCout");
  TTree * proTree = (TTree*) proFile->Get("PropTree");
  TTree * digTree = (TTree*) digFile->Get("ReconTree");
  
  if (genTree->GetEntries() != proTree->GetEntries())
    {
      cerr << "The trees in " << argv[1] << " and " << argv[2] << " have different numbers of entries."
	   << "Double check that you have the correct files.\n\n";
      exit(-2);
    }

  if (genTree->GetEntries() != digTree->GetEntries())
    {
      cerr << "The trees in " << argv[1] << " and " << argv[3] << " have different numbers of entries."
	   << "Double check that you have the correct files.\n\n";
      exit(-2);
    }

  // Load the branches
  Gen_Event * genData = NULL;
  double true_ze,true_zr,proton_t0;
  double recon_theta, recon_mom, recon_phi, recon_zr;
  int lad_plane;
  genTree->SetBranchAddress("event",&genData);
  genTree->SetBranchAddress("ze",&true_ze);
  genTree->SetBranchAddress("zr",&true_zr);
  genTree->SetBranchAddress("t0",&proton_t0);
  proTree->SetBranchAddress("lad_plane",&lad_plane);
  digTree->SetBranchAddress("mom_recon",&recon_mom);
  digTree->SetBranchAddress("theta_recon",&recon_theta);
  digTree->SetBranchAddress("phi_recon",&recon_phi);
  digTree->SetBranchAddress("z_recon",&recon_zr);

  for (int i=0 ; i<genTree->GetEntries() ; i++)
    {
      if (i%100000==0)
	cerr << "Working on event " << i << "\n";

      genTree->GetEntry(i);
      proTree->GetEntry(i);
      digTree->GetEntry(i);

      // Test if we hit LAD
      if (lad_plane <0)
	continue;

      histZRecon->Fill(recon_zr - true_ze);

      // Reconstruct some important quantities
      TVector3 electron_mom_recon = genData->particles[0].momentum;
      TVector3 proton_mom_recon(recon_mom*sin(recon_theta)*cos(recon_phi),
				recon_mom*sin(recon_theta)*sin(recon_phi),
				recon_mom*cos(recon_theta));
      TVector3 q_recon = TVector3(0.,0.,E1) - electron_mom_recon;
      
      double omega_recon = E1 - electron_mom_recon.Mag();
      double QSq_recon = 2.*E1*electron_mom_recon.Mag()*(1.-cos(electron_mom_recon.Theta()));
      double theta_qs_recon = q_recon.Angle(proton_mom_recon);
      double Es_recon = sqrt(sq(recon_mom) + sq(mP));
      double Wprime_recon = sqrt(sq(mP) - QSq_recon + 2.*omega_recon*(mD - Es_recon) + 2.*recon_mom*q_recon.Mag()*cos(theta_qs_recon));
      double xprime_recon = QSq_recon/(2.*((mD - Es_recon)*omega_recon + recon_mom*q_recon.Mag()*cos(theta_qs_recon)));

      // Restrict our events of interest
      // Cut on QSq
      if (QSq_recon < 2.) continue;
      // Cut on Wprime
      if (Wprime_recon < 1.8) continue;
      // Cut on phi
      if (fabs(recon_phi) > 17.*M_PI/180.) continue;
      // Cut on z
      if (fabs(recon_zr - true_zr) > 3.6) continue;
      // Cut on the theta_qs
      if (theta_qs_recon < 110.*M_PI/180.) continue;
      // Cut on recon_mom
      if (recon_mom < 0.275) continue;
      // Cut on LAD eDep vs timing (to be done later)

      // Fill other histograms
      histMomRes->Fill(genData->particles[1].momentum.Mag(),recon_mom- genData->particles[1].momentum.Mag(),weight);

      // Test if we are in the spectrometer acceptance
      double e_phi = genData->particles[0].momentum.Phi();
      if (e_phi < -M_PI) e_phi += 2.*M_PI;
      
      // SHMS
      if (fabs(e_phi) < shms_acc_phi)
	{
	  shmsHist->Fill(xprime_recon,electron_mom_recon.Theta()*180./M_PI,electron_mom_recon.Mag(),weight);
	}
      
      // HMS
      if (fabs(e_phi-M_PI) < hms_acc_phi)
	{
	  hmsHist->Fill(xprime_recon,electron_mom_recon.Theta()*180./M_PI,electron_mom_recon.Mag(),weight);
	}
    }

  genFile->Close();
  proFile->Close();
  digFile->Close();

  outFile->cd();
  histZRecon->Write();
  histMomRes->Write();
  shmsHist->Write();
  hmsHist->Write();
  outFile->Close();

  return 0;
}
