#include <iostream>
#include <cmath>
#include <cstdlib>

#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVector3.h"
#include "TVectorT.h"

#include "lad_helpers.h"
#include "gen_tree.h"
#include "constants.h"

using namespace std;

const double loXAng = 12.5*M_PI/180.;
const double loXMom = 4.8;
const double hiXAng = 14.5*M_PI/180.;
const double hiXMom = 5.2;

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
  const double lumi=1.E-6; // 1E36 cm^-2 s-1 in units of nb^-1 ns^-1
  const double duration=3600.*1.E9; // one hour of run time in ns
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
  TH2D * histZRecon = new TH2D("z_res","Vertex resolution;theta_e [deg.];zp_recon - ze_recon [cm];Counts",100,4.,24.,100,-10.,10.);
  TH2D * histZPRecon_byTheta = new TH2D("zp_res_by_theta","Vertex resolution;theta_p [deg.];zp_recon - ze_recon [cm];Counts",
				       100.,lad_min_theta_deg,lad_max_theta_deg,100,-10.,10.);
  TH1D * histZPRecon = new TH1D("zp_res","Proton Vertex resolution;zp_recon - ze_gen [cm];Counts",100,-10.,10.);
  TH1D * histT0True = new TH1D("t0","Events of interest;true t0 [ns];Counts",100,-20,20.);
  TH2D * histEdepRes = new TH2D("edep_res","Events hitting lad;Momentum (tof);Mom_edep - mom_tof;Counts",100,0.,1.,100.,-0.1,0.1);
  TH2D * histMomRes = new TH2D("mom_res","Recoil protons;true mom [GeV];Recon - true [GeV];Counts",100,0.,1.,100,-0.05,0.05);
  TH2D * histLAD_y = new TH2D("lad_y","Good events;LAD Bar; y_pos ; Counts",33,-0.5,32.5,84,-210.,210.);
  TH2D * histGenPhi = new TH2D("gen_phi","Good events;LAD Bar;Proton Phi [deg];Counts",33,-0.5,32.5,200,-lad_max_phi*180./M_PI,lad_max_phi*180./M_PI);
  TH3D * hmsHist = new TH3D("hms","hms_acceptance;x_prime;angle;momentum",100,0.,1.,200,4.,24.,200,0.,10.0);
  TH3D * shmsHist = new TH3D("shms","shms_acceptance;x_prime;angle;momentum",100,0.,1.,200,4.,24.,200,0.,10.0);

  // Some rate plots
  TH2D * histLoXRate = new TH2D("loX_rate","low x' settings;x';alpha_s; counts/hr",20,0.,1.,8,1.15,1.55);
  TH2D * histHiXRate = new TH2D("hiX_rate","high x' settings;x';alpha_s; counts/hr",20,0.,1.,8,1.15,1.55);

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

  // Normalize the weight by the total number of generated events.
  weight = weight/((double)genTree->GetEntries());
  cerr << "The weight for this run is " << weight << "\n";

  // Load the branches
  Gen_Event * genData = NULL;
  double true_ze,true_zr;
  double lad_x,lad_y,lad_z;
  double recon_theta, recon_mom, recon_phi, recon_zr, recon_mom_from_edep, recon_ze;
  int lad_plane, lad_bar;
  genTree->SetBranchAddress("event",&genData);
  genTree->SetBranchAddress("ze",&true_ze);
  genTree->SetBranchAddress("zr",&true_zr);
  proTree->SetBranchAddress("lad_plane",&lad_plane);
  proTree->SetBranchAddress("lad_x",&lad_x);
  proTree->SetBranchAddress("lad_y",&lad_y);
  proTree->SetBranchAddress("lad_z",&lad_z);
  digTree->SetBranchAddress("lad_bar",&lad_bar);
  digTree->SetBranchAddress("mom_recon",&recon_mom);
  digTree->SetBranchAddress("theta_recon",&recon_theta);
  digTree->SetBranchAddress("phi_recon",&recon_phi);
  digTree->SetBranchAddress("z_recon",&recon_zr);
  digTree->SetBranchAddress("ze_recon",&recon_ze);
  digTree->SetBranchAddress("mom_from_edep_recon",&recon_mom_from_edep);

  // Do preliminary loop to fit vertex position
  cerr << "Doing preliminary loop through events to establish vertex cuts...\n";
  for (int i=0 ; i<genTree->GetEntries() ; i++)
    {
      genTree->GetEntry(i);
      proTree->GetEntry(i);
      digTree->GetEntry(i);

      // Test if we hit LAD
      if (lad_plane <0)
	continue;

      double theta_r_deg = genData->particles[1].momentum.Theta()*180./M_PI;

      histZPRecon_byTheta->Fill(theta_r_deg,recon_zr - true_zr,weight);
      
      if ((theta_r_deg < 105) && ( theta_r_deg > 90))
	histZPRecon->Fill(recon_zr - true_zr,weight);
    }
  // Do a fit to establish proton vertex reconstruction resolution
  histZPRecon->Fit("gaus");
  double zr_width=histZPRecon->GetFunction("gaus")->GetParameter(2);

  // Do a second loop through, once the width of the vertex resolution has been established
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

      histEdepRes->Fill(recon_mom,recon_mom_from_edep - recon_mom,weight);

      // Reconstruct some important quantities
      double proton_t0 = genData->particles[1].t0;
      TVector3 electron_mom_recon = genData->particles[0].momentum;
      TVector3 proton_mom_recon(recon_mom*sin(recon_theta)*cos(recon_phi),
				recon_mom*sin(recon_theta)*sin(recon_phi),
				recon_mom*cos(recon_theta));
      TVector3 q_recon = TVector3(0.,0.,E1) - electron_mom_recon;

      histZRecon->Fill(electron_mom_recon.Theta()*180./M_PI,recon_zr - recon_ze,weight);
      
      double omega_recon = E1 - electron_mom_recon.Mag();
      double QSq_recon = 2.*E1*electron_mom_recon.Mag()*(1.-cos(electron_mom_recon.Theta()));
      double theta_qs_recon = q_recon.Angle(proton_mom_recon);
      double Es_recon = sqrt(sq(recon_mom) + sq(mP));
      double Wprime_recon = sqrt(sq(mD) + sq(mP) - QSq_recon + 2.*omega_recon*(mD - Es_recon) - 2.*mD*Es_recon + 2.*proton_mom_recon.Dot(q_recon));
      double xprime_recon = QSq_recon/(2.*((mD - Es_recon)*omega_recon + proton_mom_recon.Dot(q_recon)));

      // True quantities
      TVector3 electron_mom_true = genData->particles[0].momentum;
      TVector3 proton_mom_true = genData->particles[1].momentum;
      TVector3 q_true = TVector3(0.,0.,E1) - electron_mom_true;
      double omega_true = E1 - electron_mom_true.Mag();
      double QSq_true = 2.*E1*electron_mom_true.Mag()*(1.-cos(electron_mom_true.Theta()));
      double theta_qs_true = q_true.Angle(proton_mom_true);
      double Es_true = sqrt(proton_mom_true.Mag2() + sq(mP));
      double Wprime_true = sqrt(sq(mD) + sq(mP) - QSq_true + 2.*omega_true*(mD - Es_true) -2.*mD*Es_true + 2.*proton_mom_true.Dot(q_true));
      double xprime_true = QSq_true/(2.*((mD - Es_true)*omega_true + proton_mom_true.Dot(q_true)));
      double alpha_s_true = (Es_true - proton_mom_true.Dot(q_true)/q_true.Mag())/mP;

      // Restrict our events of interest
      // Cut on QSq
      if (QSq_recon < 2.) continue;
      //if (QSq_true < 2.) continue;
      // Cut on Wprime
      //if (Wprime_recon < 1.8) continue;
      if (Wprime_recon < 2.) continue;
      //if (Wprime_true < 2.) continue;
      // Cut on the theta_qs
      if (theta_qs_recon < 110.*M_PI/180.) continue;
      //if (theta_qs_true < 110.*M_PI/180.) continue;
      // Cut on recon_mom
      if (recon_mom < 0.275) continue;
      //if (proton_mom_true.Mag() < 0.275) continue;

      // Background reduction cuts
      // Cut on z
      if (fabs(recon_zr - recon_ze) > 2.*sqrt(sq(zr_width/sin(recon_theta)) + sq(0.3/sin(electron_mom_recon.Theta())))) continue;
      // Cut on LAD eDep vs timing (to be made better later)
      if (fabs(recon_mom_from_edep - recon_mom)>0.0404) continue;

      // Fill other histograms
      histT0True->Fill(proton_t0,weight);
      histMomRes->Fill(genData->particles[1].momentum.Mag(),recon_mom- genData->particles[1].momentum.Mag(),weight);
      histLAD_y->Fill(11*lad_plane + lad_bar,lad_y,weight);
      histGenPhi->Fill(11*lad_plane + lad_bar,proton_mom_true.Phi()*180/M_PI);

      // Test if we are in the spectrometer acceptance
      double e_phi = genData->particles[0].momentum.Phi();
      if (e_phi < -0.5*M_PI) e_phi += 2.*M_PI;
      
      // SHMS
      if (fabs(e_phi) < 0.5*shms_acc/(2.*shms_acc_theta*sin(electron_mom_true.Theta())))
	{
	  shmsHist->Fill(xprime_recon,electron_mom_recon.Theta()*180./M_PI,electron_mom_recon.Mag(),weight);
	  
	  // Test if we are in the nominal low x' spectrometer acceptance
	  if ((fabs(electron_mom_true.Theta() - loXAng)<shms_acc_theta) &&
	      (electron_mom_true.Mag() < loXMom*shms_acc_mom_hi) &&
	      (electron_mom_true.Mag() > loXMom/shms_acc_mom_lo))
	    {
	      histLoXRate->Fill(xprime_true,alpha_s_true,weight);
	    }
	  
	  // Test if we are in the nominal high x' spectrometer acceptance
	  if ((fabs(electron_mom_true.Theta() - hiXAng)<shms_acc_theta) &&
	      (electron_mom_true.Mag() < hiXMom*shms_acc_mom_hi) &&
	      (electron_mom_true.Mag() > hiXMom/shms_acc_mom_lo))
	    {
	      histHiXRate->Fill(xprime_true,alpha_s_true,weight);
	    }
	}
      
      // HMS
      if (fabs(e_phi-M_PI) < 0.5*hms_acc/(2.*hms_acc_theta*sin(electron_mom_true.Theta())))
	{
	  hmsHist->Fill(xprime_recon,electron_mom_recon.Theta()*180./M_PI,electron_mom_recon.Mag(),weight);

	  // Test if we are in the nominal low x' spectrometer acceptance
	  if ((fabs(electron_mom_true.Theta() - loXAng)<hms_acc_theta) &&
	      (electron_mom_true.Mag() < loXMom*hms_acc_mom) &&
	      (electron_mom_true.Mag() > loXMom/hms_acc_mom))
	    {
	      histLoXRate->Fill(xprime_true,alpha_s_true,weight);
	    }

	  // Test if we are in the nominal high x' spectrometer acceptance
	  if ((fabs(electron_mom_true.Theta() - hiXAng)<hms_acc_theta) &&
	      (electron_mom_true.Mag() < hiXMom*hms_acc_mom) &&
	      (electron_mom_true.Mag() > hiXMom/hms_acc_mom))
	    {
	      histHiXRate->Fill(xprime_true,alpha_s_true,weight);
	    }
	}

    }

  genFile->Close();
  proFile->Close();
  digFile->Close();

  outFile->cd();
  histZRecon->Write();
  histZPRecon->Write();
  histZPRecon_byTheta->Write();
  histT0True->Write();
  histEdepRes->Write();
  histMomRes->Write();
  histLAD_y->Write();
  histGenPhi->Write();
  histLoXRate->Write();
  histHiXRate->Write();
  shmsHist->Write();
  hmsHist->Write();
  outFile->Close();

  return 0;
}
