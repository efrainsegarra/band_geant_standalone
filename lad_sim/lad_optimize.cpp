#include <cmath>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TBox.h"

#include "lad_helpers.h"
#include "constants.h"

using namespace std;

const double rtd = 180./M_PI;
const int xThreshold = 51;

// Spectrometer acceptances: tunable based on which spec is being optimized
double theta_acc, mom_acc_lo, mom_acc_hi;

// Helper functions
void getBounds(double theta_c, double mom_c, double &theta_lo, double &theta_hi, double &mom_lo, double &mom_hi);
double integrateSpec(double thetaDeg_c, double mom_c, TH2D* hist);
double integrate2D(double loThDeg, double hiThDeg, double loMom, double hiMom, TH2D * hist);
void sigAndBkg(double thetaDeg_c, double mom_c, TH2D * sigHist, TH2D * bkgHist, double &sig, double &bkg);
double statError(double theta_c, double mom_c, TH2D * sig, TH2D * bkg);
void optimizeSpecPos(TH2D * sig, TH2D * bkg, double &theta_c, double &mom_c);

int main(int argc, char ** argv)
{
  if (argc != 4)
    {
      cerr << "Wrong number of arguments. Instead use: \n"
	   << "\tlad_optimize /path/to/signal/file /path/to/bkg/file /path/to/output\n";
      return -1;
    }

  TFile * sigFile = new TFile(argv[1]);
  TH3D * sigHMS = (TH3D*)sigFile->Get("hms");
  TH3D * sigSHMS = (TH3D*)sigFile->Get("shms");
  TFile * bkgFile = new TFile(argv[2]);
  TH3D * bkgHMS = (TH3D*)bkgFile->Get("hms");
  TH3D * bkgSHMS = (TH3D*)bkgFile->Get("shms");

  // Test to see which histogram has more counts
  TH3D * sigHist, *bkgHist;
  if (sigHMS->Integral() > sigSHMS->Integral())
    {
      sigHist = sigHMS;
      bkgHist = bkgHMS;
      theta_acc = hms_acc_theta*rtd;
      mom_acc_lo = hms_acc_mom;
      mom_acc_hi = hms_acc_mom;
    }
  else
    {
      sigHist = sigSHMS;
      bkgHist = bkgSHMS;
      mom_acc_lo = shms_acc_mom_lo;
      mom_acc_hi = shms_acc_mom_hi;
      theta_acc = shms_acc_theta*rtd;
    }

  TFile * outFile = new TFile(argv[3],"RECREATE");

  // Get the High-x projections
  sigHist->GetXaxis()->SetRange(xThreshold,100);
  TH2D * sigHiX = (TH2D*)sigHist->Project3D("sigHiX_zy");
  sigHiX->Write();

  bkgHist->GetXaxis()->SetRange(xThreshold,100);
  TH2D * bkgHiX = (TH2D*)bkgHist->Project3D("bkgHiX_zy");
  bkgHiX->Write();

  // Get the Low-x projections
  sigHist->GetXaxis()->SetRange(26,35);
  TH2D * sigLoX = (TH2D*)sigHist->Project3D("sigLoX_zy");
  sigLoX->Write();

  bkgHist->GetXaxis()->SetRange(26,35);
  TH2D * bkgLoX = (TH2D*)bkgHist->Project3D("bkgLoX_zy");
  bkgLoX->Write();

  // First, let's optimize the high-x settings
  double optTheta = 17;
  double optMom = 4.4;
  double loTheta,hiTheta,loMom,hiMom;
  getBounds(optTheta,optMom,loTheta,hiTheta,loMom,hiMom);
  TBox * origHiBox = new TBox(loTheta,loMom,hiTheta,hiMom);
  optimizeSpecPos(sigHiX,bkgHiX,optTheta,optMom);
  getBounds(optTheta,optMom,loTheta,hiTheta,loMom,hiMom);
  TBox * optHiBox = new TBox(loTheta,loMom,hiTheta,hiMom);

  TCanvas * hiXCan = new TCanvas("hiX");
  hiXCan->cd();
  sigHiX->Draw("COL");

  origHiBox->SetFillStyle(0);
  origHiBox->SetLineColor(1);
  origHiBox->SetLineWidth(4);

  optHiBox->SetFillStyle(0);
  optHiBox->SetLineColor(2);
  optHiBox->SetLineWidth(4);

  origHiBox->Draw("SAME");
  optHiBox->Draw("SAME");
  hiXCan->Write();

  cout << "Optimum spectrometer position for High-x kinematics: \n" 
       << optTheta << " [deg]\n"
       << optMom << " [GeV]\n";


  // Now let's update the low x settings
  TCanvas * loXCan = new TCanvas("loX");
  optTheta=13.5;
  optMom=4.4;
  getBounds(optTheta,optMom,loTheta,hiTheta,loMom,hiMom);
  TBox * origLoBox = new TBox(loTheta,loMom,hiTheta,hiMom);
  optimizeSpecPos(sigLoX,bkgLoX,optTheta,optMom);
  getBounds(optTheta,optMom,loTheta,hiTheta,loMom,hiMom);
  TBox * optLoBox = new TBox(loTheta,loMom,hiTheta,hiMom);
  loXCan->cd();
  sigLoX->Draw("COL");

  origLoBox->SetFillStyle(0);
  origLoBox->SetLineColor(1);
  origLoBox->SetLineWidth(4);

  optLoBox->SetFillStyle(0);
  optLoBox->SetLineColor(2);
  optLoBox->SetLineWidth(4);

  origLoBox->Draw("SAME");
  optLoBox->Draw("SAME");
  loXCan->Write();

  cout << "Optimum spectrometer position for Low-x kinematics: \n" 
       << optTheta << " [deg]\n"
       << optMom << " [GeV]\n";

  outFile->Close();

  return 0;
}

double integrate2D(double loThDeg, double hiThDeg, double loMom, double hiMom, TH2D * hist)
{
  if (loThDeg<4.) loThDeg=4.;
  if (hiThDeg>24.) hiThDeg=24.;
  if (loMom<0.) loMom=0.;
  if (hiMom>10.) hiMom=10.;

  int loMomBin = hist->GetYaxis()->FindBin(loMom);
  int hiMomBin = hist->GetYaxis()->FindBin(hiMom);
  int loThetaBin = hist->GetXaxis()->FindBin(loThDeg);
  int hiThetaBin = hist->GetXaxis()->FindBin(hiThDeg);
  
  // Integrate the bulk
  double integral=hist->Integral(loThetaBin+1,hiThetaBin-1,loMomBin+1,hiMomBin-1);
  
  // Integrate the borders
  double loThetaFrac = (hist->GetXaxis()->GetBinUpEdge(loThetaBin) - loThDeg)/
    (hist->GetXaxis()->GetBinWidth(loThetaBin));
  double hiThetaFrac = (hiThDeg - hist->GetXaxis()->GetBinLowEdge(hiThetaBin))/
    (hist->GetXaxis()->GetBinWidth(hiThetaBin));
  double loMomFrac = (hist->GetYaxis()->GetBinUpEdge(loMomBin) - loMom)/
    (hist->GetYaxis()->GetBinWidth(loMomBin));
  double hiMomFrac = (hiMom - hist->GetYaxis()->GetBinLowEdge(hiMomBin))/
    (hist->GetYaxis()->GetBinWidth(hiMomBin));

  // Start with the left border
  integral += loThetaFrac * hist->Integral(loThetaBin,loThetaBin,loMomBin+1,hiMomBin-1);
  
  // Lower border
  integral += loMomFrac * hist->Integral(loThetaBin+1,hiThetaBin-1,loMomBin,loMomBin);

  // Right border
  integral += hiThetaFrac * hist->Integral(hiThetaBin,hiThetaBin,loMomBin+1,hiMomBin-1);

  // Upper border
  integral += hiMomFrac * hist->Integral(loThetaBin+1,hiThetaBin-1,hiMomBin,hiMomBin);

  // Now integrate up the four corners
  integral += ( loMomFrac * loThetaFrac * hist->Integral(loThetaBin,loThetaBin,loMomBin,loMomBin)
		+ hiMomFrac * loThetaFrac * hist->Integral(loThetaBin,loThetaBin,hiMomBin,hiMomBin)
		+ loMomFrac * hiThetaFrac * hist->Integral(hiThetaBin,hiThetaBin,loMomBin,loMomBin)
		+ hiMomFrac * hiThetaFrac * hist->Integral(hiThetaBin,hiThetaBin,hiMomBin,hiMomBin));
		
  return integral;
}

void getBounds(double theta_c, double mom_c, double &theta_lo, double &theta_hi, double &mom_lo, double &mom_hi)
{
  mom_hi = mom_acc_hi*mom_c;
  mom_lo = mom_c/mom_acc_lo;
  theta_lo = theta_c - theta_acc;
  theta_hi = theta_c + theta_acc;
}

double integrateSpec(double thetaDeg_c, double mom_c, TH2D* hist)
{
  double hiMom,loMom,hiThDeg,loThDeg;
  getBounds(thetaDeg_c,mom_c,loThDeg,hiThDeg,loMom,hiMom);
  return integrate2D(loThDeg,hiThDeg,loMom,hiMom,hist);
}

void sigAndBkg(double thetaDeg_c, double mom_c, TH2D * sigHist, TH2D * bkgHist, double &sig, double &bkg)
{
  sig = integrateSpec(thetaDeg_c, mom_c, sigHist);
  bkg = integrateSpec(thetaDeg_c, mom_c, bkgHist);
}

double statError(double theta_c, double mom_c, TH2D * sig, TH2D * bkg)
{
  double sigCounts, bkgCounts;
  sigAndBkg(theta_c,mom_c,sig,bkg,sigCounts,bkgCounts);

  if (sigCounts<=0.)
    return 1.E10;

  double stat = sqrt(sigCounts+bkgCounts)/sigCounts;
  
  return stat;
}

void optimizeSpecPos(TH2D * sig, TH2D * bkg, double &theta_c, double &mom_c)
{
  double best_theta = theta_c;
  double best_mom = mom_c;
  double best_stat = statError(theta_c,mom_c,sig,bkg);

  for (double theta_guess=4.;theta_guess <24. ; theta_guess+=0.1)
    {
      //cout << "Checking angle " << theta_guess << " deg.\n";

      for (double mom_guess=0.5 ; mom_guess <9. ; mom_guess+=0.01)
	{
	  double stat = statError(theta_guess,mom_guess,sig,bkg);

	  //cout << theta_guess << " deg.    " << mom_guess << " GeV      stat: " << stat
	  //     << "   best: " << best_stat << "\n";
	  
	  if (stat < best_stat)
	    {
	      best_stat=stat;
	      best_theta = theta_guess;
	      best_mom =mom_guess;
	    }	  
	}
    }
  
  theta_c=best_theta;
  mom_c=best_mom;
}
