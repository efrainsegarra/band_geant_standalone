#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVectorT.h"

#include "constants.h"
#include "gen_tree.h"
#include "lad_helpers.h"

using namespace std;

double multScatTheta(double xOverX0, double beta, double mom);
void adjustMomVector(TVector3 &vect, double theta);
void hitPlane(TVector3 &pos, const TVector3 &mom, double &path, double &incidence, double rC, double theta);
int whichLADPlane(TVector3 &pos, const TVector3 &mom, double &path);

TRandom3 * myRand;

const double rtd = 180./M_PI;

int main(int argc, char ** argv)
{
  if (argc != 4)
    {
      cerr << "Wrong number of arguments. Instead use\n\tdet_sim /path/to/input/file /path/to/output/file [gemXOverX0]\n";
      exit(-1);
    }

  // Random number generator
  myRand = new TRandom3(0);

  // Get the files and trees
  TFile * inFile = new TFile(argv[1]);
  TFile * outFile = new TFile(argv[2],"RECREATE");
  TTree * inTree = (TTree*)inFile->Get("MCout");
  TTree * outTree = new TTree("PropTree","Propagated Tracks");

  const double gemXOverX0 = atof(argv[3]);

  // Set up the branches
  double zr, ze, gem1_x, gem1_y, gem1_z, gem2_x, gem2_y, gem2_z, lad_x, lad_y, lad_z, path_length, hit_time, t0, mom_at_lad, recon_ze;
  int lad_plane;
  Gen_Event * thisEvent = NULL;
  inTree->SetBranchAddress("ze",&ze);
  inTree->SetBranchAddress("zr",&zr);
  inTree->SetBranchAddress("event",&thisEvent);
  inTree->SetBranchAddress("t0",&t0);
  outTree->Branch("gem1_x",&gem1_x,"gem1_x/D");
  outTree->Branch("gem1_y",&gem1_y,"gem1_y/D");
  outTree->Branch("gem1_z",&gem1_z,"gem1_z/D");
  outTree->Branch("gem2_x",&gem2_x,"gem2_x/D");
  outTree->Branch("gem2_y",&gem2_y,"gem2_y/D");
  outTree->Branch("gem2_z",&gem2_z,"gem2_z/D");
  outTree->Branch("lad_x",&lad_x,"lad_x/D");
  outTree->Branch("lad_y",&lad_y,"lad_y/D");
  outTree->Branch("lad_z",&lad_z,"lad_z/D");
  outTree->Branch("t0",&t0,"t0/D");
  outTree->Branch("path_length",&path_length,"path_length/D");
  outTree->Branch("hit_time",&hit_time,"hit_time/D");
  outTree->Branch("mom_at_lad",&mom_at_lad,"mom_at_lad/D");
  outTree->Branch("recon_ze",&recon_ze,"recon_ze/D");
  outTree->Branch("lad_plane",&lad_plane,"lad_plane/I");

  // Copy the cross section
  TVectorT<double> * csVec;
  if (inFile->GetListOfKeys()->Contains("totalCS"))
    {
      csVec = (TVectorT<double>*)inFile->Get("totalCS");
      gFile=outFile;
      csVec->Write("totalCS");
    }
  else if (inFile->GetListOfKeys()->Contains("totalCSSq"))
    {
      csVec = (TVectorT<double>*)inFile->Get("totalCSSq");  
      gFile=outFile;
      csVec->Write("totalCSSq");
    }

  // Loop over the events
  for (int i=0 ; i<inTree->GetEntries(); i++)
    {
      if (i%100000==0) cerr << "Working on event " << i << "\n";

      // Get new input event
      inTree->GetEvent(i);

      // Set up the initial position and momentum
      path_length=0.;
      TVector3 protonMomVect = thisEvent->particles[1].momentum;
      TVector3 protonPosVect = TVector3(0.,0.,zr);
      double protonPhi = protonMomVect.Phi();
      double protonTheta = protonMomVect.Theta();
      double protonMom = protonMomVect.Mag();
      double protonBeta = protonMom/sqrt(sq(mP) + sq(protonMom));
      recon_ze = ze + myRand->Gaus()*0.3/sin(thisEvent->particles[0].momentum.Theta());
      
      // Adjust for multiple scattering in the target
      //cout << " Target\n";
      adjustMomVector(protonMomVect,multScatTheta((1./769.1 + 0.02/8.897)/(cos(protonPhi)*sin(protonTheta)),protonBeta, protonMom));
     
      // Multiple scattering at the scattering chamber window
      const double scatteringChamberR=30.;
      double tempA = sq(protonMomVect.X()) + sq(protonMomVect.Z());
      double tempB = -2.*sq(protonMomVect.X())*protonPosVect.Z();
      double tempC = sq(protonMomVect.X()*protonPosVect.X()) - sq(scatteringChamberR*protonMomVect.Z());
      double hitZ = (-tempB - sqrt(sq(tempB) - 4.*tempA*tempC))/(2.*tempA);
      double temp = sqrt(sq(scatteringChamberR)-sq(hitZ))/protonMomVect.X();
      path_length+= temp*protonMom;
      protonPosVect = protonPosVect + temp*protonMomVect;
      double incidenceAngle = protonMomVect.Angle(TVector3(protonPosVect.X(),0.,hitZ));
      //cout << " Scattering chamber \n";
      adjustMomVector(protonMomVect, multScatTheta(.0150/3.546/cos(incidenceAngle),protonBeta, protonMom));

      // Multiple scattering at GEM plane 1
      hitPlane(protonPosVect, protonMomVect, path_length, incidenceAngle, lad_gem_radii[0], lad_angles[1]);
      //cout << " GEM 1 \n";
      adjustMomVector(protonMomVect, multScatTheta(gemXOverX0/cos(incidenceAngle),protonBeta, protonMom));
      gem1_x = protonPosVect.X();
      gem1_y = protonPosVect.Y();
      gem1_z = protonPosVect.Z();

      // Multiple scattering at GEM plane 2
      hitPlane(protonPosVect, protonMomVect, path_length, incidenceAngle, lad_gem_radii[1], lad_angles[1]);
      //cout << " GEM 2 \n";
      adjustMomVector(protonMomVect, multScatTheta(gemXOverX0/cos(incidenceAngle),protonBeta, protonMom));
      gem2_x = protonPosVect.X();
      gem2_y = protonPosVect.Y();
      gem2_z = protonPosVect.Z();

      // Figure out what plane of LAD did we hit
      lad_plane = whichLADPlane(protonPosVect,protonMomVect, path_length);
      lad_x = protonPosVect.X();
      lad_y = protonPosVect.Y();
      lad_z = protonPosVect.Z();

      // Set the hit time
      hit_time = path_length/(cAir * protonBeta) + t0;

      // Set up the momentum
      mom_at_lad = protonMom;

      // Fill the tree
      outTree->Fill();
    }

  inFile->Close();
  outFile->cd();
  outTree->Write();
  outFile->Close();

  delete myRand;
  return 0;
}

double multScatTheta(double xOverX0, double beta, double mom)
{
  double theta_space= sqrt(2.*xOverX0)*0.0136/(beta*mom)*(1.+0.038*log(xOverX0));
  //cout << "Beta = " << beta << "   x = " << xOverX0 << "     Typical rotation is " << theta_space * rtd << " degrees\n";
  return theta_space;
}

void adjustMomVector(TVector3 &vect, double theta)
{
  TVector3 unitX = vect.Cross(TVector3(0.,0.,1.));
  unitX *= (1./unitX.Mag());
  TVector3 unitY = vect.Cross(unitX);
  unitY *= (1./unitY.Mag());

  double phi = 2.*M_PI*myRand->Rndm();
  TVector3 axis = unitX*cos(phi) + unitY*sin(phi);
  
  vect.Rotate(theta,axis);
}

void hitPlane(TVector3 &pos, const TVector3 &mom, double &path, double &incidence, double rC, double theta)
{
  double xC = rC * sin(theta);
  double zC = rC * cos(theta);

  double hitZ = (mom.Z()*zC + tan(theta)*(mom.Z()*(xC - pos.X()) + pos.Z()*mom.X()))/(mom.X()*tan(theta) + mom.Z());

  double lambda = (hitZ - pos.Z())/mom.Z();
  
  // Update position to hit position and path length
  path += lambda*mom.Mag();
  pos = pos + lambda*mom;
  
  // Calculate incidence angle
  incidence = mom.Angle(TVector3(xC,0.,zC));
}

int whichLADPlane(TVector3 &pos, const TVector3 &mom, double &path)
{
  TVector3 origPos = pos;
  double origPath = path;
  double incidence = 0.;

  hitPlane(pos,mom,path,incidence,lad_radii[1],lad_angles[1]);
  
  if ((fabs(pos.Y()) < 200.)&&(distInPlane(pos,lad_radii[1],lad_angles[1]) < 121.))
    return 1;

  pos = origPos;
  path=origPath;
  hitPlane(pos,mom,path,incidence,lad_radii[0],lad_angles[0]);
  if ((fabs(pos.Y()) < 200.)&&(distInPlane(pos,lad_radii[0],lad_angles[0]) < 121.))
    return 0;

  pos = origPos;
  path=origPath;
  hitPlane(pos,mom,path,incidence,lad_radii[2],lad_angles[2]);
  if ((fabs(pos.Y()) < 200.)&&(distInPlane(pos,lad_radii[2],lad_angles[2]) < 121.))
    return 2;

  path=0.;
  pos.SetXYZ(0.,0.,0.);
  return -1;
}
