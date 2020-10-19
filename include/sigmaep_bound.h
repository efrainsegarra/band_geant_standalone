#ifndef __SIGMAEP_BOUND_H__
#define __SIGMAEP_BOUND_H__

#include "TVector3.h"

const double mu_p=2.79;
const double mu_n=-1.91;

class sigmaep_bound
{
 public:
  sigmaep_bound();
  ~sigmaep_bound();
  double sigma_eN(double Ebeam, TVector3 k, TVector3 p, bool isProton);
  double sigmaCCn(double Ebeam, TVector3 k, TVector3 p, bool isProton, int n);
  double sigmaCC1(double Ebeam, TVector3 k, TVector3 p, bool isProton);
  double GEp(double QSq);
  double GMp(double QSq);

 private:
  static double Gdipole(double QSq);
  
};

#endif
