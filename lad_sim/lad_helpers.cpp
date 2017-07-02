// Helper functions for LAD
#include "lad_helpers.h"

#include <cmath>

inline double sq(double x)
{
  return x*x;
}

double distInPlane(const TVector3 &pos, double rC, double theta)
{
  double xC = rC * sin(theta);
  double zC = rC * cos(theta);

  double xH = pos.X();
  double zH = pos.Z();

  return sqrt(sq(xC-xH) + sq(zC-zH));
}
