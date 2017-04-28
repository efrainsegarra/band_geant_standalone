#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <cmath>

const double alpha=0.00729735253;

// Units of GeV
const double mP = 0.93827231;
const double mN = 0.93956536;
const double mD = 1.8756;
const double E1 = 10.9;

// Units of cm / ns
const double cAir = 29.9792458;

// LAD constants
const double lad_min_theta_deg=90.;
const double lad_max_theta_deg=165.;
const double lad_target_z = 25.;
const double lad_angles[3]={102.*M_PI/180.,127.*M_PI/180.,152.*M_PI/180.};
const double lad_radii[3]={564.,472.,564.}; // cm
const double lad_gem_radii[2]={100.,150.};
const double hms_acc_theta=.032;
const double hms_acc_phi=.085;
const double shms_acc_theta=.024;
const double shms_acc_phi=.040;

#endif
