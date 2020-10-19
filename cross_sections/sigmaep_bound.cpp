#include "sigmaep_bound.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "TVector3.h"
#include "constants.h"

sigmaep_bound::sigmaep_bound(){}

sigmaep_bound::~sigmaep_bound(){}

double sigmaep_bound::sigma_eN(double Ebeam,TVector3 k, TVector3 p)
{
	return sigmaCC1(Ebeam,k,p);
}

double sigmaep_bound::sigmaCC1(double Ebeam, TVector3 k, TVector3 p)
{
	return sigmaCCn(Ebeam, k, p, 1);
}

// This is the opposite from our usual (1,-1,-1,-1) convention
double dot4(double x0, TVector3 x, double y0, TVector3 y)
{
	return ((x*y)-(x0*y0));
}

double sigmaep_bound::sigmaCCn(double Ebeam, TVector3 k, TVector3 p, bool isProton, int n)
{
	TVector3 q = TVector3(0.,0.,Ebeam) - k;
	TVector3 pM = p-q;

	double omega = Ebeam - k.Mag();
	double QSq = q.Mag2() - sq(omega);
	double E = sqrt(p.Mag2() + sq(mN));
	double Ebar = sqrt(pM.Mag2() + sq(mN));
	double omegabar = E-Ebar;
	double QSqbar = q.Mag2() - sq(omegabar);

	// Calculate form factors
	double GE = GEp(QSq);
	double GM = GMp(QSq);

	double F1 = (GE + GM*QSq/(4.*sq(mN))) / (1. + QSq/(4.*sq(mN)));
	double kF2 = (GM - GE)/(1. + QSq/(4.*sq(mN)));

	double wC;
	double wT;
	double wS;
	double wI;

	if (n==1)
	{
		wC = (sq(E+Ebar)*(sq(F1) + QSqbar/(4.*mN*mN) * sq(kF2)) - q.Mag2()*sq(F1 + kF2))/(4.*E*Ebar);
		wT = QSqbar*sq(F1 + kF2)/(2.*Ebar*E);
		wS = p.Mag2() * sq(sin(p.Angle(q))) * (sq(F1) + QSqbar/(4.*mN*mN) * sq(kF2))/(E*Ebar);
		wI = -p.Mag()*sin(p.Angle(q))*(Ebar + E)*(sq(F1) + QSqbar/(4.*mN*mN) * sq(kF2))/(E*Ebar);
	}
	else if (n==2){  
		double pbarp = dot4(Ebar,pM,E,p);
		double pbarq = dot4(Ebar,pM,omega,q);
		double pq = dot4(E,p,omega,q);
		double qbarq = dot4(omegabar,q,omega,q);
		double sumq = dot4((Ebar+E),(pM+p),omega,q);

		wC = ((E*Ebar + 0.5 * (pbarp + sq(mN))) * sq(F1)
				- 0.5 * q.Mag2() * F1 * kF2
				- ((pbarq*E + pq*Ebar)*omega
					- Ebar * E * QSq
					+ pbarq * pq
					- 0.5 * (pbarp - sq(mN)) * q.Mag2())
				* sq(kF2)/(4*sq(mN)))
			/(E*Ebar);
		wT = (-(pbarp + sq(mN)) * sq(F1)
				+ qbarq * F1 * kF2
				+ (2*pbarq*pq - (pbarp - sq(mN))*QSq)
				* sq(kF2)/(4*sq(mN)))
			/(Ebar*E);
		wS = p.Mag2() * sq(sin(p.Angle(q))) * (sq(F1) + QSq * sq(kF2) / (4*sq(mN)) ) / (E*Ebar);
		wI = p.Mag() * sin(p.Angle(q)) * (-(Ebar + E) * sq(F1)
				+ (sumq * omega - (Ebar + E) * QSq)
				* sq(kF2)/(4*sq(mN)))
			/(E*Ebar);
	}

	double sigmaMott = cmSqGeVSq * k.Mag2() * sq( 2. * alpha * cos(k.Theta()/2.) / QSq);

	double cosPhi = cos(q.Cross(k).Angle( q.Cross(p) ));
	return sigmaMott *   ( sq(QSq)/sq(q.Mag2()) * wC +
			(QSq/(2.*q.Mag2()) + sq(tan(k.Theta()/2.))) * wT +
			QSq/q.Mag2() * sqrt(QSq/q.Mag2() + sq(tan(k.Theta()/2.))) * wI * cosPhi +
			(QSq/q.Mag2() * sq(cosPhi) + sq(tan(k.Theta()/2.))) * wS
			);
}


double sigmaep_bound::GEp(double QSq)
{
	return Gdipole(QSq);
}


double sigmaep_bound::GMp(double QSq)
{
	return mu_p * Gdipole(QSq);
}

double sigmaep_bound::Gdipole(double QSq){ return 1. / sq(1 + QSq/0.71); }

