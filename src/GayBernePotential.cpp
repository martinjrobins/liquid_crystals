/*
 * GayBernePotential.cpp
 *
 *  Created on: 5 Sep 2014
 *      Author: mrobins
 */

#include "GayBernePotential.h"
#include <limits>


double GayBernePotential::evaluate(const Vect3d &x1, const Vect3d &u1, const Vect3d &x2, const Vect3d &u2) const {
	Vect3d dx = x2-x1;
	double r = dx.norm();
	Vect3d rhat = dx/r;
	if (r==0) rhat = Vect3d(1,0,0);
	const double trunc = 0.000001;
	if (r/sigma_s < trunc) {
		dx = rhat*sigma_s*trunc;
		r = trunc*sigma_s;
	}

	double udotu = u1.dot(u2);
	double u1dotr = u1.dot(rhat);
	double u2dotr = u2.dot(rhat);
	double q = qfunc(udotu,u1dotr,u2dotr,r);

	double val = 4*epsilon(udotu,u1dotr,u2dotr)*(pow(q,12) - pow(q,6));
	//if (isnan(val)) std::cout << "potential = "<<val<<" at u1 = "<<u1<<" u2 = "<<u2<<" udotu = "<<udotu<<" r = "<<r<<" q = "<<q<<" sigma_s = "<<sigma_s<<" sigma = "<<sigma(udotu,u1dotr,u2dotr)<<" xi = "<<xi<<std::endl;
	return 4*epsilon(udotu,u1dotr,u2dotr)*(pow(q,12) - pow(q,6));
}

double GayBernePotential::cut_off() const {
	return 2.5*sigma_s*k;
}

GayBernePotential::~GayBernePotential() {
	// TODO Auto-generated destructor stub
}

double GayBernePotential::epsilon(const double udotu, const double u1dotr, const double u2dotr) const {
	return epsilon_0 * pow(epsilon_dash(udotu,u1dotr,u2dotr),mu) / pow(epsilon(udotu),nu);
}

double GayBernePotential::epsilon_dash(const double udotu, const double u1dotr, const double u2dotr) const {
	double first = pow(u1dotr + u2dotr,2) / (1 + xi_dash*udotu);
	double second = pow(u1dotr - u2dotr,2) / (1 - xi_dash*udotu);
	return 1 - (xi_dash/2.0)*(first + second);
}

double GayBernePotential::qfunc(const double udotu, const double u1dotr, const double u2dotr, const double r) const {
	return sigma_s /(r - sigma(udotu,u1dotr,u2dotr) + sigma_s);
}

double GayBernePotential::epsilon(const double udotu) const {
	return sqrt(1.0-pow(xi,2) * pow(udotu,2));
}

double GayBernePotential::sigma(const double udotu, const double u1dotr, const double u2dotr) const {
	double first = pow(u1dotr + u2dotr,2) / (1 + xi*udotu);
	double second = pow(u1dotr - u2dotr,2) / (1 - xi*udotu);
	//if (isnan(sigma_s*sqrt(1 - (xi/2.0) * (first + second)))) std::cout << "(xi/2.0) * (first + second) = "<<(xi/2.0) * (first + second)<<" 1+xi*udotu "<<1 + xi*udotu<< "1-xi*udotu "<<1 - xi*udotu <<std::endl;
	return sigma_s/sqrt(1 - (xi/2.0) * (first + second));
}


