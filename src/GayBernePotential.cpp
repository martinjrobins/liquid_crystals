/*
 * GayBernePotential.cpp
 *
 *  Created on: 5 Sep 2014
 *      Author: mrobins
 */

#include "GayBernePotential.h"


double GayBernePotential::evaluate(Vect3d& x1, Vect3d& u1, Vect3d& x2,
		Vect3d& u2) {
	Vect3d dx = x2-x1;
	double r = dx.norm();
	Vect3d rhat = dx/r;
	double udotu = u1.dot(u2);
	double u1dotr = u1.dot(rhat);
	double u2dotr = u2.dot(rhat);
	double q = qfunc(udotu,u1dotr,u2dotr,r);
	return 4*epsilon(udotu,u1dotr,u2dotr)(pow(q,12) - pow(q,6));
}

double GayBernePotential::cut_off() {
	return 2.5*sigma_s*k;
}

GayBernePotential::~GayBernePotential() {
	// TODO Auto-generated destructor stub
}

double GayBernePotential::epsilon(double udotu, double u1dotr, double u2dotr) {
	return epsilon_0 * pow(epsilon_dash(udotu,u1dotr,u2dotr),mu) * pow(epsilon(udotu),nu);
}

double GayBernePotential::epsilon_dash(double udotu, double u1dotr, double u2dotr) {
	double first = pow(u1dotr + u2dotr,2) / (1 + xi_dash*udotu);
	double second = pow(u1dotr - u2dotr,2) / (1 - xi_dash*udotu);
	return 1 - (xi_dash/2.0)*(first + second);
}

double GayBernePotential::qfunc(double udotu, double u1dotr, double u2dotr, double r) {
	return sigma_s /(r - sigma(udotu,u1dotr,u2dotr) + sigma_s);
}

double GayBernePotential::epsilon(double udotu) {
	return sqrt(1.0-pow(xi,2) * pow(udotu,2));
}

double GayBernePotential::sigma(double udotu, double u1dotr, double u2dotr) {
	double first = pow(u1dotr + u2dotr,2) / (1 + xi*udotu);
	double second = pow(u1dotr - u2dotr,2) / (1 - xi*udotu);
	return sigma_s*sqrt(1 - (xi/2.0) * (first + second));
}


