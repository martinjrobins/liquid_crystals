/*
 * GayBernePotential.h
 *
 *  Created on: 5 Sep 2014
 *      Author: mrobins
 */

#ifndef GAYBERNEPOTENTIAL_H_
#define GAYBERNEPOTENTIAL_H_

#include "Aboria.h"

using namespace Aboria;

class GayBernePotential {
public:
	GayBernePotential(double sigma_s, double k, double kdash, double mu, double nu, double epsilon_0):
		k(k),kdash(kdash),mu(mu),nu(nu),sigma_s(sigma_s),epsilon_0(epsilon_0) {
		xi = (pow(k,2) - 1) / (pow(k,2) + 1);
		if (mu == 1) {
			xi_dash = (1 - pow(kdash,1)) / (1 + pow(kdash,1));
		} else {
			xi_dash = (1 - pow(kdash,1.0/mu)) / (1 + pow(kdash,1.0/mu));
		}
	}
	double evaluate(const Vect3d &x1, const Vect3d &u1, const Vect3d &x2, const Vect3d &u2) const;
	double operator()(const Vect3d &x1, const Vect3d &u1, const Vect3d &x2, const Vect3d &u2) const {
		return evaluate(x1,u1,x2,u2);
	}
	double cut_off() const;
	virtual ~GayBernePotential();
private:
	double epsilon(const double udotu, const double u1dotr, const double u2dotr) const;
	double epsilon(const double udotu) const;
	double epsilon_dash(const double udotu, const double u1dotr, const double u2dotr) const;
	double qfunc(const double udotu, const double u1dotr, const double u2dotr, const double r) const;
	double sigma(const double udotu, const double u1dotr, const double u2dotr) const;
	double k,kdash,mu,nu,xi,xi_dash,sigma_s,epsilon_0;
};

#endif /* GAYBERNEPOTENTIAL_H_ */
