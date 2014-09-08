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
			xi_dash = (pow(kdash,1) - 1) / (pow(kdash,1) + 1);

		} else {
			xi_dash = (pow(kdash,1.0/mu) - 1) / (pow(kdash,1.0/mu) + 1);
		}
	}
	double evaluate(Vect3d &x1, Vect3d &u1, Vect3d &x2, Vect3d &u2);
	double cut_off();
	virtual ~GayBernePotential();
private:
	double epsilon(double udotu, double u1dotr, double u2dotr);
	double epsilon(double udotu);
	double epsilon_dash(double udotu, double u1dotr, double u2dotr);
	double qfunc(double udotu, double u1dotr, double u2dotr, double r);
	double sigma(double udotu, double u1dotr, double u2dotr);
	double k,kdash,mu,nu,xi,xi_dash,sigma_s,epsilon_0;
};

#endif /* GAYBERNEPOTENTIAL_H_ */
