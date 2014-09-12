/*
 * HGOPotential.h
 * 
 * Copyright 2014 Martin Robinson
 *
 * This file is part of liquid_crystals.
 *
 * liquid_crystals is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * liquid_crystals is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with liquid_crystals.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 11 Sep 2014
 *      Author: robinsonm
 */

#ifndef HGOPOTENTIAL_H_
#define HGOPOTENTIAL_H_

#include "Aboria.h"

using namespace Aboria;

class HGOPotential {
public:
	HGOPotential(double sigma_s, double k):
		k(k),sigma_s(sigma_s) {
		xi = (pow(k,2) - 1) / (pow(k,2) + 1);
	}
	double evaluate(const Vect3d &x1, const Vect3d &u1, const Vect3d &x2, const Vect3d &u2) const;
	double operator()(const Vect3d &x1, const Vect3d &u1, const Vect3d &x2, const Vect3d &u2) const {
		return evaluate(x1,u1,x2,u2);
	}
	double cut_off() const;
	virtual ~HGOPotential();
private:
	double sigma(const double udotu, const double u1dotr, const double u2dotr) const;
	double k,xi,sigma_s;
};

#endif /* HGOPOTENTIAL_H_ */
