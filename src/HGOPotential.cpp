/*
 * HGOPotential.cpp
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

#include "HGOPotential.h"


double HGOPotential::evaluate(const Vect3d& x1, const Vect3d& u1,
		const Vect3d& x2, const Vect3d& u2) const {
	Vect3d dx = x2-x1;
	double r = dx.norm();
	Vect3d rhat = dx/r;
	if (r <= std::numeric_limits<double>::min()) rhat = Vect3d(1,0,0);

	double udotu = u1.dot(u2);
	double u1dotr = u1.dot(rhat);
	double u2dotr = u2.dot(rhat);
	//std::cout <<"evaluating with x1 = "<<x1<<" u1 = "<<u1<<" x2 = "<<x2<<"u2 = "<<u2<<" sigma = "<<sigma(udotu,u1dotr,u2dotr)<<" u1dotr = "<<u1dotr<<" u2dotr = "<<u2dotr<<" udotu = "<<udotu<<" sigma_s = "<<sigma_s<<std::endl;
	if (r < sigma(udotu,u1dotr,u2dotr)) {
		//return 0.001*std::numeric_limits<double>::max();
		return 1.0;
	} else {
		return 0;
	}
}

double HGOPotential::cut_off() const {
	return 1.01*sigma_s*k;
}

HGOPotential::~HGOPotential() {
	// TODO Auto-generated destructor stub
}

double HGOPotential::sigma(const double udotu, const double u1dotr, const double u2dotr) const {
	const double first = pow(u1dotr + u2dotr,2) / (1.0 + xi*udotu);
	const double second = pow(u1dotr - u2dotr,2) / (1.0 - xi*udotu);
	//std::cout <<" first = "<<first<<" second = "<<second<<" xi = "<<xi<<std::endl;
	//if (isnan(sigma_s*sqrt(1 - (xi/2.0) * (first + second)))) std::cout << "(xi/2.0) * (first + second) = "<<(xi/2.0) * (first + second)<<" 1+xi*udotu "<<1 + xi*udotu<< "1-xi*udotu "<<1 - xi*udotu <<std::endl;
	return sigma_s*pow(1.0 - (xi/2.0) * (first + second),-0.5);
}

