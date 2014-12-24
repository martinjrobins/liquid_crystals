/*
 * LabwohlLasherPotential.cpp
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

#include "LabwohlLasherPotential.h"
#include <math.h>


double LabwohlLasherPotential::evaluate(const Vect3d& x1, const Vect3d& u1,
		const Vect3d& x2, const Vect3d& u2) const {
	Vect3d dx = x2-x1;
	if (dx.norm() < cut_off()) {
		return -epsilon*(1.5*pow(u1.dot(u2),2)-0.5);
		//return 1.0-pow(u1.dot(u2),2);
	} else {
		return 0;
	}
	//return epsilon*(pow(sin(acos(u1.dot(u2))),2));
}

double LabwohlLasherPotential::cut_off() const {
	return lattice_spacing*1.01;
}

LabwohlLasherPotential::~LabwohlLasherPotential() {
	// TODO Auto-generated destructor stub
}

