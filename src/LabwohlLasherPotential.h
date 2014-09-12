/*
 * LabwohlLasherPotential.h
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

#ifndef LABWOHLLASHERPOTENTIAL_H_
#define LABWOHLLASHERPOTENTIAL_H_

#include "Aboria.h"

using namespace Aboria;

class LabwohlLasherPotential {
public:
	LabwohlLasherPotential(double epsilon,double lattice_spacing):
		epsilon(epsilon),lattice_spacing(lattice_spacing) {}
	double evaluate(const Vect3d &x1, const Vect3d &u1, const Vect3d &x2, const Vect3d &u2) const;
	double operator()(const Vect3d &x1, const Vect3d &u1, const Vect3d &x2, const Vect3d &u2) const {
		return evaluate(x1,u1,x2,u2);
	}
	double cut_off() const;
	virtual ~LabwohlLasherPotential();
private:
	double lattice_spacing;
	double epsilon;
};

#endif /* LABWOHLLASHERPOTENTIAL_H_ */
