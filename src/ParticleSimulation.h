/*
 * ParticleSimulation.h
 *
 *  Created on: 6 Sep 2014
 *      Author: mrobins
 */

#ifndef PARTICLESIMULATION_H_
#define PARTICLESIMULATION_H_

#include "Types.h"
#include <random>

class ParticleSimulation {
public:
	ParticleSimulation():
		D_t(0),D_r(0),T(0),dt(0),method(0) {}
	void init_simulation();
	void step_simulation(unsigned int n);
	virtual ~ParticleSimulation();

private:
	void langevin_timestep();
	void monte_carlo_timestep();
	double Dt;
	double Dr;
	double T;
	double dt;
	unsigned int method;
	std::mt19937 generator;
    ptr<SpeciesType> particles;
};

#endif /* PARTICLESIMULATION_H_ */
