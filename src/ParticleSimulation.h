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
		Dt(0),Dr(0),T(0),dt(0),L(0) {}
	void add_particles(const unsigned int n);
	template<typename T>
	void langevin_timestep(const unsigned int n, T& potential);
	template<typename T>
	void monte_carlo_timestep(const unsigned int n, T& potential);
	virtual ~ParticleSimulation();

private:
	double Dt;
	double Dr;
	double T;
	double dt;
	double L;
	std::mt19937 generator;
    ptr<SpeciesType> particles;
};

#endif /* PARTICLESIMULATION_H_ */
