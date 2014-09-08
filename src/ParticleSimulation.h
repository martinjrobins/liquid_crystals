/*
 * ParticleSimulation.h
 *
 *  Created on: 6 Sep 2014
 *      Author: mrobins
 */

#ifndef PARTICLESIMULATION_H_
#define PARTICLESIMULATION_H_

#include "Types.h"
#include "GayBernePotential.h"
#include <random>

class ParticleSimulation {
public:
	ParticleSimulation():
		Dtrans(0),Drot(0),T(0),dt(0),L(0),particles(SpeciesType::New()) {}
	void add_particles(const unsigned int n);
	template<typename T>
	void langevin_timestep(const unsigned int n, T& potential);
	template<typename T>
	void monte_carlo_timestep(const unsigned int n, T& potential);
	virtual ~ParticleSimulation();

	double getDrot() const {
		return Drot;
	}

	void setDrot(double dr) {
		Drot = dr;
	}

	double getDt() const {
		return dt;
	}

	void setDt(double dt) {
		this->dt = dt;
	}

	double getDtrans() const {
		return Dtrans;
	}

	void setDtrans(double dtrans) {
		Dtrans = dtrans;
	}

	double getL() const {
		return L;
	}

	void setL(double l) {
		L = l;
	}

	ptr<SpeciesType> getParticles() const {
		return particles;
	}

	void setParticles(ptr<SpeciesType> particles) {
		this->particles = particles;
	}

	double getT() const {
		return T;
	}

	void setT(double t) {
		T = t;
	}

private:
	double Dtrans;
	double Drot;
	double T;
	double dt;
	double L;
	std::mt19937 generator;
	ptr<SpeciesType> particles;
};

#endif /* PARTICLESIMULATION_H_ */
