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


template <typename T>
void ParticleSimulation::langevin_timestep(const unsigned int n, T& potential) {
	const Vect3d min(0,0,0);
	const Vect3d max(L,L,L);
	const Vect3b periodic(false,false,false);
	const double diameter = potential.cut_off();
	particles->init_neighbour_search(min,max,diameter,periodic);
	for (int i = 0; i < n; ++i) {
		particles->update_positions(particles->begin(),particles->end(),[particles,dt,diameter,D,T,k_s](SpeciesType::Value& i) {
			REGISTER_SPECIES_PARTICLE(i);
			v << 0,0,0;
			U = 0;
			for (auto tpl: i.get_neighbours(particles)) {
				REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
				const double r2 = dx.squaredNorm();
				if (r2 > diameter*diameter) continue;
				if (r2 == 0) continue;

				const double r = dx.norm();
				const double overlap = diameter-r;
				if (overlap>0) {
					const Vect3d normal = dx/r;
					v += potential.grad(r,u,rj,uj);
					U += potential(r,u,rj,uj);
				}
			}
			v *= dt*D/(k_b*T);
			v += sqrt(2.0*D*dt)*Vect3d(i.rand_normal(),i.rand_normal(),i.rand_normal());
			rt += v;
			const Vect3d new_position = r+v;
			if ((new_position.array() < particles->get_low().array()).any() ||
					(new_position.array() >= particles->get_high().array()).any()) {
				exits++;
			}
			return new_position;
		});
	}
}

#endif /* PARTICLESIMULATION_H_ */
