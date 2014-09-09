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
		Dtrans(0),Drot(0),Temp(0),dt(0),L(0),particles(SpeciesType::New()) {}
	ParticleSimulation(double Dtrans, double Drot, double Temp, double dt, double L):
			Dtrans(Dtrans),Drot(Drot),Temp(Drot),dt(Drot),L(L),particles(SpeciesType::New()) {}
	void add_particles(const unsigned int n);
	template<typename T>
	void langevin_timestep(const unsigned int n, T& potential);
	template<typename T>
	void monte_carlo_timestep(const unsigned int n, T& potential);
	virtual ~ParticleSimulation();

	double getTemp() const {
		return Temp;
	}

	void setTemp(double temp) {
		Temp = temp;
	}

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

private:
	double Dtrans;
	double Drot;
	double Temp;
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
		particles->update_positions(particles->begin(),particles->end(),[&potential,this,diameter](SpeciesType::Value& i) {
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
			v *= this->dt*this->Dtrans/this->Temp;
			v += sqrt(2.0*this->Dtrans*this->dt)*Vect3d(i.rand_normal(),i.rand_normal(),i.rand_normal());
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

template<typename T>
void ParticleSimulation::monte_carlo_timestep(const unsigned int n, T& potential) {
	std::cout << "starting monte_carlo_timestep"<<std::endl;
	const Vect3d min(0,0,-L/100.0);
	const Vect3d max(L,L,L/100.0);
	const Vect3b periodic(false,false,false);
	const double diameter = potential.cut_off();
	particles->init_neighbour_search(min,max,diameter,periodic);

	std::uniform_real_distribution<double> uniformd(0,1);
	std::normal_distribution<double> normald(0,1);
	const Vect3d low = particles->get_low();
	const Vect3d high = particles->get_high();
	int index;

	for (int i = 0; i < n; ++i) {
		for (int ii = 0; ii < particles->size(); ++ii) {
			while(1) {
				/*
				 * generate new state x'
				 */
				const int index = uniformd(generator)*particles->size();
				REGISTER_SPECIES_PARTICLE(((*particles)[index]));
				const Vect3d rand_inc = sqrt(2.0*Dtrans*dt)*Vect3d(normald(generator),normald(generator),0);
				const double rand_inc2 = sqrt(2.0*Drot*dt)*normald(generator);
				const double cosinc = cos(rand_inc2);
				const double sininc = sin(rand_inc2);
				Vect3d candidate_pos = r+rand_inc;
				Vect3d candidate_u(u[0]*cosinc - u[1]*sininc, u[0]*sininc + u[1]*cosinc,u[2]);
				for (int d = 0; d < 3; ++d) {
					while (candidate_pos[d]<low[d]) {
						candidate_pos[d] += (high[d]-low[d]);
					}
					while (candidate_pos[d]>=high[d]) {
						candidate_pos[d] -= (high[d]-low[d]);
					}
				}

				double Udiff = 0;
				for (auto tpl: particles->get_neighbours(r)) {
					REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
					const double r2 = dx.squaredNorm();
					if (r2 > diameter*diameter) continue;
					if (j.get_id()==((*particles)[index]).get_id()) continue;
					const double r = dx.norm();
					const double overlap = diameter-r;
					if (overlap>0) {
						Udiff -= potential(candidate_pos,candidate_u,rj,uj);
					}
				}
				for (auto tpl: particles->get_neighbours(candidate_pos)) {
					REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
					const double r2 = dx.squaredNorm();
					if (r2 > diameter*diameter) continue;
					if (j.get_id()==((*particles)[index]).get_id()) continue;
					const double r = dx.norm();
					const double overlap = diameter-r;
					if (overlap>0) {
						Udiff += potential(candidate_pos,candidate_u,rj,uj);
					}
				}

				const double acceptance_ratio = exp(-Udiff/Temp);
				//std::cout <<"dU = "<<newU-oldU<<" acceptance_ratio = "<<acceptance_ratio<<std::endl;
				if (uniformd(generator)<acceptance_ratio) {
					//std::cout <<"accepted"<<std::endl;

					particles->update_positions_sequential(particles->begin()+index,particles->begin()+index+1,[this,&candidate_pos,&rand_inc,&candidate_u](SpeciesType::Value& i) {
						REGISTER_SPECIES_PARTICLE(i);
						rt += rand_inc;
						const Vect3d non_periodic_position = r+rand_inc;
						if ((non_periodic_position.array() < this->particles->get_low().array()).any() ||
								(non_periodic_position.array() >= this->particles->get_high().array()).any()) {
							exits++;
						}
						u = candidate_u;
						return candidate_pos;
					});
					break;
				}
			}
		}
	}
}

#endif /* PARTICLESIMULATION_H_ */
