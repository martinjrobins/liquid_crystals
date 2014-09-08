/*
 * ParticleSimulation.cpp
 *
 *  Created on: 6 Sep 2014
 *      Author: mrobins
 */

#include "ParticleSimulation.h"


ParticleSimulation::~ParticleSimulation() {
	// TODO Auto-generated destructor stub
}

void ParticleSimulation::add_particles(const unsigned int n) {
	const Vect3d min(0,0,0);
	const Vect3d max(L,L,L);
	const Vect3b periodic(true,true,true);
	particles->clear();
	std::uniform_real_distribution<double> distribution(0,L);
	auto dice = std::bind ( distribution, generator );
	particles->create_particles(n,[&dice](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);
		Vect3d canditate_position;
		r0 << dice(),dice(),dice();
		rt = r0;
		u << dice(),dice(),dice();
		u.normalize();
		v << 0,0,0;
		U = 0;
		exits = 0;
		return canditate_position;
	});

}




template<typename T>
void ParticleSimulation::monte_carlo_timestep(const unsigned int n, T& potential) {
	const Vect3d min(0,0,0);
	const Vect3d max(L,L,L);
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
				const Vect3d rand_inc = sqrt(2.0*D*dt)*Vect3d(normald(generator),normald(generator),normald(generator));
				Vect3d canditate_pos = r+rand_inc;
				for (int d = 0; d < 3; ++d) {
					while (canditate_pos[d]<low[d]) {
						canditate_pos[d] += (high[d]-low[d]);
					}
					while (canditate_pos[d]>=high[d]) {
						canditate_pos[d] -= (high[d]-low[d]);
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
						Udiff -= potential(r,u,rj,uj);
					}
				}
				for (auto tpl: particles->get_neighbours(canditate_pos)) {
					REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
					const double r2 = dx.squaredNorm();
					if (r2 > diameter*diameter) continue;
					if (j.get_id()==((*particles)[index]).get_id()) continue;
					const double r = dx.norm();
					const double overlap = diameter-r;
					if (overlap>0) {
						Udiff += potential(r,u,rj,uj);
					}
				}

				const double acceptance_ratio = exp(-Udiff/(k_b*T));
				//std::cout <<"dU = "<<newU-oldU<<" acceptance_ratio = "<<acceptance_ratio<<std::endl;
				if (uniformd(generator)<acceptance_ratio) {
					//std::cout <<"accepted"<<std::endl;

					particles->update_positions_sequential(particles->begin()+index,particles->begin()+index+1,[particles,&canditate_pos,&rand_inc](SpeciesType::Value& i) {
						REGISTER_SPECIES_PARTICLE(i);
						rt += rand_inc;
						const Vect3d non_periodic_position = r+rand_inc;
						if ((non_periodic_position.array() < particles->get_low().array()).any() ||
								(non_periodic_position.array() >= particles->get_high().array()).any()) {
							exits++;
						}
						return canditate_pos;
					});
					break;
				}
			}
		}
	}
}

//template void ParticleSimulation::monte_carlo_timestep(const unsigned int n, GayBernePotential& potential);



