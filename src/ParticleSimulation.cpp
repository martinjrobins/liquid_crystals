/*
 * ParticleSimulation.cpp
 *
 *  Created on: 6 Sep 2014
 *      Author: mrobins
 */

#include "ParticleSimulation.h"

ParticleSimulation::ParticleSimulation() {
	// TODO Auto-generated constructor stub

}

ParticleSimulation::~ParticleSimulation() {
	// TODO Auto-generated destructor stub
}

void ParticleSimulation::langevin_timestep() {


	A->update_positions(A->begin(),A->end(),[A,dt,diameter,D,T,k_s](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);
		v << 0,0,0;
		U = 0;
		for (auto tpl: i.get_neighbours(A)) {
			REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > diameter*diameter) continue;
			if (r2 == 0) continue;

			const double r = dx.norm();
			const double overlap = diameter-r;
			if (overlap>0) {
				const Vect3d normal = dx/r;
				v += (k_s*overlap)*normal;
				U += 0.5*k_s*pow(overlap,2);
			}
		}
		v *= dt*D/(k_b*T);
		v += sqrt(2.0*D*dt)*Vect3d(i.rand_normal(),i.rand_normal(),i.rand_normal());
		rt += v;
		const Vect3d new_position = r+v;
		if ((new_position.array() < A->get_low().array()).any() ||
				(new_position.array() >= A->get_high().array()).any()) {
			exits++;
		}
		return new_position;
	});

}

void ParticleSimulation::monte_carlo_timestep() {

	std::uniform_real_distribution<double> uniformd(0,1);
	std::normal_distribution<double> normald(0,1);
	const Vect3d low = A->get_low();
	const Vect3d high = A->get_high();
	int index;
	for (int ii = 0; ii < A->size(); ++ii) {
		while(1) {
			/*
			 * generate new state x'
			 */
			const int index = uniformd(generator)*A->size();
			REGISTER_SPECIES_PARTICLE(((*A)[index]));
			Vect3d canditate_pos = r+sqrt(2.0*D*dt)*Vect3d(normald(generator),normald(generator),normald(generator));
			for (int d = 0; d < 3; ++d) {
				while (canditate_pos[d]<low[d]) {
					canditate_pos[d] += (high[d]-low[d]);
				}
				while (canditate_pos[d]>=high[d]) {
					canditate_pos[d] -= (high[d]-low[d]);
				}
			}
			double Udiff = 0;
			for (auto tpl: A->get_neighbours(r)) {
				REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
				const double r2 = dx.squaredNorm();
				if (r2 > diameter*diameter) continue;
				if (j.get_id()==((*A)[index]).get_id()) continue;
				const double r = dx.norm();
				const double overlap = diameter-r;
				if (overlap>0) {
					Udiff -= k_s*pow(overlap,2);
				}
			}
			for (auto tpl: A->get_neighbours(canditate_pos)) {
				REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
				const double r2 = dx.squaredNorm();
				if (r2 > diameter*diameter) continue;
				if (j.get_id()==((*A)[index]).get_id()) continue;
				const double r = dx.norm();
				const double overlap = diameter-r;
				if (overlap>0) {
					Udiff += k_s*pow(overlap,2);
				}
			}

			const double acceptance_ratio = exp(-Udiff/(k_b*T));
			//std::cout <<"dU = "<<newU-oldU<<" acceptance_ratio = "<<acceptance_ratio<<std::endl;
			if (uniformd(generator)<acceptance_ratio) {
				//std::cout <<"accepted"<<std::endl;

				A->update_positions(A->begin()+index,A->begin()+index+1,[canditate_pos](SpeciesType::Value& i) {
					return canditate_pos;
				});
				break;
			}
		}
	}
}


