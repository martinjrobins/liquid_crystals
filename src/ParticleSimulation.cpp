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
	particles->clear();
	std::uniform_real_distribution<double> distribution(0,L);
	auto dice = std::bind ( distribution, generator );
	particles->create_particles(n,[&dice](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);
		Vect3d canditate_position(dice(),dice(),0);
		r0 << canditate_position;
		rt = r0;
		u << dice(),dice(),0;
		u.normalize();
		v << 0,0,0;
		U = 0;
		exits = 0;
		return canditate_position;
	});

}







