/*
 * Python.cpp
 *
 *  Created on: 7 Sep 2014
 *      Author: mrobins
 */


#include "Types.h"
#include "ParticleSimulation.h"
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace Aboria;
using namespace boost::python;


BOOST_PYTHON_MODULE(particle_simulation) {


	/*
	 * Particles
	 */
	class_<SpeciesType>("Particles")
	        .def(boost::python::vector_indexing_suite<SpeciesType >())
	    ;


	/*
	 * ParticleSimulation
	 */

	class_<ParticleSimulation>("ParticleSimulation")
			.add_property("add_particles",&ParticleSimulation::add_particles)
			.add_property("monte_carlo_timestep",&ParticleSimulation::monte_carlo_timestep<GayBernePotential>)
			.add_property("Dtrans", &ParticleSimulation::getDtrans, &ParticleSimulation::setDtrans)
			.add_property("Drot", &ParticleSimulation::getDrot, &ParticleSimulation::setDrot)
			.add_property("dt", &ParticleSimulation::getDt, &ParticleSimulation::setDt)
			.add_property("L", &ParticleSimulation::getL, &ParticleSimulation::setL)
			.add_property("particles", &ParticleSimulation::getParticles, &ParticleSimulation::setParticles)
			;

}
