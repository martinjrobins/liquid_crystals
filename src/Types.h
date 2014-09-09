/*
 * Types.h
 *
 *  Created on: 6 Sep 2014
 *      Author: mrobins
 */

#ifndef TYPES_H_
#define TYPES_H_

#include "Aboria.h"
using namespace Aboria;

#include <tuple>

const double k_b = 1.3806488e-23;

enum {SPECIES_VELOCITY,SPECIES_ORIENTATION,SPECIES_POTENTIAL,SPECIES_TOTAL_R,SPECIES_SAVED_R,SPECIES_SAVED_R1,SPECIES_NUM_EXITS};
typedef std::tuple<Vect3d,Vect3d,double,Vect3d,Vect3d,Vect3d,unsigned int> SpeciesTuple;
typedef Particles<SpeciesTuple> SpeciesType;
namespace Aboria {
template <>
class DataNames<SpeciesTuple> {
public:
	std::string static get(unsigned int n) {
		std::string ret[7] = {"velocity","orientation","potential","total_r","saved_r","saved_r1","num_exits"};
		return ret[n];
	}
};
}

#define GET_TUPLE(type,name,position,particle) type& name = std::get<position>(particle.get_data())
#define REGISTER_SPECIES_PARTICLE(particle) \
				const Vect3d& r = particle.get_position(); \
				const bool alive = particle.is_alive(); \
				GET_TUPLE(double,U,SPECIES_POTENTIAL,particle); \
				GET_TUPLE(unsigned int,exits,SPECIES_NUM_EXITS,particle); \
				GET_TUPLE(Vect3d,r0,SPECIES_SAVED_R,particle); \
				GET_TUPLE(Vect3d,u,SPECIES_ORIENTATION,particle); \
				GET_TUPLE(Vect3d,rt,SPECIES_TOTAL_R,particle); \
				GET_TUPLE(Vect3d,r1,SPECIES_SAVED_R1,particle); \
				GET_TUPLE(Vect3d,v,SPECIES_VELOCITY,particle);
#define GET_TUPLE_CONST(type,name,position,particle) const type& name = std::get<position>(particle.get_data_const())
#define REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tuple) \
				const Vect3d& dx = std::get<1>(tuple); \
				const SpeciesType::Value& j = std::get<0>(tuple); \
				const Vect3d& rj = j.get_position(); \
				const bool alivej = j.is_alive(); \
				GET_TUPLE_CONST(double,Uj,SPECIES_POTENTIAL,j); \
				GET_TUPLE_CONST(Vect3d,uj,SPECIES_ORIENTATION,j); \
				GET_TUPLE_CONST(Vect3d,vj,SPECIES_VELOCITY,j);


#endif /* TYPES_H_ */
