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

enum {SPECIES_VELOCITY,SPECIES_ORIENTATION,SPECIES_POTENTIAL,SPECIES_AVERAGED_ORIENTATION,SPECIES_AVERAGED_POSITION,SPECIES_NUM_MOVES};
typedef std::tuple<Vect3d,Vect3d,double,Vect3d,Vect3d,unsigned int> SpeciesTuple;
typedef Particles<SpeciesTuple> SpeciesType;
namespace Aboria {
template <>
class DataNames<SpeciesTuple> {
public:
	std::string static get(unsigned int n) {
		std::string ret[7] = {"velocity","orientation","potential","averaged_orientation","averaged_position","number_of_moves"};
		return ret[n];
	}
};
}

#define GET_TUPLE(type,name,position,particle) type& name = std::get<position>(particle.get_data())
#define REGISTER_SPECIES_PARTICLE(particle) \
				const Vect3d& r = particle.get_position(); \
				const bool alive = particle.is_alive(); \
				GET_TUPLE(double,U,SPECIES_POTENTIAL,particle); \
				GET_TUPLE(Vect3d,u,SPECIES_ORIENTATION,particle); \
				GET_TUPLE(Vect3d,ra,SPECIES_AVERAGED_POSITION,particle); \
				GET_TUPLE(Vect3d,ua,SPECIES_AVERAGED_ORIENTATION,particle); \
				GET_TUPLE(unsigned int,n,SPECIES_NUM_MOVES,particle); \
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
