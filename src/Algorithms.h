#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <map>
#include "Aboria.h"

template<typename DATA_TYPE, typename POTENTIAL>
double monte_carlo_timestep(const unsigned int Nb, ptr<Particles<DATA_TYPE> > particles, POTENTIAL& potential, std::map<std::string,double> &params) {
	std::cout << "starting monte_carlo_timestep"<<std::endl;

	//Particles<DATA_TYPE> start_copy(*particles);

	std::for_each(particles->begin(),particles->end(),[](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);
		n = 0;
		ra = r;
		thetaa = theta;
		u = Vect3d(cos(theta),sin(theta),0);
		ua = u;
	});

	const double Dtrans = params["Dtrans"];
	const double Drot = params["Drot"];
	const double Temp = params["Temp"];
	const double L = params["L"];

	std::cout <<"running with Dtrans = "<<Dtrans<<" Drot = "<<Drot<<" Temp = "<<Temp<<" L = "<<L<<std::endl;

	const double buffer = L/10.0;
	const Vect3d min(-buffer,-buffer,-2*buffer);
	const Vect3d max(L+buffer,L+buffer,2*buffer);
	const Vect3b periodic(false,false,false);
	const double diameter = potential.cut_off();
	particles->init_neighbour_search(min,max,diameter,periodic);

	std::uniform_real_distribution<double> uniformd(0,1);
	std::normal_distribution<double> normald(0,1);
	const Vect3d low = particles->get_low();
	const Vect3d high = particles->get_high();
	unsigned int rejections = 0;
	unsigned int rejects = 0;
	unsigned int accepts = 0;

	for (int j = 0; j < Nb; ++j) {
		for (int ii = 0; ii < particles->size(); ++ii) {

			/*
			 * generate new state x'
			 */
			const int index = (*particles)[0].rand_uniform()*particles->size();
			//const int index = ii;
			typename Particles<DATA_TYPE>::value_type& i = (*particles)[index];
			REGISTER_SPECIES_PARTICLE(i);
			thetaa = (theta + n*thetaa)/(n+1);
			ua = Vect3d(cos(thetaa),sin(thetaa),0);
			ra = (r + n*ra)/(n+1);
			n++;
			if (fixed) continue;

			const Vect3d rand_inc = Dtrans*Vect3d(i.rand_uniform()-0.5,i.rand_uniform()-0.5,0);
			const double rand_inc2 = Drot*(i.rand_uniform()-0.5);
			Vect3d candidate_pos = r+rand_inc;
			double candidate_theta = theta+rand_inc2;
			Vect3d current_u(cos(theta),sin(theta),0);
			Vect3d candidate_u(cos(candidate_theta),sin(candidate_theta),0);
			if ((candidate_pos[0]<0)||(candidate_pos[0]>=L)) {
				rejects++;	
				continue;
			}
			if ((candidate_pos[1]<0)||(candidate_pos[1]>=L)) {
				rejects++;	
				continue;
			}
	 

			double Udiff = 0;
			for (auto tpl: particles->get_neighbours(r)) {
				REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
				const double r2 = dx.squaredNorm();
				if (r2 > diameter*diameter) continue;
				if (j.get_id()==i.get_id()) continue;
				Udiff -= potential(r,current_u,rj,uj);
			}

			for (auto tpl: particles->get_neighbours(candidate_pos)) {
				REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
				const double r2 = dx.squaredNorm();
				if (r2 > diameter*diameter) continue;
				if (j.get_id()==i.get_id()) continue;
				Udiff += potential(candidate_pos,candidate_u,rj,uj);
			}

			const double acceptance_ratio = exp(-Udiff/Temp);
			//std::cout <<"dU = "<<Udiff<<" T = "<<Temp<<"acceptance_ratio = "<<acceptance_ratio<<std::endl;
			if ((Udiff <= 0) || ((*particles)[0].rand_uniform()<acceptance_ratio)) {
				//std::cout <<"accepted"<<std::endl;

				particles->update_positions_sequential(particles->begin()+index,particles->begin()+index+1,[&candidate_pos,&candidate_theta](SpeciesType::Value& i) {
					REGISTER_SPECIES_PARTICLE(i);
					//std::cout <<"updating position to "<<candidate_pos<<std::endl;
					theta = candidate_theta;
					u = Vect3d(cos(theta),sin(theta),0);
					return candidate_pos;
				});
				accepts++;
			} else {
				rejects++;
				//std::cout << "rejected move from r ="<<r<<" to "<<candidate_pos<<std::endl;
				//					if (rejections++ > 100) {
				//						ERROR("number of monte carlo rejections > 100");
				//					}
			}
		}
	}
	std::cout <<"finished monte carlo steps ratio of accepts to rejects is "<<double(accepts)/++rejects<<std::endl;

	particles->update_positions(particles->begin(),particles->end(),[&particles](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);
		//std::cout <<"updating position to "<<candidate_pos<<std::endl;
		const Vect3d averaged_pos = ra;
		ra = r;
		return averaged_pos;
	});

	std::for_each(particles->begin(),particles->end(),[particles,&potential,diameter,L](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);
		U = 0;
		for (auto tpl: particles->get_neighbours(r)) {
			REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > diameter*diameter) continue;
			if (j.get_id()==i.get_id()) continue;
			if (fixed && fixedj) continue;
			U += 0.5*potential(r,ua,rj,uaj);
		}

	});

	particles->update_positions(particles->begin(),particles->end(),[&particles](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);
		//std::cout <<"updating position to "<<candidate_pos<<std::endl;
		const Vect3d original_pos = ra;
		ra = r;
		return original_pos;
	});

	return std::accumulate(particles->begin(), particles->end(), 0.0, [](double U, SpeciesType::Value& i){
		GET_TUPLE(double,Ui,SPECIES_POTENTIAL,i);
		return U+Ui;
	});

}

#endif
