#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <map>
#include "Aboria.h"

template<typename DATA_TYPE, typename POTENTIAL>
double monte_carlo_timestep(const unsigned int Nb, unsigned int Na, ptr<Particles<DATA_TYPE> > particles, ptr<Particles<DATA_TYPE> > averaging_points, POTENTIAL& potential, std::map<std::string,double> &params) {
	std::cout << "starting monte_carlo_timestep"<<std::endl;

	//Particles<DATA_TYPE> start_copy(*particles);

	std::for_each(particles->begin(),particles->end(),[](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);
		n = 0;
		u = Vect3d(cos(theta),sin(theta),0);
	});

	std::for_each(averaging_points->begin(),averaging_points->end(),[](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);
		n = 0;
	});

	const double Dtrans = params["Dtrans"];
	const double Drot = params["Drot"];
	const double Temp = params["Temp"];
	const double L = params["L"];
	const double h = params["h"];

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

	const int outstep = Nb/100;
	if (Na <= 0) Na = 1;
	const int avstep = Nb/Na;
	for (int j = 0; j < Nb; ++j) {
		if ((j-1)%outstep == 0) std::cout << "." << std::flush;
		for (int ii = 0; ii < particles->size(); ++ii) {

			/*
			 * generate new state x'
			 */
			const int index = (*particles)[0].rand_uniform()*particles->size();
			if ((index < 0)||(index > particles->size()-1)) {
				std::cout << "crap" << std::endl;
				ERROR("bad index!!");
			}
			//const int index = ii;
			typename Particles<DATA_TYPE>::value_type& i = (*particles)[index];
			REGISTER_SPECIES_PARTICLE(i);
			//thetaa = (theta + n*thetaa)/(n+1);
			//ua = Vect3d(cos(thetaa),sin(thetaa),0);
			ua = (u + n*ua)/(n+1);
			ua.normalize();
			thetaa = acos(ua[0]);
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

				particles->update_positions_sequential(particles->begin()+index,particles->begin()+index+1,[&candidate_pos,&candidate_u,&candidate_theta](SpeciesType::Value& i) {
					REGISTER_SPECIES_PARTICLE(i);
					//std::cout <<"updating position to "<<candidate_pos<<std::endl;
					theta = candidate_theta;
					u = candidate_u;
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
		if (j%avstep == 0) {
			particles->reset_neighbour_search(h);
			std::for_each(averaging_points->begin(),averaging_points->end(),[particles,h](SpeciesType::Value& i) {
				REGISTER_SPECIES_PARTICLE(i);

				u = Vect3d(0,0,0);
				int count = 0;
				for (auto tpl: particles->get_neighbours(r)) {
					REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
					const double r2 = dx.squaredNorm();
					if (r2 > h*h) continue;
					u[0] += 2.0*uj[0]*uj[0] - 1.0;
					u[1] += uj[0]*uj[1];
					count++;
				}
				u[0] /= count;
				u[1] /= count;
				U = u.norm();
				const double m1 = sqrt(u[0]/(2*U)+0.5);
				const double m2 = u[1]/(2*U*m1);
				theta = acos(m1);
				ua = (u + n*ua)/(n+1);
				const double Ua = ua.norm();
				const double n1 = sqrt(ua[0]/(2*Ua)+0.5);
				const double n2 = ua[1]/(2*U*n1);
				thetaa = acos(n1);
				n++;
			});
			particles->reset_neighbour_search(diameter);
		}
	}
	std::cout <<"finished monte carlo steps ratio of accepts to rejects is "<<double(accepts)/++rejects<<std::endl;

//	particles->update_positions(particles->begin(),particles->end(),[&particles](SpeciesType::Value& i) {
//		REGISTER_SPECIES_PARTICLE(i);
//		//std::cout <<"updating position to "<<candidate_pos<<std::endl;
//		const Vect3d averaged_pos = ra;
//		ra = r;
//		return averaged_pos;
//	});

	std::for_each(particles->begin(),particles->end(),[particles,&potential,diameter,L](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);
		U = 0;
		for (auto tpl: particles->get_neighbours(r)) {
			REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > diameter*diameter) continue;
			if (j.get_id()==i.get_id()) continue;
			if (fixed && fixedj) continue;
			U += 0.5*potential(r,u,rj,uj);
		}

	});

//	particles->update_positions(particles->begin(),particles->end(),[&particles](SpeciesType::Value& i) {
//		REGISTER_SPECIES_PARTICLE(i);
//		//std::cout <<"updating position to "<<candidate_pos<<std::endl;
//		const Vect3d original_pos = ra;
//		ra = r;
//		return original_pos;
//	});

//	return std::accumulate(particles->begin(), particles->end(), 0.0, [](double U, SpeciesType::Value& i){
//		GET_TUPLE(double,Ui,SPECIES_POTENTIAL,i);
//		return U+Ui;
//	});


	Vect3d Q = std::accumulate(particles->begin(), particles->end(), Vect3d(0,0,0), [](Vect3d Q, SpeciesType::Value& i){
		REGISTER_SPECIES_PARTICLE(i);
		Q[0] += 2.0*ua[0]*ua[0] - 1.0;
		Q[1] += ua[0]*ua[1];
		return Q;
	});
	return (Q/particles->size()).norm();
}

template<typename DATA_TYPE>
void local_averaging(const double length_scale, ptr<Particles<DATA_TYPE> > averaging_points, ptr<Particles<DATA_TYPE> > source_points) {
	source_points->reset_neighbour_search(length_scale);
	std::for_each(averaging_points->begin(),averaging_points->end(),[source_points,length_scale](SpeciesType::Value& i) {
		REGISTER_SPECIES_PARTICLE(i);
		Q11 = 0;
		Q12 = 0;
		int count = 0;
		for (auto tpl: source_points->get_neighbours(r)) {
			REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > length_scale*length_scale) continue;
			if (j.get_id()==i.get_id()) continue;
			Q11 += 2.0*uj[0]*uj[0] - 1.0;
			Q12 += uj[0]*uj[1];
			count++;
		}
		Q11 /= count;
		Q12 /= count;
	});

}
#endif
