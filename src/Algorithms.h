#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <map>
#include "Aboria.h"

template<typename DATA_TYPE, typename POTENTIAL>
void monte_carlo_timestep(const unsigned int n, ptr<Particles<DATA_TYPE> > particles, POTENTIAL& potential, std::map<std::string,double> &params) {
	std::cout << "starting monte_carlo_timestep"<<std::endl;

	const double Dtrans = params["Dtrans"];
	const double Drot = params["Drot"];
	const double Temp = params["Temp"];
	const double dt = params["dt"];
	const double L = params["L"];

	std::cout <<"running with Dtrans = "<<Dtrans<<" Drot = "<<Drot<<" Temp = "<<Temp<<" dt = "<<dt<<" L = "<<L<<std::endl;

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

	for (int j = 0; j < n; ++j) {
		for (int ii = 0; ii < particles->size(); ++ii) {
			rejections = 0;
			while(1) {
				/*
				 * generate new state x'
				 */
				const int index = (*particles)[0].rand_uniform()*particles->size();
				typename Particles<DATA_TYPE>::value_type& i = (*particles)[index];
				REGISTER_SPECIES_PARTICLE(i);
				const Vect3d rand_inc = sqrt(2.0*Dtrans*dt)*Vect3d(i.rand_normal(),i.rand_normal(),0);
				const double rand_inc2 = sqrt(2.0*Drot*dt)*(i.rand_normal());
				const double cosinc = cos(rand_inc2);
				const double sininc = sin(rand_inc2);
				Vect3d candidate_pos = r+rand_inc;
				Vect3d candidate_u(u[0]*cosinc - u[1]*sininc, u[0]*sininc + u[1]*cosinc,u[2]);
				candidate_u.normalize();
				if ((candidate_pos[0]<0)||(candidate_pos[0]>=L)) continue;
				if ((candidate_pos[1]<0)||(candidate_pos[1]>=L)) continue;
//
//				for (int d = 0; d < 3; ++d) {
//					if (candidate_pos[d]<low[d]) {
//						//candidate_pos[d] += (high[d]-low[d]);
//						continue;
//					}
//					if (candidate_pos[d]>=high[d]) {
//						//candidate_pos[d] -= (high[d]-low[d]);
//						continue;
//					}
//				}

				double Udiff = 0;
				for (auto tpl: particles->get_neighbours(r)) {
					REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
					const double r2 = dx.squaredNorm();
					if (r2 > diameter*diameter) continue;
					if (j.get_id()==i.get_id()) continue;
					Udiff -= potential(r,u,rj,uj);
				}
				Udiff -= potential(r,u,Vect3d(0,r[1],0),Vect3d(0,1,0));
				Udiff -= potential(r,u,Vect3d(L,r[1],0),Vect3d(0,1,0));
				Udiff -= potential(r,u,Vect3d(r[0],0,0),Vect3d(1,0,0));
				Udiff -= potential(r,u,Vect3d(r[0],L,0),Vect3d(1,0,0));

				for (auto tpl: particles->get_neighbours(candidate_pos)) {
					REGISTER_NEIGHBOUR_SPECIES_PARTICLE(tpl);
					const double r2 = dx.squaredNorm();
					if (r2 > diameter*diameter) continue;
					if (j.get_id()==i.get_id()) continue;
					Udiff += potential(candidate_pos,candidate_u,rj,uj);
				}
				Udiff += potential(candidate_pos,candidate_u,Vect3d(0,r[1],0),Vect3d(0,1,0));
				Udiff += potential(candidate_pos,candidate_u,Vect3d(L,r[1],0),Vect3d(0,1,0));
				Udiff += potential(candidate_pos,candidate_u,Vect3d(r[0],0,0),Vect3d(1,0,0));
				Udiff += potential(candidate_pos,candidate_u,Vect3d(r[0],L,0),Vect3d(1,0,0));

				const double acceptance_ratio = exp(-Udiff/Temp);
				//std::cout <<"dU = "<<Udiff<<" T = "<<Temp<<"acceptance_ratio = "<<acceptance_ratio<<std::endl;
				if ((*particles)[0].rand_uniform()<acceptance_ratio) {
					//std::cout <<"accepted"<<std::endl;

					particles->update_positions_sequential(particles->begin()+index,particles->begin()+index+1,[&particles,&candidate_pos,&rand_inc,&candidate_u](SpeciesType::Value& i) {
						REGISTER_SPECIES_PARTICLE(i);
						rt += rand_inc;
						const Vect3d non_periodic_position = r+rand_inc;
						if ((non_periodic_position.array() < particles->get_low().array()).any() ||
								(non_periodic_position.array() >= particles->get_high().array()).any()) {
							exits++;
						}
						u = candidate_u;
						//std::cout <<"updating position to "<<candidate_pos<<std::endl;
						return candidate_pos;
					});
					accepts++;
					break;
				} else {
					rejects++;
					//std::cout << "rejected move from r ="<<r<<" to "<<candidate_pos<<std::endl;
					if (rejections++ > 100) {
						ERROR("number of monte carlo rejections > 100");
					}

				}
			}
		}
	}
	std::cout <<"finished monte carlo steps ratio of accepts to rejects is "<<double(accepts)/++rejects<<std::endl;
}

#endif
