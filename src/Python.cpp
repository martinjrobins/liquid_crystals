/*
 * Python.cpp
 *
 *  Created on: 7 Sep 2014
 *      Author: mrobins
 */


#include "Types.h"
#include "ParticleSimulation.h"
#include "GayBernePotential.h"
#include "Aboria.h"
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace Aboria;
using namespace boost::python;


template void ParticleSimulation::monte_carlo_timestep(const unsigned int n, GayBernePotential& potential);

template<typename T>
struct Vect3_from_python_list
{

	Vect3_from_python_list() {
		boost::python::converter::registry::push_back(
				&convertible,
				&construct,
				boost::python::type_id<Eigen::Matrix<T,3,1> >());
	}

	static void* convertible(PyObject* obj_ptr) {
		if (!PyList_Check(obj_ptr)) return 0;
		if (PyList_Size(obj_ptr) != 3) return 0;
		return obj_ptr;
	}

	static void construct(
			PyObject* obj_ptr,
			boost::python::converter::rvalue_from_python_stage1_data* data) {
		const int size = PyList_Size(obj_ptr);

		// Extract the character data from the python string
		const double x = extract<T>(PyList_GetItem(obj_ptr,0));
		const double y = extract<T>(PyList_GetItem(obj_ptr,1));
		const double z = extract<T>(PyList_GetItem(obj_ptr,2));


		// Grab pointer to memory into which to construct the new QString
		void* storage = (
				(boost::python::converter::rvalue_from_python_storage<Eigen::Matrix<T,3,1> >*)
				data)->storage.bytes;

		// in-place construct the new QString using the character data
		// extraced from the python object
		new (storage) Eigen::Matrix<T,3,1>(x,y,z);

		// Stash the memory chunk pointer for later use by boost.python
		data->convertible = storage;
	}

};


template<class T>
struct vtkSmartPointer_to_python {
	static PyObject *convert(const vtkSmartPointer<T> &p) {
		std::ostringstream oss;
		oss << (void*) p.GetPointer();
		std::string address_str = oss.str();

		using namespace boost::python;
		object obj = import("vtk").attr("vtkObject")(address_str);
		object obj2 = import("tvtk.api").attr("tvtk").attr("to_tvtk")(obj);
		return incref(obj2.ptr());
	}
};

#define VTK_PYTHON_CONVERSION(type) \
    /* register the to-python converter */ \
    to_python_converter<vtkSmartPointer<type>,vtkSmartPointer_to_python<type> >(); \


BOOST_PYTHON_MODULE(particleSimulation) {

	VTK_PYTHON_CONVERSION(vtkUnstructuredGrid);

	Vect3_from_python_list<double>();
	Vect3_from_python_list<int>();
	Vect3_from_python_list<bool>();


	/*
	 * Particles
	 */
	class_<SpeciesType,typename Aboria::ptr<SpeciesType> >("Particles")
	        .def(boost::python::vector_indexing_suite<SpeciesType >())
	        .def("get_grid",&SpeciesType::get_grid)
	    ;


	/*
	 * Particle
	 */
	class_<SpeciesType::value_type>("Particle")
		.add_property("id", &SpeciesType::value_type::get_id)
		.def("position", &SpeciesType::value_type::get_position,
				return_value_policy<reference_existing_object,
				with_custodian_and_ward_postcall<0,1> >())
		.def("data", &SpeciesType::value_type::get_data,
				return_value_policy<reference_existing_object,
				with_custodian_and_ward_postcall<0,1> >())
		;


	/*
	 * ParticleSimulation
	 */

	class_<ParticleSimulation>("ParticleSimulation", init<>())
			.def(init<double,double,double,double,double>(
					(arg("Dtrans"),arg("Drot"),arg("Temp"),arg("dt"),arg("L"))))
			.def("add_particles",&ParticleSimulation::add_particles)
			.def("monte_carlo_timestep",&ParticleSimulation::monte_carlo_timestep<GayBernePotential>)
			.add_property("Dtrans", &ParticleSimulation::getDtrans, &ParticleSimulation::setDtrans)
			.add_property("Drot", &ParticleSimulation::getDrot, &ParticleSimulation::setDrot)
			.add_property("Temp", &ParticleSimulation::getTemp, &ParticleSimulation::setTemp)
			.add_property("dt", &ParticleSimulation::getDt, &ParticleSimulation::setDt)
			.add_property("L", &ParticleSimulation::getL, &ParticleSimulation::setL)
			.add_property("particles", &ParticleSimulation::getParticles, &ParticleSimulation::setParticles)
			;
	/*
	 * GayBernePotential
	 */

	class_<GayBernePotential>("GayBernePotential", init<double, double, double, double, double, double>(
			(arg("sigma_s"),arg("k"),arg("kdash"),arg("mu"),arg("nu"),arg("epsilon_0"))))
			.def("evaluate",&GayBernePotential::evaluate)
			;


}
