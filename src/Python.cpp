/*
 * Python.cpp
 *
 *  Created on: 7 Sep 2014
 *      Author: mrobins
 */


#include "Types.h"
#include "ParticleSimulation.h"

BOOST_PYTHON_MODULE(particle_simulation) {


	/*
	 * vector
	 */
	class_<Aboria::Particles<> >("Vectd")
	        .def(boost::python::vector_indexing_suite<std::vector<double> >())
	    ;

	VTK_PYTHON_CONVERSION(vtkUnstructuredGrid);
	VTK_PYTHON_CONVERSION(vtkPolyData);


	Vect3_from_python_list<double>();
	Vect3_from_python_list<int>();
	Vect3_from_python_list<bool>();
	ReactionEquation_from_python_list();

	/*
	 * Species
	 */

	class_<Species,typename std::auto_ptr<Species> >("Species",boost::python::init<double>())
			.def("fill_uniform",Species_fill_uniform)
			.def("fill_uniform",Species_fill_uniform_interface)
			.def("get_concentration",Species_get_concentration1)
