/*
 * Python.cpp
 *
 *  Created on: 7 Sep 2014
 *      Author: mrobins
 */


#include "Types.h"
#include "GayBernePotential.h"
#include "HGOPotential.h"
#include "LabwohlLasherPotential.h"
#include "Algorithms.h"
#include "Aboria.h"
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>


using namespace Aboria;
using namespace boost::python;


//template void ParticleSimulation::monte_carlo_timestep(const unsigned int n, GayBernePotential& potential);

template<typename T>
struct Vect3_from_python_list
{

	Vect3_from_python_list() {
		boost::python::converter::registry::push_back(
				&convertible,
				&construct,
				boost::python::type_id<Vector<T,3> >());
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
				(boost::python::converter::rvalue_from_python_storage<Vector<T,3> >*)
				data)->storage.bytes;

		// in-place construct the new QString using the character data
		// extraced from the python object
		new (storage) Vector<T,3>(x,y,z);

		// Stash the memory chunk pointer for later use by boost.python
		data->convertible = storage;
	}

};

template<typename T>
struct Vect3_to_python
{
    static PyObject* convert(Vector<T,3> const& s)
      {
    	boost::python::list ret;
    	ret.append(s[0]);
    	ret.append(s[1]);
    	ret.append(s[2]);
        return boost::python::incref(ret.ptr());
      }
};


double getitem_Vect3d(Vect3d& v, int index) {
  if(index < 0 || index >=3)
  {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    throw boost::python::error_already_set();;
  }
  return v[index];
}

void setitem_Vect3d(Vect3d& v, int index, float val)
{
  if(index < 0 || index >=3)
  {
    PyErr_SetString(PyExc_IndexError, "index out of range");
    throw boost::python::error_already_set();;
  }
  v[index] = val;
}


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

//
// This python to C++ converter uses the fact that VTK Python objects have an
// attribute called __this__, which is a string containing the memory address
// of the VTK C++ object and its class name.
// E.g. for a vtkPoints object __this__ might be "_0000000105a64420_p_vtkPoints"
//
void* extract_vtk_wrapped_pointer(PyObject* obj)
{
    char thisStr[] = "__this__";
    //first we need to get the __this__ attribute from the Python Object
    if (!PyObject_HasAttrString(obj, thisStr))
        return NULL;

    PyObject* thisAttr = PyObject_GetAttrString(obj, thisStr);
    if (thisAttr == NULL)
        return NULL;

    char* str = PyString_AsString(thisAttr);
    if(str == 0 || strlen(str) < 1)
        return NULL;

    char hex_address[32], *pEnd;
    char *_p_ = strstr(str, "_p_vtk");
    if(_p_ == NULL) return NULL;
    char *class_name = strstr(_p_, "vtk");
    if(class_name == NULL) return NULL;
    strcpy(hex_address, str+1);
    hex_address[_p_-str-1] = '\0';

    long address = strtol(hex_address, &pEnd, 16);

    vtkObjectBase* vtk_object = (vtkObjectBase*)((void*)address);
    if(vtk_object->IsA(class_name))
    {
        return vtk_object;
    }

    return NULL;
}


#define VTK_PYTHON_CONVERSION(type) \
    /* register the to-python converter */ \
    to_python_converter<vtkSmartPointer<type>,vtkSmartPointer_to_python<type> >(); \
    /* register the from-python converter */ \
    converter::registry::insert(&extract_vtk_wrapped_pointer, type_id<type>());

template Vect3d pow<double,3,int>(Vect3d,int);
template Vect3d pow<double,3,double>(Vect3d,double);

BOOST_PYTHON_MODULE(particleSimulation) {

	VTK_PYTHON_CONVERSION(vtkUnstructuredGrid);

	Vect3_from_python_list<double>();
//	Vect3_from_python_list<int>();
//	Vect3_from_python_list<bool>();
//	to_python_converter<
//		Vect3d,
//		Vect3_to_python<double> >();
//	to_python_converter<
//		Vect3i,
//		Vect3_to_python<int> >();
//	to_python_converter<
//		Vect3b,
//		Vect3_to_python<bool> >();

//	boost::python::class_< Vect3d  >("Vect3d")
//	  .def("__getitem__", &getitem_Vect3d)
//	  .def("__setitem__", &setitem_Vect3d)
//	;

	/*
	 * map
	 */
	class_<std::map<std::string,double> >("Params")
			.def(boost::python::map_indexing_suite<std::map<std::string,double> >())
			;



	class_<Vect3d,Aboria::ptr<Vect3d> >("Vect3d",init<>())
			.def(init<const Vect3d&>())
			.def(init<double,double,double>())
			.def("norm",&Vect3d::norm)
			.def(self + self)
			.def(self - self)
			.def(self / double())
			.def(self * double())
			.def(pow(self, int()))
			.def(pow(self, double()))
			.def(self_ns::str(self))
			.def("__getitem__", &getitem_Vect3d)
			.def("__setitem__", &setitem_Vect3d)
			;


	/*
	 * Particles
	 */
	class_<SpeciesType,typename Aboria::ptr<SpeciesType> >("Particles")
	        .def(boost::python::vector_indexing_suite<SpeciesType >())
	        .def("get_grid",&SpeciesType::get_grid)
	        .def("copy_from_vtk_grid",&SpeciesType::copy_from_vtk_grid)
	    ;


	/*
	 * Particle
	 */
#define ADD_PROPERTY(NUMBER) \
		.add_property(DataNames<SpeciesTuple>::get(NUMBER).c_str(), \
					make_function(&SpeciesType::value_type::get_data_elem<NUMBER>, \
									return_value_policy<copy_non_const_reference>()), \
					&SpeciesType::value_type::set_data_elem<NUMBER>)

	class_<SpeciesType::value_type, Aboria::ptr<SpeciesType::value_type> >("Particle",init<>())
		.add_property("id", &SpeciesType::value_type::get_id)
		.add_property("position",  make_function(&SpeciesType::value_type::get_position,
									return_value_policy<copy_const_reference>()),
							&SpeciesType::value_type::set_position)
		ADD_PROPERTY(SPECIES_VELOCITY)
		ADD_PROPERTY(SPECIES_ORIENTATION)
		ADD_PROPERTY(SPECIES_POTENTIAL)
		ADD_PROPERTY(SPECIES_AVERAGED_ORIENTATION)
		ADD_PROPERTY(SPECIES_VARIANCE_ORIENTATION)
		ADD_PROPERTY(SPECIES_AVERAGED_THETA)
		ADD_PROPERTY(SPECIES_THETA)
		ADD_PROPERTY(SPECIES_AVERAGED_POSITION)
		ADD_PROPERTY(SPECIES_FIXED)
		ADD_PROPERTY(SPECIES_NUM_MOVES)
		;

//	/*
//	 * ParticleSimulation
//	 */
//
//	class_<ParticleSimulation>("ParticleSimulation", init<>())
//			.def(init<double,double,double,double,double>(
//					(arg("Dtrans"),arg("Drot"),arg("Temp"),arg("dt"),arg("L"))))
//			.def("add_particles",&ParticleSimulation::add_particles)
//			.def("monte_carlo_timestep",&ParticleSimulation::monte_carlo_timestep<GayBernePotential>)
//			.add_property("Dtrans", &ParticleSimulation::getDtrans, &ParticleSimulation::setDtrans)
//			.add_property("Drot", &ParticleSimulation::getDrot, &ParticleSimulation::setDrot)
//			.add_property("Temp", &ParticleSimulation::getTemp, &ParticleSimulation::setTemp)
//			.add_property("dt", &ParticleSimulation::getDt, &ParticleSimulation::setDt)
//			.add_property("L", &ParticleSimulation::getL, &ParticleSimulation::setL)
//			.add_property("particles", &ParticleSimulation::getParticles, &ParticleSimulation::setParticles)
//			;

	/*
	 * Algorithms
	 */
	def("monte_carlo_timestep",&monte_carlo_timestep<SpeciesTuple,GayBernePotential>);
	def("monte_carlo_timestep",&monte_carlo_timestep<SpeciesTuple,HGOPotential>);
	def("monte_carlo_timestep",&monte_carlo_timestep<SpeciesTuple,LabwohlLasherPotential>);

	def("local_averaging",&local_averaging<SpeciesTuple>);


	/*
	 * GayBernePotential
	 */

	class_<GayBernePotential>("GayBernePotential", init<double, double, double, double, double, double>(
			(arg("sigma_s"),arg("k"),arg("kdash"),arg("mu"),arg("nu"),arg("epsilon_0"))))
			.def("evaluate",&GayBernePotential::evaluate)
			;

	class_<HGOPotential>("HGOPotential", init<double, double>(
			(arg("sigma_s"),arg("k"))))
			.def("evaluate",&HGOPotential::evaluate)
			;

	class_<LabwohlLasherPotential>("LabwohlLasherPotential", init<double, double>(
			(arg("epsilon"),arg("lattice_spacing"))))
			.def("evaluate",&LabwohlLasherPotential::evaluate)
			;


}
