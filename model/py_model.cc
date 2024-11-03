#include <pybind11/pybind11.h>

#include "gsim/common/numerical_types.h"
#include "gsim/model/model.h"

PYBIND11_MODULE(py_model, m) {
  pybind11::class_<gsim::Model>(m, "Model")
    .def(pybind11::init<>())
    .def(pybind11::init<float, int, float, float>())
    .def("step", &gsim::Model::step, "Step the simulation forward by the specified amount of time")
    .def("add_particle", &gsim::Model::addParticle, "Add a particle to the simulation")
    .def("initialize", &gsim::Model::initialize, "Construct the initial state given the existing particles")
    .def("random_particles", &gsim::Model::randomParticles, "Generates a set of random particles")
    .def("write_particles", &gsim::Model::writeParticles, "Write the current state to a file")
    .def("print_stats", &gsim::Model::printStats, "Prints simulation stats to the console")
    .def("get_time", &gsim::Model::getTime, "Gets the current time of the simulation");
}
