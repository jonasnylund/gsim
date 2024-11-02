#include <pybind11/pybind11.h>

#include "gsim/model/model.h"

PYBIND11_MODULE(gsim, m) {
  pybind11::class_<gsim::Model>(m, "Model")
    .def(pybind11::init<>())
    .def("step", &gsim::Model::step);
}