#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

namespace py = pybind11;

#include "breakthrough.h"
#include "component.h"
#include "fitting.h"
#include "isotherm.h"
#include "mixture_prediction.h"
#include "multi_site_isotherm.h"

PYBIND11_MODULE(_ruptura, m)
{
  py::class_<Component>(m, "Component")
      .def(py::init<size_t, std::string, std::vector<Isotherm>, double, double, double, bool>())
      .def_readonly("GasPhaseMolFraction", &Component::Yi0)
      .def_readonly("MoleculeName", &Component::name)
      .def("__repr__", &Component::repr);
  py::class_<Isotherm>(m, "Isotherm")
      .def(py::init<std::string, std::vector<double>, size_t>())
      .def("__repr__", &Isotherm::repr);
  py::class_<MixturePrediction>(m, "MixturePrediction")
      .def(py::init<std::string, std::vector<Component>, size_t, size_t, double, double, double, size_t, size_t, size_t,
                    size_t>())
      .def("__repr__", &MixturePrediction::repr)
      .def("compute", &MixturePrediction::compute);
  py::class_<Breakthrough>(m, "Breakthrough")
      .def(py::init<std::string, std::vector<Component>, size_t, size_t, size_t, size_t, double, double, double, double,
                    double, double, double, double, size_t, bool, bool, double, const MixturePrediction>())
      .def("__repr__", &Breakthrough::repr)
      .def("compute", &Breakthrough::compute);
  py::class_<Fitting>(m, "Fitting")
      .def(py::init<std::string, std::vector<Component>, size_t>())
      .def("evaluate", &Fitting::evaluate)
      .def("compute", &Fitting::compute);
}
