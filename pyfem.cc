#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "Mesh.hh"

namespace py = pybind11;

PYBIND11_MODULE(pyfem, m) {
    using MXd  = Eigen::MatrixXd;
    using MXi  = Eigen::MatrixXi;
    using VXd  = Eigen::VectorXd;
    using VXi  = Eigen::VectorXi;
    
    py::class_<Mesh>(m, "Mesh")
        .def(py::init<const MXd, const MXi, const double, const VXi, const VXi, const VXd>(), 
        py::arg("V0"), py::arg("E"), py::arg("rho"), py::arg("Dirichlet"), py::arg("Neumann"), py::arg("externalForces"))
        .def("potential_energy", &Mesh::potential_energy, py::arg("x"))
        .def("force", &Mesh::force, py::arg("x"))
        ;
}