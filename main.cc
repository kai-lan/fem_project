#include <iostream>
#include "Simulation.hh"

using MXd = Eigen::MatrixXd;
using MXi = Eigen::MatrixXi;
using VXi = Eigen::VectorXi;

//  --- --- --- --- ---
// | / | / | / | / | / |
//  --- --- --- --- ---
// | / | / | / | / | / |
//  --- --- --- --- ---
int main(int argc, char *argv[]) {
    MXd V0(15, 2); // Vertices
    V0 << 0, 0,
        0, 1,
        0, 2,
        1, 0,
        1, 1,
        1, 2,
        2, 0,
        2, 1,
        2, 2,
        3, 0,
        3, 1,
        3, 2,
        4, 0,
        4, 1,
        4, 2;
    MXi E(16, 3); // Elements (triangles)
    E << 0, 3, 4,
        0, 4, 1,
        1, 4, 5,
        1, 5, 2,
        3, 6, 7,
        3, 7, 4,
        4, 7, 8,
        4, 8, 5,
        6, 9, 10,
        6, 10, 7,
        7, 10, 11,
        7, 11, 8,
        9, 12, 13,
        9, 13, 10,
        10, 13, 14,
        10, 14, 11;
    VXi Dirichlet(3); // Left side is clamped to wall
    Dirichlet << 0, 1, 2;
    VXi Neumann(9); // Other boundary nodes are Neumann nodes
    Neumann << 3, 5, 6, 8, 9, 11, 12, 13, 15;
    MXd externalF(9, 2);
    externalF << 0, 0,  0, 0,   0, 0,   0, 0,   0, 0,   0, 0,  0, -5,   0, -5,  0, -5;
    double rho = 1;
    Mesh M(V0, E, rho, Dirichlet, Neumann, externalF);
    double step = 0.01;
    int numSteps = 100;
    Simulation sim(M, step, numSteps);
    sim.simulate();
    for (auto& pos : sim.positions)
        std::cout << pos.transpose() << std::endl; // Output as row vectors for readability
    //std::cout << M.V0 << std::endl;
    //std::cout << M.flatten(M.V0) << std::endl;
    //std::cout << M.weight() << std::endl;
}