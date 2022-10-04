/**
 * @file Simulation.hh
 * @author Kai Lan (kai.weixian.lan@gmail.com)
 * @brief A wrapper class to compute simulations
 * 
 * @date 2022-05-30
 */
#ifndef SIMULATION
#define SIMULATION
#include "Mesh.hh"
#include "FiniteDifference.hh"

class Simulation
{
    using VXd = Eigen::VectorXd;
    using MXd = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;
private:
    /* data */
public:
    Mesh& mesh;
    double stepSize;
    int numSteps;
    std::vector<double> timestamps;
    std::vector<VXd> velocities;
    std::vector<VXd> positions;

    Simulation(Mesh& M, double step, int steps) : mesh(M), stepSize(step), numSteps(steps) {
        timestamps.push_back(0); // Start time
        velocities.emplace_back(VXd::Zero(mesh.V0.size())); // Zero initial veloticies
        positions.push_back(mesh.flatten(mesh.V0));   // Initial positions
    }
    // TODO: only support Euler's method for now
    void simulate() {
        std::function<VXd (VXd, VXd)> next_v;
        std::function<VXd (double, VXd)> next_x;
        next_v = [this](VXd v_old, VXd x_old) { // a <- M^-1 * (f(x_old) + G)
            VXd v_new = v_old + stepSize 
            * this->mesh.massMatrix().cwiseInverse().asDiagonal()
            * (this->mesh.force(this->mesh.unflatten(x_old)) + this->mesh.weight());
            this->velocities.push_back(v_new);
            return v_new;
        };
        next_x = [this, &next_v] (double t, VXd x_old) {
            return next_v(this->velocities.back(), x_old);
        };
        
        double t = 0;
        for (int i = 0; i < numSteps; ++i) {
            //std::cout <<"i: " << i << std::endl;
            t += stepSize;
            timestamps.push_back(t);
            positions.push_back(FiniteDifference::euler_step(stepSize, t, positions.back(), next_x));
        }
    }

    ~Simulation() {}
};

#endif /* SIMULATION */
