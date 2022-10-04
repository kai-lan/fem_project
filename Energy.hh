/**
 * @file Energy.hh
 * @author Kai Lan (kai.weixian.lan@gmail.com)
 * @brief
 *
 * @date 2022-05-11
 */
#ifndef ENERGY
#define ENERGY

#include <stdexcept>
#include <Eigen/Dense>
#include <cmath>

class Energy
{
    using MXd = Eigen::MatrixXd;
private:
    double mu=384.615, lambda=576.923; // Default: E = 1000, nu = 0.3
    int dim; //Dimension
public:
    Energy(int d) : dim(d) {
        if ((d != 2) && (d != 3)) throw std::runtime_error("Dimension can only be 2 or 3!");
    }
    Energy(int d, double E, double nu) : dim(d) {
        if ((d != 2) && (d != 3)) throw std::runtime_error("Dimension can only be 2 or 3!");
        mu = E / (2*(1 + nu));
        lambda = E * nu / ((1 + nu) * (1 - 2*nu));
    }
    // Compute elastic energy using corotation model: phi = mu * |F - R|^2 + lambda/2 * (det(F) - 1)^2
    double energy(const MXd &F);

    // Compute P (derivative of energy over F)
    MXd pressure(const MXd &F);

    ~Energy(){};
};

void compute_polar_decomp_2d(const Eigen::MatrixXd &F, Eigen::MatrixXd &R, Eigen::MatrixXd &S);
void compute_polar_decomp_3d(const Eigen::MatrixXd &F, Eigen::MatrixXd &R, Eigen::MatrixXd &S);
#endif /* ENERGY */
