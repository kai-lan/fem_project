#include "Energy.hh"

using MXd = Eigen::MatrixXd;
using VXd = Eigen::VectorXd;

// Compute the polar decomposition F = RS, with R rotation and S symmetric for A dimension 3 x 3
void compute_polar_decomp_2d(const MXd &F, MXd &R, MXd &S) {
    R = MXd::Zero(2, 2), S = MXd::Zero(2, 2);
    double denominator = std::sqrt(std::pow(F(1, 0) - F(0, 1), 2) + std::pow(F(0, 0) + F(1, 1), 2));
    R(0, 0) = (F(0, 0) + F(1, 1)) / denominator;
    R(1, 1) = (F(0, 0) + F(1, 1)) / denominator;
    R(0, 1) = (F(0, 1) - F(1, 0)) / denominator;
    R(1, 0) = (F(1, 0) - F(0, 1)) / denominator;
    S = R.transpose() * F;
}

// Compute the polar decomposition F = RS, with R rotation and S symmetric for A dimension 3 x 3
void compute_polar_decomp_3d(const MXd &F, MXd &R, MXd &S) {
    MXd U, V;
    VXd sigma;
    auto svd = F.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
    U = svd.matrixU();
    V = svd.matrixV();
    sigma = svd.singularValues();
    // Edit if det(R) = 1, ie, det(U) * det(V) = -1
    if (U.determinant() * V.determinant() < 0) {
        sigma(2) *= -1;
        V.col(2) *= -1;
    }
    S = V * sigma.asDiagonal() * V.transpose();
    R = U * V.transpose();
}
void compute_polar_decomp(const MXd &F, MXd &R, MXd &S, const int dim) {
    if (dim == 2) compute_polar_decomp_2d(F, R, S);
    else          compute_polar_decomp_3d(F, R, S);
}
// Compute elastic energy using corotation model: phi = mu * |F - R|^2 + lambda/2 * (det(F) - 1)^2
double Energy::energy(const MXd &F) {
    MXd R, S;
    compute_polar_decomp(F, R, S, dim);
    return mu * std::pow((F - R).norm(), 2) + lambda/2 * std::pow((F.determinant() - 1), 2);
}

// Compute P (derivative of energy over F): P = 2 * mu * (F - R) + lambda * (J - 1) * J * F^{-T}
MXd Energy::pressure(const MXd &F) {
    MXd R, S;
    compute_polar_decomp(F, R, S, dim);
    double J = F.determinant();
    return 2 * mu * (F - R) + lambda * J * (J - 1) * F.inverse().transpose();
}