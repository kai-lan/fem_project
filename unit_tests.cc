#include <iostream>
#include "Mesh.hh"

void test_polar_decomp () {
    for (int i = 0; i < 30; i++) std::cout << "-";
    std::cout << "\nTest Polar Decomposition\n";
    for (int i = 0; i < 30; i++) std::cout << "-";
    std::cout << ""<<std::endl;
    srand((unsigned int) time(0));
    int n = 3;
    Eigen::MatrixXd F(n, n), R(n, n), S(n, n), I = Eigen::Matrix3d::Identity();
    F = Eigen::MatrixXd::Random(n, n);
    compute_polar_decomp_3d(F, R, S);
    auto print = [&](int n){
        std::cout << "Dimension: " << n << std::endl; 
        std::cout << "F\n" << F << std::endl; 
        std::cout << "|F - RS| = " << (F - R*S).norm() << std::endl; 
        std::cout << "det R = " << R.determinant() << std::endl; 
        std::cout << "|R^TR - I| = " << (R.transpose() * R - I).norm() << std::endl; 
    };
    print(n);
    n = 2;
    F = Eigen::MatrixXd::Random(n, n);
    R.resize(n, n); S.resize(n, n); 
    I = Eigen::Matrix2d::Identity();
    compute_polar_decomp_2d(F, R, S);
    print(n); 
}
void test_gradient() {
    Eigen::MatrixXd V0(3, 2), V(3, 2); 
    V0 << 0, 0, 2, 0, 1, 1;
    V = V0 + Eigen::MatrixXd::Random(3, 2);
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> perturb = Eigen::MatrixXd::Random(3, 2);
    Eigen::VectorXd perturb_flat = Eigen::Map<Eigen::VectorXd>(perturb.data(), perturb.size());
    Eigen::MatrixXi E(1, 3);
    E << 0, 1, 2;
    Mesh M(V0, E, 1);
    Eigen::VectorXd eps = Eigen::VectorXd::LinSpaced(100, -8, 0).unaryExpr([](double x) { return std::pow(10, x); });
    for(int i = 0; i < eps.size(); ++i) {
        double ep = eps(i);
        double err = std::abs((M.potential_energy(V + ep/2*perturb) - M.potential_energy(V - ep/2*perturb))/ep + M.force(V).dot(perturb_flat));
        std::cout << ep << " " << err << " ";
    }
    std::cout << "" << std::endl;
}

int main(int argc, char *argv[]) {
    if (argc == 1) std::cout << "Please specify a unit test number:\n"
                "0  polar decomposition\n"
                "1  gradient of energy (force)";
    std::srand((unsigned int) time(0));
    for (int i = 1; i < argc; ++i) {
        int type = std::stoi(argv[i]);
        switch (type)
        {
            case 0:
                test_polar_decomp();
                break;
            case 1:
                test_gradient();
            default:
                break;
        }
    }   
}