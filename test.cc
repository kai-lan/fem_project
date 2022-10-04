#include <iostream>
#include <Eigen/Dense>
using MXd = Eigen::MatrixXd;
using VXd = Eigen::VectorXd;
using MXdRM = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;

VXd func(const MXdRM& M) {return VXd::Map(M.data(), M.size());}
VXd func1(const MXd& M) {return VXd::Map(M.data(), M.size());}

int main() {
    Eigen::Matrix<double, -1, -1, Eigen::ColMajor> V(2, 3);
    V << 1, 2, 3,
        4, 5, 6;
    Eigen::VectorXd v = Eigen::VectorXd::Map(V.data(), V.size());
    std::cout << "Mapped vector " << v << std::endl;
    std::cout << "From func " << func1(V) << std::endl;
    std::cout << "1st row " << V.row(0) << std::endl;
    //f.call(f.str);
}