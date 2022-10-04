/**
 * @file finiteDifference.hh
 * @author Kai Lan (kai.weixian.lan@gmail.com)
 * @brief Implemented finite difference methods for discretized time domain.
 * 
 * @date 2022-05-30
 */
#ifndef FINITEDIFFERENCE
#define FINITEDIFFERENCE
#include <functional>
#include <Eigen/Dense>

using VXd = Eigen::VectorXd;

// y'(t) = f(t, y(t)), y(t_0) = y_0
// y_{n+1} = y_n + h * f(t_n, y_n)
namespace EulerFD
{   
    double euler_step(double h, double t, double y_old, std::function<double (double, double)> f) {
        return y_old + h * f(t, y_old);
    }
    VXd euler_step(double h, double t, VXd y_old, std::function<VXd (double, VXd)> f) {
        return y_old + h * f(t, y_old);
    }
}
// y_{n+1} = y_n + h * f(t_n + h/2, y_n + h/2 * f(t_n, y_n))
namespace MidpointFD
{
    double midpoint_step(double h, double t, double y_old, std::function<double (double, double)> f) {
        return y_old + h * f(t + h/2, EulerFD::euler_step(h/2, t, y_old, f));
    }
    VXd midpoint_step(double h, double t, VXd y_old, std::function<VXd (double, VXd)> f) {
        return y_old + h * f(t + h/2, EulerFD::euler_step(h/2, t, y_old, f));
    }
}
// y'(t) = f(t, y(t)), y(t_0) = y_0
// y_{n+1} = y_n + h * f(t_{n+1}, y_{n+1})
namespace BackwardEulerFD
{
    // TODO
}

namespace NewtonFD
{
    // TODO
}

namespace FiniteDifference 
{
    using namespace EulerFD;
    using namespace MidpointFD;
    using namespace BackwardEulerFD;
    using namespace NewtonFD;
}
#endif /* FINITEDIFFERENCE */
