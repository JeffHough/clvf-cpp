#include "clvf/clvf.h"
#include <cmath>
#include "clvf/utils.h"

namespace clvf {

double CLVF::SFunction(double r, double theta) const {
    if (r < alpha_){
        return ka_ * r / alpha_ * std::sin(theta);
    }
    return ka_ * alpha_ / r * std::sin(theta);
}

double CLVF::VFunction(double r, double theta) const {
    if (std::abs(alpha_ - r) < b_){
        return kc_ * (alpha_ - r)/b_;
    }
    return kc_ * clvf::Sign(alpha_ - r);
}

double CLVF::GFunction(double r) const {
    if (r < alpha_){
        return r;
    }
    return alpha_*alpha_/r;
}

double CLVF::AccelerationBound(double w_max, double w_dot_max) const {
    double tmp0 = kc_*kc_/b_;
    double tmp1 = std::pow(ka_ + alpha_*w_max, 2.0)/alpha_;
    double tmp2 = ka_*ka_/(2.0*alpha_) + ka_*w_max + alpha_*w_dot_max;

    double tmp3 = std::pow(tmp0 + tmp1, 2.0);
    double tmp4 = std::pow(tmp2, 2.0);

    return std::sqrt(tmp3 + tmp4);
}

Eigen::Vector3d CLVF::AHatVector(
    const Eigen::Vector3d& r_vector,
    const Eigen::Vector3d& o_hat_vector
) const {
    double theta = ThetaFromTwoVectors(r_vector, o_hat_vector);

    Eigen::Vector3d r_hat_vector = r_vector.normalized();

    constexpr double kSmallNumber = 0.00001;
    if (std::abs(theta) > kSmallNumber){
        return (o_hat_vector - r_hat_vector*std::cos(theta)) / std::sin(theta);
    }

    return Eigen::Vector3d::Zero();
}

Eigen::Vector3d CLVF::DesiredVelocity(
    const Eigen::Vector3d& r_vector,
    const Eigen::Vector3d& o_hat_vector,
    const Eigen::Vector3d& omega_OI,
    const Eigen::Vector3d& d_vector_dot
) const {
    double theta = ThetaFromTwoVectors(r_vector, o_hat_vector);
    double r = r_vector.norm();
    Eigen::Vector3d r_hat_vector = r_vector.normalized();

    return VFunction(r, theta)*r_hat_vector
            + SFunction(r, theta)*AHatVector(r_vector, o_hat_vector)
             + GFunction(r)*clvf::Skew(omega_OI)*r_hat_vector
              + d_vector_dot;
}
    
} // namespace clvf