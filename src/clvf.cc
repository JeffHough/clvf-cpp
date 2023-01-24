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
    if (std::abs(r - alpha_) < b_){
        return clvf::Sign(r - alpha_)  *  (-kc_/(b_*b_)*(r-alpha_)*(r-alpha_)) - 2*kc_/b_*(r-alpha_);
    }
    return kc_ * clvf::Sign(alpha_ - r);
}

double CLVF::GFunction(double r) const {
    if (r < alpha_){
        return r;
    }
    return alpha_*alpha_/r;
}

Eigen::Vector3d CLVF::AHatVector(
    const Eigen::Vector3d& r_vector,
    const Eigen::Vector3d& o_hat_vector
) {
    double theta = ThetaFromTwoVectors(r_vector, o_hat_vector);

    Eigen::Vector3d r_hat_vector = r_vector.normalized();

    constexpr double kSmallNumber = 0.00001;
    if (std::abs(theta) > kSmallNumber){
        return (o_hat_vector - r_hat_vector*std::cos(theta)) / std::sin(theta);
    }

    return Eigen::Vector3d::Zero();
}

Eigen::Vector3d CLVF::EHatVector(
    const Eigen::Vector3d& r_vector,
    const Eigen::Vector3d& o_hat_vector
){
    auto r_hat_vector = r_vector.normalized();
    auto a_hat_vector = AHatVector(r_vector, o_hat_vector);
    return r_hat_vector.cross(a_hat_vector);
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

// Acceleration functions:
Eigen::Vector3d CLVF::RHatDot(const Eigen::Vector3d& r_vector, const Eigen::Vector3d& velocity) const {
    Eigen::Vector3d r_hat = r_vector.normalized();
    double r = r_vector.norm();

    return ((Eigen::Matrix3d::Identity() - r_hat*r_hat.transpose())/r) * velocity;
}

Eigen::Vector3d CLVF::AHatDot(
    const Eigen::Vector3d& a_hat,
    const Eigen::Vector3d& e_hat,
    const Eigen::Vector3d& o_hat,
    const Eigen::Vector3d& r_hat,
    const Eigen::Vector3d& r_hat_dot,
    const Eigen::Vector3d& omega_OI
) const {
    Eigen::Matrix3d dahat_drhat = 
        (Eigen::Matrix3d::Identity() - a_hat*a_hat.transpose())/(a_hat.dot(o_hat) + kSmallNumber)
        * (-r_hat*o_hat.transpose() - (r_hat.dot(o_hat))*Eigen::Matrix3d::Identity());

    Eigen::Matrix3d dahat_dohat = 
        (e_hat*e_hat.transpose())/(a_hat.dot(o_hat) + kSmallNumber);

    return dahat_drhat*r_hat_dot + dahat_dohat*(omega_OI.cross(o_hat));
}

double CLVF::RDot(
    const Eigen::Vector3d& r_hat,
    const Eigen::Vector3d& velocity
) const {
    return r_hat.dot(velocity);
}

double CLVF::ThetaDot(
    const Eigen::Vector3d& o_hat,
    const Eigen::Vector3d& r_hat_dot,
    const Eigen::Vector3d& e_hat,
    const Eigen::Vector3d& omega_OI,
    double theta
) const {
    return -o_hat.dot(r_hat_dot) / (std::sin(theta)+kSmallNumber) + omega_OI.dot(e_hat);
}

double CLVF::GFunctionDot(double r, double r_dot) const {
    if (r < alpha_){
    // derivative of "r":
    return r_dot;
    }

    // derivative of alpha_*alpha_/r; 
    return (-alpha_*alpha_/(r*r) * r_dot);
}

double CLVF::SFunctionDot(double theta, double r, double theta_dot, double r_dot) const {
    if (r < alpha_){
        // derivative of: ka_ * r / alpha_ * std::sin(theta);
        double dgdr = ka_/alpha_ * std::sin(theta);
        double dgdt = ka_*r/alpha_*std::cos(theta);

        return dgdr*r_dot + dgdt*theta_dot;
    }
    
    // derivative of ka_ * alpha_ / r * std::sin(theta);
    double dgdr = -ka_* alpha_/(r*r) * std::sin(theta);
    double dgdt = ka_* alpha_ /r * std::cos(theta);

    return dgdr*r_dot + dgdt*theta_dot;
}

double CLVF::VFunctionDot(double r, double r_dot) const {
    /* VFunction:
    if (std::abs(r - alpha_) < b_){
        return clvf::Sign(r - alpha_)  *  (-kc_/(b_*b_)*(r-alpha_)*(r-alpha_)) - 2*kc_/b_*(r-alpha_);
    }
    return kc_ * clvf::Sign(alpha_ - r);
    */
    if (std::abs(r - alpha_) < b_){
        return clvf::Sign(r - alpha_)  *  2*(-kc_/(b_*b_)*(r-alpha_)) - 2*kc_/b_;
    }
    return 0.0;  
}

Eigen::Vector3d CLVF::DesiredAcceleration(
    const Eigen::Vector3d& r_vector,
    const Eigen::Vector3d& velocity,
    const Eigen::Vector3d& o_hat_vector,
    const Eigen::Vector3d& omega_OI,
    const Eigen::Vector3d& omega_dot_OI,
    const Eigen::Vector3d& d_ddot
) const {
    
    auto r_hat = r_vector.normalized();
    double r = r_vector.norm();
    double r_dot = RDot(r_hat, velocity);
    double vc_dot = VFunctionDot(r, r_dot);
    double theta = ThetaFromTwoVectors(r_hat, o_hat_vector);
    double vc = VFunction(r, theta);
    auto r_hat_dot = RHatDot(r_vector, velocity);
    auto e_hat = EHatVector(r_vector, o_hat_vector);
    auto a_hat = AHatVector(r_vector, o_hat_vector);
    double sa = SFunction(r, theta);
    auto a_hat_dot = AHatDot(a_hat, e_hat, o_hat_vector, r_hat, r_hat_dot, omega_OI);
    double g = GFunction(r);
    double g_dot = GFunctionDot(r, r_dot);
    double theta_dot = ThetaDot(o_hat_vector, r_hat_dot, e_hat, omega_OI, theta);
    double sa_dot = SFunctionDot(theta, r, theta_dot, r_dot);
    
    
    return 
        vc_dot * r_hat + vc*r_hat_dot + sa_dot*a_hat + sa*a_hat_dot + 
            g_dot*Skew(omega_OI)*r_hat + g*Skew(omega_dot_OI)*r_hat + g*Skew(omega_OI)*r_hat_dot
                + d_ddot;
}

Eigen::Vector3d& CLVF::Controller(
      const Eigen::Vector3d& v_desired,
      const Eigen::Vector3d& v_actual,
      const Eigen::Vector3d& desired_acceleration
    ) const {
        return beta_*(v_desired - v_actual) + desired_acceleration;
    }
    
} // namespace clvf