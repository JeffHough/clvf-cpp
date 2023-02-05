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

double CLVF::AccelerationBound(
    double omega_max, 
    double omega_and_omega_dot_max,
    double d_ddot_max
) const {
    double tmp0 = std::pow(0.7698*kc_*kc_/b_ + ka_*ka_/alpha_, 2.0);
    double factor = 2.0/(3*std::sqrt(3));
    double tmp1 = std::pow(ka_*ka_/(2*alpha_) + factor*ka_*omega_max, 2.0);

    double tmp2 = std::sqrt(tmp0 + tmp1);
    double tmp3 = 2.385*ka_*omega_max + omega_and_omega_dot_max*alpha_;

    return tmp2 + tmp3;
}

bool CLVF::InSwitchRange(
    const Eigen::Vector3d& r_vector,
    const Eigen::Vector3d& o_vector_hat
) const {
    double r = r_vector.norm();
    double theta = clvf::ThetaFromTwoVectors(r_vector, o_vector_hat);

    return ((std::abs(r - alpha_) < radius_error_before_changing_to_LVF_) && (std::abs(theta) < theta_error_before_changing_to_LVF_));
}

std::pair<Eigen::Vector3d, double> CLVF::OHatAndAlphaFromLVFValues(
    const Eigen::Vector3d& o_hat_B_LVF,
    double alhpa_LVF,
    const Eigen::Vector3d& d_vector_LVF
) {
    // Get the full vector:
    auto alpha_times_o_hat = d_vector_LVF + o_hat_B_LVF.normalized()*alhpa_LVF;

    // Extract the size and direction:
    return std::make_pair(alpha_times_o_hat.normalized(), alpha_times_o_hat.norm());
}

// LVF Library:
double LVF::ThetaN(double theta) const {
    return theta < theta_docking_cone_ ? (theta/theta_docking_cone_)*M_PI/2 : M_PI/2;
}

double LVF::RN(double r) const {
    return r >= alpha_prime_ ? final_angle_ : r/alpha_prime_ * final_angle_;
}

Eigen::Vector3d LVF::AHatVector(
    const Eigen::Vector3d& r_vector,
    const Eigen::Vector3d& o_hat_vector
) {
    Eigen::Vector3d r_hat_vector = r_vector.normalized();
    double theta = clvf::ThetaFromTwoVectors(r_hat_vector, o_hat_vector);
    return r_hat_vector.cross(o_hat_vector.cross(r_hat_vector))/(std::sin(theta) + clvf::kSmallNumber);
}

Eigen::Vector3d LVF::DesiredVelocity(
    const Eigen::Vector3d& r_vector,
    const Eigen::Vector3d& o_hat_vector,
    const Eigen::Vector3d& omega_OI,
    const Eigen::Vector3d& d_dot
) const {
    
    double r = r_vector.norm();
    Eigen::Vector3d r_hat = r_vector/r;
    double theta = clvf::ThetaFromTwoVectors(r_vector, o_hat_vector);
    
    Eigen::Vector3d a_hat = AHatVector(r_vector, o_hat_vector);

    auto theta_N = ThetaN(theta);
    auto r_N = RN(r);
    auto v_rel = VRel(r_N);
    auto vc = VFunction(v_rel, theta_N);
    auto sa = SFunction(v_rel, theta_N);
    auto g = GFunction(r);

    return vc*r_hat + sa*a_hat + g*omega_OI.cross(r_hat) + d_dot;

}

Eigen::Vector3d LVF::RHatDot(
    const Eigen::Vector3d& r_vector,
    const Eigen::Vector3d& velocity
) const {
    auto r_hat = r_vector.normalized();
    auto r = r_vector.norm();

    return (Eigen::Matrix3d::Identity() - r_hat*r_hat.transpose())/r * velocity;
}

Eigen::Vector3d LVF::AHatDot(
    const Eigen::Vector3d& a_hat,
    const Eigen::Vector3d& e_hat,
    const Eigen::Vector3d& o_hat,
    const Eigen::Vector3d& r_hat,
    const Eigen::Vector3d& r_hat_dot,
    const Eigen::Vector3d& omega_OI
) const {
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    double theta = clvf::ThetaFromTwoVectors(r_hat, o_hat);
    auto dahat_drhat = (I - a_hat*a_hat.transpose())/(std::sin(theta) + clvf::kSmallNumber)
        * (2*o_hat*r_hat.transpose()-r_hat*o_hat.transpose() - (r_hat.transpose()*o_hat)*I);

    auto dahat_dohat = (I - r_hat*r_hat.transpose())*(I - a_hat*a_hat.transpose())/(std::sin(theta) + clvf::kSmallNumber);

    return dahat_drhat*r_hat_dot + dahat_dohat * omega_OI.cross(o_hat);
}

double LVF::ThetaDot(
    const Eigen::Vector3d& o_hat,
    const Eigen::Vector3d& r_hat_dot,
    const Eigen::Vector3d& e_hat,
    const Eigen::Vector3d& omega_OI,
    const Eigen::Vector3d& r_hat,
    double theta
) const {
    return (-o_hat.dot(r_hat_dot) + (-r_hat.transpose()*(omega_OI.cross(o_hat))))/(std::sin(theta) + clvf::kSmallNumber);
}

double LVF::VFunctionDot(
    double theta,
    double r,
    double theta_dot, 
    double r_dot
) const {
    double dvreldr;
    auto r_N = RN(r);
    if (r <= alpha_prime_){
        dvreldr = v_max_*std::cos(r_N)*final_angle_/alpha_prime_;
    } else {
        dvreldr = 0;
    }
        
    auto theta_N = ThetaN(theta);
    auto dvcdr = -dvreldr*std::cos(theta_N);
    
    auto v_rel = VRel(r_N);

    double dvcdtheta;

    if (theta < theta_docking_cone_) {
        dvcdtheta = v_rel*std::sin(theta_N) * M_PI/(2*theta_docking_cone_);
    }
    else {
        dvcdtheta = 0;
    }
                    
    return dvcdr*r_dot + dvcdtheta*theta_dot;
}

double LVF::SFunctionDot(double theta, double r, double theta_dot, double r_dot) const {
    
    double dvreldr;
    auto r_N = RN(r);
    if (r <= alpha_prime_){
        dvreldr = v_max_*std::cos(r_N)*final_angle_/alpha_prime_;
    } else {
        dvreldr = 0;
    }
        
    auto theta_N = ThetaN(theta);
    auto dvcdr = -dvreldr*std::cos(theta_N);
    
    auto v_rel = VRel(r_N);
    
    
    auto dsadr = dvreldr*std::sin(theta_N);
    
    double dsadtheta;
    if (theta < theta_docking_cone_){
        dsadtheta = v_rel*std::cos(theta_N) * M_PI/(2*theta_docking_cone_);
    } else {
        dsadtheta = 0;
    }

    return dsadr*r_dot + dsadtheta*theta_dot;
}

Eigen::Vector3d LVF::DesiredAcceleration(
    const Eigen::Vector3d& r_vector,
    const Eigen::Vector3d& velocity,
    const Eigen::Vector3d& o_hat_vector,
    const Eigen::Vector3d& omega_OI,
    const Eigen::Vector3d& omega_dot_OI,
    const Eigen::Vector3d& d_ddot
) const {

    // Get all of the intermediate values:
    auto r = r_vector.norm();
    auto r_hat = r_vector.normalized();
    auto theta = clvf::ThetaFromTwoVectors(r_vector, o_hat_vector);
    auto a_hat = AHatVector(r_vector, o_hat_vector);

    auto r_dot = RDot(r_hat, velocity);
    auto r_hat_dot = RHatDot(r_vector, velocity);
    auto e_hat = EHatVector(r_vector, o_hat_vector);
    auto theta_dot = ThetaDot(o_hat_vector, r_hat_dot, e_hat, omega_OI, r_hat, theta);

    auto r_N = RN(r);
    auto v_rel = VRel(r_N);
    auto theta_N = ThetaN(theta);

    auto vc = VFunction(v_rel, theta_N);
    double sa = SFunction(v_rel, theta_N);
    auto vc_dot = VFunctionDot(theta, r, theta_dot, r_dot);

    auto sa_dot = SFunctionDot(theta, r, theta_dot, r_dot);
    auto a_hat_dot = AHatDot(a_hat, e_hat, o_hat_vector, r_hat, r_hat_dot, omega_OI);
    auto g_dot = GFunctionDot(r_dot);
    auto g = GFunction(r);

    // Compute the final value:
    return vc_dot*r_hat + vc*r_hat_dot + sa_dot*a_hat + sa*a_hat_dot + 
        g_dot*omega_OI.cross(r_hat) + g*omega_dot_OI.cross(r_hat) + 
            g*omega_OI.cross(r_hat_dot) + d_ddot;
}

double LVF::AccelerationBound(
    double omega_max, 
    double omega_and_omega_dot_max,
    double d_ddot_max
) const {



//   N1: dockingPortNorm
//   N2: rotNorm
//   N3: w_max for Coriolis stuff

    return 0.724612*(1 + M_PI/(2*theta_docking_cone_))*(v_max_*v_max_*final_angle_/alpha_prime_)+ 2*omega_max*v_max_ + d_ddot_max + omega_and_omega_dot_max*alpha_prime_;

}
    
} // namespace clvf