#include "clvf/spacecraft.h"
#include "clvf/utils.h"

namespace clvf {
Spacecraft::Spacecraft(
  const Eigen::Matrix3d& inertia, 
  double mass, 
  double beta,
  const Eigen::Vector3d& o_hat_B,
  const Eigen::Vector3d& docking_port_B,
  double angle_of_acceptance
)
: inertia_{inertia}, 
mass_{mass}, 
inertia_inverse_{inertia.inverse()}, 
beta_{beta},
o_hat_B_{o_hat_B.normalized()},
docking_port_B_{docking_port_B},
angle_of_acceptance_{angle_of_acceptance}{};

Eigen::Vector3d Spacecraft::OmegaDot(
  const Eigen::Vector3d& omega, 
  const Eigen::Vector3d& tau_external) const {
    // Euler dynamics:
    return inertia_inverse_*(tau_external - omega.cross(inertia_*omega));
}

Eigen::Vector3d Spacecraft::RInertialDDot(
  const Eigen::Vector3d& r_inertial, 
  const Eigen::Vector3d& u_inertial,
  double mu) const {
    double r_magnitude_squared = r_inertial.dot(r_inertial);

    // Acceleration due to gravity and external force:
    return (-mu * r_inertial.normalized()/r_magnitude_squared + u_inertial/mass_);
}

Eigen::Matrix3d Spacecraft::RotationMatrixDot(
    const Eigen::Matrix3d& C,
    const Eigen::Vector3d& w
) const {
  return -clvf::Skew(w)*C;
}

Eigen::Vector3d Spacecraft::Control(
  const Eigen::Vector3d& desired_speed,
  const Eigen::Vector3d& actual_speed,
  const Eigen::Vector3d& desired_acceleration
) const {
  return mass_*(beta_ * (desired_speed - actual_speed) + desired_acceleration);
}

// Convert orbital elements to initial position and velocity:
std::pair<Eigen::Vector3d, Eigen::Vector3d> Spacecraft::OrbitalElementsToPosVel(
  double semi_major_axis,
  double eccentricity,
  double inclination,
  double RAAN,
  double argument_of_perigee,
  double true_anomaly
) {
    // First, getting distance "r" at this instant:
    double r = semi_major_axis*(1.0 - eccentricity*eccentricity)/(1+eccentricity*std::cos(true_anomaly));

    // Next, get the perifocal frame coordinates:
    Eigen::Vector3d rP({r*std::cos(true_anomaly), r*std::sin(true_anomaly), 0});
              
    // Next, get the rotation matrix from perifocal into ECI:
    Eigen::Matrix3d C_PI = clvf::C3(argument_of_perigee)*clvf::C1(inclination)*clvf::C3(RAAN);
    Eigen::Matrix3d C_IP = C_PI.transpose();
        
    // Get rI
    auto rI = C_IP*rP;
        
    // Get the semi-latus rectum (p)
    auto p = semi_major_axis*(1-eccentricity*eccentricity);
        
    // Get perifocal frame velocity:
    Eigen::Vector3d vP ({-sqrt(clvf::kMu/p)*std::sin(true_anomaly), sqrt(clvf::kMu/p)*(eccentricity + std::cos(true_anomaly)), 0 });

    // Get the inertial velocity:
    auto vI = C_IP*vP;

    // Return the position velocity pair:
    return std::make_pair(rI, vI);
}

}