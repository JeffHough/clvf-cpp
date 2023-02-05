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
    return (-mu * r_inertial/r_magnitude_squared + u_inertial/mass_);
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

}