#include "clvf/spacecraft.h"

namespace clvf {
Spacecraft::Spacecraft(const Eigen::Matrix3d& inertia, double mass)
: inertia_{inertia}, mass_{mass}, inertia_inverse_{inertia.inverse()}{};

Eigen::Vector3d Spacecraft::OmegaDot(
  const Eigen::Vector3d& omega, 
  const Eigen::Vector3d& tau_external) const {
    // Euler dynamics:
    return inertia_inverse_*(tau_external - omega.cross3(inertia_*omega));
}

Eigen::Vector3d Spacecraft::RInertialDDot(
  const Eigen::Vector3d& r_inertial, 
  const Eigen::Vector3d& u_inertial,
  double mu) const {
    double r_magnitude_squared = r_inertial.dot(r_inertial);

    // Acceleration due to gravity and external force:
    return (-mu * r_inertial/r_magnitude_squared + u_inertial/mass_);
}


}