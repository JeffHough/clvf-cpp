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

Eigen::Matrix3d Spacecraft::RotationMatrixDot(
    const Eigen::Matrix3d& C,
    const Eigen::Vector3d& w
) const {
  return -Skew(w)*C;
}
  
Eigen::Matrix3d Spacecraft::Skew(
  const Eigen::Vector3d& v
) const {
  Eigen::Matrix3d skew_matrix;
  skew_matrix << 
    0,    -v(2),  v(1),
    v(2),   0,    v(0),
    v(1),  -v(0),   0;

  return skew_matrix;
}

}