#ifndef CLVF_SPACECRAFT_H_
#define CLVF_SPACECRAFT_H_

#include <Eigen/Dense>

namespace clvf {

class Spacecraft {
  private:
    // Inertia matrix of the spacecraft:
    const Eigen::Matrix3d inertia_;
    const Eigen::Matrix3d inertia_inverse_;
    const double mass_;

    // The control parameter for the spacecraft control:
    const double beta_;

  public:
    Spacecraft() = delete;
    Spacecraft(const Eigen::Matrix3d& inertia, double mass, double beta);

    // Euler dynamics of the spacecraft:
    Eigen::Vector3d OmegaDot(
      const Eigen::Vector3d& omega, 
      const Eigen::Vector3d& tau_external) const;

    // Translational dynamics of the spacecraft:
    Eigen::Vector3d RInertialDDot(
      const Eigen::Vector3d& r_inertial, 
      const Eigen::Vector3d& u_external, 
      double mu) const;

    // The rotational kinematics:
    Eigen::Matrix3d RotationMatrixDot(
      const Eigen::Matrix3d& C,
      const Eigen::Vector3d& w
    ) const;

    // The controller:
    Eigen::Vector3d Control(
      const Eigen::Vector3d& desired_speed,
      const Eigen::Vector3d& actual_speed,
      const Eigen::Vector3d& desired_acceleration
    ) const;
  };

}

#endif