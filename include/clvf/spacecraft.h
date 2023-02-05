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

    // Assuming this is a target spacecraft, it needs a docking port,
    // an orientation vector, and a docking cone half-angle of acceptance.
    const Eigen::Vector3d o_hat_B_;
    const Eigen::Vector3d docking_port_B_;
    const double angle_of_acceptance_;

    // The control parameter for the spacecraft control:
    const double beta_;

    // Where is the docking port in the body frame?
    const Eigen::Vector3d d_vector_;

  public:
    Spacecraft() = delete;
    Spacecraft(
      const Eigen::Matrix3d& inertia, 
      double mass, 
      double beta,
      const Eigen::Vector3d& o_hat_B,
      const Eigen::Vector3d& docking_port_B,
      double angle_of_acceptance);

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

    // Getter function for the D-Vector:
    const Eigen::Vector3d& DVectorB() const {return d_vector_;}

    // Getter for the O-hat-B vector:
    const Eigen::Vector3d& OHatB() const {return o_hat_B_;}

    // Getter for the angle of acceptance:
    double AngleOfAcceptance() const {return angle_of_acceptance_;}

    // Getter for the mass:
    double Mass() const {return mass_;}

    // Convert orbital elements to initial position and velocity:
    static std::pair<Eigen::Vector3d, Eigen::Vector3d> OrbitalElementsToPosVel(
      double semi_major_axis,
      double eccentricity,
      double inclination,
      double RAAN,
      double argument_of_latitude /*?*/,
      double true_anomaly
    );
  };

}

#endif