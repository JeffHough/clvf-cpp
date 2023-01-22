#ifndef CLVF_CLVF_H_
#define CLVF_CLVF_H_

#include <Eigen/Dense>

namespace clvf {

// Guidance law based on CLVFs:
class CLVF {
  private:
    // Parameters of the vector field:
    const double kc_;
    const double ka_;
    const double b_;
    const double alpha_;

    // The circulation component:
    double SFunction(
      double r,
      double theta
    ) const;

    // The contraction component:
    double VFunction(
      double r,
      double theta
    ) const;

    // The g-function to compensate for target rotation:
    double GFunction(
      double r
    ) const;

    // Computation for alignment direction:
    static Eigen::Vector3d AHatVector(
      const Eigen::Vector3d& r_vector,
      const Eigen::Vector3d& o_hat_vector
    );

    // Computation for the e-vector:
    static Eigen::Vector3d EHatVector(
      const Eigen::Vector3d& r_vector,
      const Eigen::Vector3d& o_hat_vector
    );

    // BELOW: several intermediate derivatives to compute the desired acceleration:
    // double ThetaDot(
    //   const Eigen::Vector3d& r_vector,
    //   const Eigen::Vector3d& o_hat_vector,
    //   const Eigen::Vector3d& omega_OI
    // ) const;

  public:
    // No default constructor:
    CLVF() = delete;
    CLVF(
      double kc,
      double ka,
      double b,
      double alpha
    ) : kc_{kc}, ka_{ka}, b_{b}, alpha_{alpha}{};

    Eigen::Vector3d DesiredVelocity(
      const Eigen::Vector3d& r_vector,
      const Eigen::Vector3d& o_hat_vector,
      const Eigen::Vector3d& omega_OI,
      const Eigen::Vector3d& d_vector_dot
    ) const;

    Eigen::Vector3d DesiredAcceleration(
      const Eigen::Vector3d& r_vector
    ) const;

    // TODO - THIS IS NOT CORRECT! NEED THE 3D VERSION.
    // double AccelerationBound(double w_max, double w_dot_max) const;
};

// TODO - CREATE THE DESIRED ACCELERATION FUNCTION.

// Controller based on CLVFs:
class CLVFController {

};
}

#endif