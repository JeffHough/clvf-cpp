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
    Eigen::Vector3d AHatVector(
      const Eigen::Vector3d& r_vector,
      const Eigen::Vector3d& o_hat_vector
    ) const;

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

    double AccelerationBound(double w_max, double w_dot_max) const;
};

// Controller based on CLVFs:
class CLVFController {

};
}

#endif