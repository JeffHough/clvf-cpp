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

    Eigen::Vector3d CirculationComponent(
      const Eigen::Vector3d& position
    ) const;

    Eigen::Vector3d ContractionComponent(
      const Eigen::Vector3d& position
    ) const;

  public:
    CLVF(
      double kc,
      double ka,
      double b
    ) : kc_{kc}, ka_{ka}, b_{b}{};

    Eigen::Vector3d DesiredVelocity(
      const Eigen::Vector3d& position
    ) const;

    Eigen::Vector3d DesiredAcceleration(
      const Eigen::Vector3d& position
    ) const;

};

// Controller based on CLVFs:
class CLVFController {

};
}

#endif