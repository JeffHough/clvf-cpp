#ifndef CLVF_UTILS_H_
#define CLVF_UTILS_H_

#include <Eigen/Dense>
#include <cmath>

namespace clvf {

// Some constants:
constexpr double kMu = 3.986*1e5; // [km^3/s^2]
constexpr double kEarthRadius = 6378.137; // [km]

constexpr double kMetersToKm = 0.001;
constexpr double kKmToMeters = 1000;

constexpr double kD2R = M_PI/180.0;
constexpr double kR2D = 1.0/kD2R;

// Skew symmetric matrix of a vector:
inline Eigen::Matrix3d Skew(
    const Eigen::Vector3d& v
){
  Eigen::Matrix3d skew_matrix;
  skew_matrix << 
    0,    -v(2),  v(1),
    v(2),   0,    v(0),
    v(1),  -v(0),   0;

  return skew_matrix;
};

// Utility funciton, because there isn't a std::sign function in c++:
inline double Sign(double input){
    if (input < 0.0){
        return -1.0;
    }

    if (input > 0.0) {
        return 1.0;
    }

    return 0.0;
}

// Find the angle between two vectors:
inline double ThetaFromTwoVectors(
    const Eigen::Vector3d& a, 
    const Eigen::Vector3d& b
){
    // NOTE: Assumes that neither a or b are the null vector {0,0,0}.
    return std::acos(a.normalized().dot(b.normalized()));
}

// An integrator:
template <int N, int M>
Eigen::Matrix<double, N, M> EulerIntegrate(
    const Eigen::Matrix<double, N, M>& derivative, 
    const Eigen::Matrix<double, N, M>& previous_value,
    double dt
){
    // Just does Euler interation (for now)
    return previous_value + derivative*dt;
}

}

#endif
