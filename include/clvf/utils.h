#ifndef CLVF_UTILS_H_
#define CLVF_UTILS_H_

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>

namespace clvf {

constexpr double kSmallNumber = 1e-7;

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
    return std::acos(a.dot(b)/(a.norm()*b.norm() + clvf::kSmallNumber));
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

// An integrator for quaternions::
inline Eigen::Quaterniond EulerIntegrate(
    const Eigen::Quaterniond& derivative, 
    const Eigen::Quaterniond& previous_value,
    double dt
){
    // Just does Euler interation (for now)
    Eigen::Quaterniond result;
    result.x() = previous_value.x() + derivative.x()*dt;
    result.y() = previous_value.y() + derivative.y()*dt;
    result.z() = previous_value.z() + derivative.z()*dt;
    result.w() = previous_value.w() + derivative.w()*dt;

    result.normalize();
    return result;
}

// The Rotation matrix for the x-direction:
inline Eigen::Matrix3d C1(double x){
    double cx = std::cos(x);
    double sx = std::sin(x);
    return Eigen::Matrix3d({
        {1, 0, 0},
        {0, cx, sx},
        {0, -sx, cx}
    });
}

// The Rotation matrix for the y-direction:
inline Eigen::Matrix3d C2(double y){
    double cy = std::cos(y);
    double sy = std::sin(y);
    return Eigen::Matrix3d({
        {cy, 0, -sy},
        {0, 1, 0},
        {sy, 0, cy}
    });
}

// The Rotation matrix for the z-direction:
inline Eigen::Matrix3d C3(double z){
    double cz = std::cos(z);
    double sz = std::sin(z);
    return Eigen::Matrix3d({
        {cz, sz, 0},
        {-sz, cz, 0},
        {0, 0, 1}
    });
}

inline Eigen::Matrix3d NormalizeRotationMatrix(const Eigen::Matrix3d& C){
    Eigen::Quaterniond q_temp(C);
    q_temp.normalize();
    return q_temp.matrix();
}

inline Eigen::Vector3d RotateVectorByQuaternion(const Eigen::Vector3d& v, const Eigen::Quaterniond& quaternion){
    Eigen::Quaterniond tmp;
    tmp.x() = v(0);
    tmp.y() = v(1);
    tmp.z() = v(2);

    //auto rotated_v = quaternion*tmp*quaternion.conjugate();
    auto rotated_v = quaternion.conjugate()*tmp*quaternion;
    return Eigen::Vector3d({rotated_v.x(), rotated_v.y(), rotated_v.z()});
}

inline Eigen::Quaterniond QuaternionDerivative(const Eigen::Vector3d& omega, const Eigen::Quaterniond& quaternion){
    Eigen::Quaterniond tmp;
    tmp.x() = 0.5 *omega(0);
    tmp.y() = 0.5 *omega(1);
    tmp.z() = 0.5 *omega(2);
    tmp.w() = 0.0;

    return (tmp * quaternion);
}

}

#endif
