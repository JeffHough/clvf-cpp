#include "clvf/utils.h"
#include <iostream>
#include "clvf/spacecraft.h"

int main(){
  // Define some vector:
  Eigen::Vector3d v = {1.0, -2.0, 3.7};

  // Define some quaternion:
  const double theta = 0.5;
  const Eigen::Vector3d rotation_direction = {1.0, -3.0, 5.5};
  auto unit_rotation_vector = rotation_direction.normalized();
  double w = std::cos(theta/2);
  auto rotation_vector = std::sin(theta/2) * unit_rotation_vector;

  // The quaternion:
  Eigen::Quaterniond q(w, rotation_vector(0), rotation_vector(1), rotation_vector(2));

  // Rotate the vector:
  auto rotated_vector_1 = clvf::RotateVectorByQuaternion(v, q);

  // Can I now integrate the rotation matrix/quaternion??
  Eigen::Vector3d omega = {-1.0, 0.5, 0.9};
  Eigen::Matrix3d C_dot = -clvf::Skew(omega)*q.matrix();

  double dt = 0.001;
  Eigen::Matrix3d new_C = clvf::NormalizeRotationMatrix( clvf::EulerIntegrate(C_dot, q.matrix(), dt) );
  Eigen::Quaterniond q_new_1(new_C);

  Eigen::Quaterniond q_dot = clvf::QuaternionDerivative(omega, q);
  Eigen::Quaterniond q_new_2 = clvf::EulerIntegrate(q_dot, q, dt);

  std::cout << "q_update_1:\n" << q_new_1 << "\nq_update_2:\n" << q_new_2 << std::endl;
  
}