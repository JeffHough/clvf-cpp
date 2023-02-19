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

  double dt = 0.1;
  Eigen::Matrix3d new_C = clvf::NormalizeRotationMatrix( clvf::EulerIntegrate(C_dot, q.matrix(), dt) );
  Eigen::Quaterniond q_new_1(new_C);

  Eigen::Quaterniond q_dot = clvf::QuaternionDerivative(omega, q);
  Eigen::Quaterniond q_new_2 = clvf::EulerIntegrate(q_dot, q, dt);

  std::cout << "q initial was:\n" << q << std::endl;
  std::cout << "q_update_1:\n" << q_new_1 << "\nq_update_2:\n" << q_new_2 << std::endl;

  // Spin a particular vector?
  Eigen::Vector3d v_B = {1.0, 0.0, 0.0};
  Eigen::Vector3d w_BI = {0.0, 1.0, 0.0};

  // Quaternion to put the vector into the INERTIAL frame:
  Eigen::Quaterniond q_BI;

  // A 90 degree rotation about y-axis:
  q_BI.x() = 0.0;
  q_BI.y() = 0.7071;
  q_BI.z() = 0.0;
  q_BI.w() = 0.7071;
  
  // Integrate v_B by pure rotation:
  Eigen::Vector3d v_I = clvf::RotateVectorByQuaternion(v_B, q_BI.conjugate());
  Eigen::Vector3d v_I_dot = w_BI.cross(v_I);

  Eigen::Vector3d v_I2_by_integration = v_I + dt*v_I_dot;
  auto q_BI_dot = clvf::QuaternionDerivative(w_BI, q_BI);
  Eigen::Quaterniond q_BI_new = clvf::EulerIntegrate(q_BI_dot, q_BI, dt);
  Eigen::Vector3d v_I2_by_quaterion_integration = clvf::RotateVectorByQuaternion(v_B, q_BI_new.conjugate());

  std::cout << "\noriginal vector (inertial): \n" << v_I << 
  "\nRotating by w_BI = \n" << w_BI <<
  "\nBy pure integration:\n" << v_I2_by_integration <<
  "\nBy quaternion integration:\n" << v_I2_by_quaterion_integration <<"\n";

}