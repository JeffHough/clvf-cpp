#include "clvf/clvf.h"
#include "clvf/spacecraft.h"

int main(){
  // Set up a control parameter:
  double beta = 5.0; // [rads/s]
  // Set up a target spacecraft:
  Eigen::Matrix3d target_inertia = Eigen::Vector3d({3.0, 5.0, 7.0}).asDiagonal();
  double target_mass = 1.0;
  clvf::Spacecraft target(target_inertia, target_mass, beta);

  // Set up a chaser spacecraft:
  Eigen::Matrix3d chaser_inertia = Eigen::Matrix3d::Identity();
  double chaser_mass = 1.0;
  clvf::Spacecraft chaser(chaser_inertia, chaser_mass, beta);

  // Set up the CLVF and LVF parameters:

}