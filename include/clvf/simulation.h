#ifndef CLVF_SIMULATION_H_
#define CLVF_SIMULATION_H_

#include "clvf/spacecraft.h"
#include "clvf/clvf.h"

namespace clvf {

struct SimulationData {
  // Dynamics data:
  Eigen::Vector3d target_orbital_position;
  Eigen::Vector3d target_orbital_velocity;
  Eigen::Vector3d target_omega;
  Eigen::Matrix3d target_rotation_matrix;
  Eigen::Vector3d chaser_relative_position;
  Eigen::Vector3d chaser_relative_velocity;

  // Guidance data:
  Eigen::Vector3d desired_speed;
  Eigen::Vector3d desired_acceleration;

  // Control data:
  Eigen::Vector3d control_vector;

  // Simulation complete:
  bool simulation_complete;
};


class Simulation {
  private:
    // elements of the simulation:
    CLVF clvf_;
    LVF lvf_;
    Spacecraft target_spacecraft_;
    Spacecraft chaser_spacecraft_;

    // an integer to keep track of how many steps were spent inside the CLVF switch region:
    int steps_in_switch_region_{0};

    // an integer to keep track of how many steps were spent in the end region of the LVF:
    int steps_in_end_region_{0};

    // a maximum number of steps which should be spent inside of the switch region:
    const int max_steps_in_switch_region_;

    // a maximum number of steps which is spent inside of the "end" region before completing.
    const int max_steps_in_end_region_;

  public:
    Simulation(
      const CLVF& clvf,
      const LVF& lvf,
      const Spacecraft& target_spacecraft,
      const Spacecraft& chaser_spacecraft,
      int max_steps_in_switch_region,
      int max_steps_in_end_region
    )
    :
    clvf_{clvf},
    lvf_{lvf},
    target_spacecraft_{target_spacecraft},
    chaser_spacecraft_{chaser_spacecraft},
    max_steps_in_switch_region_{max_steps_in_switch_region},
    max_steps_in_end_region_{max_steps_in_end_region}
    {};

    // Single step in the simulation:
    SimulationData Step(const SimulationData& sim_data_k);

    // Running the simulation:
    bool Run(const SimulationData& sim_data_initial);
};
}

#endif