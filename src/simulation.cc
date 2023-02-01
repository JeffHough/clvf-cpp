#include "clvf/simulation.h"

namespace clvf {

SimulationData Simulation::Step(
  const SimulationData& sim_data_k
) {
  // The data we will give OUT:
  SimulationData data_out;

  // Check if we should run the CLVF or LVF:
  auto should_run_CLVF = (steps_in_switch_region_ < max_steps_in_switch_region_);

  // Check if the simulation should end:
  auto should_end = (steps_in_end_region_ >= max_steps_in_end_region_);

  // Compute the desired speed of the chaser spacecraft:
  // if (should_run_CLVF){
  //   data_out.desired_speed = clvf_.DesiredVelocity(
  //     sim_data_k.chaser_relative_position,
  //     sim_data_k.o_hat_vector_CLVF,
  //     sim_data_k.target_omega,
  //     sim_data_k.d_vector_dot
  //   );

  //   data_out.desired_acceleration = clvf_.DesiredAcceleration(
  //     sim_data_k.chaser_relative_position,
  //     sim_data_k.chaser_relative_velocity,
  //     sim_data_k.o_hat_vector_CLVF,
  //     sim_data_k.target_omega,
  //     sim_data_k.target_omega_dot,
  //     sim_data_k.d_ddot
  //   );
  // } else {

  // }

  
  // Get the Orbital position of the chaser:
  auto chaser_orbital_position = sim_data_k.chaser_relative_position + sim_data_k.target_orbital_position;

  // Run
}

}