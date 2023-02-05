#include "clvf/simulation.h"
#include "clvf/utils.h"
#include <fstream>

namespace clvf {

SimulationData Simulation::Step(
  const SimulationData& sim_data_k
) const {
  // The data we will give OUT:
  SimulationData data_out;

  // Check if we should run the CLVF:
  auto should_run_CLVF = (sim_data_k.steps_in_switch_region < max_steps_in_switch_region_);

  // convenience:
  const auto& C_BI = sim_data_k.target_C_BI;

  // Switch case 1: whether we should run the CLVF?
  if (should_run_CLVF) {

    // Record that we are in the CLVF:
    data_out.in_CLVF = true;

    // Compute the desired speed of the chaser:
    data_out.desired_speed = clvf_.DesiredVelocity(
      sim_data_k.chaser_relative_position,
      sim_data_k.target_o_hat_vector_I_CLVF,
      sim_data_k.target_omega,
      sim_data_k.target_d_vector_I
    );

    // Compute the desired acceleration:
    data_out.desired_acceleration = clvf_.DesiredAcceleration(
      sim_data_k.chaser_relative_position,
      sim_data_k.chaser_relative_velocity,
      sim_data_k.target_o_hat_vector_I_CLVF,
      sim_data_k.target_omega,
      sim_data_k.target_omega_dot_OI,
      sim_data_k.target_d_ddot_I
    );

    // Check if we increment/reset the switch value of the CLVF:
    if (clvf_.InSwitchRange(sim_data_k.chaser_relative_position, sim_data_k.target_o_hat_vector_I_CLVF)){
      data_out.steps_in_switch_region = sim_data_k.steps_in_switch_region + 1;
    } else {
      data_out.steps_in_switch_region = 0;
    }

  } else {
    // Record that we are in the LVF:
    data_out.in_CLVF = false;

    // Compute the desired speed and acceleration according to the LVF:
    data_out.desired_speed = lvf_.DesiredVelocity(
      sim_data_k.chaser_relative_position,
      sim_data_k.target_o_hat_vector_I_LVF,
      sim_data_k.target_omega,
      sim_data_k.target_d_dot_I
    );

    data_out.desired_acceleration = lvf_.DesiredAcceleration(
      sim_data_k.chaser_relative_position,
      sim_data_k.chaser_relative_velocity,
      sim_data_k.target_o_hat_vector_I_LVF,
      sim_data_k.target_omega,
      sim_data_k.target_omega_dot_OI,
      sim_data_k.target_d_ddot_I
    );

    // Check if we need to increment/reset the "end simulation" counter:
    if (lvf_.InEndRange(sim_data_k.chaser_relative_position, sim_data_k.target_d_vector_I)){
      data_out.steps_in_end_region = sim_data_k.steps_in_end_region + 1;
    } else {
      data_out.steps_in_end_region = 0;
    }

  }

  // Update chaser orbit velocity and position:
  data_out.chaser_orbital_velocity = clvf::EulerIntegrate<3,1>(sim_data_k.chaser_orbital_acceleration, sim_data_k.chaser_orbital_velocity, dt_);
  data_out.chaser_orbital_position = clvf::EulerIntegrate<3,1>(sim_data_k.chaser_orbital_velocity, sim_data_k.chaser_orbital_position, dt_);

  // Update for the target orbit dynamics:
  data_out.target_orbital_velocity = clvf::EulerIntegrate<3,1>(sim_data_k.target_orbital_acceleration, sim_data_k.target_orbital_velocity, dt_);
  data_out.target_orbital_position = clvf::EulerIntegrate<3,1>(sim_data_k.target_orbital_velocity, sim_data_k.target_orbital_position, dt_);

  // Integrate to get the new target omega:
  data_out.target_omega = clvf::EulerIntegrate(sim_data_k.target_omega_dot_OI, sim_data_k.target_omega, dt_);

  // Update the target rotation matrix:
  // TODO - SHOULD DO THIS IN TERMS OF QUATERNIONS
  data_out.target_C_BI = clvf::EulerIntegrate<3,3>(sim_data_k.target_C_BI_dot, sim_data_k.target_C_BI, dt_);

  // Check if the simulation should end:
  auto should_end = ((data_out.steps_in_end_region >= max_steps_in_end_region_) || (data_out.time >= max_time_));

  // Switch case 0: simulation is already over:
  if (should_end){
    data_out.simulation_complete = true;
  }

  // Increment the time:
  data_out.time = sim_data_k.time + dt_;

  // this will compute all of the control vectors, accelerations, relative positions, etc.
  data_out.ComputeDependentData(
    clvf_.OHatB(),
    chaser_spacecraft_,
    target_spacecraft_
  );

  // Finally, end the step:
  return data_out;
}

// Run a full simulation:
void Simulation::Run(
  const SimulationData& sim_data_initial,
  const std::string& data_file_name
) const {
  // Create the sime data:
  SimulationData sim_data = sim_data_initial;

  // Make sure the dependent data is already calculated:
  sim_data.ComputeDependentData(clvf_.OHatB(), chaser_spacecraft_, target_spacecraft_);

  // Create a data file to log the simulation:
  std::ofstream data_stream;
  data_stream.open(data_file_name);

  LogHeaders(data_stream);

  while(!sim_data.simulation_complete){
    LogData(data_stream, sim_data);
    sim_data = Step(sim_data);
  }
  LogData(data_stream, sim_data);
  data_stream.close();
}

}