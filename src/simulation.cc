#include "clvf/simulation.h"
#include "clvf/utils.h"

namespace clvf {

SimulationData Simulation::Step(
  const SimulationData& sim_data_k
) {
  // The data we will give OUT:
  SimulationData data_out;

  // Extract the C_BI from previous sim data:
  auto C_BI = sim_data_k.target_rotation_matrix;

  // Extract the docking port in the inertial-frame:
  auto target_d_vector_I = C_BI.transpose() * target_spacecraft_.DVectorB();

  // Undisturbed target Euler dynamics, put back into the inertial frame:
  auto target_omega_dot_OI = C_BI.transpose() * target_spacecraft_.OmegaDot(
    C_BI*sim_data_k.target_omega, 
    Eigen::Vector3d::Zero());

  // The velocity of the docking port as viewed in the inertial frame:
  auto target_d_dot_I = 
    sim_data_k.target_omega.cross(target_d_vector_I);

  // The acceleration of the docking port as viewed in the inertial frame:
  auto target_d_ddot_I = 
    target_omega_dot_OI.cross(target_d_vector_I) + 
      sim_data_k.target_omega.cross(sim_data_k.target_omega.cross(target_d_vector_I));

  // Check if we should run the CLVF:
  auto should_run_CLVF = (steps_in_switch_region_ < max_steps_in_switch_region_);

  // Switch case 1: whether we should run the CLVF?
  if (should_run_CLVF) {

    // Record that we are in the CLVF:
    data_out.in_CLVF = true;

    auto target_o_hat_vector_I = C_BI.transpose() * clvf_.OHatB();

    // Compute the desired speed of the chaser:
    data_out.desired_speed = clvf_.DesiredVelocity(
      sim_data_k.chaser_relative_position,
      target_o_hat_vector_I,
      sim_data_k.target_omega,
      target_d_vector_I
    );

    // Compute the desired acceleration:
    data_out.desired_acceleration = clvf_.DesiredAcceleration(
      sim_data_k.chaser_relative_position,
      sim_data_k.chaser_relative_velocity,
      target_o_hat_vector_I,
      sim_data_k.target_omega,
      target_omega_dot_OI,
      target_d_ddot_I
    );

    // Check if we increment/reset the switch value of the CLVF:
    if (clvf_.InSwitchRange(sim_data_k.chaser_relative_position, target_o_hat_vector_I)){
      steps_in_switch_region_++;
    } else {
      steps_in_switch_region_ = 0;
    }

  } else {
    // Record that we are in the LVF:
    data_out.in_CLVF = false;

    // Get the o_hat_vector from the LVF:
    auto target_o_hat_vector_I = C_BI.transpose() * target_spacecraft_.OHatB();

    // Compute the desired speed and acceleration according to the LVF:
    data_out.desired_speed = lvf_.DesiredVelocity(
      sim_data_k.chaser_relative_position,
      target_o_hat_vector_I,
      sim_data_k.target_omega,
      target_d_dot_I
    );

    data_out.desired_acceleration = lvf_.DesiredAcceleration(
      sim_data_k.chaser_relative_position,
      sim_data_k.chaser_relative_velocity,
      target_o_hat_vector_I,
      sim_data_k.target_omega,
      target_omega_dot_OI,
      target_d_ddot_I
    );

    // Check if we need to increment/reset the "end simulation" counter:
    if (lvf_.InEndRange(sim_data_k.chaser_relative_position, target_d_vector_I)){
      steps_in_end_region_++;
    } else {
      steps_in_end_region_ = 0;
    }

  }

  // Compute the remaining control and dynamics:
  data_out.control_vector = chaser_spacecraft_.Control(
    data_out.desired_speed, 
    sim_data_k.chaser_relative_velocity,
    data_out.desired_acceleration);

  // compute the change in orbital position:
  Eigen::Vector3d chaser_orbital_position = sim_data_k.target_orbital_position + sim_data_k.chaser_relative_position * clvf::kMetersToKm;
  Eigen::Vector3d chaser_orbit_velocity = sim_data_k.target_orbital_velocity + sim_data_k.chaser_relative_velocity * clvf::kMetersToKm;
  
  // Compute the chaser orbital acceleration:
  auto chaser_r_ddot = chaser_spacecraft_.RInertialDDot(
    chaser_orbital_position,
    data_out.control_vector * clvf::kMetersToKm,
    clvf::kMu
  );

  // Update chaser:
  auto chaser_orbital_velocity = clvf::EulerIntegrate<3,1>(chaser_r_ddot, chaser_orbit_velocity, dt_);
  chaser_orbital_position = clvf::EulerIntegrate<3,1>(chaser_orbital_velocity, chaser_orbital_position, dt_);

  // Update for the target:
  auto target_r_ddot = target_spacecraft_.RInertialDDot(sim_data_k.target_orbital_position, Eigen::Vector3d::Zero(), clvf::kMu);
  data_out.target_orbital_velocity = clvf::EulerIntegrate<3,1>(target_r_ddot, data_out.target_orbital_velocity, dt_);
  data_out.target_orbital_position = clvf::EulerIntegrate<3,1>(data_out.target_orbital_velocity, data_out.target_orbital_position, dt_);

  // record new chaser values in [m] of relative position:
  data_out.chaser_relative_position = (chaser_orbital_position - data_out.target_orbital_position) * clvf::kKmToMeters;
  data_out.chaser_relative_velocity = (chaser_orbital_velocity - data_out.target_orbital_velocity) * clvf::kKmToMeters;

  // Update the target rotation matrix:
  // TODO - SHOULD DO THIS IN TERMS OF QUATERNIONS
  data_out.target_rotation_matrix = clvf::EulerIntegrate<3,3>(
    target_spacecraft_.RotationMatrixDot(sim_data_k.target_rotation_matrix, sim_data_k.target_omega), sim_data_k.target_rotation_matrix, dt_
  );

  // Check if the simulation should end:
  auto should_end = (steps_in_end_region_ >= max_steps_in_end_region_);

  // Switch case 0: simulation is already over:
  if (should_end){
    data_out.simulation_complete = true;
  }

  // Increment the time:
  data_out.time = sim_data_k.time + dt_;

  // Finally, end the step:
  return data_out;
}

// TODO - WOULD IT BE EASIER TO RECORD MORE DATA?? I THINK YES.

}