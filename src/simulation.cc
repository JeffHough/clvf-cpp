#include "clvf/simulation.h"
#include "clvf/utils.h"
#include <fstream>
#include <vector>
#include <string>

namespace clvf {

SimulationData Simulation::Step(
  const SimulationData& sim_data_k,
  std::shared_ptr<IntegratorBase<3,1>> chaser_orbital_velocity_integrator,
  std::shared_ptr<IntegratorBase<3,1>> chaser_orbital_position_integrator,
  std::shared_ptr<IntegratorBase<3,1>> target_orbital_velocity_integrator,
  std::shared_ptr<IntegratorBase<3,1>> target_orbital_position_integrator,
  std::shared_ptr<IntegratorBase<3,1>> target_omega_integrator,
  std::shared_ptr<QuaternionIntegratorBase> target_q_BI_integrator
) const {
  // Hold some memory to prevent resetting back to CLVF:
  SimulationData data_out = sim_data_k;

  // Check if we should run the CLVF:
  auto should_run_CLVF = (sim_data_k.steps_in_switch_region < max_steps_in_switch_region_);

  // convenience:
  const auto& q_BI = sim_data_k.target_q_BI;

  // Switch case 1: whether we should run the CLVF?
  if (should_run_CLVF) {

    // Record that we are in the CLVF:
    data_out.in_CLVF = true;

    // Compute the desired speed of the chaser:
    data_out.desired_speed = clvf_.DesiredVelocity(
      sim_data_k.chaser_relative_position,
      sim_data_k.target_o_hat_vector_I_CLVF,
      sim_data_k.target_omega,
      sim_data_k.target_d_dot_I
    );

    // Compute the desired acceleration:
    data_out.desired_acceleration = (data_out.desired_speed - sim_data_k.desired_speed)/dt_;
    
    /*clvf_.DesiredAcceleration(
      sim_data_k.chaser_relative_position,
      sim_data_k.chaser_relative_velocity,
      sim_data_k.target_o_hat_vector_I_CLVF,
      sim_data_k.target_omega,
      sim_data_k.target_omega_dot_OI,
      sim_data_k.target_d_ddot_I
    );*/

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
      sim_data_k.chaser_relative_position-sim_data_k.target_d_vector_I,
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
  data_out.chaser_orbital_velocity = chaser_orbital_velocity_integrator->Integrate(sim_data_k.chaser_orbital_acceleration);
  data_out.chaser_orbital_position = chaser_orbital_position_integrator->Integrate(data_out.chaser_orbital_velocity);

  // Update for the target orbit dynamics:
  data_out.target_orbital_velocity = target_orbital_velocity_integrator->Integrate(sim_data_k.target_orbital_acceleration);
  data_out.target_orbital_position = target_orbital_position_integrator->Integrate(data_out.target_orbital_velocity);

  // Integrate to get the new target omega:
  data_out.target_omega = target_omega_integrator->Integrate(sim_data_k.target_omega_dot_OI);

  // Update quaternion:
  data_out.target_q_BI = target_q_BI_integrator->Integrate(sim_data_k.target_q_BI_dot);

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
    clvf_.Alpha(),
    chaser_spacecraft_,
    target_spacecraft_,
    sim_data_k.chaser_relative_velocity
  );

  // Finally, end the step:
  return data_out;
}

SimulationData Simulation::KinematicStep(
  const SimulationData& sim_data_k,
  std::shared_ptr<IntegratorBase<3,1>> chaser_relative_position_integrator,
  std::shared_ptr<IntegratorBase<3,1>> target_omega_integrator,
  std::shared_ptr<QuaternionIntegratorBase> target_q_BI_integrator
) const {
  // Hold some memory to prevent resetting back to CLVF:
  SimulationData data_out = sim_data_k;

  // Check if we should run the CLVF:
  auto should_run_CLVF = (sim_data_k.steps_in_switch_region < max_steps_in_switch_region_);

  // convenience:
  const auto& q_BI = sim_data_k.target_q_BI;

  // Switch case 1: whether we should run the CLVF?
  if (should_run_CLVF) {

    // Record that we are in the CLVF:
    data_out.in_CLVF = true;

    // Compute the desired speed of the chaser:
    data_out.desired_speed = clvf_.DesiredVelocity(
      sim_data_k.chaser_relative_position,
      sim_data_k.target_o_hat_vector_I_CLVF,
      sim_data_k.target_omega,
      sim_data_k.target_d_dot_I
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
      sim_data_k.chaser_relative_position-sim_data_k.target_d_vector_I,
      sim_data_k.target_o_hat_vector_I_LVF,
      sim_data_k.target_omega,
      sim_data_k.target_d_dot_I
    );

    // Check if we need to increment/reset the "end simulation" counter:
    if (lvf_.InEndRange(sim_data_k.chaser_relative_position, sim_data_k.target_d_vector_I)){
      data_out.steps_in_end_region = sim_data_k.steps_in_end_region + 1;
    } else {
      data_out.steps_in_end_region = 0;
    }

  }

  // Integrate to get the new target omega:
  data_out.target_omega = target_omega_integrator->Integrate(sim_data_k.target_omega_dot_OI);

  // Update the target rotation matrix:
  data_out.target_q_BI = target_q_BI_integrator->Integrate(sim_data_k.target_q_BI_dot);

  // Check if the simulation should end:
  auto should_end = ((data_out.steps_in_end_region >= max_steps_in_end_region_) || (data_out.time >= max_time_));

  // Update chaser relative position and velocity in meters
  data_out.chaser_relative_velocity = data_out.desired_speed;
  data_out.chaser_relative_position = chaser_relative_position_integrator->Integrate(data_out.desired_speed);

  // Switch case 0: simulation is already over:
  if (should_end){
    data_out.simulation_complete = true;
  }

  // Increment the time:
  data_out.time = sim_data_k.time + dt_;

  // this will compute all of the control vectors, accelerations, relative positions, etc.
  data_out.ComputeDependentKinematicData(
    clvf_.OHatB(),
    clvf_.Alpha(),
    chaser_spacecraft_,
    target_spacecraft_,
    sim_data_k.desired_speed
  );

  // Finally, end the step:
  return data_out;
}

// Run a full simulation:
SimulationResult Simulation::Run(
  const SimulationData& sim_data_initial,
  const std::string& data_file_name
) const {
  // Create the sim data:
  SimulationData sim_data = sim_data_initial;

  // Make sure the dependent data is already calculated:
  sim_data.ComputeDependentData(
    clvf_.OHatB(), 
    clvf_.Alpha(), 
    chaser_spacecraft_, 
    target_spacecraft_,
    sim_data.chaser_relative_velocity);

  // Create all of the relavent integrators:
  auto chaser_orbital_velocity_integrator = std::make_shared<EulerIntegrator<3,1>>(sim_data.chaser_orbital_velocity, dt_);
  auto chaser_orbital_position_integrator = std::make_shared<EulerIntegrator<3,1>>(sim_data.chaser_orbital_position, dt_);
  auto target_orbital_velocity_integrator = std::make_shared<EulerIntegrator<3,1>>(sim_data.target_orbital_velocity, dt_);
  auto target_orbital_position_integrator = std::make_shared<EulerIntegrator<3,1>>(sim_data.target_orbital_position, dt_);
  auto target_omega_integrator = std::make_shared<EulerIntegrator<3,1>>(sim_data.target_omega, dt_);
  auto target_q_BI_integrator = std::make_shared<QuaternionEulerIntegrator>(sim_data.target_q_BI, dt_);

  // Create a data file to log the simulation,
  // only if a file name is given.
  std::ofstream data_stream;
  if (data_file_name != ""){
    data_stream.open(data_file_name);
    LogHeaders(data_stream);
  }

  while(!sim_data.simulation_complete){
    if (data_file_name != ""){
      LogData(data_stream, sim_data);
    }
    sim_data = Step(
      sim_data,
      chaser_orbital_velocity_integrator,
      chaser_orbital_position_integrator,
      target_orbital_velocity_integrator,
      target_orbital_position_integrator,
      target_omega_integrator,
      target_q_BI_integrator);
  }

  if (data_file_name != ""){
    LogData(data_stream, sim_data);
    data_stream.close();
  }

  SimulationResult result{};
  result.total_delta_v = sim_data.accumulated_delta_v;
  result.total_time = sim_data.time;

  return result;
}

// Run a full simulation:
SimulationResult Simulation::KinematicRun(
  const SimulationData& sim_data_initial,
  const std::string& data_file_name
) const {
  // Create the sime data:
  SimulationData sim_data = sim_data_initial;

  // Make sure the dependent data is already calculated:
  sim_data.ComputeDependentKinematicData(
    clvf_.OHatB(), 
    clvf_.Alpha(), 
    chaser_spacecraft_, 
    target_spacecraft_,
    Eigen::Vector3d::Zero());

  // Create the integrators:
  auto chaser_relative_position_integrator = std::make_shared<EulerIntegrator<3,1>>(sim_data.chaser_relative_position, dt_);
  auto target_omega_integrator = std::make_shared<EulerIntegrator<3,1>>(sim_data.target_omega, dt_);
  auto target_q_BI_integrator = std::make_shared<QuaternionEulerIntegrator>(sim_data.target_q_BI, dt_);

  // Create a data file to log the simulation:
  std::ofstream data_stream;
  if (data_file_name != ""){
    data_stream.open(data_file_name);
    LogHeaders(data_stream);
  }

  while(!sim_data.simulation_complete){
    if (data_file_name != ""){
      LogData(data_stream, sim_data);
    }
    sim_data = KinematicStep(
      sim_data,
      chaser_relative_position_integrator,
      target_omega_integrator,
      target_q_BI_integrator);
  }
  if (data_file_name != ""){
    LogData(data_stream, sim_data);
    data_stream.close();
  }

  SimulationResult result{};
  result.total_delta_v = sim_data.accumulated_delta_v;
  result.total_time = sim_data.time;

  return result;
}

void Simulation::LogHeaders(std::ofstream& data_stream){
  // TIME:
  data_stream << "time,";

  // Target oribtal position:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "target_orbital_position_" << i <<",";
  }

  // Target oribtal velocity:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "target_orbital_velocity_" << i <<",";
  }

  // Target orbital acceleration:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "target_orbital_acceleration_" << i <<",";
  }

  // target omega:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "target_omega_OI_" << i <<",";
  }  
  
  // target_omega_dot:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "target_omega_OI_dot_" << i <<",";
  }  

  // target q_BI:
  std::vector<std::string> names = {"x", "y", "z", "w"};
  for (int i = 0 ; i < 4 ; ++i){
    data_stream << "target_q_BI_" << names[i] << ",";
  }

  // target q_BI_dot:
  for (int i = 0 ; i < 4 ; ++i){
    data_stream << "target_q_BI_dot_" << names[i] << ",";
  }

  // Chaser oribtal position:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "chaser_orbital_position_" << i <<",";
  }

  // Chaser oribtal velocity:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "chaser_orbital_velocity_" << i <<",";
  }

  // Chaser orbital acceleration:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "chaser_orbital_acceleration_" << i <<",";
  }

  // Chaser relative position:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "chaser_relative_position_" << i <<",";
  }

  // Chaser relative position body fixed::
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "chaser_relative_position_B_" << i <<",";
  }

  // Chaser relative velocity:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "chaser_relative_velocity_" << i <<",";
  }

  // target_d_vector_I:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "target_d_vector_I_" << i <<",";
  }

  // target_d_dot_vector_I:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "target_d_dot_vector_I_" << i <<",";
  }

  // target_d_ddot_vector_I:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "target_d_ddot_vector_I_" << i <<",";
  }

  // target_o_hat_vector_I_CLVF:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "target_o_hat_vector_I_CLVF_" << i <<",";
  }
  
  // target_o_hat_vector_I_LVF:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "target_o_hat_vector_I_LVF_" << i <<",";
  }

  // desired_speed
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "desired_speed_" << i <<",";
  }

  // desired_acceleration:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "desired_acceleration_" << i <<",";
  }

  // control_vector:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "control_vector_" << i <<",";
  }

  // target_d_vector_B:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "target_d_vector_B_" << i <<",";
  }

  // target_o_hat_vector_B_LVF:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "target_o_hat_vector_B_LVF_" << i <<",";
  }

  // target_o_hat_vector_B_CLVF:
  for (int i = 0 ; i < 3 ; ++i){
    data_stream << "target_o_hat_vector_B_CLVF_" << i <<",";
  }

  data_stream << "alpha_CLVF,";

  data_stream << "r_CLVF,";
  data_stream << "theta_CLVF,";
  data_stream << "r_LVF,";
  data_stream << "theta_CLVF,";

  data_stream << "simulation_complete,";
  data_stream << "steps_in_switch_region,";
  data_stream << "steps_in_end_region,";
  data_stream << "in_CLVF,";
  data_stream << "delta_v,";
  data_stream << "accumulated_delta_v\n";
}

void Simulation::LogData(std::ofstream& data_stream, const SimulationData& sim_data){
  // TIME:
  data_stream << sim_data.time << ",";

  // Target oribtal position:
  for (const auto& num : sim_data.target_orbital_position){
    data_stream << num <<",";
  }

  // Target oribtal velocity:
  for (const auto& num : sim_data.target_orbital_velocity){
    data_stream << num <<",";
  }

  // Target orbital acceleration:
  for (const auto& num : sim_data.target_orbital_acceleration){
    data_stream << num <<",";
  }

  // target omega:
  for (const auto& num : sim_data.target_omega){
    data_stream << num <<",";
  }  
  
  // target_omega_dot:
  for (const auto& num : sim_data.target_omega_dot_OI){
    data_stream << num <<",";
  }

  // target C_BI:
  data_stream << sim_data.target_q_BI.x() << ",";
  data_stream << sim_data.target_q_BI.y() << ",";
  data_stream << sim_data.target_q_BI.z() << ",";
  data_stream << sim_data.target_q_BI.w() << ",";

  // target C_BI_dot:
  data_stream << sim_data.target_q_BI.x() << ",";
  data_stream << sim_data.target_q_BI.y() << ",";
  data_stream << sim_data.target_q_BI.z() << ",";
  data_stream << sim_data.target_q_BI.w() << ",";

  // Chaser oribtal position:
  for (const auto& num : sim_data.chaser_orbital_position){
    data_stream << num <<",";
  }

  // Chaser oribtal velocity:
  for (const auto& num : sim_data.chaser_orbital_velocity){
    data_stream << num <<",";
  }

  // Chaser orbital acceleration:
  for (const auto& num : sim_data.chaser_orbital_acceleration){
    data_stream << num <<",";
  }

  // Chaser relative position:
  for (const auto& num : sim_data.chaser_relative_position){
    data_stream << num <<",";
  }

  // Chaser relative position body fixed::
  for (const auto& num : sim_data.chaser_relative_position_B){
    data_stream << num <<",";
  }

  // Chaser relative velocity:
  for (const auto& num : sim_data.chaser_relative_velocity){
    data_stream << num <<",";
  }

  // target_d_vector_I:
  for (const auto& num : sim_data.target_d_vector_I){
    data_stream << num <<",";
  }

  // target_d_dot_vector_I:
  for (const auto& num : sim_data.target_d_dot_I){
    data_stream << num <<",";
  }

  // target_d_ddot_vector_I:
  for (const auto& num : sim_data.target_d_ddot_I){
    data_stream << num <<",";
  }

  // target_o_hat_vector_I_CLVF:
  for (const auto& num : sim_data.target_o_hat_vector_I_CLVF){
    data_stream << num <<",";
  }
  
  // target_o_hat_vector_I_LVF:
  for (const auto& num : sim_data.target_o_hat_vector_I_LVF){
    data_stream << num <<",";
  }

  // desired_speed
  for (const auto& num : sim_data.desired_speed){
    data_stream << num <<",";
  }

  // desired_acceleration:
  for (const auto& num : sim_data.desired_acceleration){
    data_stream << num <<",";
  }

  // control_vector:
  for (const auto& num : sim_data.control_vector){
    data_stream << num <<",";
  }

  // target_d_vector_B:
  for (const auto& num : sim_data.target_d_vector_B){
    data_stream << num <<",";
  }

  // target_o_hat_vector_B_LVF:
  for (const auto& num : sim_data.target_o_hat_vector_B_LVF){
    data_stream << num <<",";
  }

  // target_o_hat_vector_B_CLVF:
  for (const auto& num : sim_data.target_o_hat_vector_B_CLVF){
    data_stream << num <<",";
  }

  data_stream << sim_data.alpha_CLVF << ",";

  data_stream << sim_data.r_CLVF << ",";
  data_stream << sim_data.theta_CLVF << ",";
  data_stream << sim_data.r_LVF << ",";
  data_stream << sim_data.theta_CLVF << ",";

  data_stream << sim_data.simulation_complete << ",";
  data_stream << sim_data.steps_in_switch_region << ",";
  data_stream << sim_data.steps_in_end_region << ",";
  data_stream << sim_data.in_CLVF << ",";
  data_stream << sim_data.delta_v << ",";
  data_stream << sim_data.accumulated_delta_v << "\n";
}

}