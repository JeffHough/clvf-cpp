#include "clvf/clvf.h"
#include "clvf/spacecraft.h"
#include "clvf/utils.h"
#include "clvf/simulation.h"
#include <iostream>

int main(){
  // Set up a control parameter:
  double beta = 5.0; // [rads/s]

  // Set up a target spacecraft:
  Eigen::Matrix3d target_inertia = Eigen::Vector3d({3.0, 5.0, 7.0}).asDiagonal();
  double target_mass = 1.0;
  Eigen::Vector3d o_hat_B = {1.0, 0., 0.};
  Eigen::Vector3d d_vector_B = {1., 0., 0.};
  double angle_of_acceptance = 30 * clvf::kD2R;
  clvf::Spacecraft target(target_inertia, target_mass, beta, o_hat_B, d_vector_B, angle_of_acceptance);

  // For now, just make the chaser identical. They do different things anyways 
  // (chaser is pure translation, target is mostly rotation).
  clvf::Spacecraft chaser(target_inertia, target_mass, beta, o_hat_B, d_vector_B, angle_of_acceptance);

  // Set up the LVF parameters:
  double v_max = 1.0;
  double frac = 0.99;
  double alpha_prime = 5.0;
  double end_region_radius = 0.01;
  
  clvf::LVF lvf_guidance(
    v_max, 
    frac, 
    alpha_prime, 
    target.AngleOfAcceptance(), 
    end_region_radius
  );

  // Set up the CLVF parameters:
  double kc = 5.0;
  double ka = 1.0;
  double b = 15.0;
  double radius_error_before_changing_to_LVF = 0.01;
  double theta_error_before_changing_to_LVF = 0.01;

  clvf::CLVF clvf_guidance(
    kc, 
    ka, 
    b, 
    alpha_prime, 
    target.OHatB(), 
    target.DVectorB(), 
    radius_error_before_changing_to_LVF,
    theta_error_before_changing_to_LVF
  );

  // Instantiate the simulation:
  double dt = 0.01;
  int max_steps_in_switch_region = static_cast<int>(5.0 / dt);
  int max_steps_in_end_region = static_cast<int>(5.0 / dt);
  double max_time = 5000.0; // seconds of simulation time.

  clvf::Simulation sim(
    clvf_guidance, 
    lvf_guidance,
    target,
    chaser,
    dt,
    max_steps_in_switch_region,
    max_steps_in_end_region,
    max_time
  );

  // Run a simulation with particular initial conditions:
  clvf::SimulationData initial_data{};

  // Set all of the independent simulation data:
  double semi_major_axis = 8000; //km
  double eccentricity = 0.01;
  double inclination = 0;
  double argument_of_perigee = 0;
  double RAAN = 0.0;
  double true_anomaly = 0.0;

  auto pos_vel = target.OrbitalElementsToPosVel(
    semi_major_axis,
    eccentricity,
    inclination,
    RAAN,
    argument_of_perigee,
    true_anomaly
  );

  std::cout << "LVF target-vector\n" << target.OHatB() << 
              "\nCLVF o-hat-vector\n"<<clvf_guidance.OHatB() << 
              "\nCLVF target-vector"<< clvf_guidance.OHatB()*clvf_guidance.Alpha() << 
              "\nLVF D-vector_B\n" << target.DVectorB() << "\n";

  initial_data.target_orbital_position = pos_vel.first;
  initial_data.target_orbital_velocity = pos_vel.second;
  initial_data.target_omega = {-1.5*clvf::kD2R, 3.5*clvf::kD2R, 2.0*clvf::kD2R};
  initial_data.target_q_BI = Eigen::Quaterniond::Identity();  

  // TODO: should the relative position be the "non-dependent" all the time?
  // initial_data.chaser_orbital_position = initial_data.target_orbital_position + Eigen::Vector3d({0.05, -0.03, 0.01}); // km;
  // initial_data.chaser_orbital_velocity = initial_data.target_orbital_velocity;
  initial_data.chaser_relative_position = {50, -30, 100};
  initial_data.simulation_complete = false;
  initial_data.in_CLVF = true;

  // sim.Run(initial_data, "test.csv");
  // Test with JUST the kinematics (follow exactly the desired speed).
  auto sim_result = sim.KinematicRun(initial_data);//, "test.csv");

  std::cout << "\nThe total delta-v: " << sim_result.total_delta_v 
            << "\nThe total simulated time: " << sim_result.total_time
            << std::endl;

  return 0;
}