#ifndef CLVF_SIMULATION_H_
#define CLVF_SIMULATION_H_

#include "clvf/spacecraft.h"
#include "clvf/clvf.h"
#include <fstream>
#include "clvf/utils.h"

namespace clvf {

struct SimulationData {
  // Time
  double time;

  // Dynamics data:
  Eigen::Vector3d target_orbital_position;
  Eigen::Vector3d target_orbital_velocity;
  Eigen::Vector3d target_orbital_acceleration; // dependent
  Eigen::Vector3d target_omega;
  Eigen::Vector3d target_omega_dot_OI; // dependent
  Eigen::Quaterniond target_q_BI;
  Eigen::Quaterniond target_q_BI_dot; // dependent
  
  Eigen::Vector3d chaser_orbital_position;
  Eigen::Vector3d chaser_orbital_velocity;
  Eigen::Vector3d chaser_orbital_acceleration; // dependent
  Eigen::Vector3d chaser_relative_position; //dependent
  Eigen::Vector3d chaser_relative_velocity; //dependent
  Eigen::Vector3d chaser_relative_position_B; //dependent

  // Some geometry:
  Eigen::Vector3d target_d_vector_I; // dependent
  Eigen::Vector3d target_d_dot_I; // dependent
  Eigen::Vector3d target_d_ddot_I; // dependent
  Eigen::Vector3d target_o_hat_vector_I_CLVF; // dependent
  Eigen::Vector3d target_o_hat_vector_I_LVF; // dependent

  // Guidance data:
  Eigen::Vector3d desired_speed;
  Eigen::Vector3d desired_acceleration;

  // Control data:
  Eigen::Vector3d control_vector; // dependent

  // Body-fixed things:
  Eigen::Vector3d target_d_vector_B;          // dependent
  Eigen::Vector3d target_o_hat_vector_B_CLVF; // dependent
  Eigen::Vector3d target_o_hat_vector_B_LVF;  // dependent
  double alpha_CLVF;                          // dependent

  // CLVF core values:
  double r_CLVF;        // dependent
  double theta_CLVF;    // dependent

  // LVF core values:
  double r_LVF;         // dependent
  double theta_LVF;     // dependent

  // Simulation complete:
  bool simulation_complete;

  // Number of steps completed in the "switch" and "end" regions:
  int steps_in_switch_region;
  int steps_in_end_region;

  // Are we in the CLVF?
  bool in_CLVF;

  // A helper function to compute all of the "dependent data" given the true initial values:
  void ComputeDependentData(
    const Eigen::Vector3d& o_hat_B_CLVF,
    double alpha_CLVF_in,
    const Spacecraft& chaser_spacecraft,
    const Spacecraft& target_spacecraft
  ) {

    // grab the body fixed values again:
    target_d_vector_B = target_spacecraft.DVectorB();
    target_o_hat_vector_B_LVF = target_spacecraft.OHatB();
    target_o_hat_vector_B_CLVF = o_hat_B_CLVF;

    alpha_CLVF = alpha_CLVF_in;

    // Extract the C_BI from previous sim data:
    const auto& q_BI = target_q_BI;
    Eigen::Quaterniond q_IB = q_BI.conjugate();

    // Angular dynamics:
    target_omega_dot_OI = 
      RotateVectorByQuaternion( target_spacecraft.OmegaDot(
        RotateVectorByQuaternion(target_omega, q_BI), 
        Eigen::Vector3d::Zero()
      ), q_IB);


    target_q_BI_dot = QuaternionDerivative(target_omega, target_q_BI);
    // target_C_BI_dot = target_spacecraft.RotationMatrixDot(C_BI, target_omega);

    // update compute the docking port vector:
    target_d_vector_I =  RotateVectorByQuaternion(target_spacecraft.DVectorB(), q_IB);

    // The velocity of the docking port as viewed in the inertial frame:
    target_d_dot_I = 
      target_omega.cross(target_d_vector_I);

    // The acceleration of the docking port as viewed in the inertial frame:
    target_d_ddot_I = 
      target_omega_dot_OI.cross(target_d_vector_I) + 
        target_omega.cross(target_omega.cross(target_d_vector_I));

    // The O_hat_vector of the CVLF:
    target_o_hat_vector_I_CLVF = RotateVectorByQuaternion(o_hat_B_CLVF, q_IB);

    // And for the LVF:
    target_o_hat_vector_I_LVF = RotateVectorByQuaternion(target_spacecraft.OHatB(), q_IB);

    // Chaser relative position and velocity in meters:
    chaser_relative_position = (chaser_orbital_position - target_orbital_position) * clvf::kKmToMeters;
    chaser_relative_velocity = (chaser_orbital_velocity - target_orbital_velocity) * clvf::kKmToMeters;

    chaser_relative_position_B = RotateVectorByQuaternion(chaser_relative_position, q_BI);

    // Update the control vector:
    control_vector = chaser_spacecraft.Control(
      desired_speed, 
      chaser_relative_velocity,
      desired_acceleration
    );

    // Compute the chaser orbital acceleration:
    chaser_orbital_acceleration = chaser_spacecraft.RInertialDDot(
      chaser_orbital_position,
      control_vector * clvf::kMetersToKm,
      clvf::kMu
    );

    // compute the target orbital acceleration:
    target_orbital_acceleration = target_spacecraft.RInertialDDot(
      target_orbital_position, 
      Eigen::Vector3d::Zero(), 
      clvf::kMu
    );

    // Compute the r and theta values for the CLVF:
    r_CLVF = chaser_relative_position.norm();
    theta_CLVF = ThetaFromTwoVectors(chaser_relative_position_B, o_hat_B_CLVF);

    // Compute for the LVF:
    r_LVF = (chaser_relative_position - target_d_vector_I).norm();
    theta_CLVF = ThetaFromTwoVectors((chaser_relative_position - target_d_vector_I), target_o_hat_vector_I_LVF);
  }

  // A helper function to compute all of the "dependent data" given the true initial values:
  void ComputeDependentKinematicData(
    const Eigen::Vector3d& o_hat_B_CLVF,
    double alpha_CLVF_in,
    const Spacecraft& chaser_spacecraft,
    const Spacecraft& target_spacecraft
  ) {

    // grab the body fixed values again:
    target_d_vector_B = target_spacecraft.DVectorB();
    target_o_hat_vector_B_LVF = target_spacecraft.OHatB();
    target_o_hat_vector_B_CLVF = o_hat_B_CLVF;

    alpha_CLVF = alpha_CLVF_in;

    // Extract the C_BI from previous sim data:
    const auto& q_BI = target_q_BI;
    Eigen::Quaterniond q_IB = q_BI.conjugate();

    // Angular dynamics:
    target_omega_dot_OI = 
      RotateVectorByQuaternion( target_spacecraft.OmegaDot(
        RotateVectorByQuaternion(target_omega, q_BI), 
        Eigen::Vector3d::Zero()
      ), q_IB);


    target_q_BI_dot = QuaternionDerivative(target_omega, target_q_BI);

    // update compute the docking port vector:
    target_d_vector_I =  RotateVectorByQuaternion(target_spacecraft.DVectorB(), q_IB);

    // The velocity of the docking port as viewed in the inertial frame:
    target_d_dot_I = 
      target_omega.cross(target_d_vector_I);

    // The acceleration of the docking port as viewed in the inertial frame:
    target_d_ddot_I = 
      target_omega_dot_OI.cross(target_d_vector_I) + 
        target_omega.cross(target_omega.cross(target_d_vector_I));

    // The O_hat_vector of the CVLF:
    target_o_hat_vector_I_CLVF = RotateVectorByQuaternion(o_hat_B_CLVF, q_IB);

    // And for the LVF:
    target_o_hat_vector_I_LVF = RotateVectorByQuaternion(target_spacecraft.OHatB(), q_IB);
    chaser_relative_position_B = RotateVectorByQuaternion(chaser_relative_position, q_BI);

    // Compute the r and theta values for the CLVF:
    r_CLVF = chaser_relative_position.norm();
    theta_CLVF = ThetaFromTwoVectors(chaser_relative_position_B, o_hat_B_CLVF);

    // Compute for the LVF:
    r_LVF = (chaser_relative_position - target_d_vector_I).norm();
    theta_CLVF = ThetaFromTwoVectors((chaser_relative_position - target_d_vector_I), target_o_hat_vector_I_LVF);
  }

};

class Simulation {
  private:
    // elements of the simulation:
    CLVF clvf_;
    LVF lvf_;
    Spacecraft target_spacecraft_;
    Spacecraft chaser_spacecraft_;
    const double dt_;

    // a maximum number of steps which should be spent inside of the switch region:
    const int max_steps_in_switch_region_;

    // a maximum number of steps which is spent inside of the "end" region before completing.
    const int max_steps_in_end_region_;

    // a maximum amount of time, before the simulation is forced to end:
    const double max_time_;

    // For Logging the datafile headers and data:
    static void LogHeaders(std::ofstream& data_stream);
    static void LogData(std::ofstream& data_stream, const SimulationData& sim_data);

  public:
    Simulation(
      const CLVF& clvf,
      const LVF& lvf,
      const Spacecraft& target_spacecraft,
      const Spacecraft& chaser_spacecraft,
      double dt,
      int max_steps_in_switch_region,
      int max_steps_in_end_region,
      double max_time
    )
    :
    clvf_{clvf},
    lvf_{lvf},
    target_spacecraft_{target_spacecraft},
    chaser_spacecraft_{chaser_spacecraft},
    dt_{dt},
    max_steps_in_switch_region_{max_steps_in_switch_region},
    max_steps_in_end_region_{max_steps_in_end_region},
    max_time_{max_time}
    {};

    // Single step in the simulation:
    SimulationData Step(const SimulationData& sim_data_k) const;

    // Running the simulation:
    void Run(
      const SimulationData& sim_data_initial,
      const std::string& data_file_name
    ) const;

    // Pure kinematic step:
    SimulationData KinematicStep(const SimulationData& sim_data_k) const;

    // Run a pure kinematic simulation:
    void KinematicRun(
      const SimulationData& sim_data_initial,
      const std::string& data_file_name
    ) const;};
}

#endif