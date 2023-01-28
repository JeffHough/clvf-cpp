#ifndef CLVF_CLVF_H_
#define CLVF_CLVF_H_

#include <Eigen/Dense>

namespace clvf {

constexpr double kSmallNumber = 1e-7;

// Guidance law based on CLVFs:
class CLVF {
  private:

    // Parameters of the vector field:
    const double kc_;
    const double ka_;
    const double b_;
    const double alpha_;

    // Parameter of the controller:
    const double beta_;

    // The circulation component:
    double SFunction(
      double r,
      double theta
    ) const;

    // The contraction component:
    double VFunction(
      double r,
      double theta
    ) const;

    // The g-function to compensate for target rotation:
    double GFunction(
      double r
    ) const;

    // Computation for alignment direction:
    static Eigen::Vector3d AHatVector(
      const Eigen::Vector3d& r_vector,
      const Eigen::Vector3d& o_hat_vector
    );

    // Computation for the e-vector:
    static Eigen::Vector3d EHatVector(
      const Eigen::Vector3d& r_vector,
      const Eigen::Vector3d& o_hat_vector
    );

    // BELOW: several intermediate derivatives to compute the desired acceleration:
    Eigen::Vector3d RHatDot(
      const Eigen::Vector3d& r_vector,
      const Eigen::Vector3d& velocity
    ) const;

    Eigen::Vector3d AHatDot(
      const Eigen::Vector3d& a_hat,
      const Eigen::Vector3d& e_hat,
      const Eigen::Vector3d& o_hat,
      const Eigen::Vector3d& r_hat,
      const Eigen::Vector3d& r_hat_dot,
      const Eigen::Vector3d& omega_OI
    ) const;

    double RDot(
      const Eigen::Vector3d& r_hat,
      const Eigen::Vector3d& velocity
    ) const;

    double ThetaDot(
      const Eigen::Vector3d& o_hat,
      const Eigen::Vector3d& r_hat_dot,
      const Eigen::Vector3d& e_hat,
      const Eigen::Vector3d& omega_OI,
      double theta
    ) const;

    double GFunctionDot(double r, double r_dot) const;
    double SFunctionDot(double theta, double r, double theta_dot, double r_dot) const;
    double VFunctionDot(double r, double r_dot) const;

  public:
    // No default constructor:
    CLVF() = delete;
    CLVF(
      double kc,
      double ka,
      double b,
      double alpha,
      double beta
    ) : kc_{kc}, ka_{ka}, b_{b}, alpha_{alpha}, beta_{beta}{};

    Eigen::Vector3d DesiredVelocity(
      const Eigen::Vector3d& r_vector,
      const Eigen::Vector3d& o_hat_vector,
      const Eigen::Vector3d& omega_OI,
      const Eigen::Vector3d& d_vector_dot
    ) const;

    Eigen::Vector3d DesiredAcceleration(
      const Eigen::Vector3d& r_vector,
      const Eigen::Vector3d& velocity,
      const Eigen::Vector3d& o_hat_vector,
      const Eigen::Vector3d& omega_OI,
      const Eigen::Vector3d& omega_dot_OI,
      const Eigen::Vector3d& d_ddot
    ) const;

    double AccelerationBound(
      double omega_max, 
      double omega_and_omega_dot_max,
      double d_ddot_max
    ) const;
};

// Class for standard Lyapunov Vector Field:
class LVF {
  private:
    const double v_max_;
    const double final_angle_;
    const double alpha_prime_;
    const double theta_docking_cone_;

    // Parameter of the controller:
    const double beta_;

    // Various helper functions for the LVF:
    double ThetaN(
      double theta
    ) const;

    double RN(
      double r
    ) const;

    double VRel(
      double r_N
    ) const {
      return v_max_ * std::sin(r_N);
    }

    // The circulation component:
    double SFunction(
      double v_rel,
      double theta_N
    ) const {
      return v_rel * std::sin(theta_N);
    }

    // The contraction component:
    double VFunction(
      double v_rel,
      double theta_N
    ) const {
      return -v_rel * std::cos(theta_N);
    }

    // The g-function to compensate for target rotation:
    inline double GFunction(
      double r
    ) const {
      return r;
    }

    // Computation for alignment direction:
    static Eigen::Vector3d AHatVector(
      const Eigen::Vector3d& r_vector,
      const Eigen::Vector3d& o_hat_vector
    );

    // Computation for the e-vector:
    static Eigen::Vector3d EHatVector(
      const Eigen::Vector3d& r_vector,
      const Eigen::Vector3d& o_hat_vector
    ) {
      return (r_vector.normalized()).cross(AHatVector(r_vector, o_hat_vector));
    };

    // BELOW: several intermediate derivatives to compute the desired acceleration:
    Eigen::Vector3d RHatDot(
      const Eigen::Vector3d& r_vector,
      const Eigen::Vector3d& velocity
    ) const ;

    Eigen::Vector3d AHatDot(
      const Eigen::Vector3d& a_hat,
      const Eigen::Vector3d& e_hat,
      const Eigen::Vector3d& o_hat,
      const Eigen::Vector3d& r_hat,
      const Eigen::Vector3d& r_hat_dot,
      const Eigen::Vector3d& omega_OI
    ) const;

    double RDot(
      const Eigen::Vector3d& r_hat,
      const Eigen::Vector3d& velocity
    ) const {
      return r_hat.dot(velocity);
    }

    double ThetaDot(
      const Eigen::Vector3d& o_hat,
      const Eigen::Vector3d& r_hat_dot,
      const Eigen::Vector3d& e_hat,
      const Eigen::Vector3d& omega_OI,
      const Eigen::Vector3d& r_hat,
      double theta
    ) const;

    double GFunctionDot(double r_dot) const {return r_dot;}
    double SFunctionDot(double theta, double r, double theta_dot, double r_dot) const;
    
    double VFunctionDot(    
      double theta,
      double r,
      double theta_dot, 
      double r_dot
    ) const;

  public:
    LVF() = delete;
    LVF(
      double v_max,
      double frac,
      double alpha_prime,
      double theta_docking_cone,
      double beta
    ):
    v_max_{v_max}, 
    final_angle_{frac*M_PI}, 
    alpha_prime_{alpha_prime}, 
    theta_docking_cone_{theta_docking_cone},
    beta_{beta}{};

    Eigen::Vector3d DesiredVelocity(
      const Eigen::Vector3d& r_vector,
      const Eigen::Vector3d& o_hat_vector,
      const Eigen::Vector3d& omega_OI,
      const Eigen::Vector3d& d_dot
    ) const;

    Eigen::Vector3d DesiredAcceleration(
      const Eigen::Vector3d& r_vector,
      const Eigen::Vector3d& velocity,
      const Eigen::Vector3d& o_hat_vector,
      const Eigen::Vector3d& omega_OI,
      const Eigen::Vector3d& omega_dot_OI,
      const Eigen::Vector3d& d_ddot
    ) const;

    double AccelerationBound(
      double omega_max, 
      double omega_and_omega_dot_max) const;
};

}

#endif