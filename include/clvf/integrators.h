#ifndef CLVF_INTEGRATORS_H_
#define CLVF_INTEGRATORS_H_

#include <Eigen/Dense>

namespace clvf {

// Integrtor BASE class:
template <int N, int M>
class IntegratorBase {
    public:
        virtual Eigen::Matrix<double,N,M> Integrate(const Eigen::Matrix<double,N,M>& derivative) = 0;
};

class QuaternionIntegratorBase {
    public:
        virtual Eigen::Quaterniond Integrate(const Eigen::Quaterniond& derivative) = 0;
};

// EULER INTEGRATORS:
template <int N, int M>
class EulerIntegrator : public IntegratorBase<N,M> {

    private:
            Eigen::Matrix<double, N, M> previous_value_;
            const double dt_;

    public:
        // Constructor:
        EulerIntegrator() = delete;
        EulerIntegrator(
            const Eigen::Matrix<double,N,M>& initial_value,
            const double dt
        ) : previous_value_{initial_value}, dt_{dt}{};

        // Integrator function:
        Eigen::Matrix<double,N,M> Integrate(
            const Eigen::Matrix<double, N, M>& derivative
        ) override {
            previous_value_ = previous_value_ + derivative*dt_;
            return previous_value_;
        }
};

// EULER INTEGRATOR FOR QUATERNIONS:
class QuaternionEulerIntegrator : public QuaternionIntegratorBase {
    private:
        Eigen::Quaterniond previous_quaternion_;
        const double dt_;

    public:

        // Constructors:
        QuaternionEulerIntegrator()=delete;
        QuaternionEulerIntegrator(
            const Eigen::Quaterniond& initial_quaternion,
            const double dt
        ) : previous_quaternion_{initial_quaternion.normalized()}, dt_{dt}{};

        // Integration function:
        Eigen::Quaterniond Integrate(
            const Eigen::Quaterniond& derivative){
            Eigen::Quaterniond result;
            result.x() = previous_quaternion_.x() + derivative.x()*dt_;
            result.y() = previous_quaternion_.y() + derivative.y()*dt_;
            result.z() = previous_quaternion_.z() + derivative.z()*dt_;
            result.w() = previous_quaternion_.w() + derivative.w()*dt_;
            result.normalize();

            previous_quaternion_ = result;
            return previous_quaternion_;
}
};

}

#endif