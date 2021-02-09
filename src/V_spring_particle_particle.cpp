#include <V_spring_particle_particle.h>

void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3, 3);
    Eigen::MatrixXd neg_I = -1 * Eigen::MatrixXd::Identity(3, 3);
    Eigen::MatrixXd B(I.rows(), I.cols()+I.cols());
    B << neg_I, I;
    // cat q0 and q1 vertically
    Eigen::VectorXd q(6,1);
    q << q0, q1;
    double temp = (q.transpose() * B.transpose() * B * q).array().sqrt()(0);

    V = 0.5 * stiffness * pow((temp - l0), 2);
    
}