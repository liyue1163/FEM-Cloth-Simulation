#include <V_membrane_corotational.h>
#include <igl/svd3x3.h>
#include "iostream"

//Allowed to use libigl SVD or Eigen SVD for this part
void V_membrane_corotational(double &energy, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    Eigen::Vector3d q0, q1, q2;
    q0 = q.segment(3 * element(0), 3);
    q1 = q.segment(3 * element(1), 3);
    q2 = q.segment(3 * element(2), 3);

    Eigen::Vector3d X0, X1, X2;
    X0 = V.row(element(0));
    X1 = V.row(element(1));
    X2 = V.row(element(2));

    // deformed gradient matrix
    Eigen::Matrix3d F;
    // deformed space unit normal
    Eigen::Vector3d n = (q1 - q0).cross(q2 - q0).normalized();

    // ** UN ** deformed space unit normal
    Eigen::Vector3d N = (X1 - X0).cross(X2 - X0).normalized();
    // world space position matrix
    Eigen::Matrix34d x_world;
    x_world.col(0) = q0;
    x_world.col(1) = q1;
    x_world.col(2) = q2;
    x_world.col(3) = n;
    // ** UN **deformed space position matrix
    Eigen::Matrix43d X_undeform;
    X_undeform.block<3, 3>(0, 0) = dX;
    X_undeform.block<1, 3>(3, 0) = N.transpose();
    X_undeform(3,0) = N[0];
    X_undeform(3,1) = N[1];
    X_undeform(3,2) = N[2];


    F = x_world * X_undeform;

    // compute svd
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector3d singularVals = svd.singularValues();
    double psi = mu * (pow(singularVals[0] - 1.0, 2) + pow(singularVals[1] - 1.0, 2) + pow(singularVals[2] - 1.0, 2)) +
                 lambda/2.0 * pow(singularVals.sum() - 3.0, 2);

    energy = area * psi;

}
