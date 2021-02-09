#include <dphi_cloth_triangle_dX.h>

//compute 3x3 deformation gradient 
void dphi_cloth_triangle_dX(Eigen::Matrix3d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {

    Eigen::Vector3d X0, X1, X2;
    X0 = V.row(element(0));
    X1 = V.row(element(1));
    X2 = V.row(element(2));

    Eigen::MatrixXd T;
    T.resize(3, 2);
    T.col(0) = X1 - X0;
    T.col(1) = X2 - X0;

    Eigen::Vector2d Ones = Eigen::Vector2d::Ones();
    dphi.block<1, 3>(0, 0) = -1.0 * Ones.transpose() * (T.transpose() * T).inverse() * T.transpose();
    dphi.block<2, 3>(1, 0) = (T.transpose() * T).inverse() * T.transpose();

}