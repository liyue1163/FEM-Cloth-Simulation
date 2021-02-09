#include <assemble_forces.h>
#include <iostream>
#include <dphi_cloth_triangle_dX.h>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0,
                     double mu, double lambda) {
    f = Eigen::VectorXd::Zero(q.rows());

    // iterate through each triangle
    for (int i=0; i<F.rows(); i++) {

        Eigen::Vector9d dV;
        Eigen::RowVectorXi element = F.row(i);

        // convert 1x9 vector to 3x3 matrix
        Eigen::Matrix<double, 1, 9> tmp_row;
        tmp_row = dX.row(i); //ei is the triangle index.
        Eigen::Matrix3d dphi_dX = Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::ColMajor>>(tmp_row.data());

        // triangle area
        double area = a0(i);
        dV_membrane_corotational_dq(dV, q, dphi_dX, V, element, area, mu, lambda);


        int tri_idx1 = element(0);
        int tri_idx2 = element(1);
        int tri_idx3 = element(2);

        // directly update relavent elements from f,
        f[3 * tri_idx1]     -= dV[0];
        f[3 * tri_idx1 + 1] -= dV[1];
        f[3 * tri_idx1 + 2] -= dV[2];

        f[3 * tri_idx2]     -= dV[3];
        f[3 * tri_idx2 + 1] -= dV[4];
        f[3 * tri_idx2 + 2] -= dV[5];

        f[3 * tri_idx3]     -= dV[6];
        f[3 * tri_idx3 + 1] -= dV[7];
        f[3 * tri_idx3 + 2] -= dV[8];
    }
};
