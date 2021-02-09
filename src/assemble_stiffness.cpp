#include <assemble_stiffness.h>
#include <dphi_cloth_triangle_dX.h>
#include <iostream>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0,
                     double mu, double lambda) {

    K.resize(q.rows(), q.rows());
    K.setZero();
    std::vector<Eigen::Triplet<double>> tripletList;

    // iterate through each triangle
    for (int tri=0; tri<F.rows(); tri++) {

        Eigen::Matrix99d K_tri;
        Eigen::RowVectorXi element = F.row(tri);

        // convert 1x9 vector to 3x3 matrix
        Eigen::Matrix<double, 1, 9> tmp_row;
        tmp_row = dX.row(tri); //ei is the triangle index.
        Eigen::Matrix3d dphi_dX = Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::ColMajor>>(tmp_row.data());
        double area_i = a0(tri);
        d2V_membrane_corotational_dq2(K_tri, q, dphi_dX, V, element, area_i, mu, lambda);

        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {

                Eigen::Matrix3d K_ij = -1.0 * K_tri.block<3,3>(i * 3,j * 3);
                int row_idx = element(i);
                int col_idx = element(j);

                for (int row_offset=0; row_offset<3; row_offset++) {
                    for (int col_offset=0; col_offset<3; col_offset++) {
                        tripletList.push_back({row_idx * 3 + row_offset, col_idx * 3 + col_offset, K_ij(row_offset, col_offset)});
                    }
                }
            }
        }
    }
    K.setFromTriplets(tripletList.begin(), tripletList.end());

};
