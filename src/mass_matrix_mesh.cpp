#include <mass_matrix_mesh.h>


// helper function
void triangle_mass_matrix(Eigen::Matrix99d &M, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double area) {
    Eigen::Matrix99d M_partial;
    M_partial.setZero();
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            if (i == j) {
                M_partial.block<3, 3>(i*3, j*3) = 1.0/12.0 * Eigen::Matrix3d::Identity();
            }
            else {
                M_partial.block<3, 3>(i*3, j*3) = 1.0/24.0 * Eigen::Matrix3d::Identity();
            }
        }
    }

    M = 2.0 * density * area * M_partial;
}


void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, 
                         Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F,
                         double density, Eigen::Ref<const Eigen::VectorXd> areas) {

    M.resize(q.rows(), q.rows());
    std::vector<Eigen::Triplet<double>> tripletList;

    for (int tri=0; tri<F.rows(); tri++) {
        Eigen::Matrix99d M_j;
        Eigen::RowVectorXi element_j = F.row(tri);
        double area_tri = areas(tri);
        triangle_mass_matrix(M_j, element_j, density, area_tri);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                int row_idx = element_j(i);
                int col_idx = element_j(j);

                for (int k=0; k<3; k++) {
                    tripletList.push_back({row_idx * 3 + k, col_idx * 3 + k, M_j(i * 3 + k, j * 3 + k)});
                }
            }
        }
    }
    M.setFromTriplets(tripletList.begin(), tripletList.end());

}
 
