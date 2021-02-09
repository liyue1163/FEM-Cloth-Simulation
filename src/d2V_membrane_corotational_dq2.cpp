#include <d2V_membrane_corotational_dq2.h>
#include <iostream>

// helper function
void get_dn_dq_j( Eigen::MatrixXd &dn_dq_j, Eigen::Vector3d n, Eigen::Vector3d q0, Eigen::Vector3d q1, Eigen::Vector3d q2) {

    Eigen::Vector3d dx1 = q1 - q0;
    Eigen::Vector3d dx2 = q2 - q0;

    Eigen::Matrix3d skew_x1;
    skew_x1 << 0, -1.0*dx1[2], dx1[1],
            dx1[2], 0, -1.0*dx1[0],
            -1.0*dx1[1], dx1[0], 0;

    Eigen::Matrix3d skew_x2;
    skew_x2 << 0, -1.0*dx2[2], dx2[1],
            dx2[2], 0, -1.0*dx2[0],
            -1.0*dx2[1], dx2[0], 0;

    Eigen::MatrixXd concat_Identity1;
    concat_Identity1.resize(3, 9);
    concat_Identity1 << -1.0 * Eigen::MatrixXd::Identity(3, 3),
            Eigen::MatrixXd::Zero(3, 3),
            Eigen::MatrixXd::Identity(3, 3);

    Eigen::MatrixXd concat_Identity2;
    concat_Identity2.resize(3, 9);
    concat_Identity2 << -1.0 * Eigen::MatrixXd::Identity(3, 3),
            Eigen::MatrixXd::Identity(3, 3),
            Eigen::MatrixXd::Zero(3, 3);

    dn_dq_j = (1.0 / (dx1.cross(dx2)).norm()) * (Eigen::MatrixXd::Identity(3, 3) - n*n.transpose()) * (skew_x1 * concat_Identity1 - skew_x2 * concat_Identity2);

}


void d2V_membrane_corotational_dq2(Eigen::Matrix99d &H, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX,
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area,
                          double mu, double lambda) {


    //SVD = USW^T
    Eigen::Matrix3d U;
    Eigen::Vector3d S;
    Eigen::Matrix3d W;
    Eigen::Matrix3d F; //deformation gradient

    double tol = 1e-5;

    //Compute SVD of F here

    Eigen::Vector3d q0, q1, q2;
    q0 = q.segment(3 * element(0), 3);
    q1 = q.segment(3 * element(1), 3);
    q2 = q.segment(3 * element(2), 3);
    Eigen::Vector9d q_j;
    q_j << q0, q1, q2;

    Eigen::Vector3d X0, X1, X2;
    X0 = V.row(element(0));
    X1 = V.row(element(1));
    X2 = V.row(element(2));

    // deformed space unit normal vec
    Eigen::Vector3d n = (q1 - q0).cross(q2 - q0).normalized();
    // ** UN **deformed space unit normal vec
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

    F = x_world * X_undeform;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);

    U = svd.matrixU();
    S = svd.singularValues();
    W = svd.matrixV();

    //deal with singularity in the svd gradient
    if(std::fabs(S[0] - S[1]) < tol || std::fabs(S[1] - S[2]) < tol || std::fabs(S[0] - S[2]) < tol) {

        F += Eigen::Matrix3d::Random()*tol;
        Eigen::JacobiSVD<Eigen::Matrix3d> svd2(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        U = svd2.matrixU();
        W = svd2.matrixV();
        S = svd2.singularValues();
    }

    //Fix for inverted elements (thanks to Danny Kaufman)
    double det = S[0]*S[1];

     if(det <= -1e-10)
    {
        if(S[0] < 0) S[0] *= -1;
        if(S[1] < 0) S[1] *= -1;
        if(S[2] < 0) S[2] *= -1;
    }

    if(U.determinant() <= 0)
    {
        U(0, 2) *= -1;
        U(1, 2) *= -1;
        U(2, 2) *= -1;
    }

    if(W.determinant() <= 0)
    {
        W(0, 2) *= -1;
        W(1, 2) *= -1;
        W(2, 2) *= -1;
    }

    //TODO: compute H, the hessian of the corotational energy
    Eigen::Tensor3333d dU;
    Eigen::Tensor333d dS;
    Eigen::Tensor3333d dV;
    // compute dU, dS, dV
    dsvd(dU, dS, dV, F);

   Eigen::DiagonalMatrix<double, 3> ds;
   Eigen::Vector3d ds_diag;
   for (int i=0; i < ds_diag.size(); i++) {
       ds_diag[i] = 2.0 * mu * (S[i] - 1.0) + lambda * (S.sum() - 3.0) * 1.0;
   }
   ds = ds_diag.asDiagonal();


    // compute ds_2
    Eigen::Matrix3d ds_2 = 1.0 * lambda * Eigen::Matrix3d::Ones(3, 3);

    for (int k = 0; k < ds_2.rows(); k++) {
        ds_2(k, k) += 2.0 * mu;
    }

    Eigen::Matrix99d d2psi_dF2;
    d2psi_dF2.setZero();
    Eigen::Matrix3d d2psi_dF2_ij;
    for (int i=0; i < d2psi_dF2_ij.rows(); i++) {
        for (int j=0; j < d2psi_dF2_ij.cols(); j++) {


            Eigen::Matrix3d dU_dF_ij = dU[i][j];
            Eigen::Matrix3d dV_dF_ij = dV[i][j];

            Eigen::Vector3d diag_vec = ds_2 * dS[i][j];
            Eigen::Matrix3d diag = Eigen::Matrix3d::Zero();
            diag(0,0) = diag_vec[0];
            diag(1,1) = diag_vec[1];
            diag(2,2) = diag_vec[2];

            d2psi_dF2_ij = dU_dF_ij * ds * W.transpose() + U * diag * W.transpose() + U * ds * dV_dF_ij.transpose();
            // row-major flatten
            Eigen::Vector9d vec_d2psi_dF2_ij;
            for (int k=0; k < d2psi_dF2_ij.cols(); k++) {
                vec_d2psi_dF2_ij.segment(3*k, 3) = d2psi_dF2_ij.row(k);
            }
            d2psi_dF2.row(3 * i + j) = vec_d2psi_dF2_ij;
        }
    }


    // compute dF_dq_j
    Eigen::Matrix99d dF_dq_j;
    // require compute dn_dq_j
    Eigen::MatrixXd dn_dq_j;
    dn_dq_j.resize(3, 9);
    get_dn_dq_j(dn_dq_j, n, q0, q1, q2);
    Eigen::Matrix99d B_j;
    B_j <<  dX(0, 0), 0, 0, dX(1, 0), 0, 0, dX(2, 0), 0, 0,
            dX(0, 1), 0, 0, dX(1, 1), 0, 0, dX(2, 1), 0, 0,
            dX(0, 2), 0, 0, dX(1, 2), 0, 0, dX(2, 2), 0, 0,
            0, dX(0, 0), 0, 0, dX(1, 0), 0, 0, dX(2, 0), 0,
            0, dX(0, 1), 0, 0, dX(1, 1), 0, 0, dX(2, 1), 0,
            0, dX(0, 2), 0, 0, dX(1, 2), 0, 0, dX(2, 2), 0,
            0, 0, dX(0, 0), 0, 0, dX(1, 0), 0, 0, dX(2, 0),
            0, 0, dX(0, 1), 0, 0, dX(1, 1), 0, 0, dX(2, 1),
            0, 0, dX(0, 2), 0, 0, dX(1, 2), 0, 0, dX(2, 2);

    Eigen::MatrixXd N_j = Eigen::MatrixXd::Zero(9, 3);
    N_j.block<3,1>(0, 0) = N;
    N_j.block<3,1>(3, 1) = N;
    N_j.block<3,1>(6, 2) = N;

    dF_dq_j = B_j + N_j * dn_dq_j;

    H = area * dF_dq_j.transpose() * d2psi_dF2 * dF_dq_j;


    //fix errant eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix99d> es(H);

    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();

    for (int i = 0; i < 9; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }

    H = Evec * DiagEval * Evec.transpose();

}
