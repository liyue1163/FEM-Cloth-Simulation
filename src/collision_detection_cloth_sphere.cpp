#include <collision_detection_cloth_sphere.h>
#include <iostream>
void collision_detection_cloth_sphere(std::vector<unsigned int> &cloth_index, std::vector<Eigen::Vector3d> &normals, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Vector3d> center, double radius) {

    cloth_index.clear();
    normals.clear();

    int V_num = q.size() / 3;
    for (int v = 0; v < V_num; v++) {
        Eigen::Vector3d vert = q.segment(3 * v, 3);
        if ( !((vert - center).squaredNorm() > pow(radius, 2)) ) {
            Eigen::Vector3d n = (vert - center).normalized();

            cloth_index.push_back( (unsigned int) v);
            normals.push_back(n);
        }
    }
}
