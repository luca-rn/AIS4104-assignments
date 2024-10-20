#include <iostream>
#include <functional>

#include "math/math.h"

void print_pose(const Eigen::Matrix4d &tf, const std::string &label = std::string()) {
    Eigen::Vector3d pos = tf.block<3, 1>(0, 3);
    Eigen::Vector3d euler = math::euler_zyx_from_rotation(math::rotation(tf));
    if (!label.empty()) {
        std::cout << label << " ";
    }
    std::cout << "pos " << pos.transpose() << "  ZYX: " << euler.transpose() * math::RAD_TO_DEG << std::endl;
}

void ur3e_test_fk() {
    std::cout << "Forward kinematics tests" << std::endl;
    std::vector<double> v1{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    print_pose(math::ur3e_space_fk(math::std_vector_to_eigen(v1) * math::DEG_TO_RAD));
    print_pose(math::ur3e_body_fk(math::std_vector_to_eigen(v1) * math::DEG_TO_RAD));
    std::cout << std::endl;

    std::vector<double> v2{0.0, 0.0, 0.0, -90.0, 0.0, 0.0};
    print_pose(math::ur3e_space_fk(math::std_vector_to_eigen(v2) * math::DEG_TO_RAD));
    print_pose(math::ur3e_body_fk(math::std_vector_to_eigen(v2) * math::DEG_TO_RAD));
    std::cout << std::endl;

    std::vector<double> v3{0.0, 0.0, -180.0, 0.0, 0.0, 0.0};
    print_pose(math::ur3e_space_fk(math::std_vector_to_eigen(v3) * math::DEG_TO_RAD));
    print_pose(math::ur3e_body_fk(math::std_vector_to_eigen(v3) * math::DEG_TO_RAD));
    std::cout << std::endl;

    std::vector<double> v4{0.0, 0.0, -90.0, 0.0, 0.0, 0.0};
    print_pose(math::ur3e_space_fk(math::std_vector_to_eigen(v4) * math::DEG_TO_RAD));
    print_pose(math::ur3e_body_fk(math::std_vector_to_eigen(v4) * math::DEG_TO_RAD));
}

void test_newton_raphson_root_find(const std::function<double(double)> &f, double x0) {
    auto [iterations, x_hat] = math::newton_raphson_root_find(f, x0);
    std::cout << "NR root f, x0=" << x0 << " -> it=" << iterations << " x=" << x_hat << " f(x)=" << f(x_hat) <<
            std::endl;
}

void test_gradient_descent_root_find(const std::function<double(double)> &f, double x0) {
    auto [iterations, x_hat] = math::gradient_descent_root_find(f, x0);
    std::cout << "GD root f, x0=" << x0 << " -> it=" << iterations << " x=" << x_hat << " f(x)=" << f(x_hat) <<
            std::endl;
}

void test_root_find() {
    std::cout << "Root finding tests" << std::endl;
    auto f1 = [](double x) {
        return (x - 3.0) * (x - 3.0) - 1.0;
    };
    test_newton_raphson_root_find(f1, -20);
    test_gradient_descent_root_find(f1, -20);
}

void ur3e_test_jacobian(const Eigen::VectorXd &joint_positions) {
    Eigen::Matrix4d tsb = math::ur3e_body_fk(joint_positions);
    auto [m, space_screws] = math::ur3e_space_chain();
    Eigen::MatrixXd jb = math::ur3e_body_jacobian(joint_positions);
    Eigen::MatrixXd js = math::ur3e_space_jacobian(joint_positions);
    Eigen::MatrixXd ad_tsb = math::adjoint_matrix(tsb);
    Eigen::MatrixXd ad_tbs = math::adjoint_matrix(tsb.inverse());
    std::cout << "Jb: " << std::endl << jb << std::endl << "Ad_tbs*Js:" << std::endl << ad_tbs * js << std::endl <<
            std::endl;
    std::cout << "Js: " << std::endl << js << std::endl << "Ad_tsb*Jb:" << std::endl << ad_tsb * jb << std::endl <<
            std::endl;
    std::cout << "d Jb: " << std::endl << jb - ad_tbs * js << std::endl << std::endl;
    std::cout << "d Js: " << std::endl << js - ad_tsb * jb << std::endl << std::endl;
}

void ur3e_test_jacobian() {
    std::cout << "Jacobian matrix tests" << std::endl;
    std::vector<double> v1{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    ur3e_test_jacobian(math::std_vector_to_eigen(v1));
    std::vector<double> v2{45.0, -20.0, 10.0, 2.5, 30.0, -50.0};
    ur3e_test_jacobian(math::std_vector_to_eigen(v2));
}

void ur3e_ik_test_pose(const Eigen::Vector3d &pos, const Eigen::Vector3d &zyx, const Eigen::VectorXd &j0)
{
    std::cout << "Test from pose" << std::endl;
    Eigen::Matrix4d t_sd = math::transformation_matrix(math::rotation_matrix_from_euler_zyx(zyx), pos);
    auto [iterations, j_ik] = math::ur3e_ik_body(t_sd, j0);
    Eigen::Matrix4d t_ik = math::ur3e_body_fk(j_ik);
    print_pose(t_ik, " IK pose");
    print_pose(t_sd, "Desired pose");
    std::cout << "Converged after " << iterations << " iterations" << std::endl;
    std::cout << "J_0: " << j0.transpose() * math::RAD_TO_DEG << std::endl;
    std::cout << "J_ik: " << j_ik.transpose() * math::RAD_TO_DEG << std::endl << std::endl;
}

void ur3e_ik_test_configuration(const Eigen::VectorXd &joint_positions, const Eigen::VectorXd &j0)
{
    std::cout << "Test from configuration" << std::endl;
    Eigen::Matrix4d t_sd = math::ur3e_space_fk(joint_positions);
    auto [iterations, j_ik] = math::ur3e_ik_body(t_sd, j0);
    Eigen::Matrix4d t_ik = math::ur3e_body_fk(j_ik);
    print_pose(t_ik, " IK pose");
    print_pose(t_sd, "Desired pose");
    std::cout << "Converged after " << iterations << " iterations" << std::endl;
    std::cout << "J_0: " << j0.transpose() << std::endl;
    std::cout << "J_d: " << joint_positions.transpose() << std::endl;
    std::cout << "J_ik: " << j_ik.transpose() * math::RAD_TO_DEG << std::endl << std::endl;
}

void ur3e_ik_test()
{
    std::vector<double> v1 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    Eigen::VectorXd j_t0 = math::std_vector_to_eigen(v1);
    std::vector<double> v2 {0.0, 0.0, -89.0, 0.0, 0.0, 0.0};
    Eigen::VectorXd j_t1 = math::std_vector_to_eigen(v2);
    ur3e_ik_test_pose(Eigen::Vector3d{0.3289, 0.22315, 0.36505}, Eigen::Vector3d{0.0, 90.0, -90.0}*math::RAD_TO_DEG, j_t0);
    ur3e_ik_test_pose(Eigen::Vector3d{0.3289, 0.22315, 0.36505}, Eigen::Vector3d{0.0, 90.0, -90.0}*math::RAD_TO_DEG, j_t1);
    std::vector<double> v3{50.0, -30.0, 20, 0.0, -30.0, 50.0};
    Eigen::VectorXd j_t2 = math::std_vector_to_eigen(v3);
    std::vector<double> v4{45.0, -20.0, 10.0, 2.5, 30.0,-50.0};
    Eigen::VectorXd j_d1 = math::std_vector_to_eigen(v4);
    ur3e_ik_test_configuration(j_d1, j_t0);
    ur3e_ik_test_configuration(j_d1, j_t2);
}

int main() {
    ur3e_test_fk();
    test_root_find();
    ur3e_test_jacobian();
    ur3e_ik_test();
    return 0;
}
