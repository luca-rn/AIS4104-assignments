//
// Created by Admin on 5/09/2024.
//
#ifndef MATH_H
#define MATH_H
#include <Eigen/Dense>

namespace math {
    const double DEG_TO_RAD = 0.0174532925;
    const double RAD_TO_DEG = 57.2957795;

    bool floatEquals(double a, double b);

    Eigen::Matrix3d skew_symmetric(Eigen::Vector3d v);

    Eigen::Matrix3d rotation(const Eigen::Matrix4d &tf);

    Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p);

    Eigen::Matrix3d rotate_x(double radians);

    Eigen::Matrix3d rotate_y(double radians);

    Eigen::Matrix3d rotate_z(double radians);

    Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e);

    Eigen::Vector3d euler_zyx_from_rotation(const Eigen::Matrix3d &r);

    Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v);

    Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h);

    Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf);

    double cot(double x);

    Eigen::Matrix3d matrix_exponential(const Eigen::Vector3d &w, double theta);

    std::pair<Eigen::Vector3d, double> matrix_logarithm(const Eigen::Matrix3d &r);

    Eigen::Matrix4d matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta);

    std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix4d &t);

    Eigen::Matrix4d planar_3r_fk_transform(const std::vector<double> &joint_positions);

    Eigen::Matrix4d planar_3r_fk_screw(const std::vector<double> &joint_positions);

    Eigen::Matrix4d ur3e_fk_screw(const std::vector<double> &joint_positions);

    Eigen::Matrix4d ur3e_fk_transform(const std::vector<double> &joint_positions);

    Eigen::VectorXd std_vector_to_eigen(std::vector<double> &v);

    Eigen::Matrix4d matrix_exponential(Eigen::VectorXd &screw, double theta);

    bool is_average_below_eps(const std::vector<double> &values, double eps, const uint8_t n_values);

    std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd> > ur3e_space_chain();

    Eigen::Matrix4d ur3e_space_fk(const Eigen::VectorXd &joint_positions);

    std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd> > ur3e_body_chain();

    Eigen::Matrix4d ur3e_body_fk(const Eigen::VectorXd &joint_positions);

    std::pair<uint32_t, double> newton_raphson_root_find(const std::function<double(double)> &f, double x_0,
                                                         double dx_0 = 0.5, double eps = 10e-7);

    std::pair<uint32_t, double> gradient_descent_root_find(const std::function<double(double)> &f, double x_0,
                                                           double gamma = 0.1, double dx_0 = 0.5, double eps = 10e-7);

    Eigen::MatrixXd ur3e_space_jacobian(const Eigen::VectorXd &current_joint_positions);
    Eigen::Matrix4d matrix_negexponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta);
    Eigen::MatrixXd ur3e_body_jacobian(const Eigen::VectorXd &current_joint_positions);
    std::pair<size_t, Eigen::VectorXd> ur3e_ik_body(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd &current_joint_positions, double gamma = 1e-2, double v_e = 4e-3, double w_e = 4e-3);
}

#endif //MATH_H
