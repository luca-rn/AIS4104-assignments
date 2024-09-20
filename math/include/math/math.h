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
};

#endif //MATH_H
