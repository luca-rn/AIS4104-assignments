//
// Created by Admin on 5/09/2024.
//
#ifndef MATH_H
#define MATH_H
#include <Eigen/Dense>

namespace math {
    bool floatEquals(double a, double b);
    Eigen::Matrix3d rotate_x(double radians);
    Eigen::Matrix3d rotate_y(double radians);
    Eigen::Matrix3d rotate_z(double radians);
    Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e);
    Eigen::Vector3d euler_xyz_from_rotation(const Eigen::Matrix3d &r);
};

#endif //MATH_H
