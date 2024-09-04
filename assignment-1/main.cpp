#include <iostream>

#include <Eigen/Dense>

double rad_to_deg(double radians) {
    return radians * 57.2957795;
}

double deg_to_rad(double degrees) {
    return degrees * 0.0174532925;
}

Eigen::Matrix3d skew_symmetric(Eigen::Vector3d v) {
    Eigen::Matrix3d matrix;
    matrix <<
        0.0, -v.z(), v.y(),
        v.z(), 0.0, -v.x(),
        -v.y(), v.x(), 0.0;
    return matrix;
}

void skew_symmetric_test()
{
Eigen::Matrix3d skew_matrix = skew_symmetric(Eigen::Vector3d{0.5, 0.5, 0.707107});
std::cout << "Skew-symmetric matrix: " << std::endl;
std::cout << skew_matrix << std::endl;
std::cout << "Skew-symmetric matrix transposition: " << std::endl;
std::cout << -skew_matrix.transpose() << std::endl;
}

Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x,
const Eigen::Vector3d &y,
const Eigen::Vector3d &z)
{
    Eigen::Matrix3d matrix;
    // implement the necessary equations and functionality.
    matrix <<
        x,y,z;
    return matrix;
}

Eigen::Matrix3d rotate_x(double degrees)
{
    double radians = deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    // implement the necessary equations and functionality.
    matrix <<
        1.0, 0.0, 0.0,
        0.0, std::cos(radians), -std::sin(radians),
        0.0, std::sin(radians), std::cos(radians);
    return matrix;
}

Eigen::Matrix3d rotate_y(double degrees)
{
    double radians = deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    // implement the necessary equations and functionality.
    matrix <<
        std::cos(radians), 0.0, std::sin(radians),
        0.0, 1.0, 0.0,
        -std::sin(radians), 0.0, std::cos(radians);

    return matrix;
}

Eigen::Matrix3d rotate_z(double degrees)
{
    double radians = deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    // implement the necessary equations and functionality.
    matrix <<
        std::cos(radians), -std::sin(radians), 0.0,
        std::sin(radians), std::cos(radians), 0.0,
        0.0, 0.0, 1.0;

    return matrix;
}

Eigen::CommaInitializer<Eigen::Matrix<double, 3, 3>> operator^(const Eigen::CommaInitializer<Eigen::Matrix<double, 3, 3>>& lhs, double rhs);

Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double degrees)
{
    double radians = deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    // implement the necessary equations and functionality.
    matrix <<
        std::cos(radians)+axis.x()*axis.x()*(1-std::cos(radians)),
        axis.x()*axis.y()*(1-std::cos(radians))-axis.z()*std::sin(radians),
        axis.x()*axis.z()*(1-std::cos(radians))+axis.y()*std::sin(radians),
        axis.x()*axis.y()*(1-std::cos(radians))+axis.z()*std::sin(radians),
        std::cos(radians)+axis.y()*axis.y()*(1-std::cos(radians)),
        axis.y()*axis.z()*(1-std::cos(radians))-axis.x()*std::sin(radians),
        axis.x()*axis.z()*(1-std::cos(radians))-axis.y()*std::sin(radians),
        axis.y()*axis.z()*(1-std::cos(radians))+axis.x()*std::sin(radians),
        std::cos(radians)+axis.z()*axis.z()*(1-std::cos(radians));
    return matrix;
}

// e = [α β γ]
Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e)
{
    Eigen::Matrix3d matrix = rotate_z(e(0))*rotate_y(e(1))*rotate_x(e(2));
    return matrix;
}

void rotation_matrix_test()
{
    Eigen::Matrix3d rot =
    rotation_matrix_from_euler_zyx(Eigen::Vector3d{45.0, -45.0, 90.0});
    Eigen::Matrix3d rot_aa =
    rotation_matrix_from_axis_angle(Eigen::Vector3d{0.8164966, 0.0, 0.5773503}, 120.0);
    Eigen::Matrix3d rot_fa =
    rotation_matrix_from_frame_axes(Eigen::Vector3d{0.5, 0.5, 0.707107},
    Eigen::Vector3d{-0.5, -0.5, 0.707107},
    Eigen::Vector3d{0.707107, -0.707107, 0.0});
    std::cout << "Rotation matrix from Euler: " << std::endl;
    std::cout << rot << std::endl << std::endl;
    std::cout << "Rotation matrix from axis-angle pair: " << std::endl;
    std::cout << rot_aa << std::endl << std::endl;
    std::cout << "Rotation matrix from frame axes: " << std::endl;
    std::cout << rot_fa << std::endl << std::endl;
}

Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
{
    Eigen::Matrix4d matrix;
    // implement the necessary equations and functionality.
    matrix <<
        r(0,0), r(0,1), r(0,2), p(0),
        r(1,0), r(1,1), r(1,2), p(1),
        r(2,0), r(2,1), r(2,2), p(2),
        0, 0, 0, 1;
return matrix;
}

void transformation_matrix_test()
{
    Eigen::Matrix3d r = rotation_matrix_from_euler_zyx(Eigen::Vector3d{45, -45.0, 90.0});
    Eigen::Vector3d v{1.0, -2.0, 3.0};
    std::cout << "transformation_matrix: " << std::endl;
    std::cout << transformation_matrix(r, v) << std::endl;
}

void transform_vector() {
    Eigen::Vector4d v_a_4d = {2.5, 3.0, -10.0, 1.0};
    Eigen::Vector3d translation = {0.0, 0.0, 10.0};
    Eigen::Vector3d eulerZYX = {60.0, 45.0, 0.0};
    Eigen::Matrix3d rotation = rotation_matrix_from_euler_zyx(eulerZYX);
    Eigen::Matrix4d t_a_w = transformation_matrix(rotation, translation);

    Eigen::Vector4d v_w_4d = t_a_w * v_a_4d;
    Eigen::Vector3d v_w = {v_w_4d(0), v_w_4d(1), v_w_4d(2)};

    std::cout << "Vw" << std::endl;
    std::cout << v_w << std::endl;
}
void example(double constant)
{
    Eigen::Matrix3d identity;
    identity <<
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;
    std::cout << "I: " << std::endl << identity << std::endl << std::endl;
    std::cout << constant <<"*I: " << std::endl << constant * identity << std::endl << std::endl;
}

int main()
{
    skew_symmetric_test();
    rotation_matrix_test();
    transformation_matrix_test();
    transform_vector();
    return 0;
}
