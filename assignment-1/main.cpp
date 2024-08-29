#include <iostream>

#include <Eigen/Dense>

double rad_to_deg(double radians) {
    return radians * 57.2957795;
}

double deg_to_rad(double degrees) {
    return degrees * 0.0174532925;
}

Eigen::Matrix3d rotate_x(double degrees) {
    double radians = deg_to_rad(degrees);
    Eigen::Matrix3d matrix{
        1.0, 0.0, 0.0,
        0.0, std::cos(radians), -std::sin(radians),
        0.0, std::sin(radians), std::cos(radians)
    };
    return matrix;
}

Eigen::Matrix3d skew_symmetric(Eigen::Vector3d v) {
    Eigen::Matrix3d matrix{
        0.0, -v.z(), v.y(),
        v.z(), 0.0, -v.x(),
        -v.y(), v.x(), 0.0
    };
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
    example(2.0);
    return 0;
}
