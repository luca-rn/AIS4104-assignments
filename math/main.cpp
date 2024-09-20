//
// Created by Admin on 5/09/2024.
//
#include <iostream>
// #include <Eigen/Dense>
#include "include/math/math.h"

bool math::floatEquals(double a, double b) {
    return std::abs(a - b) < 1e-6;
}

Eigen::Matrix3d math::skew_symmetric(Eigen::Vector3d v) {
    Eigen::Matrix3d matrix;
    matrix <<
            0.0, -v.z(), v.y(),
            v.z(), 0.0, -v.x(),
            -v.y(), v.x(), 0.0;
    return matrix;
}

Eigen::Matrix4d math::transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p) {
    Eigen::Matrix4d matrix;
    // implement the necessary equations and functionality.
    matrix <<
            r(0, 0), r(0, 1), r(0, 2), p(0),
            r(1, 0), r(1, 1), r(1, 2), p(1),
            r(2, 0), r(2, 1), r(2, 2), p(2),
            0, 0, 0, 1;
    return matrix;
}

Eigen::Matrix3d math::rotate_x(double radians) {
    Eigen::Matrix3d matrix;
    // implement the necessary equations and functionality.
    matrix <<
            1.0, 0.0, 0.0,
            0.0, std::cos(radians), -std::sin(radians),
            0.0, std::sin(radians), std::cos(radians);
    return matrix;
}

Eigen::Matrix3d math::rotate_y(double radians) {
    Eigen::Matrix3d matrix;
    // implement the necessary equations and functionality.
    matrix <<
            std::cos(radians), 0.0, std::sin(radians),
            0.0, 1.0, 0.0,
            -std::sin(radians), 0.0, std::cos(radians);

    return matrix;
}

Eigen::Matrix3d math::rotate_z(double radians) {
    Eigen::Matrix3d matrix;
    // implement the necessary equations and functionality.
    matrix <<
            std::cos(radians), -std::sin(radians), 0.0,
            std::sin(radians), std::cos(radians), 0.0,
            0.0, 0.0, 1.0;

    return matrix;
}

Eigen::Matrix3d math::rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e) {
    Eigen::Matrix3d matrix = rotate_z(e(0)) * rotate_y(e(1)) * rotate_x(e(2));
    return matrix;
}

Eigen::Vector3d math::euler_zyx_from_rotation(const Eigen::Matrix3d &r) {
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    if (floatEquals(r(2, 0), -1.0)) {
        b = EIGEN_PI / 2.0;
        a = 0.0;
        c = std::atan2(r(0, 1), r(1, 1));
    } else if (floatEquals(r(2, 0), 1.0)) {
        b = -(EIGEN_PI / 2.0);
        a = 0.0;
        c = -std::atan2(r(0, 1), r(1, 1));
    } else {
        b = std::atan2(-r(2, 0), std::sqrt(r(0, 0) * r(0, 0) + r(1, 0) * r(1, 0)));
        a = std::atan2(r(1, 0), r(0, 0));
        c = std::atan2(r(2, 1), r(2, 2));
    }
    return Eigen::Vector3d{a, b, c};
}

Eigen::VectorXd math::twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v) {
    Eigen::VectorXd twist(6);
    twist << w(0), w(1), w(2), v(0), v(1), v(2);
    return twist;
}

Eigen::VectorXd math::screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h) {
    Eigen::VectorXd screw_axis(6);
    screw_axis << s, skew_symmetric(s) * q + h * s;
    return screw_axis;
}

Eigen::MatrixXd math::adjoint_matrix(const Eigen::Matrix4d &tf) {
    Eigen::Matrix3d r = tf.block(0, 0, 3, 3);
    Eigen::Vector3d p = tf.block(0, 3, 3, 1);
    Eigen::MatrixXd adjoint(6, 6);
    adjoint << r, Eigen::Matrix3d::Zero(), skew_symmetric(p) * r, r;
    /*
    Eigen::MatrixXd adjoint(6,6);
    adjoint << tf(0,0), tf(0,1), tf(0,2), 0, 0, 0,
    tf(1,0), tf(1,1), tf(1,2), 0, 0, 0,
    tf(2,0), tf(2,1), tf(2,2), 0, 0, 0,
    -tf(2,3)*tf(1,0)+tf(1,3)*tf(2,0),
    -tf(2,3)*tf(1,1)+tf(1,3)*tf(2,1),
    -tf(2,3)*tf(1,2)+tf(1,3)*tf(2,2), tf(0,0), tf(0,1), tf(0,2),
    tf(2,3)*tf(0,0)-tf(0,3)*tf(2,0),
    tf(2,3)*tf(0,1)-tf(0,3)*tf(2,1),
    tf(2,3)*tf(0,2)-tf(0,3)*tf(2,2), tf(1,0), tf(1,1), tf(1,2),
    -tf(1,3)*tf(0,0)+tf(0,3)*tf(1,0),
    -tf(1,3)*tf(0,1)+tf(0,3)*tf(1,1),
    -tf(1,3)*tf(0,2)+tf(0,3)*tf(1,2), tf(2,0), tf(2,1), tf(2,2);
    */
    return adjoint;
}

double math::cot(double x) {
    return std::cos(x * DEG_TO_RAD) / std::sin(x * DEG_TO_RAD);
}

Eigen::Matrix3d math::matrix_exponential(const Eigen::Vector3d &w, double theta) {
    return Eigen::Matrix3d::Identity() + std::sin(theta * DEG_TO_RAD) * skew_symmetric(w) + (
               1 - std::cos(theta * DEG_TO_RAD)) * skew_symmetric(w) * skew_symmetric(w);
}

std::pair<Eigen::Vector3d, double> math::matrix_logarithm(const Eigen::Matrix3d &r) {
    Eigen::Vector3d w = Eigen::Vector3d::Zero();
    double theta = 0.0;
    if (floatEquals(r.trace(), -1.0)) {
        theta = EIGEN_PI;
        if (!floatEquals(r(2, 2), -1.0)) {
            w = 1 / (2 * (1 + r(2, 2))) * Eigen::Vector3d{r(0, 2), r(1, 2), 1 + r(2, 2)};
        } else if (!floatEquals(r(1, 1), -1.0)) {
            w = 1 / (2 * (1 + r(1, 1))) * Eigen::Vector3d{r(0, 1), 1 + r(1, 1), r(2, 1)};
        } else {
            w = 1 / (2 * (1 + r(0, 0))) * Eigen::Vector3d{1 + r(0, 0), r(1, 0), r(2, 0)};
        }
    } else if (!r.isIdentity()) {
        theta = 1 / std::cos(0.5 * (r.trace() - 1));
        Eigen::Matrix3d w_skew = 1 / (2 * std::sin(theta)) * (r - r.transpose());
        w = {w_skew(2, 1), w_skew(1, 0), w_skew(0, 2)};
    }
    return std::make_pair(w, theta);
}

Eigen::Matrix4d math::matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta) {
    Eigen::Matrix4d t = Eigen::Matrix4d::Zero();
    if (floatEquals(w.norm(), 1.0)) {
        Eigen::Matrix3d e_w = matrix_exponential(w, theta);
        Eigen::Vector3d g_v = (Eigen::Matrix3d::Identity() * theta * DEG_TO_RAD + (
                                   1 - std::cos(theta * DEG_TO_RAD)) * skew_symmetric(w) + (
                                   theta * DEG_TO_RAD - std::sin(theta * DEG_TO_RAD)) * skew_symmetric(w) *
                               skew_symmetric(w)) * v;
        t = transformation_matrix(e_w, g_v);
    } else if (floatEquals(v.norm(), 1.0) && floatEquals(w.norm(), 0.0)) {
        t = transformation_matrix(Eigen::Matrix3d::Identity(), v * theta * DEG_TO_RAD);
    }
    return t;
}

std::pair<Eigen::VectorXd, double> math::matrix_logarithm(const Eigen::Matrix4d &t) {
    double theta = 0.0;
    Eigen::VectorXd s(6);
    Eigen::Vector3d w = Eigen::Vector3d::Zero();
    Eigen::Vector3d v = Eigen::Vector3d::Zero();
    Eigen::Matrix3d r = t.block(0, 0, 3, 3);
    Eigen::Vector3d p = t.block(0, 3, 3, 1);
    if (!r.isIdentity()) {
        std::pair<Eigen::Vector3d, double> x = math::matrix_logarithm(r);
        w = x.first;
        theta = x.second;
        v = (1 / (theta * DEG_TO_RAD) * Eigen::Matrix3d::Identity() - 0.5 * skew_symmetric(w) + (
                 1 / (theta * DEG_TO_RAD) - 0.5 * cot(theta) / 2) * skew_symmetric(w) * skew_symmetric(
                 w)) * p;
    } else {
        v = p / p.norm();
        theta = p.norm();
    }
    s << w, v;
    return std::make_pair(s, theta);
}

Eigen::Matrix4d math::planar_3r_fk_transform(const std::vector<double> &joint_positions) {
    double l1 = 10.0, l2 = 10.0, l3 = 10.0;
    double theta1 = joint_positions.at(0) * DEG_TO_RAD, theta2 = joint_positions.at(1) * DEG_TO_RAD;
    double theta3 = joint_positions.at(2) * DEG_TO_RAD;
    Eigen::Matrix4d t01 = transformation_matrix(rotate_z(theta1), {0, 0, 0});
    Eigen::Matrix4d t12 = transformation_matrix(rotate_z(theta2), {l1, 0, 0});
    Eigen::Matrix4d t23 = transformation_matrix(rotate_z(theta3), {l2, 0, 0});
    Eigen::Matrix4d t34 = transformation_matrix(Eigen::Matrix3d::Identity(), {l3, 0, 0});

    /*
    Eigen::Matrix4d t01 = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d t12 = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d t23 = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d t34 = Eigen::Matrix4d::Identity();
    t01.block(0, 0, 2, 2) <<
            std::cos(theta1), -std::sin(theta1),
            std::sin(theta1), std::cos(theta1);
    t12.block(0, 0, 2, 2) <<
        std::cos(theta2), -std::sin(theta2),
            std::sin(theta2), std::cos(theta2);
    t23.block(0, 0, 2, 2) <<
        std::cos(theta3), -std::sin(theta3),
            std::sin(theta3), std::cos(theta3);
    t12(0,3) = l1;
    t23(0,3) = l2;
    t34(0,3) = l3;
    */
    return t01 * t12 * t23 * t34;
}

Eigen::Matrix4d math::planar_3r_fk_screw(const std::vector<double> &joint_positions) {
    double l1 = 10.0, l2 = 10.0, l3 = 10.0;
    double theta1 = joint_positions.at(0), theta2 = joint_positions.at(1);
    double theta3 = joint_positions.at(2);
    Eigen::Vector3d w = {0.0, 0.0, 1.0};
    Eigen::Vector3d q1 = {0.0, 0.0, 0.0};
    Eigen::Vector3d q2 = {l1, 0.0, 0.0};
    Eigen::Vector3d q3 = {l1 + l2, 0.0, 0.0};
    Eigen::VectorXd s1 = screw_axis(w, q1, 0);
    Eigen::VectorXd s2 = screw_axis(w, q2, 0);
    Eigen::VectorXd s3 = screw_axis(w, q3, 0);

    Eigen::Vector3d v1 = {s1(3), s1(4), s1(5)};
    Eigen::Vector3d v2 = {s2(3), s2(4), s2(5)};
    Eigen::Vector3d v3 = {s3(3), s3(4), s3(5)};

    Eigen::Matrix4d m = Eigen::Matrix4d::Identity();
    m(0, 3) = l1 + l2 + l3;
    Eigen::Matrix4d t04 = matrix_exponential(w, v1, theta1) * matrix_exponential(w, v2, theta2) *
                          matrix_exponential(w, v3, theta3) * m;
    return t04;
}

Eigen::Matrix4d math::ur3e_fk_screw(const std::vector<double> &joint_positions) {
    /*
    double h1 = 0.089, h2 = 0.095;
    double l1 = 0.425, l2 = 0.392;
    double w1 =  0.109, w2 = 0.082;
    */
    double h1 = 0.15185, h2 = 0.08535;
    double l1 = 0.24355, l2 = 0.2132;
    double w1 = 0.13105, w2 = 0.0921;

    Eigen::Vector3d w_1 = {0.0, 0.0, 1.0};
    Eigen::Vector3d w_2 = {0.0, 1.0, 0.0};
    Eigen::Vector3d w_3 = {0.0, 1.0, 0.0};
    Eigen::Vector3d w_4 = {0.0, 1.0, 0.0};
    Eigen::Vector3d w_5 = {0.0, 0.0, -1.0};
    Eigen::Vector3d w_6 = {0.0, 1.0, 0.0};

    Eigen::Vector3d v_1 = {0.0, 0.0, 0.0};
    Eigen::Vector3d v_2 = {-h1, 0.0, 0.0};
    Eigen::Vector3d v_3 = {-h1, 0.0, l1};
    Eigen::Vector3d v_4 = {-h1, 0.0, l1 + l2};
    Eigen::Vector3d v_5 = {-w1, l1 + l2, 0.0};
    Eigen::Vector3d v_6 = {h2 - h1, 0.0, l1 + l2};

    double theta1 = joint_positions.at(0);
    double theta2 = joint_positions.at(1);
    double theta3 = joint_positions.at(2);
    double theta4 = joint_positions.at(3);
    double theta5 = joint_positions.at(4);
    double theta6 = joint_positions.at(5);

    Eigen::Matrix4d m;
    m <<
            -1.0, 0.0, 0.0, l1 + l2,
            0.0, 0.0, 1.0, w1 + w2,
            0.0, 1.0, 0.0, h1 - h2,
            0.0, 0.0, 0.0, 1.0;

    Eigen::Matrix4d t01 = matrix_exponential(w_1, v_1, theta1);
    Eigen::Matrix4d t12 = matrix_exponential(w_2, v_2, theta2);
    Eigen::Matrix4d t23 = matrix_exponential(w_3, v_3, theta3);
    Eigen::Matrix4d t34 = matrix_exponential(w_4, v_4, theta4);
    Eigen::Matrix4d t45 = matrix_exponential(w_5, v_5, theta5);
    Eigen::Matrix4d t56 = matrix_exponential(w_6, v_6, theta6);
    Eigen::Matrix4d t06 = t01 * t12 * t23 * t34 * t45 * t56 * m;
    return t06;
}

Eigen::Matrix4d math::ur3e_fk_transform(const std::vector<double> &joint_positions) {
    /*
    double h1 = 0.089, h2 = 0.095;
    double l1 = 0.425, l2 = 0.392;
    double w1 =  0.109, w2 = 0.082;
    */

    double h1 = 0.15185, h2 = 0.08535;
    double l1 = 0.24355, l2 = 0.2132;
    double w1 = 0.13105, w2 = 0.0921;

    double theta1 = joint_positions.at(0) * DEG_TO_RAD, theta2 = joint_positions.at(1) * DEG_TO_RAD;
    double theta3 = joint_positions.at(2) * DEG_TO_RAD, theta4 = joint_positions.at(3) * DEG_TO_RAD;
    double theta5 = joint_positions.at(4) * DEG_TO_RAD, theta6 = joint_positions.at(5) * DEG_TO_RAD;

    Eigen::Matrix3d r_sb;
    r_sb <<
            -1.0, 0.0, 0.0,
            0.0, 0.0, 1.0,
            0.0, 1.0, 0.0;

    Eigen::Matrix4d t01 = transformation_matrix(rotate_z(theta1), {0.0, 0.0, 0.0});
    Eigen::Matrix4d t12 = transformation_matrix(rotate_y(theta2), {0.0, 0.0, h1});
    Eigen::Matrix4d t23 = transformation_matrix(rotate_y(theta3), {l1, 0.0, 0.0});
    Eigen::Matrix4d t34 = transformation_matrix(rotate_y(theta4), {l2, 0.0, 0.0});
    Eigen::Matrix4d t45 = transformation_matrix(rotate_z(-theta5), {0.0, w1, 0.0});
    Eigen::Matrix4d t56 = transformation_matrix(rotate_y(theta6), {0.0, w2, -h2});
    Eigen::Matrix4d t67 = transformation_matrix(r_sb, {0.0, 0.0, 0.0});
    Eigen::Matrix4d t07 = t01 * t12 * t23 * t34 * t45 * t56 * t67;

    /*
    Eigen::Matrix4d t01 = transformation_matrix(rotate_z(theta1),{0.0,0.0,h1});
    Eigen::Matrix4d t12 = transformation_matrix(rotate_y(theta2),{l1,0.0,0.0});
    Eigen::Matrix4d t23 = transformation_matrix(rotate_y(theta3), {l2, 0.0, 0.0});
    Eigen::Matrix4d t34 = transformation_matrix(rotate_y(theta4), {0.0, w1, 0.0});
    Eigen::Matrix4d t45 = transformation_matrix(rotate_z(-theta5), {0.0, 0.0, -h2});
    Eigen::Matrix4d t56 = transformation_matrix(rotate_y(-theta6), {0.0, w2, 0.0});
    Eigen::Matrix4d t06 = t01*t12*t23*t34*t45*t56;
    */
    return t07;
}
