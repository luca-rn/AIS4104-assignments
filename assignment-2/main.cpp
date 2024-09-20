#include <iostream>

#include "math/math.h"

void wrench_in_frames() {
    Eigen::Vector3d f_w{-30, 0 , 0};
    Eigen::Vector3d m_s{0, 0, 2};
    Eigen::Vector3d e_ws{60, -60, 0};

    Eigen::Matrix3d r_ws = math::rotate_y(e_ws(0)*math::DEG_TO_RAD)*math::rotate_z(e_ws(1)*math::DEG_TO_RAD)*math::rotate_x(e_ws(2)*math::DEG_TO_RAD);
    Eigen::Vector3d f_s = r_ws.transpose()*f_w;
    Eigen::Vector3d m_w = r_ws*m_s;

    std::cout << "f_w: " << f_w.transpose() << std::endl;
    std::cout << "m_w: " << m_w.transpose() << std::endl << std::endl;
    std::cout << "f_s: " << f_s.transpose() << std::endl;
    std::cout << "m_s: " << m_s.transpose() << std::endl << std::endl;
}

void sum_of_wrenches_different_frames() {
    double m_a = 0.1, g = 10.0, m_h = 0.5;
    Eigen::VectorXd f_h(6);
    f_h << 0, 0, 0, 0, -5, 0;
    Eigen::VectorXd f_a(6);
    f_a << 0, 0, 0, 0, 0, 1;

    Eigen::Matrix4d t_hf;
    t_hf << 1, 0, 0, -0.1,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1;
    Eigen::Matrix4d t_af;
    t_af << 1, 0, 0, -0.25,
    0, 0, 1, 0,
    0, -1, 0, 0,
    0, 0, 0, 1;

    /*
    std::cout << "AdjThf " << std::endl <<adjoint_matrix(t_hf).transpose() << std::endl;
    std::cout << "first term: " << std::endl <<adjoint_matrix(t_hf.transpose())*f_h << std::endl;
    */
    Eigen::VectorXd f_f(6);
    f_f = math::adjoint_matrix(t_hf).transpose()*f_h + math::adjoint_matrix(t_af).transpose()*f_a;

    std::cout << "f_f: " << f_f.transpose() << std::endl;
}

void print_pose(const std::string &label, const Eigen::Matrix4d &tf) {
    Eigen::Matrix3d r = tf.block(0, 0, 3, 3);
    Eigen::Vector3d p = tf.block(0, 3, 3, 1);
    Eigen::Vector3d zyx = math::euler_zyx_from_rotation(r);

    std::cout << label << std::endl;
    std::cout << "Euler ZYX Angles: " << zyx.transpose()*math::RAD_TO_DEG << std::endl;
    std::cout << "Linear Position: " << p.transpose() << std::endl << std::endl;
}

void test_planar_3r_fk_transform() {
    std::vector<double> j1{0.0, 0.0, 0.0};
    std::string l1 = "j1 test_planar_3r_fk_transform";
    std::vector<double> j2{90.0, 0.0, 0.0};
    std::string l2 = "j2 test_planar_3r_fk_transform";
    std::vector<double> j3{0.0, 90.0, 0.0};
    std::string l3 = "j3 test_planar_3r_fk_transform";
    std::vector<double> j4{0.0, 0.0, 90.0};
    std::string l4 = "j4 test_planar_3r_fk_transform";
    std::vector<double> j5{10.0, -15.0, 2.75};
    std::string l5 = "j5 test_planar_3r_fk_transform";

    print_pose(l1,math::planar_3r_fk_transform(j1));
    print_pose(l2,math::planar_3r_fk_transform(j2));
    print_pose(l3,math::planar_3r_fk_transform(j3));
    print_pose(l4,math::planar_3r_fk_transform(j4));
    print_pose(l5,math::planar_3r_fk_transform(j5));
}

void test_planar_3r_fk_screw() {
    std::vector<double> j1{0.0, 0.0, 0.0};
    std::string l1 = "j1 test_planar_3r_fk_screw";
    std::vector<double> j2{90.0, 0.0, 0.0};
    std::string l2 = "j2 test_planar_3r_fk_screw";
    std::vector<double> j3{0.0, 90.0, 0.0};
    std::string l3 = "j3 test_planar_3r_fk_screw";
    std::vector<double> j4{0.0, 0.0, 90.0};
    std::string l4 = "j4 test_planar_3r_fk_screw";
    std::vector<double> j5{10.0, -15.0, 2.75};
    std::string l5 = "j5 test_planar_3r_fk_screw";

    print_pose(l1,math::planar_3r_fk_screw(j1));
    print_pose(l2,math::planar_3r_fk_screw(j2));
    print_pose(l3,math::planar_3r_fk_screw(j3));
    print_pose(l4,math::planar_3r_fk_screw(j4));
    print_pose(l5,math::planar_3r_fk_screw(j5));
}

void test_ur3e_fk_screw() {
    std::vector<double> j1{0.0, 0.0, 0.0, -90.0, 0.0, 0.0};
    std::string l1 = "j1 test_ur3e_fk_screw";
    std::vector<double> j2{0.0, -180.0, 0.0, 0.0, 0.0, 0.0};
    std::string l2 = "j2 test_ur3e_fk_screw";
    std::vector<double> j3{0.0, -90.0, 0.0, 0.0, 0.0, 0.0};
    std::string l3 = "j3 test_ur3e_fk_screw";

    print_pose(l1,math::ur3e_fk_screw(j1));
    print_pose(l2,math::ur3e_fk_screw(j2));
    print_pose(l3,math::ur3e_fk_screw(j3));
}

void test_ur3e_fk_transform() {
    std::vector<double> j1{0.0, 0.0, 0.0, -90.0, 0.0, 0.0};
    std::string l1 = "j1 test_ur3e_fk_transform";
    std::vector<double> j2{0.0, -180.0, 0.0, 0.0, 0.0, 0.0};
    std::string l2 = "j2 test_ur3e_fk_transform";
    std::vector<double> j3{0.0, -90.0, 0.0, 0.0, 0.0, 0.0};
    std::string l3 = "j3 test_ur3e_fk_transform";

    print_pose(l1,math::ur3e_fk_transform(j1));
    print_pose(l2,math::ur3e_fk_transform(j2));
    print_pose(l3,math::ur3e_fk_transform(j3));
}

int main()
{
    wrench_in_frames();
    sum_of_wrenches_different_frames();
    test_planar_3r_fk_transform();
    test_planar_3r_fk_screw();
    test_ur3e_fk_screw();
    test_ur3e_fk_transform();
    return 0;
}
