//
// Created by yao on 19-12-6.
//

#include "ros/ros.h"
#include "acado_toolkit.hpp"
USING_NAMESPACE_ACADO

template <typename T>
void print(std::string name, T vari, int N) {
    std::cout<< name <<" = [";
    for (int j = 0; j < N; ++j) {
        std::cout<<"("<<vari(j, 0) << ", " << vari(j, 1) << "), ";
    }
    std::cout << "\b\b]" << std::endl;
}

int main(int argc, char **argv) {
    std::vector<int> time;
    ros::init(argc, argv, "topic_publisher");
    ros::NodeHandle nh;
    ros::Rate loop_rate(10);

    // — state variables (acadoVariables.x)—
    DifferentialState x;
    DifferentialState y;
    DifferentialState phi;
    DifferentialState v;
    DifferentialState w;
// — control inputs —
    Control a;
    Control w_dot;
// —- differential equations —-
    double g = 9.81;      // vehicle wheel base
    DifferentialEquation f;
    f << dot(x) == v * cos(phi);
    f << dot(y) == v * sin(phi);
    f << dot(phi) == w;
    f << dot(v) == a;
    f << dot(w) == w_dot;

    Function h;
    h << x << y << phi << v << a << w_dot;

    // LSQ coefficient matrix
    DMatrix Q(6, 6);
    Q(0, 0) = 1.0;
    Q(1, 1) = 1.0;
    Q(2, 2) = 1.0;
    Q(3, 3) = 1.0;
    Q(4, 4) = 1.0;
    Q(5, 5) = 1.0;

    double N = 12;
    double dt = 0.1;

    while (ros::ok()) {
        clock_t t_start = clock();
        // Reference
        Grid time_grid(0., (N - 1) * dt, N);
        VariablesGrid reference(6, time_grid);

        reference(0, 0) = 0.;
        reference(0, 1) = 1.;
        reference(0, 2) = - M_PI / 6;
        reference(0, 3) = 5.;
        reference(0, 4) = 0.;
        reference(0, 5) = 0.;

        for (int i = 1; i < N; ++i) {
            reference(i, 0) = reference(i - 1, 0) + reference(i - 1, 3) * cos(reference(i - 1, 2)) * dt;
            reference(i, 1) = reference(i - 1, 1) + reference(i - 1, 3) * sin(reference(i - 1, 2)) * dt;
            reference(i, 2) = reference(i - 1, 2) + (double(rand()) / RAND_MAX - 0.5) / 10;
            reference(i, 3) = reference(i - 1, 3) + double(rand()) / RAND_MAX - 0.5;
            reference(i, 4) = 0.;
            reference(i, 5) = 0.;
        }
        print("refer_point", reference, N);

        OCP ocp(reference.getTimePoints());
        ocp.minimizeLSQ(Q, h, reference);

        ocp.subjectTo(f);
        ocp.subjectTo(-0.4 * g <= v * w <= 0.4 * g);
        ocp.subjectTo(-2 <= a <= 2);
        ocp.subjectTo(-2 <= w_dot <= 2);
        ocp.subjectTo(AT_START, x == 0.);
        ocp.subjectTo(AT_START, y == 0.);
        ocp.subjectTo(AT_START, v == 1.);
        ocp.subjectTo(AT_START, phi == -M_PI/6);
        ocp.subjectTo(AT_START, w == 0.);

        OptimizationAlgorithm algorithm(ocp);
        algorithm.set(PRINTLEVEL, NONE);
        algorithm.set(PRINT_COPYRIGHT, NONE);
        algorithm.solve();

        VariablesGrid states, controls;

        algorithm.getDifferentialStates(states);
        algorithm.getControls(controls);

        print("control_list", controls, N);

        time.emplace_back(int(1000 * (clock() - t_start) / CLOCKS_PER_SEC) + 1);
        ROS_INFO_STREAM("cost time: " << time.back());
        loop_rate.sleep();
    }
    double mean_time = std::accumulate(time.begin(), time.end(), 0) / time.size();
    std::cout<<time.size()<<std::endl;
    std::cout<<mean_time<<std::endl;
    return 0;
}