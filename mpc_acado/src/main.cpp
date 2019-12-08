//
// Created by yao on 19-12-6.
//

#include"ros/ros.h"
#include "acado_toolkit.hpp"
USING_NAMESPACE_ACADO

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
    DifferentialState delta;
// — control inputs —
    Control a;
    Control delta_rate;
// —- differential equations —-
    double L = 1.32;      // vehicle wheel base
    DifferentialEquation f;
    f << dot(x) == v * cos(phi);
    f << dot(y) == v * sin(phi);
    f << dot(phi) == v * tan(delta) / L;
    f << dot(v) == a;
    f << dot(delta) == delta_rate;

    Function h;
    h << x << y << phi << v << a << delta_rate;

    // LSQ coefficient matrix
    DMatrix Q(6, 6);
    Q(0, 0) = 1.0;
    Q(1, 1) = 1.0;
    Q(2, 2) = 10.0;
    Q(3, 3) = 1.0;
    Q(4, 4) = 10.0;
    Q(5, 5) = 10.0;

    double N = 12;
    double dt = 0.15;

    while (ros::ok()) {
        clock_t t_start = clock();
        // Reference
        Grid time_grid(0., N * dt, N);
        VariablesGrid reference(6, time_grid);

        reference(0, 0) = 0.;
        reference(0, 1) = 1.;
        reference(0, 2) = M_PI / 6;
        reference(0, 3) = 5.;
        reference(0, 4) = 0.;
        reference(0, 5) = 0.;

        for (int i = 1; i < N; ++i) {
            reference(i, 0) =
                reference(i - 1, 0) + reference(i - 1, 3) * cos(reference(i - 1, 3)) * dt + double(rand()) / RAND_MAX;
            reference(i, 1) =
                reference(i - 1, 1) + reference(i - 1, 3) * sin(reference(i - 1, 3)) * dt + double(rand()) / RAND_MAX;
            reference(i, 2) = reference(i - 1, 2) + double(rand()) / RAND_MAX / 10;
            reference(i, 3) = reference(i - 1, 3) + double(rand()) / RAND_MAX;
            reference(i, 4) = 0.;
            reference(i, 5) = 0.;
        }
//        reference.print();

        OCP ocp(reference.getTimePoints());
        ocp.minimizeLSQ(Q, h, reference);

        ocp.subjectTo(f);
        ocp.subjectTo(-M_PI / 6 <= delta <= M_PI / 6);
        ocp.subjectTo(-1 <= a <= 1);
        ocp.subjectTo(-M_PI / 18 <= delta_rate <= M_PI / 18);
        ocp.subjectTo(AT_START, x == 0.);
        ocp.subjectTo(AT_START, y == 0.);
        ocp.subjectTo(AT_START, v == 1. + double(rand()) / RAND_MAX);
        ocp.subjectTo(AT_START, phi == 0.);

        OptimizationAlgorithm algorithm(ocp);
        algorithm.set(PRINTLEVEL, NONE);
        algorithm.set(PRINT_COPYRIGHT, NONE);
        algorithm.solve();

        VariablesGrid states, controls;

        algorithm.getDifferentialStates(states);
        algorithm.getControls(controls);
//        controls.print();
        time.emplace_back(int(1000 * (clock() - t_start) / CLOCKS_PER_SEC) + 1);
        ROS_INFO_STREAM("cost time: " << time.back());
        loop_rate.sleep();
    }
    double mean_time = std::accumulate(time.begin(), time.end(), 0) / time.size();
    std::cout<<time.size()<<std::endl;
    std::cout<<mean_time<<std::endl;
    return 0;
}