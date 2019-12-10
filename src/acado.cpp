//
// Created by yao on 19-12-6.
//

#include "ros/ros.h"
#include "acado_toolkit.hpp"
#include "params.h"
USING_NAMESPACE_ACADO

template<typename T>
void print(std::string name, T vari, int N) {
    std::cout << name << " = [";
    for (int j = 0; j < N; ++j) {
        std::cout << "(" << vari(j, 0) << ", " << vari(j, 1) << "), ";
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
    Q(0, 0) = weight_x;
    Q(1, 1) = weight_y;
    Q(2, 2) = weight_phi;
    Q(3, 3) = weight_v;
    Q(4, 4) = weight_a;
    Q(5, 5) = weight_wd;

    int N = step_N;

    while (ros::ok()) {
        clock_t t_start = clock();
        // Reference
        Grid time_grid(0., (N - 1) * dt, N);
        VariablesGrid reference(6, time_grid);

        reference(0, 0) = 0.;
        reference(0, 1) = 1.;
        reference(0, 2) = -M_PI / 6;
        reference(0, 3) = 5.;
        reference(0, 4) = 0.;
        reference(0, 5) = 0.;

        for (int i = 1; i < N; ++i) {
            reference(i, 0) = reference(i - 1, 0) + reference(i - 1, 3) * cos(reference(i - 1, 2)) * dt;
            reference(i, 1) = reference(i - 1, 1) + reference(i - 1, 3) * sin(reference(i - 1, 2)) * dt;
            reference(i, 2) = reference(i - 1, 2);// + (double(rand()) / RAND_MAX - 0.5) / 10;
            reference(i, 3) = reference(i - 1, 3);// + double(rand()) / RAND_MAX - 0.5;
            reference(i, 4) = 0.;
            reference(i, 5) = 0.;
        }
        print("refer_point", reference, N);

        OCP ocp(reference.getTimePoints());
        ocp.minimizeLSQ(Q, h, reference);

        ocp.subjectTo(f);
        ocp.subjectTo(0. <= v <= 10.);
        ocp.subjectTo(-10. <= w <= 10.);
        ocp.subjectTo(-2. <= a <= 2.);
        ocp.subjectTo(-2. <= w_dot <= 2.);
        ocp.subjectTo(-0.4 * g <= v * w <= 0.4 * g);
        ocp.subjectTo(AT_START, x == x_init);
        ocp.subjectTo(AT_START, y == y_init);
        ocp.subjectTo(AT_START, phi == phi_init);
        ocp.subjectTo(AT_START, v == v_init);
        ocp.subjectTo(AT_START, w == w_init);

        OptimizationAlgorithm algorithm(ocp);
        algorithm.set(PRINTLEVEL, NONE);
        algorithm.set(PRINT_COPYRIGHT, NONE);
        algorithm.solve();

        VariablesGrid states, controls;

        algorithm.getDifferentialStates(states);
        algorithm.getControls(controls);

        print("solver_point", states, N);
        print("control_list", controls, N);

        time.emplace_back(int(1000 * (clock() - t_start) / CLOCKS_PER_SEC) + 1);
        ROS_INFO_STREAM("cost time: " << time.back());
        loop_rate.sleep();
    }
    double mean_time = std::accumulate(time.begin(), time.end(), 0) / time.size();
    std::cout << time.size() << std::endl;
    std::cout << mean_time << std::endl;
    return 0;
}