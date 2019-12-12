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
    int fail = 0;

    std::vector<std::vector<double>> *refer = nullptr;

    ros::init(argc, argv, "acado_node");
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

    while (ros::ok()) {
        clock_t t_start = clock();

        double x_init, y_init, phi_init, v_init, w_init;
        generateInitState(x_init, y_init, phi_init, v_init, w_init, random_state);
        generateReferPoint(refer, random_state);

        Grid time_grid(0., (step_N - 1) * dt, step_N);
        VariablesGrid reference(6, time_grid);

        for (int i = 0; i < step_N; ++i) {
            for (int j = 0; j < 4; ++j)
                reference(i, j) = refer->at(i)[j];
            reference(i, 4) = 0.;
            reference(i, 5) = 0.;
        }

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
        returnValueType result = algorithm.solve().getType();

        VariablesGrid states, controls;
        algorithm.getDifferentialStates(states);
        algorithm.getControls(controls);

        if (result != 0) {
            fail++;
            printInitState(x_init, y_init, phi_init, v_init, w_init);
            print("refer_point", reference, step_N);
            print("solver_point", states, step_N);
            print("control_list", controls, step_N);
        }

        time.emplace_back(int(1000 * (clock() - t_start) / CLOCKS_PER_SEC) + 1);
        std::cout << "cost time: " << time.back() << std::endl;
        loop_rate.sleep();
    }

    delete refer;

    std::cout << "fail: " << fail << std::endl;
    printTime(time);
    return 0;
}