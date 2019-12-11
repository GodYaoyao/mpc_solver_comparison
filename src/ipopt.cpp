//
// Created by yao on 19-12-9.
//
#include <iostream>
#include "ros/ros.h"
#include "params.h"
#include "cppad/cppad.hpp"
#include "cppad/ipopt/solve.hpp"
using CppAD::AD;

class FgEval {
public:
    FgEval(std::vector<std::vector<double>> *refer) : refer(refer) {}

public:
    std::vector<std::vector<double>> *refer;

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

    void operator()(ADvector &fg, const ADvector &vars) {
//        ROS_INFO_STREAM("fg: " << fg.size());
//        ROS_INFO_STREAM("vars: " << vars.size());

        fg[0] = 0.0;
        for (int i = 0; i < step_N; ++i) {
            fg[0] += weight_x * pow(vars[x_begin + i] - refer->at(i)[0], 2);
            fg[0] += weight_y * pow(vars[y_begin + i] - refer->at(i)[1], 2);
            fg[0] += weight_phi * pow(vars[phi_begin + i] - refer->at(i)[2], 2);
            fg[0] += weight_v * pow(vars[v_begin + i] - refer->at(i)[3], 2);
        }
        for (int i = 0; i < step_N - 1; ++i) {
            fg[0] += weight_a * pow(vars[a_begin + i], 2);
            fg[0] += weight_wd * pow(vars[wd_begin + i], 2);
        }

        // Initial constraints
        // We add 1 to each of the starting indices due to cost being located at
        // index 0 of `fg`. This bumps up the position of all the other values.
        fg[1 + x_begin] = vars[x_begin];
        fg[1 + y_begin] = vars[y_begin];
        fg[1 + phi_begin] = vars[phi_begin];
        fg[1 + v_begin] = vars[v_begin];
        fg[1 + w_begin] = vars[w_begin];

        // The rest of the constraints
        for (int i = 0; i < step_N - 1; i++) {
            // The state at time t+1 .
            AD<double> x1 = vars[x_begin + i + 1];
            AD<double> y1 = vars[y_begin + i + 1];
            AD<double> phi1 = vars[phi_begin + i + 1];
            AD<double> v1 = vars[v_begin + i + 1];
            AD<double> w1 = vars[w_begin + i + 1];
            // The state at time t.
            AD<double> x0 = vars[x_begin + i];
            AD<double> y0 = vars[y_begin + i];
            AD<double> phi0 = vars[phi_begin + i];
            AD<double> v0 = vars[v_begin + i];
            AD<double> w0 = vars[w_begin + i];
            // Only consider the actuation at time t.
            AD<double> a0 = vars[a_begin + i];
            AD<double> wd0 = vars[wd_begin + i];

            fg[2 + x_begin + i] = x1 - (x0 + (v0 * dt + a0 * dt * dt / 2) * CppAD::cos(phi0));
            fg[2 + y_begin + i] = y1 - (y0 + (v0 * dt + a0 * dt * dt / 2) * CppAD::sin(phi0));
            fg[2 + phi_begin + i] = phi1 - (phi0 + w0 * dt + wd0 * dt * dt / 2);
            fg[2 + v_begin + i] = v1 - (v0 + a0 * dt);
            fg[2 + w_begin + i] = w1 - (w0 + wd0 * dt);
        }
    }
};

int main(int argc, char **argv) {

    std::vector<int> time;
    int fail = 0;

    ros::init(argc, argv, "ipopt_node");
    ros::NodeHandle nh;
    ros::Rate loop_rate(10);

    while (ros::ok()) {
        clock_t t_start = clock();
        double x_init = 0.;
        double y_init = 0.;
        double phi_init = -M_PI / 6 + (double(rand()) / RAND_MAX - 0.5) / 10;
        double v_init = 1.;
        double w_init = 0.;

        ROS_INFO_STREAM("phi_init" << phi_init);

        std::vector<std::vector<double>> *refer = new std::vector<std::vector<double>>(step_N, std::vector<double>(4));
        refer->at(0)[0] = 0.;
        refer->at(0)[1] = 1.;
        refer->at(0)[2] = -M_PI / 6;
        refer->at(0)[3] = 5.;
        for (int i = 1; i < step_N; ++i) {
            refer->at(i)[0] = refer->at(i - 1)[0] + refer->at(i - 1)[3] * cos(refer->at(i - 1)[2]) * dt;
            refer->at(i)[1] = refer->at(i - 1)[1] + refer->at(i - 1)[3] * sin(refer->at(i - 1)[2]) * dt;
            refer->at(i)[2] = refer->at(i - 1)[2] + (double(rand()) / RAND_MAX - 0.5) / 10;
            refer->at(i)[3] = refer->at(i - 1)[3] + double(rand()) / RAND_MAX - 0.5;
        }

        bool ok = true;
        typedef CPPAD_TESTVECTOR(double) Dvector;

        Dvector vars(n_vars);
        for (unsigned int i = 0; i < n_vars; i++) {
            vars[i] = 0;
        }
        vars[x_begin] = x_init;
        vars[y_begin] = y_init;
        vars[phi_begin] = phi_init;
        vars[v_begin] = v_init;
        vars[w_begin] = w_init;
        vars[a_begin] = 0.;
        vars[wd_begin] = 0.;

        Dvector vars_lowerbound(n_vars);
        Dvector vars_upperbound(n_vars);
        for (int i = 0; i < v_begin; i++) {
            vars_lowerbound[i] = -1.0e19;
            vars_upperbound[i] = 1.0e19;
        }
        for (int i = v_begin; i < w_begin; i++) {
            vars_lowerbound[i] = 0.;
            vars_upperbound[i] = 10.;
        }
        for (int i = w_begin; i < a_begin; i++) {
            vars_lowerbound[i] = -10.;
            vars_upperbound[i] = 10.;
        }
        for (int i = a_begin; i < wd_begin; i++) {
            vars_lowerbound[i] = -2.;
            vars_upperbound[i] = 2.;
        }
        for (int i = wd_begin; i < n_vars; i++) {
            vars_lowerbound[i] = -2.;
            vars_upperbound[i] = 2.;
        }

        Dvector constraints_lowerbound(n_constrains);
        Dvector constraints_upperbound(n_constrains);
        for (unsigned int i = 0; i < n_constrains; i++) {
            constraints_lowerbound[i] = 0.0;
            constraints_upperbound[i] = 0.0;
        }
        constraints_lowerbound[x_begin] = x_init;
        constraints_upperbound[x_begin] = x_init;

        constraints_lowerbound[y_begin] = y_init;
        constraints_upperbound[y_begin] = y_init;

        constraints_lowerbound[phi_begin] = phi_init;
        constraints_upperbound[phi_begin] = phi_init;

        constraints_lowerbound[v_begin] = v_init;
        constraints_upperbound[v_begin] = v_init;

        constraints_lowerbound[w_begin] = w_init;
        constraints_upperbound[w_begin] = w_init;

        // options for IPOPT solver
        std::string options;
        // Uncomment this if you'd like more print information
        options += "Integer print_level  0\n";
        // NOTE: Setting sparse to true allows the solver to take advantage
        // of sparse routines, this makes the computation MUCH FASTER. If you
        // can uncomment 1 of these and see if it makes a difference or not but
        // if you uncomment both the computation time should go up in orders of
        // magnitude.
        options += "Sparse  true        forward\n";
        options += "Sparse  true        reverse\n";
        // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
        // Change this as you see fit.
        options += "Numeric max_cpu_time          0.5\n";

        // place to return solution
        CppAD::ipopt::solve_result<Dvector> solution;
        // object that computes objective and constraints
        FgEval fg_eval(refer);
        // solve the problem
        CppAD::ipopt::solve<Dvector, FgEval>(options, vars,
                                             vars_lowerbound, vars_upperbound,
                                             constraints_lowerbound, constraints_upperbound,
                                             fg_eval, solution);

        // Check some of the solution values
        ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

        if (!ok) {
            fail++;
            printInitState(x_init, y_init, phi_init, v_init, w_init);
            printReferPoint(refer);
            printSolutionResult(solution.x, "solver_point", x_begin, y_begin, step_N);
            printSolutionResult(solution.x, "control_list", a_begin, wd_begin, step_N - 1);
        }

        delete refer;
        time.emplace_back(int(1000 * (clock() - t_start) / CLOCKS_PER_SEC) + 1);
        ROS_INFO_STREAM("cost time: " << time.back());
        loop_rate.sleep();
    }
    std::cout << "fail: " << fail << std::endl;
    printTime(time);
    return 0;
}