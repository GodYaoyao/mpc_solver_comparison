//
// Created by yao on 19-12-9.
//
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
//#include "ros/ros.h"
#include "snopt/snopt.h"
#include "snopt/snoptProblem.hpp"
#include "params.h"

std::vector<std::vector<double>> *refer = nullptr;

void fgFunction(int *Status, int *n, double x[],
                int *needF, int *neF, double F[],
                int *needG, int *neG, double G[],
                char *cu, int *lencu,
                int iu[], int *leniu,
                double ru[], int *lenru) {

    F[0] = 0.;
    for (int i = 0; i < step_N; ++i) {
        F[0] += weight_x * pow(x[x_begin + i] - refer->at(i)[0], 2);
        F[0] += weight_y * pow(x[y_begin + i] - refer->at(i)[1], 2);
        F[0] += weight_phi * pow(x[phi_begin + i] - refer->at(i)[2], 2);
        F[0] += weight_v * pow(x[v_begin + i] - refer->at(i)[3], 2);
    }
    for (int i = 0; i < step_N - 1; ++i) {
        F[0] += weight_a * pow(x[a_begin + i], 2);
        F[0] += weight_wd * pow(x[wd_begin + i], 2);
    }

    F[1 + x_begin] = x[x_begin];
    F[1 + y_begin] = x[y_begin];
    F[1 + phi_begin] = x[phi_begin];
    F[1 + v_begin] = x[v_begin];
    F[1 + w_begin] = x[w_begin];

    for (int j = 0; j < step_N - 1; ++j) {

        double x1 = x[x_begin + j + 1];
        double y1 = x[y_begin + j + 1];
        double phi1 = x[phi_begin + j + 1];
        double v1 = x[v_begin + j + 1];
        double w1 = x[w_begin + j + 1];

        double x0 = x[x_begin + j];
        double y0 = x[y_begin + j];
        double phi0 = x[phi_begin + j];
        double v0 = x[v_begin + j];
        double w0 = x[w_begin + j];

        double a0 = x[a_begin + j];
        double wd0 = x[wd_begin + j];

        F[2 + x_begin + j] = x1 - (x0 + (v0 * dt + a0 * dt * dt / 2) * cos(phi0));
        F[2 + y_begin + j] = y1 - (y0 + (v0 * dt + a0 * dt * dt / 2) * sin(phi0));
        F[2 + phi_begin + j] = phi1 - (phi0 + (w0 * dt + wd0 * dt * dt / 2));
        F[2 + v_begin + j] = v1 - (v0 + a0 * dt);
        F[2 + w_begin + j] = w1 - (w0 + wd0 * dt);
    }
}

int main(int argc, char **argv) {

    std::vector<int> time;
    int no_one = 0;

    double *x = new double[n_vars]{0};
    double *xlow = new double[n_vars]{0};
    double *xupp = new double[n_vars]{0};
    double *xmul = new double[n_vars]{0};
    int *xstate = new int[n_vars]{0};

    double *F = new double[n_constrains]{0};
    double *Flow = new double[n_constrains]{0};
    double *Fupp = new double[n_constrains]{0};
    double *Fmul = new double[n_constrains]{0};
    int *Fstate = new int[n_constrains]{0};

    for (int K = 0; K < 100; ++K) {

        clock_t t_start = clock();
        double x_init = 0.;
        double y_init = 0.;
        double phi_init = -M_PI / 6 + (double(rand()) / RAND_MAX - 0.5) / 10;
        double v_init = 1.;
        double w_init = 0.;

        refer = new std::vector<std::vector<double>>(step_N, std::vector<double>(4));
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

        // Allocate and initialize;
        int nS = 0, nInf;
        double sInf;

        // Set the upper and lower bounds.
        for (int i = 0; i < v_begin; ++i) {
            xlow[i] = -1e20;
            xupp[i] = 1e20;
            xstate[i] = 0;
        }
        for (int i = v_begin; i < w_begin; ++i) {
            xlow[i] = 0;
            xupp[i] = 10;
            xstate[i] = 0;
        }
        for (int i = w_begin; i < a_begin; ++i) {
            xlow[i] = -10.;
            xupp[i] = 10;
            xstate[i] = 0;
        }
        for (int i = a_begin; i < n_vars; ++i) {
            xlow[i] = -2.;
            xupp[i] = 2.;
            xstate[i] = 0;
        }
        x[x_begin] = x_init;
        x[y_begin] = y_init;
        x[phi_begin] = phi_init;
        x[v_begin] = v_init;
        x[w_begin] = w_init;
        x[a_begin] = 0.;
        x[wd_begin] = 0.;

        for (int j = 0; j < n_constrains; ++j) {
            Flow[j] = 0.;
            Fupp[j] = 0.;
        }

        Flow[0] = -1e20;
        Fupp[0] = 1e20;

        Flow[1 + x_begin] = x_init;
        Fupp[1 + x_begin] = x_init;

        Flow[1 + y_begin] = y_init;
        Fupp[1 + y_begin] = y_init;

        Flow[1 + phi_begin] = phi_init;
        Fupp[1 + phi_begin] = phi_init;

        Flow[1 + v_begin] = v_init;
        Fupp[1 + v_begin] = v_init;

        Flow[1 + w_begin] = w_init;
        Fupp[1 + w_begin] = w_init;

        int ObjRow = 0;
        double ObjAdd = 0;

        int Cold = 0, Basis = 1, Warm = 2;

        snoptProblemA ToyProb;
        ToyProb.initialize("", 0);
        // snopta will compute the Jacobian by finite-differences.
        // The user has the option of calling  snJac  to define the
        // coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
        ToyProb.setIntParameter("Derivative option", 0);
        ToyProb.setIntParameter("Verify level ", 3);

        // Solve the problem.
        // snJac is called implicitly in this case to compute the Jacobian.
        int solution = ToyProb.solve(Cold, n_constrains, n_vars,
                                     ObjAdd, ObjRow, fgFunction,
                                     xlow, xupp, Flow, Fupp,
                                     x, xstate, xmul,
                                     F, Fstate, Fmul,
                                     nS, nInf, sInf);

        if (solution != 1) {
            std::cout << "solution: " << solution << std::endl;
            no_one++;
            printInitState(x_init, y_init, phi_init, v_init, w_init);
            printReferPoint(refer);
            printSolutionResult(x, "solver_point", x_begin, y_begin, step_N);
            printSolutionResult(x, "control_list", a_begin, wd_begin, step_N - 1);
        }

        time.emplace_back(int(1000 * (clock() - t_start) / CLOCKS_PER_SEC) + 1);
    }

    delete refer;

    delete[]x;
    delete[]xlow;
    delete[]xupp;
    delete[]xmul;
    delete[]xstate;

    delete[]F;
    delete[]Flow;
    delete[]Fupp;
    delete[]Fmul;
    delete[]Fstate;

    std::cout << "fail: " << no_one << std::endl;
    printTime(time);
    return 0;
}
