//
// Created by yao on 19-12-9.
//
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include "ros/ros.h"
#include "snopt/snopt.h"
#include "snopt/snoptProblem.hpp"
#include "params.h"

std::vector<std::vector<double>> *refer = nullptr;

void fFunction(int *Status, int *n, double x[],
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

void fgFunction(int *Status, int *n, double x[],
                int *needF, int *neF, double F[],
                int *needG, int *neG, double G[],
                char *cu, int *lencu,
                int iu[], int *leniu,
                double ru[], int *lenru) {
    if (*needF > 0) {
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

        for (int i = 0; i < step_N - 1; ++i) {

            double x1 = x[x_begin + i + 1];
            double y1 = x[y_begin + i + 1];
            double phi1 = x[phi_begin + i + 1];
            double v1 = x[v_begin + i + 1];
            double w1 = x[w_begin + i + 1];

            double x0 = x[x_begin + i];
            double y0 = x[y_begin + i];
            double phi0 = x[phi_begin + i];
            double v0 = x[v_begin + i];
            double w0 = x[w_begin + i];

            double a0 = x[a_begin + i];
            double wd0 = x[wd_begin + i];

            F[2 + x_begin + i] = x1 - (x0 + (v0 * dt + a0 * dt * dt / 2) * cos(phi0));
            F[2 + y_begin + i] = y1 - (y0 + (v0 * dt + a0 * dt * dt / 2) * sin(phi0));
            F[2 + phi_begin + i] = phi1 - (phi0 + (w0 * dt + wd0 * dt * dt / 2));
            F[2 + v_begin + i] = v1 - (v0 + a0 * dt);
            F[2 + w_begin + i] = w1 - (w0 + wd0 * dt);
        }
    }

    if (*needG > 0) {
        for (int i = 0; i < step_N; ++i) {
            G[i] = 2 * weight_x * (x[x_begin + i] - refer->at(i)[0]);
            G[12 + i] = 2 * weight_y * (x[y_begin + i] - refer->at(i)[1]);
        } // G[23], G.size() = 24

        for (int i = 0; i < step_N; ++i) {
            G[24 + 3 * i] = 2 * weight_phi * (x[phi_begin + i] - refer->at(i)[2]);
        }
        for (int i = 0; i < step_N - 1; ++i) {
            G[24 + 3 * i + 1] = (x[v_begin + i] * dt + x[a_begin + i] * dt * dt / 2) * sin(x[phi_begin + i]);
            G[24 + 3 * i + 2] = -(x[v_begin + i] * dt + x[a_begin + i] * dt * dt / 2) * cos(x[phi_begin + i]);
        } // G[57], G.size() = 58

        for (int i = 0; i < step_N; ++i) {
            G[58 + 3 * i] = 2 * weight_v * (x[v_begin + i] - refer->at(i)[3]);
        }
        for (int i = 0; i < step_N - 1; ++i) {
            G[58 + 3 * i + 1] = -dt * cos(x[phi_begin + i]);
            G[58 + 3 * i + 2] = -dt * sin(x[phi_begin + i]);
        } // G[91], G.size() = 92

        for (int i = 0; i < step_N - 1; ++i) {
            G[92 + 3 * i] = 2 * weight_a * x[a_begin + i];
            G[92 + 3 * i + 1] = -dt * dt / 2 * cos(x[phi_begin + i]);
            G[92 + 3 * i + 2] = -dt * dt / 2 * sin(x[phi_begin + i]);
        } // G[124], G.size() = 125

        for (int i = 0; i < step_N - 1; ++i) {
            G[125 + i] = 2 * weight_wd * x[wd_begin + i];
        } // G[135], G.size() = 136

        for (int i = 0; i < step_N - 1; ++i) {
            G[136 + 2 * i] = -1;
            G[136 + 2 * i + 1] = 1;
        } // G[157], G.size() = 158

        for (int i = 0; i < step_N - 1; ++i) {
            G[158 + 2 * i] = -1;
            G[158 + 2 * i + 1] = 1;
        } // G[179], G.size() = 180
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

    int lenA = 200;
    int *iAfun = new int[lenA]{0};
    int *jAvar = new int[lenA]{0};
    double *A = new double[lenA]{0};

    int lenG = 200;
    int *iGfun = new int[lenG]{0};
    int *jGvar = new int[lenG]{0};

    int neA, neG; // neA and neG must be defined when providing dervatives

    bool jacobi = false;
    bool gradient = true;

    ros::init(argc, argv, "snopt_node");
    ros::NodeHandle nh;
    ros::Rate loop_rate(10);

    while (ros::ok()) {
        clock_t t_start = clock();

        double x_init, y_init, phi_init, v_init, w_init;
        generateInitState(x_init, y_init, phi_init, v_init, w_init, random_state);
        generateReferPoint(refer, random_state);

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

        int solution = 0;
        snoptProblemA ToyProb;
        ToyProb.initialize("", 0);
        // snopta will compute the Jacobian by finite-differences.
        // The user has the option of calling  snJac  to define the
        // coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).

        if (!gradient) {
            ToyProb.setIntParameter("Derivative option", 0);
            ToyProb.setIntParameter("Verify level ", 3);

            // Solve the problem.
            // snJac is called implicitly in this case to compute the Jacobian.
            solution = ToyProb.solve(Cold, n_constrains, n_vars,
                                     ObjAdd, ObjRow, fFunction,
                                     xlow, xupp, Flow, Fupp,
                                     x, xstate, xmul,
                                     F, Fstate, Fmul,
                                     nS, nInf, sInf);
        } else {
            if (!jacobi) {
                ToyProb.computeJac(n_constrains, n_vars, fgFunction,
                                   x, xlow, xupp,
                                   iAfun, jAvar, A, neA,
                                   iGfun, jGvar, neG);
                for (int j = 0; j < lenA; ++j) {
                    if (*(iAfun + j) > 0) {
                        --*(iAfun + j);
                        --*(jAvar + j);
                    }
                }
                for (int j = 0; j < lenG; ++j) {
                    if (*(iGfun + j) > 0) {
                        --*(iGfun + j);
                        --*(jGvar + j);
                    }
                }
                if (neA != 115 || neG != 180)
                    jacobi = false;
                else
                    jacobi = true;
            }

            ToyProb.setIntParameter("Derivative option", 1);
            ToyProb.setIntParameter("Major Iteration limit", 250);
            ToyProb.setIntParameter("Verify level ", 3);
            solution = ToyProb.solve(Cold, n_constrains, n_vars, ObjAdd,
                                     ObjRow, fgFunction,
                                     iAfun, jAvar, A, neA, iGfun, jGvar, neG,
                                     xlow, xupp, Flow, Fupp,
                                     x, xstate, xmul,
                                     F, Fstate, Fmul,
                                     nS, nInf, sInf);
        }

        if (solution != 1) {
            std::cout << "solution: " << solution << std::endl;
            no_one++;
            printInitState(x_init, y_init, phi_init, v_init, w_init);
            printReferPoint(refer);
            printSolutionResult(x, "solver_point", x_begin, y_begin, step_N);
            printSolutionResult(x, "control_list", a_begin, wd_begin, step_N - 1);
        }

        time.emplace_back(int(1000 * (clock() - t_start) / CLOCKS_PER_SEC) + 1);
        std::cout << "cost time: " << time.back() << std::endl;
        loop_rate.sleep();
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

    delete[]iAfun;
    delete[]jAvar;
    delete[]A;
    delete[]iGfun;
    delete[]jGvar;

    std::cout << "fail: " << no_one << std::endl;
    printTime(time);
    return 0;
}
