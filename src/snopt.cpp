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

std::vector<std::vector<double>> *refer = nullptr;

double dt = 0.1;
int N = 12;
int state = 5;
int control = 2;

int n_vars = state * N + control * (N - 1);
int n_constrains = 1 + state * N;
// State
int x_begin = 0;
int y_begin = x_begin + N;
int phi_begin = y_begin + N;
int v_begin = phi_begin + N;
int w_begin = v_begin + N;
// Control
int a_begin = w_begin + N;
int wd_begin = a_begin + N - 1;

void fgFunction(int *Status, int *n, double x[],
                int *needF, int *neF, double F[],
                int *needG, int *neG, double G[],
                char *cu, int *lencu,
                int iu[], int *leniu,
                double ru[], int *lenru) {

    F[0] = 0.;
    for (int i = 0; i < N; ++i) {
        F[0] += pow(x[x_begin + i] - refer->at(i)[0], 2);
        F[0] += pow(x[y_begin + i] - refer->at(i)[1], 2);
        F[0] += pow(x[phi_begin + i] - refer->at(i)[2], 2);
        F[0] += pow(x[v_begin + i] - refer->at(i)[3], 2);
    }
    for (int i = 0; i < N - 1; ++i) {
        F[0] += pow(x[a_begin + i], 2);
        F[0] += pow(x[wd_begin + i], 2);
    }

    F[1 + x_begin] = x[x_begin];
    F[1 + y_begin] = x[y_begin];
    F[1 + phi_begin] = x[phi_begin];
    F[1 + v_begin] = x[v_begin];
    F[1 + w_begin] = x[w_begin];

    for (int j = 0; j < N - 1; ++j) {

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

//    std::vector<std::vector<double>> A1(N, std::vector<double>(N, 0));
//    for (int i = 0; i < N; ++i) {
//        A1[i][i] = 1;
//    }
//    for (int i = 0; i < N - 1; ++i) {
//        A1[i + 1][i] = -1;
//    }
//
//    std::vector<std::vector<double>> A2(N - 1, std::vector<double>(N - 1, 0));
//    for (int i = 0; i < N - 1; ++i) {
//        A2[i][i] = -dt;
//    }
//
//    for (int i = 0; i < A2.size(); ++i) {
//        for (int j = 0; j < A2[0].size(); ++j) {
//            std::cout << A2[i][j] << "  ";
//        }
//        std::cout << std::endl;
//    }

//    int lenA = state + 2 * (N - 1) + 2 * (N - 1) + 4 * (N - 1) + 3 * (N - 1) + 3 * (N - 1);
//    int *iAfun = new int[lenA]{0};
//    int *jAvar = new int[lenA]{0};
//    double *A = new double[lenA]{0};
//    int neA = lenA;
//
//    for (int i = 0; i < state; ++i) {
//        iAfun[i] = i*N+1;
//        jAvar[i] = i*N;
//        A[i] = 1;
//    }
//
//    for (int k = 0; k < state * (N-1); ++k) {
//        iAfun[state + 2 * k + 1] = k + 2;
//        jAvar[state + 2 * k + 1] = k;
//        A[state + 2 * k + 1] = -1;
//
//        iAfun[state + 2 * k] = k + 2;
//        jAvar[state + 2 * k] = k + 1;
//        A[state + 2 * k] = 1;
//    }
//
//    for (int k = 0; k < 3*N; ++k) {
//        iAfun[2*state * N - 1+k] = 25 + k;
//        jAvar[2*state * N - 1+k] = 47 + k;
//        A[2*state * N - 1+k] = -dt;
//    }
//
//    for (int l = 0; l < lenA; ++l) {
//        std::cout <<l<<" : " << "(" << *(iAfun + l) << ", " << *(jAvar + l) << ") : " << *(A + l) << std::endl;
//    }

    for (int K = 0; K < 100; ++K) {

        clock_t t_start = clock();
        double x_init = 0.;
        double y_init = 0.;
        double phi_init = -M_PI / 6;// + (double(rand()) / RAND_MAX - 0.5) / 10;
        double v_init = 1.;
        double w_init = 0.;

        refer = new std::vector<std::vector<double>>(N, std::vector<double>(4));
        refer->at(0)[0] = 0.;
        refer->at(0)[1] = 1.;
        refer->at(0)[2] = -M_PI / 6;
        refer->at(0)[3] = 5.;

        for (int i = 1; i < N; ++i) {
            refer->at(i)[0] = refer->at(i - 1)[0] + refer->at(i - 1)[3] * cos(refer->at(i - 1)[2]) * dt;
            refer->at(i)[1] = refer->at(i - 1)[1] + refer->at(i - 1)[3] * sin(refer->at(i - 1)[2]) * dt;
            refer->at(i)[2] = refer->at(i - 1)[2];// + (double(rand()) / RAND_MAX - 0.5) / 10;
            refer->at(i)[3] = refer->at(i - 1)[3];// + double(rand()) / RAND_MAX - 0.5;
        }

        // Allocate and initialize;
        int nS = 0, nInf;
        double sInf;

        double *x = new double[n_vars];
        double *xlow = new double[n_vars];
        double *xupp = new double[n_vars];
        double *xmul = new double[n_vars];
        int *xstate = new int[n_vars];

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

        double *F = new double[n_constrains];
        double *Flow = new double[n_constrains];
        double *Fupp = new double[n_constrains];
        double *Fmul = new double[n_constrains];
        int *Fstate = new int[n_constrains];

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
            no_one++;

            std::cout << "solution: " << solution << std::endl;
            std::cout << "x = " << x_init << std::endl;
            std::cout << "y = " << y_init << std::endl;
            std::cout << "phi = " << phi_init << std::endl;
            std::cout << "v = " << v_init << std::endl;
            std::cout << "w = " << w_init << std::endl;

            std::cout << "refer_point" << " = [";
            for (int i = 0; i < N; ++i) {
                std::cout << "(" << refer->at(i)[0] << ", " << refer->at(i)[1] << "), ";
            }
            std::cout << "\b\b]" << std::endl;

            std::cout << "solver_point" << " = [";
            for (int i = 0; i < N; ++i) {
                std::cout << "(" << x[x_begin + i] << ", " << x[y_begin + i] << "), ";
            }
            std::cout << "\b\b]" << std::endl;

            std::cout << "control_list" << " = [";
            for (int i = 0; i < N - 1; ++i) {
                std::cout << "(" << x[a_begin + i] << ", " << x[wd_begin + i] << "), ";
            }
            std::cout << "\b\b]" << std::endl;
        }

        delete refer;
        refer = nullptr;

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

        time.emplace_back(int(1000 * (clock() - t_start) / CLOCKS_PER_SEC) + 1);
    }
    std::cout << "no_one: " << no_one << std::endl;
    double mean_time = std::accumulate(time.begin(), time.end(), 0) / time.size();
    std::cout << time.size() << std::endl;
    std::cout << mean_time << std::endl;

    return 0;
}
