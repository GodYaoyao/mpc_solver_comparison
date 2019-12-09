//
// Created by yao on 19-12-9.
//
#include <iostream>
#include <cmath>
#include <vector>
#include "snopt/snopt.h"
#include "snopt/snoptProblem.hpp"
using namespace std;

//bool snoptSolve(const State &start, const State &end, double *param,
//                const double *lower_bounds,
//                const double *upper_bounds, ParamKnot *solution) {
//    snoptProblemA nlp;
//    // Allocate and initialize;
//    int neF = 4;
//    int n = 3;
//
//    int nS = 0, nInf;
//    double sInf;
//
//    double *x = new double[n];
//    double *xlow = new double[n];
//    double *xupp = new double[n];
//    double *xmul = new double[n];
//    int *xstate = new int[n];
//
//    double *F = new double[neF];
//    double *Flow = new double[neF];
//    double *Fupp = new double[neF];
//    double *Fmul = new double[neF];
//    int *Fstate = new int[neF];
//
//    int ObjRow = 0;
//    double ObjAdd = 0;
//
//    int Cold = 0, Basis = 1, Warm = 2;
//
//    int lenG = 12;
//    int *iGfun = new int[lenG];
//    int *jGvar = new int[lenG];
//    int neG = lenG;
//    // Set the upper and lower bounds.
//    for (int i = 0; i < n; i++) {
//        xlow[i] = lower_bounds[i];
//        xupp[i] = upper_bounds[i];
//    }
//
//    x[0] = param[0];
//    x[1] = param[1];
//    x[2] = param[2];
//    xstate[0] = 0;
//    xstate[1] = 0;
//    xstate[2] = 0;
//
//    Flow[0] = -1e20;
//    Flow[1] = end.x;
//    Flow[2] = end.y;
//    Flow[3] = end.z;
//    Fupp[0] = 1e20;
//    Fupp[1] = end.x;
//    Fupp[2] = end.y;
//    Fupp[3] = end.z;
//    Fmul[0] = 0;
//    Fmul[1] = 0;
//    Fmul[2] = 0;
//    Fmul[3] = 0;
//
//    // snopta will compute the Jacobian by finite-differences.
//    // snJac will be called  to define the
//    // coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
//    // Provide the elements of the Jacobian explicitly.
//    iGfun[0] = 0;
//    jGvar[0] = 0;
//
//    iGfun[1] = 0;
//    jGvar[1] = 1;
//
//    iGfun[2] = 0;
//    jGvar[2] = 2;
//
//    iGfun[3] = 1;
//    jGvar[3] = 0;
//
//    iGfun[4] = 1;
//    jGvar[4] = 1;
//
//    iGfun[5] = 1;
//    jGvar[5] = 2;
//
//    iGfun[6] = 2;
//    jGvar[6] = 0;
//
//    iGfun[7] = 2;
//    jGvar[7] = 1;
//
//    iGfun[8] = 2;
//    jGvar[8] = 2;
//
//    iGfun[9] = 3;
//    jGvar[9] = 0;
//
//    iGfun[10] = 3;
//    jGvar[10] = 1;
//
//    iGfun[11] = 3;
//    jGvar[11] = 2;
//
//    int lenA = 1;
//    int *iAfun = new int[lenA];
//    int *jAvar = new int[lenA];
//    double *A = new double[lenA];
//    int neA = 0;
//
//    // Load the data for ToyProb ...
//    nlp.initialize("", 0);  // no print file; summary on
//    nlp.setIntParameter("Derivative option", 1);
//    nlp.setIntParameter("Verify level ", 3);
//    nlp.setProbName("NLP");
//    nlp.setIntParameter("Major Iteration limit", 100);
//    nlp.setIntParameter("Verify level ", 3);
//    nlp.solve(Cold, neF, n, ObjAdd, ObjRow, usrfun, iAfun, jAvar, A, neA, iGfun,
//              jGvar, neG, xlow, xupp, Flow, Fupp, x, xstate, xmul, F, Fstate,
//              Fmul, nS, nInf, sInf);
//
//    solution->p0 = 0;
//    solution->p1 = x[0];
//    solution->p2 = x[1];
//    solution->p3 = 0;
//    solution->sf = x[2];
//
//    delete[] iAfun;
//    delete[] jAvar;
//    delete[] A;
//    delete[] iGfun;
//    delete[] jGvar;
//
//    delete[] x;
//    delete[] xlow;
//    delete[] xupp;
//    delete[] xmul;
//    delete[] xstate;
//
//    delete[] F;
//    delete[] Flow;
//    delete[] Fupp;
//    delete[] Fmul;
//    delete[] Fstate;
//
//    return true;
//}

void toyusrf_(int    *Status, int *n,    double x[],
              int    *needF,  int *neF,  double F[],
              int    *needG,  int *neG,  double G[],
              char      *cu,  int *lencu,
              int    iu[],    int *leniu,
              double ru[],    int *lenru )
{
    //==================================================================
    // Computes the nonlinear objective and constraint terms for the toy
    // problem featured in the SnoptA users guide.
    // neF = 3, n = 2.
    //
    //   Minimize     x(2)
    //
    //   subject to   x(1)**2      + 4 x(2)**2  <= 4,
    //               (x(1) - 2)**2 +   x(2)**2  <= 5,
    //                x(1) >= 0.
    //
    //==================================================================
    F[0] =  pow(x[1] - 1, 2);
    F[1] =  x[0]*x[0] + 4*x[1]*x[1];
    F[2] = (x[0] - 2)*(x[0] - 2) + x[1]*x[1];
}

void toyusrfg_( int    *Status, int *n,    double x[],
                int    *needF,  int *neF,  double F[],
                int    *needG,  int *neG,  double G[],
                char      *cu,  int *lencu,
                int    iu[],    int *leniu,
                double ru[],    int *lenru )
{
    //==================================================================
    // Computes the nonlinear objective and constraint terms for the toy
    // problem featured in the SnoptA users guide.
    // neF = 3, n = 2.
    //
    //   Minimize     x(2)
    //
    //   subject to   x(1)**2      + 4 x(2)**2  <= 4,
    //               (x(1) - 2)**2 +   x(2)**2  <= 5,
    //                x(1) >= 0.
    //
    // The triples (g(k),iGfun(k),jGvar(k)), k = 1:neG, define
    // the sparsity pattern and values of the nonlinear elements
    // of the Jacobian.
    //==================================================================

    if ( *needF > 0 ) {
        F[0] =  pow(x[1] - 1, 2); //  Objective row
        F[1] =  x[0]*x[0] + 4*x[1]*x[1];
        F[2] = (x[0] - 2)*(x[0] - 2) + x[1]*x[1];
    }


    if ( *needG > 0 ) {
        // iGfun[0] = 1
        // jGvar[0] = 0
        G[0] = 2*x[0];

        // iGfun[1] = 1
        // jGvar[1] = 1
        G[1] = 8*x[1];

        // iGfun[2] = 2
        // jGvar[2] = 0
        G[2] = 2*(x[0] - 2);

        // iGfun[3] = 2
        // jGvar[3] = 1
        G[3] = 2*x[1];
    }
}

void example(){
    snoptProblemA ToyProb;

    // Allocate and initialize;
    int n     =  2;
    int neF   =  3;

    int nS = 0, nInf;
    double sInf;

    double *x      = new double[n];
    double *xlow   = new double[n];
    double *xupp   = new double[n];
    double *xmul   = new double[n];
    int    *xstate = new    int[n];

    double *F      = new double[neF];
    double *Flow   = new double[neF];
    double *Fupp   = new double[neF];
    double *Fmul   = new double[neF];
    int    *Fstate = new int[neF];

    int    ObjRow  = 0;
    double ObjAdd  = 0;

    int Cold = 0, Basis = 1, Warm = 2;

    // Set the upper and lower bounds.
    xlow[0]   =  0.0;  xlow[1]   = -1e20;
    xupp[0]   = 1e20;  xupp[1]   =  1e20;
    xstate[0] =    0;  xstate[1] =  0;

    Flow[0] = -1e20; Flow[1] = -1e20; Flow[2] = -1e20;
    Fupp[0] =  1e20; Fupp[1] =   4.0; Fupp[2] =  5.0;
    x[0]    = 1.0;
    x[1]    = 1.0;

    // Load the data for ToyProb ...
    ToyProb.initialize("", 0);
//    ToyProb.setProbName   ("Toy0");
//    ToyProb.setPrintFile  ( "Toy0.out" );

    // snopta will compute the Jacobian by finite-differences.
    // The user has the option of calling  snJac  to define the
    // coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
    ToyProb.setIntParameter( "Derivative option", 0 );
    ToyProb.setIntParameter( "Verify level ", 3 );

    // Solve the problem.
    // snJac is called implicitly in this case to compute the Jacobian.
    ToyProb.solve( Cold, neF, n,
                   ObjAdd, ObjRow, toyusrf_,
                   xlow, xupp, Flow, Fupp,
                   x, xstate, xmul,
                   F, Fstate, Fmul,
                   nS, nInf, sInf);

    cout<<x[0]<<endl;
    cout<<x[1]<<endl;

    cout<<xstate[0]<<endl;
    cout<<xstate[1]<<endl;

    printf("\nSolving toy1 problem using derivatives...\n");

    // Reset the variables and solve ...

    int lenA   = 6;
    int *iAfun = new int[lenA];
    int *jAvar = new int[lenA];
    double *A  = new double[lenA];

    int lenG   = 6;
    int *iGfun = new int[lenG];
    int *jGvar = new int[lenG];

    int neA, neG; // neA and neG must be defined when providing dervatives

    xstate[0] =   0;  xstate[1] = 0;
    Fmul[0]   =   0;  Fmul[0]   = 0; Fmul[0] =    0;
    x[0]      = 1.0;
    x[1]      = 1.0;


    // Provide the elements of the Jacobian explicitly.
    iGfun[0] = 1;
    jGvar[0] = 0;

    iGfun[1] = 1;
    jGvar[1] = 1;

    iGfun[2] = 2;
    jGvar[2] = 0;

    iGfun[3] = 2;
    jGvar[3] = 1;
    neG      = 4;

    iAfun[0] = 0;
    jAvar[0] = 1;
    A[0]     = 1.0;
    neA      = 1;

    ToyProb.initialize("", 0);
//    ToyProb.setProbName    ( "Toy1" );         // Give the problem a new name for Snopt.
//    ToyProb.setPrintFile   ( "Toy1.out" );
//    ToyProb.setSpecsFile   ( "sntoya.spc" );
    ToyProb.setIntParameter( "Derivative option", 1 );
    ToyProb.setIntParameter( "Major Iteration limit", 250 );
    ToyProb.setIntParameter( "Verify level ", 3 );
    ToyProb.solve          ( Cold, neF, n, ObjAdd,
                             ObjRow, toyusrfg_,
                             iAfun, jAvar, A, neA, iGfun, jGvar, neG,
                             xlow, xupp, Flow, Fupp,
                             x, xstate, xmul,
                             F, Fstate, Fmul,
                             nS, nInf, sInf);

    for (int i = 0; i < n; i++ ){
        cout << "x = " << x[i] << " xstate = " << xstate[i] << endl;
    }
    for (int i = 0; i < neF; i++ ){
        cout << "F = " << F[i] << " Fstate = " << Fstate[i] << endl;
    }

    ToyProb.solve          ( Warm, neF, n, ObjAdd,
                             ObjRow, toyusrfg_,
                             iAfun, jAvar, A, neA, iGfun, jGvar, neG,
                             xlow, xupp, Flow, Fupp,
                             x, xstate, xmul,
                             F, Fstate, Fmul,
                             nS, nInf, sInf);

    for (int i = 0; i < n; i++ ){
        cout << "x = " << x[i] << " xstate = " << xstate[i] << endl;
    }
    for (int i = 0; i < neF; i++ ){
        cout << "F = " << F[i] << " Fstate = " << Fstate[i] << endl;
    }


    delete []iAfun;  delete []jAvar;  delete []A;
    delete []iGfun;  delete []jGvar;

    delete []x;      delete []xlow;   delete []xupp;
    delete []xmul;   delete []xstate;

    delete []F;      delete []Flow;   delete []Fupp;
    delete []Fmul;   delete []Fstate;
}


int counter=0;
std::vector<std::vector<double>> *refer = nullptr;

double dt = 0.1;
int N = 12;
int state = 5;
int control = 2;

int n_vars =  state * N + control *(N-1);
int n_constrains =  1 + state * N;
// State
int x_begin = 0;
int y_begin = x_begin + N;
int phi_begin = y_begin + N;
int v_begin = phi_begin + N;
int w_begin = v_begin + N;
// Control
int a_begin = w_begin + N;
int wd_begin = a_begin + N - 1;

void fgFunction(int    *Status, int *n,    double x[],
                int    *needF,  int *neF,  double F[],
                int    *needG,  int *neG,  double G[],
                char      *cu,  int *lencu,
                int    iu[],    int *leniu,
                double ru[],    int *lenru )  {

    counter++;

    F[0] = 0.;
    for (int i = 0; i < N; ++i) {
        F[0] += 10 * pow(x[x_begin + i] - refer->at(i)[0], 2);
        F[0] += 10 * pow(x[y_begin + i] - refer->at(i)[1], 2);
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

        F[2 + x_begin + j] = x1 - (x0 + (v0 * dt + a0 * dt * dt/2) * cos(phi0));
        F[2 + y_begin + j] = y1 - (y0 + (v0 * dt + a0 * dt * dt/2) * sin(phi0));
        F[2 + phi_begin + j] = phi1 - (phi0 + (w0 * dt + wd0 * dt * dt/2));
        F[2 + v_begin + j] = v1 - (v0 + a0 * dt);
        F[2 + w_begin + j] = w1 - (w0 + wd0 * dt);
    }
}

void test(){

    clock_t t_start = clock();
    double x_init = 0.;
    double y_init = 0.;
    double phi_init = -M_PI / 6 + (double(rand()) / RAND_MAX - 0.5) / 10;
    double v_init = 1.;
    double w_init = 0.;

    refer = new std::vector<std::vector<double>>(N, std::vector<double>(4));
    refer->at(0)[0] = 0.;
    refer->at(0)[1] = 1.;
    refer->at(0)[2] = - M_PI / 6;
    refer->at(0)[3] = 5.;

    for (int i = 1; i < N; ++i) {
        refer->at(i)[0] = refer->at(i - 1)[0] + refer->at(i - 1)[3] * cos(refer->at(i - 1)[2]) * dt;
        refer->at(i)[1] = refer->at(i - 1)[1] + refer->at(i - 1)[3] * sin(refer->at(i - 1)[2]) * dt;
        refer->at(i)[2] = refer->at(i - 1)[2] + (double(rand()) / RAND_MAX - 0.5) / 10;
        refer->at(i)[3] = refer->at(i - 1)[3] + double(rand()) / RAND_MAX - 0.5;
    }

    // Allocate and initialize;
    int nS = 0, nInf;
    double sInf;

    double *x      = new double[n_vars];
    double *xlow   = new double[n_vars];
    double *xupp   = new double[n_vars];
    double *xmul   = new double[n_vars];
    int    *xstate = new    int[n_vars];

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

    double *F      = new double[n_constrains];
    double *Flow   = new double[n_constrains];
    double *Fupp   = new double[n_constrains];
    double *Fmul   = new double[n_constrains];
    int    *Fstate = new int[n_constrains];

    for (int j = 0; j < n_constrains; ++j) {
        Flow[j] = 0.;
        Fupp[j] = 0.;
    }

    Flow[0] = -1e20;
    Fupp[0] =  1e20;

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

    int    ObjRow  = 0;
    double ObjAdd  = 0;

    int Cold = 0, Basis = 1, Warm = 2;

    snoptProblemA ToyProb;
    ToyProb.initialize("", 0);
    // snopta will compute the Jacobian by finite-differences.
    // The user has the option of calling  snJac  to define the
    // coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
    ToyProb.setIntParameter( "Derivative option", 0 );
    ToyProb.setIntParameter( "Verify level ", 3 );

    // Solve the problem.
    // snJac is called implicitly in this case to compute the Jacobian.
    int solution = ToyProb.solve( Cold, n_constrains, n_vars,
                                  ObjAdd, ObjRow, fgFunction,
                                  xlow, xupp, Flow, Fupp,
                                  x, xstate, xmul,
                                  F, Fstate, Fmul,
                                  nS, nInf, sInf);
    cout<<solution<<endl;

    cout<<counter<<endl;

//    std::cout<< "refer_point" <<" = [";
//    for (int i = 0; i < N; ++i) {
//        std::cout<<"("<<refer->at(i)[0] << ", " << refer->at(i)[1] << "), ";
//    }
//    std::cout << "\b\b]" << std::endl;
//
//    std::cout<< "calculate_list" <<" = [";
//    for (int i = 0; i < N; ++i) {
//        std::cout<<"("<< x[x_begin+i] << ", " << x[y_begin+i] << "," << x[phi_begin+i] << "), ";
//    }
//    std::cout << "\b\b]" << std::endl;
//
//    std::cout<< "control_list" <<" = [";
//    for (int i = 0; i < N-1; ++i) {
//        std::cout<<"("<< x[a_begin+i] << ", " << x[wd_begin+i] << "), ";
//    }
//    std::cout << "\b\b]" << std::endl;

    delete refer;
    refer = nullptr;

    delete []x;      delete []xlow;   delete []xupp;
    delete []xmul;   delete []xstate;

    delete []F;      delete []Flow;   delete []Fupp;
    delete []Fmul;   delete []Fstate;

    cout<<int(1000 * (clock() - t_start) / CLOCKS_PER_SEC) + 1 << "ms" <<endl;
}

int main(int argc, char **argv) {

//    example();
    for (int i = 0; i < 10; ++i) {
        test();
    }

    return 0;
}