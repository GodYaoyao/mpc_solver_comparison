//
// Created by yao on 19-12-10.
//
#include <iostream>
#include <cmath>
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

void toyusrf_(int *Status, int *n, double x[],
              int *needF, int *neF, double F[],
              int *needG, int *neG, double G[],
              char *cu, int *lencu,
              int iu[], int *leniu,
              double ru[], int *lenru) {
    //==================================================================
    // Computes the nonlinear objective and constraint terms for the toy
    // problem featured in the SnoptA users guide.
    // neF = 3, n = 2.
    //
    //   Minimize     pow(x[1] - 1, 2)
    //
    //   subject to   x(0)**2      + 4 x(1)**2  <= 4,
    //               (x(0) - 2)**2 +   x(1)**2  <= 5,
    //                x(0) >= 0.
    //
    //==================================================================
//    cout << *n << endl;
//    cout << *neF << endl;
    F[0] = pow(x[1] - 1, 2);
    F[1] = x[0] * x[0] + 4 * x[1] * x[1];
    F[2] = (x[0] - 2) * (x[0] - 2) + x[1] * x[1];
}

void toyusrfg_(int *Status, int *n, double x[],
               int *needF, int *neF, double F[],
               int *needG, int *neG, double G[],
               char *cu, int *lencu,
               int iu[], int *leniu,
               double ru[], int *lenru) {
    //==================================================================
    // Computes the nonlinear objective and constraint terms for the toy
    // problem featured in the SnoptA users guide.
    // neF = 3, n = 2.
    //
    //   Minimize     pow(x[1] - 1, 2)
    //
    //   subject to   x(0)**2      + 4 x(1)**2  <= 4,
    //               (x(0) - 2)**2 +   x(0) * x(1)  <= 5,
    //                x(0) >= 0.
    //
    // The triples (g(k),iGfun(k),jGvar(k)), k = 1:neG, define
    // the sparsity pattern and values of the nonlinear elements
    // of the Jacobian.
    //==================================================================

    if (*needF > 0) {
        F[0] = pow(x[1] - 1, 2); //  Objective row
        F[1] = x[0] * x[0] + 4 * x[1] * x[1];
        F[2] = (x[0] - 2) * (x[0] - 2) + x[0] * x[1];
    }

    if (*needG > 0) {

        G[0] = 2 * (x[1] - 1);

        G[1] = 2 * x[0];

        G[2] = 8 * x[1];

        G[3] = 2 * (x[0] - 2) + x[1];

        G[4] = x[0];
    }
}

int main(int argc, char **argv) {

    snoptProblemA ToyProb;

    // Allocate and initialize;
    int n = 2;
    int neF = 3;

    int nS = 0, nInf;
    double sInf;

    double *x = new double[n];
    double *xlow = new double[n];
    double *xupp = new double[n];
    double *xmul = new double[n];
    int *xstate = new int[n];

    double *F = new double[neF];
    double *Flow = new double[neF];
    double *Fupp = new double[neF];
    double *Fmul = new double[neF];
    int *Fstate = new int[neF];

    int ObjRow = 0;
    double ObjAdd = 0;

    int Cold = 0, Basis = 1, Warm = 2;

    // Set the upper and lower bounds.
    xlow[0] = 0.0;
    xlow[1] = -1e20;
    xupp[0] = 1e20;
    xupp[1] = 1e20;
    xstate[0] = 0;
    xstate[1] = 0;

    Flow[0] = -1e20;
    Flow[1] = -1e20;
    Flow[2] = -1e20;
    Fupp[0] = 1e20;
    Fupp[1] = 4.0;
    Fupp[2] = 5.0;
    x[0] = 1.0;
    x[1] = 1.0;

    // Load the data for ToyProb ...
    ToyProb.initialize("", 0);
//    ToyProb.setProbName   ("Toy0");
//    ToyProb.setPrintFile  ( "Toy0.out" );

    // snopta will compute the Jacobian by finite-differences.
    // The user has the option of calling  snJac  to define the
    // coordinate arrays (iAfun,jAvar,A) and (iGfun, jGvar).
    ToyProb.setIntParameter("Derivative option", 0);
    ToyProb.setIntParameter("Verify level ", 3);

    // Solve the problem.
    // snJac is called implicitly in this case to compute the Jacobian.
    ToyProb.solve(Cold, neF, n,
                  ObjAdd, ObjRow, toyusrf_,
                  xlow, xupp, Flow, Fupp,
                  x, xstate, xmul,
                  F, Fstate, Fmul,
                  nS, nInf, sInf);

    for (int i = 0; i < n; i++) {
        cout << "x = " << x[i] << " xstate = " << xstate[i] << endl;
    }
    for (int i = 0; i < neF; i++) {
        cout << "F = " << F[i] << " Fstate = " << Fstate[i] << endl;
    }

    printf("\nSolving toy1 problem using derivatives...\n");

    // Reset the variables and solve ...

    int lenA = 6;
    int *iAfun = new int[lenA];
    int *jAvar = new int[lenA];
    double *A = new double[lenA];

    int lenG = 6;
    int *iGfun = new int[lenG];
    int *jGvar = new int[lenG];

    int neA, neG; // neA and neG must be defined when providing dervatives

    xstate[0] = 0;
    xstate[1] = 0;
    Fmul[0] = 0;
    Fmul[0] = 0;
    Fmul[0] = 0;
    x[0] = 1.0;
    x[1] = 1.0;


    // Provide the elements of the Jacobian explicitly.
    iGfun[0] = 0;
    jGvar[0] = 1;

    iGfun[1] = 1;
    jGvar[1] = 0;

    iGfun[2] = 1;
    jGvar[2] = 1;

    iGfun[3] = 2;
    jGvar[3] = 0;

    iGfun[4] = 2;
    jGvar[4] = 1;

    neG = 5;

//    iAfun[0] = 0;
//    jAvar[0] = 0;
//    A[0] = 1.0;
    neA = 0;

    ToyProb.initialize("", 0);
//    ToyProb.setProbName    ( "Toy1" );         // Give the problem a new name for Snopt.
//    ToyProb.setPrintFile   ( "Toy1.out" );
//    ToyProb.setSpecsFile   ( "sntoya.spc" );
    ToyProb.setIntParameter("Derivative option", 1);
    ToyProb.setIntParameter("Major Iteration limit", 250);
    ToyProb.setIntParameter("Verify level ", 3);
    ToyProb.solve(Cold, neF, n, ObjAdd,
                  ObjRow, toyusrfg_,
                  nullptr, nullptr, nullptr, 0, iGfun, jGvar, neG,
                  xlow, xupp, Flow, Fupp,
                  x, xstate, xmul,
                  F, Fstate, Fmul,
                  nS, nInf, sInf);

    for (int i = 0; i < n; i++) {
        cout << "x = " << x[i] << " xstate = " << xstate[i] << endl;
    }
    for (int i = 0; i < neF; i++) {
        cout << "F = " << F[i] << " Fstate = " << Fstate[i] << endl;
    }

    ToyProb.solve(Warm, neF, n, ObjAdd,
                  ObjRow, toyusrfg_,
                  nullptr, nullptr, nullptr, 0, iGfun, jGvar, neG,
                  xlow, xupp, Flow, Fupp,
                  x, xstate, xmul,
                  F, Fstate, Fmul,
                  nS, nInf, sInf);

    for (int i = 0; i < n; i++) {
        cout << "x = " << x[i] << " xstate = " << xstate[i] << endl;
    }
    for (int i = 0; i < neF; i++) {
        cout << "F = " << F[i] << " Fstate = " << Fstate[i] << endl;
    }

    delete[]iAfun;
    delete[]jAvar;
    delete[]A;
    delete[]iGfun;
    delete[]jGvar;

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

    return 0;
}
