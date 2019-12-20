//
// Created by yao on 19-12-19.
//
#include <iostream>
#include "ceres/ceres.h"

using namespace ceres;

struct CostFunctor {
    template<typename T>
    bool operator()(const T *const x, T *residual) const {
        residual[0] = (T(5.0) - x[0]);
        residual[1] = (T(-1.0) - x[1]);
        return true;
    }
};

int main(int argc, char **argv) {
    // The variable to solve for with its initial value.
    double initial_x = 0.5;
    double initial_y = 5.0;
    double x[2];
    x[0] = initial_x;
    x[1] = initial_y;

    // Build the problem.
    Problem problem;

    // Set up the only cost function (also known as residual). This uses
    // auto-differentiation to obtain the derivative (jacobian).
    CostFunction *cost_function =
        new AutoDiffCostFunction<CostFunctor, 2, 2>(new CostFunctor);
    problem.AddResidualBlock(cost_function, NULL, x);
    for (int i = 0; i < 2; ++i) {
        problem.SetParameterLowerBound(x, i, -2.);
        problem.SetParameterUpperBound(x, i, 6.);
    }

    // Run the solver!
    Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = false;
    Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    std::cout << summary.BriefReport() << "\n";
    std::cout << "x : " << initial_x << " -> " << x[0] << "\n";
    std::cout << "y : " << initial_y << " -> " << x[1] << "\n";

    return 0;
}