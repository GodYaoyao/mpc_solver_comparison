//
// Created by yao on 19-12-10.
//

#ifndef SOLVER_COMPARISON_PARAMS_H
#define SOLVER_COMPARISON_PARAMS_H
const double g = 9.81;
const int step_N = 12;
const double dt = 0.1;
const int weight_x = 1;
const int weight_y = 1;
const int weight_phi = 1;
const int weight_v = 1;
const int weight_a = 1;
const int weight_wd = 1;

const int state = 5;
const int control = 2;

const int n_vars = state * step_N + control * (step_N - 1);
const int n_constrains = 1 + state * step_N;

// State
const int x_begin = 0;
const int y_begin = x_begin + step_N;
const int phi_begin = y_begin + step_N;
const int v_begin = phi_begin + step_N;
const int w_begin = v_begin + step_N;
// Control
const int a_begin = w_begin + step_N;
const int wd_begin = a_begin + step_N - 1;

void printInitState(double x_init, double y_init, double phi_init, double v_init, double w_init) {
    std::cout << "x = " << x_init << std::endl;
    std::cout << "y = " << y_init << std::endl;
    std::cout << "phi = " << phi_init << std::endl;
    std::cout << "v = " << v_init << std::endl;
    std::cout << "w = " << w_init << std::endl;
}

void printReferPoint(std::vector<std::vector<double>> *refer) {
    std::cout << "refer_point" << " = [";
    for (int i = 0; i < refer->size(); ++i) {
        std::cout << "(" << refer->at(i)[0] << ", " << refer->at(i)[1] << "), ";
    }
    std::cout << "\b\b]" << std::endl;
}

template<typename T>
void printSolutionResult(const T &sol, std::string name, int x1_begin, int x2_begin, int N) {
    std::cout << name << " = [";
    for (int i = 0; i < N; ++i) {
        std::cout << "(" << sol[x1_begin + i] << ", " << sol[x2_begin + i] << "), ";
    }
    std::cout << "\b\b]" << std::endl;
}

template<typename T>
void printTime(const std::vector<T> &time) {
    if (time.empty())
        return;
    T a = 0.;
    double mean_time = std::accumulate(time.begin(), time.end(), a) / time.size();
    std::cout << time.size() << std::endl;
    std::cout << mean_time << std::endl;
}
#endif //SOLVER_COMPARISON_PARAMS_H