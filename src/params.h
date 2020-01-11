//
// Created by yao on 19-12-10.
//

#ifndef SOLVER_COMPARISON_PARAMS_H
#define SOLVER_COMPARISON_PARAMS_H
const bool random_state = true;
const double g = 9.81;
const int step_N = 15;
const double dt = 0.1;
const int weight_x = 10;
const int weight_y = 10;
const int weight_phi = 10;
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

void generateInitState(double &x_init,
                       double &y_init,
                       double &phi_init,
                       double &v_init,
                       double &w_init,
                       bool random = true) {
    x_init = (random ? (double(rand()) / RAND_MAX) : 0.);
    y_init = (random ? (double(rand()) / RAND_MAX) : 0.);
    phi_init = -M_PI / 6 + (random ? (double(rand()) / RAND_MAX - 0.5) / 10 : 0.);
    v_init = 1.;
    w_init = 0.;
}

void generateReferPoint(std::vector<std::vector<double>> *&refer, bool random = true) {
    refer->at(0)[0] = 1.;
    refer->at(0)[1] = 0.;
    refer->at(0)[2] = -M_PI / 6;
    refer->at(0)[3] = 5.;
    for (int i = 1; i < step_N; ++i) {
        refer->at(i)[0] = refer->at(i - 1)[0] + refer->at(i - 1)[3] * cos(refer->at(i - 1)[2]) * dt;
        refer->at(i)[1] = refer->at(i - 1)[1] + refer->at(i - 1)[3] * sin(refer->at(i - 1)[2]) * dt;
        refer->at(i)[2] = refer->at(i - 1)[2] + (random ? (double(rand()) / RAND_MAX - 0.5) / 10 : 0.);
        refer->at(i)[3] = refer->at(i - 1)[3] + (random ? double(rand()) / RAND_MAX - 0.5 : 0.);
    }
}

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
    int mean_time = std::accumulate(time.begin(), time.end(), a) / time.size() + 1;
    std::cout << time.size() << std::endl;
    std::cout << mean_time << std::endl;
}
#endif //SOLVER_COMPARISON_PARAMS_H
