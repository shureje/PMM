#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <corecrt_math_defines.h>


class NonlinearHeatSolver {
private:
    double L, T, beta;
    int N, M;
    double h, tau;
    std::vector<double> x, u_curr, u_next;
    
public:
    NonlinearHeatSolver(double L_, double T_, int N_, int M_, double beta_) 
        : L(L_), T(T_), N(N_), M(M_), beta(beta_) {
        h = L / (N - 1);
        tau = T / M;
        x.resize(N);
        u_curr.resize(N);
        u_next.resize(N);
        for (int i = 0; i < N; i++) {
            x[i] = i * h;
        }
    }
    
    double k(double u) {
        double kappa0 = 1.0;
        double beta = 2.0;
        return kappa0 * pow(fabs(u) + 1e-6, beta);
    }
    
    double source(double x, double t, double u) {
        return 0.0;
    }
    
    void setInitialConditions() {
        double xi0 = 0.5;
        for (int i = 0; i < N; i++) {
            if (x[i] < xi0) {
                u_curr[i] = pow(xi0 - x[i], 1.0/beta);  
            } else {
                u_curr[i] = 0.0;
            }
        }
    }
    
    void setBoundaryConditions(double t) {
        u_next[0] = u_next[1];
        u_next[N-1] = u_next[N-2];
    }
    
    void timeStep(double t) {
        for (int i = 1; i < N-1; i++) {
            double k_right = k((u_curr[i] + u_curr[i+1]) / 2.0);
            double k_left = k((u_curr[i-1] + u_curr[i]) / 2.0);
            double diff_term = (k_right * (u_curr[i+1] - u_curr[i]) - 
                               k_left * (u_curr[i] - u_curr[i-1])) / (h * h);
            double src = source(x[i], t, u_curr[i]);
            u_next[i] = u_curr[i] + tau * (diff_term + src);
        }
        setBoundaryConditions(t + tau);
        u_curr = u_next;
    }
    
    bool checkStability() {
        double max_k = 0.0;
        for (int i = 0; i < N; i++) {
            max_k = std::max(max_k, k(u_curr[i]));
        }
        return tau <= h * h / (2.0 * max_k);
    }
    
    void solve() {
        setInitialConditions();
        if (!checkStability()) {
            std::cout << "Warning: stability condition may be violated!" << std::endl;
            return;
        }
        std::ofstream file("solution.dat");
        for (int i = 0; i < N; i++) {
            file << x[i] << " " << 0.0 << " " << u_curr[i] << std::endl;
        }
        for (int n = 0; n < M; n++) {
            double t = n * tau;
            timeStep(t);
            if (n % (M/10) == 0) {
                for (int i = 0; i < N; i++) {
                    file << x[i] << " " << t + tau << " " << u_curr[i] << std::endl;
                }
            }
        }
        file.close();
        std::cout << "Solution saved to solution.dat" << std::endl;
    }
};

int main() {
    double L = 1.0;
    double T = 1.0;
    double beta = 2.0;
    int N = 101;
    int M = 50000;
    NonlinearHeatSolver solver(L, T, N, M, beta);
    solver.solve();
    return 0;
}