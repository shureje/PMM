#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

class NonlinearHeatSolver {
private:
    double L, dt, T_max, k0, beta;
    int N, M;
    double h;
    vector<double> x, T_old, T_new;
    

    
public:
    NonlinearHeatSolver(double L_, int N_, int M_, double T_max_, double k0_, double beta_) 
        : L(L_), N(N_), M(M_), T_max(T_max_), k0(k0_), beta(beta_) {
        h = L / (N - 1);
        dt = T_max / M;
        x.resize(N);
        T_old.resize(N);
        T_new.resize(N);
        for (int i = 0; i < N; i++) {
            x[i] = i * h;
        }
    }

    double k(double T) {
        if (T <= 0) return 1e-10;
        double k_val = k0 * pow(T, beta);
        return min(k_val, 2.0);  // ограничение коэф. теплопроводности
    }
    
    void solve() {
        // Начальное условие T(x,0) = 0
        for (int i = 0; i < N; i++) {
            T_old[i] = 0.1;
        }
        
        double max_k = k0 * pow(1.0, beta);
        double stability = max_k * dt / (h * h);
        cout << "Число устойчивости: " << stability << " (должно быть < 0.5)" << endl;
        
        if (stability >= 0.5) {
            cout << "ВНИМАНИЕ: Схема неустойчива! Уменьшите dt или увеличьте h" << endl;
            return;
        }
        
        ofstream file("solution.dat");
        
        for (int n = 0; n < M; n++) {
            double t = n * dt;
            
            for (int i = 1; i < N-1; i++) {
                double k_left = 0.5 * (k(T_old[i]) + k(T_old[i-1]));
                double k_right = 0.5 * (k(T_old[i]) + k(T_old[i+1]));
                
                double flux_left = k_left * (T_old[i] - T_old[i-1]) / h;
                double flux_right = k_right * (T_old[i+1] - T_old[i]) / h;
                
                T_new[i] = T_old[i] + dt * (flux_right - flux_left) / h;
                if (T_new[i] < 0) T_new[i] = 0.01;
            }
            
            // Граничное условие T_x(0,t) = -2
            T_new[0] = T_new[1] + 2 * h;
            T_new[N-1] = 1.0;
            
            for (int i = 0; i < N; i++) {
                file << x[i] << " " << t << " " << T_new[i] << endl;
            }
            
            T_old = T_new;
        }
        
        file.close();
        cout << "Solution saved to solution.dat" << endl;
    }
};

int main() {
    double L = 1.0;
    double T_max = 2;
    double k0 = 1.0;
    double beta = 1.2; // должно быть > 1
    int N = 50;
    int M = 20000;

    setlocale(LC_ALL, "ru.utf-8");
    
    NonlinearHeatSolver solver(L, N, M, T_max, k0, beta);
    solver.solve();
    return 0;
}