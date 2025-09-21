#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>


using namespace std;

int main() {
    const int N = 100;
    const int M = 200;
    const double L = 1.0;
    const double T = 1.0;
    const double dx = L / (N - 1);
    const double dt = T / (M - 1);
    const double c = 1.0;
    const double r = c * dt / dx;
    
    vector<vector<double>> u(M, vector<double>(N, 0.0));
    
    // НУ: u(x,0) = 0, ut(x,0) = sin(2πx)
    for (int i = 0; i < N; i++) {
        double x = i * dx;
        u[0][i] = 0.0;
        if (i > 0 && i < N-1) {
            u[1][i] = u[0][i] + dt * sin(2 * M_PI * x);
        }
    }
    
    // ГУ: u(0,t) = u(1,t) = 0
    for (int j = 0; j < M; j++) {
        u[j][0] = 0.0;
        u[j][N-1] = 0.0;
    }
    
    // Схема 
    for (int j = 1; j < M-1; j++) {
        for (int i = 1; i < N-1; i++) {
            u[j+1][i] = 2*u[j][i] - u[j-1][i] + r*r*(u[j][i+1] - 2*u[j][i] + u[j][i-1]);
        }
    }
    
    // Вывод
    ofstream file("wave_data.txt");
    for (int j = 0; j < M; j += 5) {
        double t = j * dt;
        file << "t=" << t << "\n";
        for (int i = 0; i < N; i++) {
            double x = i * dx;
            file << x << " " << u[j][i] << "\n";
        }
        file << "\n";
    }
    file.close();
    
    cout << "Данные записаны в wave_data.txt" << endl;
    system("python plot.py");
    
    return 0;
}