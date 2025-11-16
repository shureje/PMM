#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

const int NX = 50;
const int NY = 50;
const double T = 1.5;
const double a = 2.0;
const double b = 3.0;
const double dx = a / (NX - 1);
const double dy = b / (NY - 1);
const double dt = 0.25 * min(dx*dx, dy*dy);
const int TIME_STEPS = T / dt;

double f(double x, double y) {
    return -exp(-x*y*y);
}

void solve_heat_equation() {
    vector<vector<double>> u(NX, vector<double>(NY, 0.0));
    ofstream file("heat_animation.txt");
    
    for (int step = 0; step < TIME_STEPS; step++) {
        vector<vector<double>> u_new = u;
        
        for (int i = 0; i < NX; i++) {
            double x = i * dx;
            u_new[i][0] = x*x - 2*x;
            u_new[i][NY-1] = x*x*x - 4*x;
        }
        for (int j = 0; j < NY; j++) {
            double y = j * dy;
            u_new[0][j] = y*y - 3*y;
            u_new[NX-1][j] = y*y*y - 9*y;
        }
        
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                double x = i * dx;
                double y = j * dy;
                double laplacian = (u[i-1][j] - 2*u[i][j] + u[i+1][j])/(dx*dx) +
                                 (u[i][j-1] - 2*u[i][j] + u[i][j+1])/(dy*dy);
                u_new[i][j] = u[i][j] + dt * (laplacian + f(x, y));
            }
        }
        
        u = u_new;
        
        if (step % 10 == 0) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    file << i*dx << " " << j*dy << " " << step << " " << u[i][j] << endl;
                }
            }
        }
    }
    
    file.close();
}

int main() {
    solve_heat_equation();
    return 0;
}