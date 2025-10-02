#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>

using namespace std;

class Transport2D {
private:
    int nx, ny;
    double dx, dy, dt;
    double D, alpha, Q;
    vector<vector<double>> rho, rho_new;
    vector<vector<double>> u, v;
    ofstream output_file;

public:
    Transport2D(int nx, int ny, double Lx, double Ly, double D, double alpha, double Q) 
        : nx(nx), ny(ny), D(D), alpha(alpha), Q(Q) {
        dx = Lx / (nx - 1);
        dy = Ly / (ny - 1);
        
        rho.resize(nx, vector<double>(ny, 0.0));
        rho_new.resize(nx, vector<double>(ny, 0.0));
        u.resize(nx, vector<double>(ny, 0.0));
        v.resize(nx, vector<double>(ny, 0.0));
        
        dt = 0.25 *  min(dx*dx/(2*D), dy*dy/(2*D));

        output_file.open("transport_data.dat");
    }
    
    ~Transport2D() {
        if(output_file.is_open()) {
            output_file.close();
        }
    }
    
    void setVelocity(double u_val, double v_val) {
        for(int i = 0; i < nx; i++) {
            for(int j = 0; j < ny; j++) {
                u[i][j] = u_val;
                v[i][j] = v_val;
            }
        }
    }
    
    void setInitialCondition() {
        double x0 = nx/2.0, y0 = ny/2.0;
        double sigma = 5.0;
        
        for(int i = 0; i < nx; i++) {
            for(int j = 0; j < ny; j++) {
                double r2 = (i-x0)*(i-x0) + (j-y0)*(j-y0);
                rho[i][j] = exp(-r2/(2*sigma*sigma));
            }
        }
    }
    
    void timeStep() {

        double cfl_adv = max(abs(u[0][0]) * dt / dx, abs(v[0][0]) * dt / dy);
        if(cfl_adv > 0.5) {
            cout << "Warning: CFL = " << cfl_adv << " > 0.5\n";
        }
        
        for(int i = 1; i < nx-1; i++) {
            for(int j = 1; j < ny-1; j++) {
                double adv_x, adv_y;
                if(u[i][j] > 0) {
                    adv_x = -u[i][j] * (rho[i][j] - rho[i-1][j]) / dx;
                } else {
                    adv_x = -u[i][j] * (rho[i+1][j] - rho[i][j]) / dx;
                }
                
                if(v[i][j] > 0) {
                    adv_y = -v[i][j] * (rho[i][j] - rho[i][j-1]) / dy;
                } else {
                    adv_y = -v[i][j] * (rho[i][j+1] - rho[i][j]) / dy;
                }
                
                double diff_x = D * (rho[i+1][j] - 2*rho[i][j] + rho[i-1][j]) / (dx*dx);
                double diff_y = D * (rho[i][j+1] - 2*rho[i][j] + rho[i][j-1]) / (dy*dy);
                
                rho_new[i][j] = rho[i][j] + dt * (adv_x + adv_y + diff_x + diff_y + Q - alpha*rho[i][j]);
            }
        }
        

        for(int i = 0; i < nx; i++) {
            rho_new[i][0] = rho_new[i][1];
            rho_new[i][ny-1] = rho_new[i][ny-2];
        }
        for(int j = 0; j < ny; j++) {
            rho_new[0][j] = rho_new[1][j];
            rho_new[nx-1][j] = rho_new[nx-2][j];
        }
        
        rho = rho_new;
    }
    
    void saveData(int step) {
        output_file << "# Step " << step << "\n";
        for(int i = 0; i < nx; i++) {
            for(int j = 0; j < ny; j++) {
                output_file << i*dx << " " << j*dy << " " << rho[i][j] << "\n";
            }
        }
        output_file << "\n\n";
    }
};

int main() {
    int nx = 100, ny = 100;
    double Lx = 10.0, Ly = 10.0;
    double D = 0.1, alpha = 0.01, Q = 0.0;
    
    Transport2D solver(nx, ny, Lx, Ly, D, alpha, Q);
    solver.setVelocity(1.0, 0.5);
    solver.setInitialCondition();
    
    for(int step = 0; step <= 100; step++) {
        solver.saveData(step);
        if(step % 1 == 0) {
            cout << "Step " << step << " saved\n";
        }
        solver.timeStep();
    }
    
    return 0;
}
