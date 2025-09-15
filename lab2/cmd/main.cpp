#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

class HeatEquation
{
private:
    double L;
    double T;
    int N;
    int M;
    double h;
    double tau;
    double a;
    double r;

    vector<vector<double>> u;
    
public:
    HeatEquation(double L, double T, int N, int M, double a): L(L), T(T), N(N), M(M), a(a) {
        h = L/N;
        tau = T/M;
        r = a*tau/(h*h);

        if (r > 0.5) {
            cout << "Внимание: r = " << r << " > 0.5. Условие устойчивости не выполнено." << endl;
            return;
        }

        u.resize(N, vector<double>(M + 1));
    }


    void SetInitialCondition() {
        for (int i = 0; i < N; i++) {
            u[i][0] = 0;
        }
        u[0][0] = u[1][0] + 2*h;
        u[N-1][0] = 1;
    }

    void setBoundaryCondition(int n) {
       
        u[0][n+1] = u[1][n+1] + 2*h;
        u[N-1][n+1] = 1;

    }

    void solve() {
        SetInitialCondition();

        for (int n=0; n < M; n++) {
            
            for (int i = 1; i < N-1; i++) {
                u[i][n+1] = u[i][n] + r*(u[i+1][n] - 2*u[i][n] + u[i-1][n]);
            }
            setBoundaryCondition(n);
        }

       
    }

    void saveToFile() {
        ofstream file("result.txt");
        if (file.is_open()) {
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < M; j+=50) {
                    file << u[i][j] << " ";
                }
                file << endl;
            }
            file.close();
        } else {
            cout << "Не удалось открыть файл для записи." << endl;
        }
    }
};


int main() {
    setlocale(LC_ALL, "ru.utf-8");
    double L = 1.0;
    double T = 1.0;
    int N = 50;
    int M = 20000;
  
    HeatEquation heatEquation(L, T, N, M, 1.0);

    heatEquation.solve();
    heatEquation.saveToFile();


    return 0;
}