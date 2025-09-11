#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

void tridiagonal(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, vector<double>& x, int n) {
    for (int i = 1; i < n; i++) {
        double m = a[i] / b[i-1];
        b[i] -= m * c[i-1];
        d[i] -= m * d[i-1];
    }
    x[n-1] = d[n-1] / b[n-1];
    for (int i = n-2; i >= 0; i--) {
        x[i] = (d[i] - c[i] * x[i+1]) / b[i];
    }
}

int main() {
    setlocale(LC_ALL, "ru.utf-8");

    double L = 1.0;
    double TimeMax = 0.1;
    double alpha = 1.0;
    int N = 101;
    int M = 1000;
    
    double h = L / (N - 1);
    double tau = TimeMax / (M - 1);
    double sigma = alpha * tau / (h * h);
    
    vector<vector<double>> T(M, vector<double>(N));
    
    // Начальные условия T(x,0) = 0
    for (int i = 0; i < N; i++) {
        T[0][i] = 0.0;
    }
    
    // Граничные условия
    for (int j = 0; j < M; j++) {
        T[j][N-1] = 1.0; // T(1,t) = 1
    }
    
    vector<double> a(N), b(N), c(N), d(N), x(N);
    
    for (int j = 1; j < M; j++) {
        // Заполнение трехдиагональной матрицы
        for (int i = 1; i < N-1; i++) {
            a[i] = -sigma;
            b[i] = 1 + 2*sigma;
            c[i] = -sigma;
            d[i] = T[j-1][i];
        }
        
        // Граничные условия в системе
        // При x=0: dT/dx = -2, используем T[1] - T[0] = -2*h
        b[0] = -1; c[0] = 1; d[0] = -2.0 * h;
        a[N-1] = 0; b[N-1] = 1; d[N-1] = 1.0;
        
        tridiagonal(a, b, c, d, x, N);
        
        for (int i = 0; i < N; i++) {
            T[j][i] = x[i];
        }
    }
    
    // Сохранение результатов
    ofstream file("../temperature_data.txt");
    for (int j = 0; j < M; j += 10) { // каждый 10-й временной слой
        for (int i = 0; i < N; i++) {
            file << T[j][i] << " ";
        }
        file << "\n";
    }
    file.close();
    
    cout << "Расчет завершен. Данные сохранены в temperature_data.txt" << endl;
}