#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// Константы
const int nx = 200; // Количество узлов по x
const int ny = 200; // Количество узлов по y
const double hx = 1.0; // Шаг по x
const double hy = 1.0; // Шаг по y
const double alpha = 0.5; // Коэффициент для граничных условий
const double boundary_beta = 1.0;  // Свободный член для граничных условий
const double tol = 1e-14;  // Точность сходимости
const int max_iter = 10000; // Максимальное количество итераций

// Функция для разметки области π
void mark_pi_shape(vector<vector<bool>>& is_boundary) {
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            // Горизонтальная крыша
            if (j == ny - 1) {
                if (i >= nx / 4 && i <= 3 * nx / 4) {
                    is_boundary[i][j] = true;
                }
            }
            // Левая ножка
            if (i == nx / 4) {
                if (j >= 0 && j <= ny - 1) {
                    is_boundary[i][j] = true;
                }
            }
            // Правая ножка
            if (i == 3 * nx / 4) {
                if (j >= 0 && j <= ny - 1) {
                    is_boundary[i][j] = true;
                }
            }
        }
    }
}

// Логирование границ 
/*
void log_boundary(const vector<vector<bool>>& is_boundary) {
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            cout << (is_boundary[i][j] ? "1 " : "0 ");
        }
        cout << endl;
    }
}
*/

// Функция для заполнения матрицы системы
void fill_matrix(const vector<vector<bool>>& is_boundary, vector<double>& Diag, vector<double>& Lower,
                 vector<double>& Upper, vector<double>& Right) {
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            int idx = i * ny + j;

            if (is_boundary[i][j]) {
                // Граничные условия
                Diag[idx] = 1.0;
                Right[idx] = boundary_beta / alpha;
            } else {
                // Внутренние узлы
                Diag[idx] = -2.0 / (hx * hx) - 2.0 / (hy * hy);
                if (i > 0) Lower[idx - ny] = 1.0 / (hx * hx);
                if (i < nx - 1) Upper[idx] = 1.0 / (hx * hx);
                if (j > 0) Lower[idx - 1] = 1.0 / (hy * hy);
                if (j < ny - 1) Upper[idx] = 1.0 / (hy * hy);
                Right[idx] = 0.0; // Свободный член для внутренних узлов
            }
        }
    }
}

// Метод Гаусса-Зейделя
void gauss_seidel(const vector<double>& Diag, const vector<double>& Lower, const vector<double>& Upper,
                  const vector<double>& Right, vector<double>& U) {
    double error = 1.0;
    int iter = 0;

    while (error > tol && iter < max_iter) {
        error = 0.0;

        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                int idx = i * ny + j;

                double old_U = U[idx];
                double new_U = Right[idx];

                if (i > 0) new_U -= Lower[idx - ny] * U[idx - ny];
                if (i < nx - 1) new_U -= Upper[idx] * U[idx + ny];
                if (j > 0) new_U -= Lower[idx - 1] * U[idx - 1];
                if (j < ny - 1) new_U -= Upper[idx] * U[idx + 1];

                new_U /= Diag[idx];
                U[idx] = new_U;

                error += fabs(new_U - old_U);
            }
        }

        error /= (nx * ny);
        //cout << "Iteration " << iter << ": Error = " << error << endl;
        iter++;
    }

    //cout << "Converged in " << iter << " iterations" << endl;
}

// Сохранение результата в файл
void save_result(const vector<double>& U, const string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return;
    }

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            file << U[i * ny + j] << " ";
        }
        file << endl;
    }

    file.close();
}

int main() {
    // Инициализация
    vector<vector<bool>> is_boundary(nx, vector<bool>(ny, false));
    vector<double> Diag(nx * ny, 0.0);
    vector<double> Lower(nx * ny, 0.0);
    vector<double> Upper(nx * ny, 0.0);
    vector<double> Right(nx * ny, 0.0);
    vector<double> U(nx * ny, 0.0);

    // Разметка области π
    mark_pi_shape(is_boundary);

    // Логирование границ
    // log_boundary(is_boundary);

    // Заполнение матрицы системы
    fill_matrix(is_boundary, Diag, Lower, Upper, Right);

    // Решение методом Гаусса-Зейделя
    gauss_seidel(Diag, Lower, Upper, Right, U);

    // Сохранение результата
    save_result(U, "result.txt");

    return 0;
}