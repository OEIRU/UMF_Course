#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib> // Для вызова system
using namespace std;

const int nx = 100;
const int ny = 100;
const double hx = 1.0;
const double hy = 1.0;
const double tol = 1e-14;
const int max_iter = 10000;

// Аналитическое решение
double analytical_solution(double x, double y) {
    return x + y;
}

// Разметка области π
void mark_pi_shape(vector<vector<int>> &is_boundary) {
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            is_boundary[i][j] = 0;  // Сброс значений
            if (j == ny - 1 && i >= nx/4 && i <= 3*nx/4) {
                is_boundary[i][j] = 1;
            }
            if ((i == nx/4 || i == 3*nx/4) && j >= 0 && j <= ny-1) {
                is_boundary[i][j] = 1;
            }
        }
    }
}

// Метод Гаусса-Зейделя
void gauss_seidel(vector<double> &U, const vector<vector<int>> &is_boundary, 
                 const vector<double> &X, const vector<double> &Y) {
    double error = 1.0;
    int iter = 0;
    const double hx2 = hx * hx;
    const double hy2 = hy * hy;
    const double denominator = -2.0/hx2 - 2.0/hy2;

    while (error > tol && iter < max_iter) {
        error = 0.0;
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                if (is_boundary[i][j]) continue;
                
                double sum = 0.0;
                if (i > 0)    sum += U[(i-1)*ny + j]/hx2;
                if (i < nx-1) sum += U[(i+1)*ny + j]/hx2;
                if (j > 0)    sum += U[i*ny + (j-1)]/hy2;
                if (j < ny-1) sum += U[i*ny + (j+1)]/hy2;
                
                double new_val = sum / denominator;
                error += abs(new_val - U[i*ny + j]);
                U[i*ny + j] = new_val;
            }
        }
        error /= (nx * ny);
        iter++;
    }
    cout << "Converged in " << iter << " iterations (error: " << error << ")\n";
}

// Сохранение данных в файл
void save_to_file(const vector<double>& U, const string& filename) {
    ofstream out(filename);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            out << U[i*ny + j] << " ";
        }
        out << "\n";
    }
    cout << "Data saved to " << filename << endl;
}

// Визуализация с помощью GNUPLOT
void plot_with_gnuplot(const string& data_filename, const string& title) {
    string gnuplot_command = R"(
        set pm3d map
        set palette rgbformulae 33,13,10
        set xrange [0:]
        set yrange [0:]
        set title ')" + title + R"('
        splot ')" + data_filename + R"(' matrix with pm3d notitle
    )";

    string gnuplot_script = "gnuplot_script.tmp";
    ofstream script_file(gnuplot_script);
    if (!script_file.is_open()) {
        cerr << "Error: Cannot create GNUPLOT script file" << endl;
        exit(1);
    }

    script_file << gnuplot_command;
    script_file.close();

    cout << "Plotting data using GNUPLOT..." << endl;
    system(("gnuplot -persist " + gnuplot_script).c_str());

    // Удаление временного файла скрипта
    remove(gnuplot_script.c_str());
}

// Сравнение численного и аналитического решений
void compare_solutions(const vector<double>& U, const vector<double>& X, const vector<double>& Y) {
    double max_diff = 0.0;
    double avg_diff = 0.0;
    int count = 0;

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double analytical = analytical_solution(X[i], Y[j]);
            double numerical = U[i*ny + j];
            double diff = abs(analytical - numerical);

            if (diff > max_diff) max_diff = diff;
            avg_diff += diff;
            count++;
        }
    }

    avg_diff /= count;

    cout << "Comparison of analytical and numerical solutions:\n";
    cout << "Max difference: " << max_diff << "\n";
    cout << "Average difference: " << avg_diff << "\n";
}

int main() {
    vector<vector<int>> is_boundary(nx, vector<int>(ny, 0));
    vector<double> U(nx * ny, 0.0);
    
    // Инициализация координат
    vector<double> X(nx), Y(ny);
    for (int i = 0; i < nx; ++i) X[i] = i * hx;
    for (int j = 0; j < ny; ++j) Y[j] = j * hy;

    mark_pi_shape(is_boundary);

    // Инициализация граничных условий
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (is_boundary[i][j]) {
                U[i*ny + j] = analytical_solution(X[i], Y[j]);
            }
        }
    }

    // Решение методом Гаусса-Зейделя
    gauss_seidel(U, is_boundary, X, Y);

    // Сохранение численного решения
    string numerical_filename = "numerical_solution.txt";
    save_to_file(U, numerical_filename);

    // Вычисление аналитического решения
    vector<double> analytical_U(nx * ny, 0.0);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            analytical_U[i*ny + j] = analytical_solution(X[i], Y[j]);
        }
    }

    // Сохранение аналитического решения
    string analytical_filename = "analytical_solution.txt";
    save_to_file(analytical_U, analytical_filename);

    // Сравнение решений
    compare_solutions(U, X, Y);

    // Визуализация численного решения
    plot_with_gnuplot(numerical_filename, "Numerical Solution");

    // Визуализация аналитического решения
    plot_with_gnuplot(analytical_filename, "Analytical Solution");

    return 0;
}