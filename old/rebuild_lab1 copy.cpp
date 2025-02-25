#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib> // Для вызова system
#include <iomanip>

using namespace std;

// Константы
const int nx = 100; // Количество узлов по x
const int ny = 100; // Количество узлов по y
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
        iter++;
    }
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

// Функция для чтения данных из файла result.txt
vector<vector<double>> read_data(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }
    vector<vector<double>> data;
    string line;
    while (getline(file, line)) {
        vector<double> row;
        stringstream ss(line);
        double value;
        while (ss >> value) {
            row.push_back(value);
        }
        if (!row.empty()) {
            data.push_back(row);
        }
    }
    file.close();
    return data;
}

// Функция для создания временного файла данных для GNUPLOT
void create_temp_data_file(const vector<vector<double>>& data, const string& temp_filename) {
    ofstream temp_file(temp_filename);
    if (!temp_file.is_open()) {
        cerr << "Error: Cannot create temporary file " << temp_filename << endl;
        exit(1);
    }
    for (size_t i = 0; i < data.size(); ++i) {
        for (size_t j = 0; j < data[i].size(); ++j) {
            temp_file << i << " " << j << " " << data[i][j] << endl;
        }
        temp_file << endl; // Пустая строка между строками данных
    }
    temp_file.close();
}

// Главная функция
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

    // Заполнение матрицы системы
    fill_matrix(is_boundary, Diag, Lower, Upper, Right);

    // Решение методом Гаусса-Зейделя
    gauss_seidel(Diag, Lower, Upper, Right, U);

    // Сохранение результата в файл
    string result_filename = "result.txt";
    save_result(U, result_filename);

    // Чтение данных из файла
    vector<vector<double>> data = read_data(result_filename);

    // Создание временного файла данных для GNUPLOT
    string temp_data_filename = "temp_data.txt";
    create_temp_data_file(data, temp_data_filename);

    // Команда для GNUPLOT
    string gnuplot_command = R"(
        set pm3d map
        set palette rgbformulae 33,13,10
        set xrange [0:]
        set yrange [0:]
        splot ')" + temp_data_filename + R"(' using 1:2:3 with pm3d notitle
    )";

    // Запись команды в скрипт GNUPLOT
    string gnuplot_script = "gnuplot_script.tmp";
    ofstream script_file(gnuplot_script);
    if (!script_file.is_open()) {
        cerr << "Error: Cannot create GNUPLOT script file" << endl;
        exit(1);
    }
    script_file << gnuplot_command;
    script_file.close();

    // Вызов GNUPLOT
    cout << "Plotting data using GNUPLOT..." << endl;
    system(("gnuplot -persist " + gnuplot_script).c_str());

    // Удаление временных файлов
    remove(temp_data_filename.c_str());
    remove(gnuplot_script.c_str());

    return 0;
}