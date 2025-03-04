#include <iostream>
#include <fstream> 
#include <vector>
#include <cmath>
#include <iomanip>
#include <cstdlib>

using namespace std;

enum BoundaryType {
    DIRICHLET, // Условие первого рода
    ROBIN      // Условие третьего рода
};

// Базовый класс для решения СЛАУ методом блочной релаксации
class SLAU {
protected:
    int N; // Размерность системы
    int nx, ny; // Размеры сетки
    vector<vector<double>> A; // Матрица коэффициентов (в диагональном формате)
    vector<double> b; // Вектор правой части
    vector<double> x; // Вектор неизвестных
    int maxiter; // Максимальное число итераций
    double tol; // Требуемая точность

public:
    SLAU(int N, int nx, int ny, int maxiter, double tol) 
         : N(N), nx(nx), ny(ny), maxiter(maxiter), tol(tol), x(N, 0.0) {}

    // Метод блочной релаксации для решения СЛАУ
    void solve() {
        for (int k = 0; k < maxiter; k++) {
            double residual = 0.0;
            for (int i = 0; i < N; i++) {
                double S = A[i][2] * x[i];
                if (i > 0) S += A[i][1] * x[i - 1]; // Левый сосед
                if (i < N - 1) S += A[i][3] * x[i + 1]; // Правый сосед
                if (i >= nx) S += A[i][0] * x[i - nx]; // Верхний сосед
                if (i < N - nx) S += A[i][4] * x[i + nx]; // Нижний сосед

                double new_x = x[i] + (b[i] - S) / A[i][2];
                residual += abs(b[i] - S);
                x[i] = new_x;
            }
            if (residual < tol) {
                cout << "Converged in " << k << " iterations." << endl;
                break;
            }
        }
    }

    // Установка матрицы A и вектора b
    void set_system(const vector<vector<double>> &A_, const vector<double> &b_) {
        A = A_;
        b = b_;
    }

    // Получение решения
    const vector<double>& get_solution() const { return x; }
};

// Класс для решения задачи теплопроводности методом конечных разностей
class MKR {
private:
    int nx, ny; // Размеры сетки
    double hx, hy; // Шаги сетки
    double x0, y0, x1, y1; // Границы области
    double xc1, yc1, xc2, yc2; // Границы выреза
    double k, q; // Коэффициенты уравнения
    vector<double> u_exact; // Аналитическое решение
    vector<double> x; // Решение
    BoundaryType boundary_type_left = DIRICHLET;    // Тип КУ на левой границе
    BoundaryType boundary_type_right = DIRICHLET;  // Тип КУ на правой границе
    BoundaryType boundary_type_bottom = DIRICHLET; // Тип КУ на нижней границе
    BoundaryType boundary_type_top = DIRICHLET;    // Тип КУ на верхней границе

public:
    MKR(double x0, double y0, double x1, double y1, double xc1, double yc1, double xc2, double yc2,
        int nx, int ny, double k, double q)
        : x0(x0), y0(y0), x1(x1), y1(y1), xc1(xc1), yc1(yc1), xc2(xc2), yc2(yc2),
          nx(nx), ny(ny), k(k), q(q) {
        hx = (x1 - x0) / nx;
        hy = (y1 - y0) / ny;
    }

    // Создание сетки
    void create_grid() {
        vector<double> X(nx + 2), Y(ny + 2); // +2 для включения границ
        for (int i = 0; i <= nx; i++) {
            X[i] = x0 + i * hx;
        }
        for (int j = 0; j <= ny; j++) {
            Y[j] = y0 + j * hy;
        }

        int active_nodes = 0;
        for (int j = 1; j <= ny; j++) {
            for (int i = 1; i <= nx; i++) {
                double x_coord = X[i];
                double y_coord = Y[j];

                // Исключаем узлы внутри "выреза" и ниже него
                if (x_coord >= xc1 && x_coord <= xc2 && y_coord >= yc1 && y_coord <= yc2) {

                    continue;
                }

                active_nodes++;
            }
        }

        x.resize(active_nodes);
        u_exact.resize(active_nodes);
    }

    // Создание СЛАУ
    void create_SLAE() {
        int N = x.size();
        vector<vector<double>> A(N, vector<double>(5, 0.0)); // Диагональный формат
        vector<double> b(N, 0.0);

        int idx = 0;
        for (int j = 1; j <= ny; j++) {
            for (int i = 1; i <= nx; i++) {
                double x_coord = x0 + i * hx;
                double y_coord = y0 + j * hy;

                // Исключаем узлы внутри "выреза" и ниже него
                if (x_coord >= xc1 && x_coord <= xc2 && y_coord >= yc1 && y_coord <= yc2) {

                    continue;
                }

                // Заполняем матрицу A и вектор b
                A[idx][0] = -k / (hy * hy); // Верхний сосед
                A[idx][1] = -k / (hx * hx); // Левый сосед
                A[idx][2] = 2 * k * (1 / (hx * hx) + 1 / (hy * hy)) + q; // Центральный элемент
                A[idx][3] = -k / (hx * hx); // Правый сосед
                A[idx][4] = -k / (hy * hy); // Нижний сосед

                b[idx] = f_function(x_coord, y_coord);

                // Учет граничных условий
                if (x_coord == x0 || x_coord == x1 || y_coord == y0 || y_coord == y1) {
                    double bc_value = boundary_condition(x_coord, y_coord);
                    if (x_coord == x0) {
                        b[idx] += k / (hx * hx) * bc_value; // Левая граница
                    }
                    if (x_coord == x1) {
                        b[idx] += k / (hx * hx) * bc_value; // Правая граница
                    }
                    if (y_coord == y0) {
                        b[idx] += k / (hy * hy) * bc_value; // Нижняя граница
                    }
                    if (y_coord == y1) {
                        b[idx] += k / (hy * hy) * bc_value; // Верхняя граница
                    }
                }

                idx++;
            }
        }

        SLAU slau(N, nx, ny, 1000, 1e-12);
        slau.set_system(A, b);
        slau.solve();
        x = slau.get_solution();
    }

    void apply_boundary_conditions(vector<vector<double>> &A, vector<double> &b) {
        double alpha = 1.0; // Коэффициент в условии Робена
    
        // Левая граница (x = x0)
        for (int j = 0; j < ny; j++) {
            int idx = j * nx;
            if (boundary_type_left == DIRICHLET) {
                // Условия первого рода (Dirichlet)
                A[idx][0] = A[idx][1] = A[idx][3] = A[idx][4] = 0.0;
                A[idx][2] = 1.0;
                b[idx] = boundary_condition(x0, y0 + j * hy);
            } else if (boundary_type_left == ROBIN) {
                // Условия третьего рода (Robin)
                double x_coord = x0;
                double y_coord = y0 + j * hy;
                double g_value = k * 2 * x_coord + alpha * (x_coord * x_coord + y_coord * y_coord);
                A[idx][0] = 0.0; // Верхний сосед
                A[idx][1] = -k / (2 * hx); // Левый сосед (центральная разность)
                A[idx][2] = k / (2 * hx) + alpha; // Центральный элемент
                A[idx][3] = 0.0; // Правый сосед
                A[idx][4] = 0.0; // Нижний сосед
                b[idx] = g_value;
            }
        }
    
        // Правая граница (x = x1)
        for (int j = 0; j < ny; j++) {
            int idx = j * nx + (nx - 1);
            if (boundary_type_right == DIRICHLET) {
                // Условия первого рода (Dirichlet)
                A[idx][0] = A[idx][1] = A[idx][3] = A[idx][4] = 0.0;
                A[idx][2] = 1.0;
                b[idx] = boundary_condition(x1, y0 + j * hy);
            } else if (boundary_type_right == ROBIN) {
                // Условия третьего рода (Robin)
                double x_coord = x1;
                double y_coord = y0 + j * hy;
                double g_value = k * 2 * x_coord + alpha * (x_coord * x_coord + y_coord * y_coord);
                A[idx][0] = 0.0; // Верхний сосед
                A[idx][1] = 0.0; // Левый сосед
                A[idx][2] = k / (2 * hx) + alpha; // Центральный элемент
                A[idx][3] = -k / (2 * hx); // Правый сосед (центральная разность)
                A[idx][4] = 0.0; // Нижний сосед
                b[idx] = g_value;
            }
        }
    
        // Нижняя граница (y = y0)
        for (int i = 0; i < nx; i++) {
            int idx = i;
            if (boundary_type_bottom == DIRICHLET) {
                // Условия первого рода (Dirichlet)
                A[idx][0] = A[idx][1] = A[idx][3] = A[idx][4] = 0.0;
                A[idx][2] = 1.0;
                b[idx] = boundary_condition(x0 + i * hx, y0);
            } else if (boundary_type_bottom == ROBIN) {
                // Условия третьего рода (Robin)
                double x_coord = x0 + i * hx;
                double y_coord = y0;
                double g_value = k * 2 * y_coord + alpha * (x_coord * x_coord + y_coord * y_coord);
                A[idx][0] = -k / (2 * hy); // Верхний сосед (центральная разность)
                A[idx][1] = 0.0; // Левый сосед
                A[idx][2] = k / (2 * hy) + alpha; // Центральный элемент
                A[idx][3] = 0.0; // Правый сосед
                A[idx][4] = 0.0; // Нижний сосед
                b[idx] = g_value;
            }
        }
    
        // Верхняя граница (y = y1)
        for (int i = 0; i < nx; i++) {
            int idx = (ny - 1) * nx + i;
            if (boundary_type_top == DIRICHLET) {
                // Условия первого рода (Dirichlet)
                A[idx][0] = A[idx][1] = A[idx][3] = A[idx][4] = 0.0;
                A[idx][2] = 1.0;
                b[idx] = boundary_condition(x0 + i * hx, y1);
            } else if (boundary_type_top == ROBIN) {
                // Условия третьего рода (Robin)
                double x_coord = x0 + i * hx;
                double y_coord = y1;
                double g_value = k * 2 * y_coord + alpha * (x_coord * x_coord + y_coord * y_coord);
                A[idx][0] = 0.0; // Верхний сосед
                A[idx][1] = 0.0; // Левый сосед
                A[idx][2] = k / (2 * hy) + alpha; // Центральный элемент
                A[idx][3] = 0.0; // Правый сосед
                A[idx][4] = -k / (2 * hy); // Нижний сосед (центральная разность)
                b[idx] = g_value;
            }
        }
    }

    // Функция для правой части уравнения
    double f_function(double x, double y) {
        return -4.0 * k + q * (x * x + y * y); // Для u(x, y) = x^2 + y^2
    }

    // Функция для граничных условий первого рода
    double boundary_condition(double x, double y) {
        return x * x + y * y; // Для u(x, y) = x^2 + y^2
    }

    // Тестирование на полиномах
    void test_polynomial() {
        // Вычисление точного решения
        int idx = 0;
        for (int j = 1; j <= ny; j++) {
            for (int i = 1; i <= nx; i++) {
                double x_coord = x0 + i * hx;
                double y_coord = y0 + j * hy;

                // Исключаем узлы внутри "выреза" и ниже него
                if (x_coord >= xc1 && x_coord <= xc2 && y_coord >= yc1 && y_coord <= yc2) {

                    continue;
                }

                u_exact[idx++] = x_coord * x_coord + y_coord * y_coord;
            }
        }

        // Вычисление ошибки
        double error = 0.0;
        for (int i = 0; i < x.size(); i++) {
            error += pow(x[i] - u_exact[i], 2);
        }
        error = sqrt(error / x.size());
        cout << "Error: " << error << endl;
    }

    // Исследование порядка аппроксимации
    void study_approximation_order() {
        vector<int> n_values = {10, 20, 40, 80};
        vector<double> errors;

        for (auto n : n_values) {
            MKR mkr(x0, y0, x1, y1, xc1, yc1, xc2, yc2, n, n, k, q);
            mkr.create_grid();
            mkr.create_SLAE();
            mkr.test_polynomial();
            errors.push_back(mkr.get_error());
        }

        // Вычисление порядка аппроксимации
        for (size_t i = 1; i < errors.size(); i++) {
            double order = log(errors[i - 1] / errors[i]) / log(2.0);
            cout << "Approximation order for n=" << n_values[i] << ": " << order << endl;
        }
    }

    // Сохранение данных для gnuplot
    void save_to_file(const string& filename) {
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: Unable to open file for writing." << endl;
            return;
        }

        int idx = 0;
        for (int j = 1; j <= ny; j++) {
            for (int i = 1; i <= nx; i++) {
                double x_coord = x0 + i * hx;
                double y_coord = y0 + j * hy;

                // Исключаем узлы внутри "выреза" и ниже него
                if (x_coord >= xc1 && x_coord <= xc2 && y_coord >= yc1 && y_coord <= yc2) {

                    file << x_coord << " " << y_coord << " " << 0.0 << "\n"; // Внутри выреза значение 0
                } else {
                    file << x_coord << " " << y_coord << " " << x[idx++] << "\n";
                }
            }
            file << "\n"; // Разделитель между строками для gnuplot
        }

        file.close();
        cout << "Data saved to " << filename << endl;
        
    }

    // Получение ошибки
    double get_error() const {
        double error = 0.0;
        for (int i = 0; i < x.size(); i++) {
            error += pow(x[i] - u_exact[i], 2);
        }
        return sqrt(error / x.size());
    }
};

int plot_output() {
    ofstream script("plot_script.gp");
    if (!script) {
        cerr << "Ошибка создания скрипта." << endl;
        return 1;
    }
    script << "set pm3d map\n";
    script << "splot 'solution_data.txt'\n"; 
    script << "pause -1\n";
    script.close();

    system("gnuplot plot_script.gp");
    remove("plot_script.gp");
    return 0;

}
// Главная функция
int main() {
    // Параметры задачи
    double x0 = 0.0, y0 = 0.0, x1 = 1.0, y1 = 1.0;
    double xc1 = 0.25, yc1 = 0.0, xc2 = 0.75, yc2 = 0.75; // Границы выреза
    int nx = 100, ny = 100;
    double k = 1.0, q = 0.0;

    // Создание объекта MKR
    MKR mkr(x0, y0, x1, y1, xc1, yc1, xc2, yc2, nx, ny, k, q);

    // Создание сетки
    mkr.create_grid();

    // Создание и решение СЛАУ
    mkr.create_SLAE();

    // Тестирование на полиномах
    mkr.test_polynomial();

    // Исследование порядка аппроксимации
    mkr.study_approximation_order();

    // Сохранение данных для gnuplot
    mkr.save_to_file("solution_data.txt");

    plot_output();

    return 0;
}