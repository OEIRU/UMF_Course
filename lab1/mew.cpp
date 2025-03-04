#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

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
    double k, q; // Коэффициенты уравнения
    vector<double> f; // Правая часть уравнения
    vector<double> u_exact; // Аналитическое решение
    vector<double> x; // Решение

public:
    MKR(double x0, double y0, double x1, double y1, int nx, int ny, double k, double q)
        : x0(x0), y0(y0), x1(x1), y1(y1), nx(nx), ny(ny), k(k), q(q), x(nx * ny, 0.0) {
        hx = (x1 - x0) / nx;
        hy = (y1 - y0) / ny;
    }

    // Создание сетки
    void create_grid() {
        // Генерация координат узлов
        vector<double> X(nx + 2), Y(ny + 2); // +2 для включения границ
        for (int i = 0; i <= nx + 1; i++) {
            X[i] = x0 + i * hx; // Координаты по x
        }
        for (int j = 0; j <= ny + 1; j++) {
            Y[j] = y0 + j * hy; // Координаты по y
        }
    
        // Вывод сетки
        cout << "Grid created:" << endl;
        for (int j = 0; j <= ny + 1; j++) {
            for (int i = 0; i <= nx + 1; i++) {
                cout << "(" << X[i] << ", " << Y[j] << ") ";
            }
            cout << endl;
        }
    }

    // Конечноразностная аппроксимация
    void create_SLAE() {
        int N = nx * ny;
        vector<vector<double>> A(N, vector<double>(5, 0.0)); // Диагональный формат
        vector<double> b(N, 0.0);

        for (int j = 1; j <= ny; j++) {
            for (int i = 1; i <= nx; i++) {
                int idx = (j - 1) * nx + (i - 1);

                // Коэффициенты матрицы A
                A[idx][0] = -k / (hy * hy); // Верхний сосед
                A[idx][1] = -k / (hx * hx); // Левый сосед
                A[idx][2] = 2.0 * k / (hx * hx) + 2.0 * k / (hy * hy) + q; // Центральный элемент
                A[idx][3] = -k / (hx * hx); // Правый сосед
                A[idx][4] = -k / (hy * hy); // Нижний сосед

                // Правая часть
                double x_coord = x0 + (i - 0.5) * hx;
                double y_coord = y0 + (j - 0.5) * hy;
                b[idx] = f_function(x_coord, y_coord);
            }
        }

        // Учет граничных условий
        apply_boundary_conditions(A, b);

        // Передача системы в решатель СЛАУ
        SLAU slau(N, nx, ny, 10000, 1e-8);
        slau.set_system(A, b);
        slau.solve();

        // Сохранение решения
        x = slau.get_solution();
    }

    // Применение граничных условий
    enum BoundaryType { DIRICHLET, ROBIN };

    void apply_boundary_conditions(vector<vector<double>> &A, vector<double> &b) {
        double alpha = 1.0; // Коэффициент в условии Robin
    
        // Левая граница (x = x0)
        for (int j = 0; j < ny; j++) {
            int idx = j * nx;
            if (boundary_type_left == DIRICHLET) {
                // Условия первого рода (Dirichlet)
                A[idx][0] = A[idx][1] = A[idx][3] = A[idx][4] = 0.0;
                A[idx][2] = 1.0;
                b[idx] = boundary_condition(x0, y0 + (j + 0.5) * hy);
            } else if (boundary_type_left == ROBIN) {
                // Условия третьего рода (Robin)
                double x_coord = x0;
                double y_coord = y0 + (j + 0.5) * hy;
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
                b[idx] = boundary_condition(x1, y0 + (j + 0.5) * hy);
            } else if (boundary_type_right == ROBIN) {
                // Условия третьего рода (Robin)
                double x_coord = x1;
                double y_coord = y0 + (j + 0.5) * hy;
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
                b[idx] = boundary_condition(x0 + (i + 0.5) * hx, y0);
            } else if (boundary_type_bottom == ROBIN) {
                // Условия третьего рода (Robin)
                double x_coord = x0 + (i + 0.5) * hx;
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
                b[idx] = boundary_condition(x0 + (i + 0.5) * hx, y1);
            } else if (boundary_type_top == ROBIN) {
                // Условия третьего рода (Robin)
                double x_coord = x0 + (i + 0.5) * hx;
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
        u_exact.resize(nx * ny);
        for (int j = 1; j <= ny; j++) {
            for (int i = 1; i <= nx; i++) {
                int idx = (j - 1) * nx + (i - 1);
                double x_coord = x0 + (i - 0.5) * hx;
                double y_coord = y0 + (j - 0.5) * hy;
                u_exact[idx] = x_coord * x_coord + y_coord * y_coord;
            }
        }

        // Вычисление ошибки
        double error = 0.0;
        for (int i = 0; i < nx * ny; i++) {
            error += pow(x[i] - u_exact[i], 2);
        }
        error = sqrt(error / (nx * ny));
        cout << "Polynomial test error: " << error << endl;
    }

    // Исследование порядка аппроксимации
    void study_approximation_order() {
        vector<int> n_values = {10, 20, 40, 80};
        vector<double> errors;

        for (auto n : n_values) {
            MKR mkr(x0, y0, x1, y1, n, n, k, q);
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

    // Получение ошибки
    double get_error() const {
        double error = 0.0;
        for (int i = 0; i < nx * ny; i++) {
            error += pow(x[i] - u_exact[i], 2);
        }
        return sqrt(error / (nx * ny));
    }

private:
    BoundaryType boundary_type_left = DIRICHLET;   // Тип КУ на левой границе
    BoundaryType boundary_type_right = DIRICHLET;  // Тип КУ на правой границе
    BoundaryType boundary_type_bottom = DIRICHLET; // Тип КУ на нижней границе
    BoundaryType boundary_type_top = DIRICHLET;    // Тип КУ на верхней границе
};

// Главная функция
int main() {
    // Параметры задачи
    double x0 = 0.0, y0 = 0.0, x1 = 1.0, y1 = 1.0;
    int nx = 20, ny = 20;
    double k = 1.0, q = 0.0;

    // Создание объекта MKR
    MKR mkr(x0, y0, x1, y1, nx, ny, k, q);

    // Создание сетки
    mkr.create_grid();

    // Создание и решение СЛАУ
    mkr.create_SLAE();

    // Тестирование на полиномах
    mkr.test_polynomial();

    // Исследование порядка аппроксимации
    mkr.study_approximation_order();

    return 0;
}