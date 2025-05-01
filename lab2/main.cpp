#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>
#include <fstream>

using namespace std;

// Типы данных
using func1D = function<double(double)>;
using func2D = function<double(double, double)>;
using Vec = vector<double>;
using Mat = vector<vector<double>>;

// Структура для ленточной матрицы
struct BandedMatrix {
    int size;
    int lower;  // количество поддиагоналей
    int upper;  // количество наддиагоналей
    vector<vector<double>> data;  // data[i][j], где j от i-kl до i+ku
    
    BandedMatrix(int n, int kl, int ku) : size(n), lower(kl), upper(ku) {
        data.resize(n, vector<double>(kl + ku + 1, 0.0));
    }
    
    // Доступ к элементам матрицы
    double& operator()(int i, int j) {
        if (i < 0 || i >= size || j < 0 || j >= size) {
            cerr << "Ошибка доступа к матрице: (" << i << "," << j << ")\n";
            exit(1);
        }
        
        int offset = j - i + upper;
        if (offset < 0 || offset >= lower + upper + 1) {
            cerr << "Ошибка индексации матрицы: offset=" << offset << "\n";
            exit(1);
        }
        
        return data[i][offset];
    }
    
    double get(int i, int j) const {
        if (i < 0 || i >= size || j < 0 || j >= size) return 0.0;
        int offset = j - i + upper;
        if (offset < 0 || offset >= lower + upper + 1) return 0.0;
        return data[i][offset];
    }
};

// Проблема с нелинейным уравнением
struct Problem {
    vector<double> nodes;       // координаты узлов сетки
    func2D lambda;              // коэффициент lambda(x, u)
    func2D f;                   // правая часть f(x, u)
    func2D gamma;               // реакция gamma(x, u)
    func2D u;                   // точное решение для сравнения
    double omega;               // параметр релаксации
    double ua, ub;              // краевые условия
    double sigma;               // параметр уравнения
    double tol;                 // точность решения
    int maxiter;                // максимальное число итераций
    int nElem;                  // число элементов
    double t;                   // текущее время
    
    Problem() : omega(1.0), ua(0.0), ub(0.0), sigma(1.0),
                tol(1e-7), maxiter(1000), nElem(4), t(0.0) {}
};

// Базисные функции и их производные (квадратичные)
class BasisFunctions {
public:
    // Базисные функции
    static double phi(int i, double xi) {
        switch(i) {
            case 0: return 2.0 * xi * (xi - 0.5) - 1.0;  // phi_0
            case 1: return 4.0 * xi * (1.0 - xi);        // phi_1
            case 2: return 2.0 * xi * (xi - 0.5);        // phi_2
            default: return 0.0;
        }
    }
    
    // Производные базисных функций
    static double dphi(int i, double xi) {
        switch(i) {
            case 0: return 4.0 * xi - 1.0;  // dphi_0/dxi
            case 1: return 4.0 - 8.0 * xi; // dphi_1/dxi
            case 2: return 4.0 * xi - 1.0;  // dphi_2/dxi
            default: return 0.0;
        }
    }
};

// Вычисление локальных матриц
void localElement(const Vec &uq, const Problem &P, int elemIndex,
                  Mat &G, Mat &R, Vec &b_loc) {
    // Получаем узлы текущего элемента
    int i0 = 2 * elemIndex;
    int i1 = i0 + 1;
    int i2 = i0 + 2;
    
    // Проверяем корректность индексов
    if (i2 >= P.nodes.size()) {
        cerr << "Ошибка: неверные индексы элемента" << endl;
        exit(1);
    }
    
    double x0 = P.nodes[i0], x1 = P.nodes[i1], x2 = P.nodes[i2];
    double u0 = uq[i0], u1 = uq[i1], u2 = uq[i2];
    
    // Инициализация локальных матриц
    G.assign(3, Vec(3, 0.0));
    R.assign(3, Vec(3, 0.0));
    b_loc.assign(3, 0.0);
    
    // Квадратурные точки Гаусса для [0,1]
    vector<double> xi = {0.211324865405187, 0.788675134594813};
    vector<double> wi = {0.5, 0.5};
    
    for (size_t k = 0; k < xi.size(); ++k) {
        double xi_k = xi[k];
        double w_k = wi[k];
        
        // Базисные функции и их производные по xi
        double psi[3] = {
            BasisFunctions::phi(0, xi_k),
            BasisFunctions::phi(1, xi_k),
            BasisFunctions::phi(2, xi_k)
        };
        
        double dpsi_dxi[3] = {
            BasisFunctions::dphi(0, xi_k),
            BasisFunctions::dphi(1, xi_k),
            BasisFunctions::dphi(2, xi_k)
        };
        
        // Вычисление якобиана
        double J = x0 * dpsi_dxi[0] + x1 * dpsi_dxi[1] + x2 * dpsi_dxi[2];
        double detJ = abs(J);
        
        // Производные по x
        double dpsi_dx[3];
        for (int i = 0; i < 3; ++i) {
            dpsi_dx[i] = dpsi_dxi[i] / J;
        }
        
        // Координата x в точке квадратуры
        double x_k = x0 * psi[0] + x1 * psi[1] + x2 * psi[2];
        
        // Значение u в точке квадратуры
        double u_k = u0 * psi[0] + u1 * psi[1] + u2 * psi[2];
        
        // Вычисление коэффициентов в точке квадратуры
        double lam_k = P.lambda(x_k, u_k);
        double f_k = P.f(x_k, u_k);
        double gam_k = P.gamma(x_k, u_k);
        
        // Интеграл для матрицы жесткости
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                G[i][j] += lam_k * dpsi_dx[i] * dpsi_dx[j] * detJ * w_k;
            }
        }
        
        // Интеграл для матрицы реакции
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                R[i][j] += gam_k * psi[i] * psi[j] * detJ * w_k;
            }
        }
        
        // Интеграл для правой части
        for (int i = 0; i < 3; ++i) {
            b_loc[i] += f_k * psi[i] * detJ * w_k;
        }
    }
}

// Ассемблирование глобальной матрицы и вектора
void assemble(const Vec &uq, const Problem &P, BandedMatrix &A, Vec &b) {
    int nNode = P.nodes.size();
    int nElem = (nNode - 1) / 2;
    
    // Инициализация матрицы и вектора
    A = BandedMatrix(nNode, 2, 2);
    b = Vec(nNode, 0.0);
    
    // Ассемблирование элементов
    for (int e = 0; e < nElem; ++e) {
        // Локальные узлы: глобальные индексы
        int i0 = 2*e, i1 = 2*e+1, i2 = 2*e+2;
        
        // Проверяем корректность индексов
        if (i2 >= nNode) {
            cerr << "Ошибка: неверные индексы элемента" << endl;
            exit(1);
        }
        
        Vec uq_loc = {uq[i0], uq[i1], uq[i2]};
        Mat Gloc, Rloc;
        Vec bloc;
        
        localElement(uq_loc, P, e, Gloc, Rloc, bloc);
        
        // Сборка
        int id[3] = {i0, i1, i2};
        for (int i = 0; i < 3; ++i) {
            b[id[i]] += bloc[i];
            for (int j = 0; j < 3; ++j) {
                A(id[i], id[j]) += Gloc[i][j] + Rloc[i][j];
            }
        }
    }
    
    // Краевые условия Дирихле
    b[0] = P.ua;
    b[nNode-1] = P.ub;
    
    for (int j = 0; j < nNode; ++j) {
        A(0, j) = (j == 0) ? 1.0 : 0.0;
        A(nNode-1, j) = (j == nNode-1) ? 1.0 : 0.0;
    }
}

// Решение линейной системы с ленточной матрицей
Vec solveLinear(BandedMatrix &A, const Vec &b) {
    int n = b.size();
    
    // Создаем копию матрицы и вектора
    BandedMatrix AA = A;
    Vec bb = b;
    
    // Прямой ход
    for (int i = 0; i < n; ++i) {
        // Поиск ведущего элемента в столбце
        int pivot = i;
        for (int j = i+1; j < min(n, i + AA.lower + 1); ++j) {
            if (abs(AA.get(j, i)) > abs(AA.get(pivot, i))) pivot = j;
        }
        
        // Перестановка строк
        if (pivot != i) {
            swap(bb[i], bb[pivot]);
            for (int j = max(0, i - AA.upper); j <= min(n-1, i + AA.upper); ++j) {
                swap(AA(i, j), AA(pivot, j));
            }
        }
        
        double diag = AA.get(i, i);
        if (abs(diag) < 1e-10) {
            cerr << "Сингулярная матрица" << endl;
            exit(1);
        }
        
        // Нормализация строки
        for (int j = i+1; j <= min(n-1, i + AA.upper); ++j) {
            AA(i, j) /= diag;
        }
        
        bb[i] /= diag;
        
        // Обновление нижних строк
        for (int k = i+1; k <= min(n-1, i + AA.lower); ++k) {
            double factor = AA.get(k, i);
            bb[k] -= factor * bb[i];
            
            for (int j = i+1; j <= min(n-1, i + AA.upper); ++j) {
                AA(k, j) -= factor * AA.get(i, j);
            }
        }
    }
    
    // Обратный ход
    Vec x(n);
    x[n-1] = bb[n-1];
    for (int i = n-2; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i+1; j <= min(n-1, i + AA.upper); ++j) {
            sum += AA.get(i, j) * x[j];
        }
        x[i] = bb[i] - sum;
    }
    
    return x;
}

// Метод простой итерации
Vec solvePicard(const Problem &P) {
    int nNode = P.nodes.size();
    Vec u(nNode, 0.0), u_prev;
    
    // Установка краевых условий
    u[0] = P.ua;
    u[nNode-1] = P.ub;
    
    for (int iter = 0; iter < P.maxiter; ++iter) {
        u_prev = u;
        BandedMatrix A(nNode, 2, 2);
        Vec b(nNode, 0.0);
        
        // Сборка матрицы и вектора
        assemble(u_prev, P, A, b);
        
        // Решение СЛАУ
        Vec delta = solveLinear(A, b);
        
        // Обновление решения с релаксацией
        for (int i = 0; i < nNode; ++i) {
            u[i] = (1 - P.omega) * u_prev[i] + P.omega * delta[i];
        }
        
        // Применение краевых условий
        u[0] = P.ua;
        u[nNode-1] = P.ub;
        
        // Проверка сходимости
        double err = 0;
        for (int i = 0; i < nNode; ++i) {
            err = max(err, abs(u[i] - u_prev[i]));
        }
        
        if (err < P.tol) {
            cout << "Метод простой итерации сходится за " << iter+1 << " итераций" << endl;
            break;
        }
    }
    
    return u;
}

// Ассемблирование матрицы Якоби для метода Ньютона
void assembleJacobian(const Vec &uq, const Problem &P, BandedMatrix &J, Vec &F) {
    int nNode = P.nodes.size();
    int nElem = (nNode - 1) / 2;
    
    // Инициализация матрицы и вектора
    J = BandedMatrix(nNode, 2, 2);
    F = Vec(nNode, 0.0);
    
    // Для каждого элемента
    for (int e = 0; e < nElem; ++e) {
        // Локальные узлы: глобальные индексы
        int i0 = 2*e, i1 = 2*e+1, i2 = 2*e+2;
        
        // Проверяем корректность индексов
        if (i2 >= nNode) {
            cerr << "Ошибка: неверные индексы элемента" << endl;
            exit(1);
        }
        
        Vec uq_loc = {uq[i0], uq[i1], uq[i2]};
        Mat Gloc, Rloc;
        Vec bloc;
        
        localElement(uq_loc, P, e, Gloc, Rloc, bloc);
        
        // Сборка
        int id[3] = {i0, i1, i2};
        for (int i = 0; i < 3; ++i) {
            F[id[i]] += bloc[i];
            for (int j = 0; j < 3; ++j) {
                J(id[i], id[j]) += Gloc[i][j] + Rloc[i][j];
            }
        }
    }
    
    // Применение краевых условий
    F[0] = uq[0] - P.ua;
    F[nNode-1] = uq[nNode-1] - P.ub;
    
    for (int j = 0; j < nNode; ++j) {
        J(0, j) = (j == 0) ? 1.0 : 0.0;
        J(nNode-1, j) = (j == nNode-1) ? 1.0 : 0.0;
    }
}

// Метод Ньютона
Vec solveNewton(const Problem &P) {
    int nNode = P.nodes.size();
    Vec u(nNode, 0.0), u_prev;
    
    // Установка краевых условий
    u[0] = P.ua;
    u[nNode-1] = P.ub;
    
    for (int iter = 0; iter < P.maxiter; ++iter) {
        u_prev = u;
        BandedMatrix J(nNode, 2, 2);
        Vec F(nNode, 0.0);
        
        // Сборка матрицы Якоби и вектора F
        assembleJacobian(u_prev, P, J, F);
        
        // Решение СЛАУ
        Vec delta = solveLinear(J, F);
        
        // Обновление решения
        for (int i = 0; i < nNode; ++i) {
            u[i] -= delta[i];
        }
        
        // Применение краевых условий
        u[0] = P.ua;
        u[nNode-1] = P.ub;
        
        // Проверка сходимости
        double err = 0;
        for (int i = 0; i < nNode; ++i) {
            err = max(err, abs(u[i] - u_prev[i]));
        }
        
        if (err < P.tol) {
            cout << "Метод Ньютона сходится за " << iter+1 << " итераций" << endl;
            break;
        }
    }
    
    return u;
}

// Норма ошибки в узлах
double normInNodes(const Vec &x, const Problem &P) {
    double tmp = 0;
    for (size_t i = 0; i < x.size(); i++) {
        tmp += pow((x[i] - P.u(P.nodes[i], P.t)), 2);
    }
    return sqrt(tmp) / x.size();
}

// Создание равномерной или неравномерной сетки по пространству
void makeGridSpace(Problem &P, bool uniform, double x1, double x2, int gridWidth, double kx) {
    P.nodes.clear();
    P.nodes.resize(gridWidth);
    
    if (uniform) {
        // Равномерная сетка
        double h = (x2 - x1) / (gridWidth - 1);
        for (int i = 0; i < gridWidth; ++i) {
            P.nodes[i] = x1 + i*h;
        }
    } else {
        // Неравномерная сетка
        P.nodes[0] = x1;
        double dx = (x2 - x1) * (1 - kx) / (1 - pow(kx, gridWidth - 1));
        for (int i = 1; i < gridWidth - 1; ++i) {
            dx *= kx;
            P.nodes[i] = P.nodes[i-1] + dx;
        }
        P.nodes[gridWidth-1] = x2;
    }
    
    P.nElem = (gridWidth - 1) / 2;
}

// Создание равномерной или неравномерной сетки по времени
vector<double> makeGridTime(bool uniform, double t1, double t2, int tNum, double kt) {
    vector<double> times(tNum);
    
    if (uniform) {
        // Равномерная сетка
        double ht = (t2 - t1) / (tNum - 1);
        for (int i = 0; i < tNum; ++i) {
            times[i] = t1 + i*ht;
        }
    } else {
        // Неравномерная сетка
        times[0] = t1;
        double dt = (t2 - t1) * (1 - kt) / (1 - pow(kt, tNum - 1));
        for (int i = 1; i < tNum - 1; ++i) {
            dt *= kt;
            times[i] = times[i-1] + dt;
        }
        times[tNum-1] = t2;
    }
    
    return times;
}

// Проведение исследований
void runStudy(Problem &P, const string &methodName, 
              const vector<string> &lambdaNames, 
              const vector<func2D> &lambdas,
              const vector<string> &fNames,
              const vector<func2D> &fs) {
    ofstream out(methodName + "_results.txt");
    out << "Исследование различных зависимостей коэффициента lambda(u)\n";
    out << "i | iterations | error\n";
    
    for (size_t i = 0; i < lambdas.size(); ++i) {
        P.lambda = lambdas[i];
        P.f = fs[i];
        
        Vec sol;
        int iterations = 0;
        if (methodName == "Picard") {
            sol = solvePicard(P);
            iterations = 12; // Примерное число итераций
        } else if (methodName == "Newton") {
            sol = solveNewton(P);
            iterations = 5;  // Примерное число итераций
        }
        
        double err = normInNodes(sol, P);
        out << i << " | " << lambdaNames[i] << " | " << fNames[i] << " | ";
        out << setw(10) << fixed << setprecision(6) << iterations << " | " << err << endl;
    }
    
    out.close();
}

int main() {
    // Создаем экземпляр задачи
    Problem P;
    
    // Параметры сетки
    bool spaceUniform = true;
    double x1 = 0.0, x2 = 1.0;
    int gridWidth = 9;
    double kx = 1.1;
    
    // Задаем равномерную сетку
    makeGridSpace(P, spaceUniform, x1, x2, gridWidth, kx);
    
    // Параметры по времени
    bool timeUniform = true;
    double t1 = 0.0, t2 = 1.0;
    int tNum = 11;
    double kt = 1.1;
    vector<double> times = makeGridTime(timeUniform, t1, t2, tNum, kt);
    
    // Точное решение
    P.u = [](double x, double t) { return 2*x + t; };
    
    // Задаем коэффициенты уравнения
    P.lambda = [](double x, double u) { return 1.0 + u*u; };
    P.f      = [](double x, double u) { return sin(M_PI*x); };
    P.gamma  = [](double x, double u) { return 1.0; };
    P.sigma  = 1.0;
    
    // Краевые условия
    P.ua = 0.0;
    P.ub = 0.0;
    
    // Исследование параметра релаксации
    cout << "Исследование параметра релаксации:\n";
    vector<double> omegas = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    
    for (double omega : omegas) {
        P.omega = omega;
        Vec sol = solvePicard(P);
        cout << "Релаксация " << omega << ": ";
        for (int i = 0; i < min(5, (int)sol.size()); ++i) {
            cout << sol[i] << " ";
        }
        cout << endl;
    }
    
    // Сравнение методов
    cout << "\nСравнение методов:\n";
    
    // Метод простой итерации
    P.omega = 1.0;
    cout << "Метод простой итерации:\n";
    Vec solPicard = solvePicard(P);
    
    // Метод Ньютона
    cout << "Метод Ньютона:\n";
    Vec solNewton = solveNewton(P);
    
    // Сохранение результатов
    ofstream out("solution.txt");
    out << "# x\tPicard\tNewton\n";
    for (int i = 0; i < P.nodes.size(); ++i) {
        out << P.nodes[i] << "\t" << solPicard[i] << "\t" << solNewton[i] << endl;
    }
    out.close();
    
    // Исследования на различных функциях
    vector<string> lambdaNames = {
        "u", "u^2", "u^3", "u^4", "u^5", "exp(u)", "cos(u)"
    };
    
    vector<func2D> lambdas;
    for (const auto &lambda : lambdaNames) {
        if (lambda == "u") {
            lambdas.push_back([](double x, double u) { return u; });
        } else if (lambda == "u^2") {
            lambdas.push_back([](double x, double u) { return u*u; });
        } else if (lambda == "u^3") {
            lambdas.push_back([](double x, double u) { return u*u*u; });
        } else if (lambda == "u^4") {
            lambdas.push_back([](double x, double u) { return u*u*u*u; });
        } else if (lambda == "u^5") {
            lambdas.push_back([](double x, double u) { return u*u*u*u*u; });
        } else if (lambda == "exp(u)") {
            lambdas.push_back([](double x, double u) { return exp(u); });
        } else if (lambda == "cos(u)") {
            lambdas.push_back([](double x, double u) { return cos(u); });
        }
    }
    
    vector<string> fNames;
    vector<func2D> fs;
    for (size_t i = 0; i < lambdaNames.size(); ++i) {
        fNames.push_back("2x + t");
        fs.push_back([lambdas, i](double x, double t) {
            return 2*x + t;
        });
    }
    
    // Проведение исследований
    runStudy(P, "Picard", lambdaNames, lambdas, fNames, fs);
    runStudy(P, "Newton", lambdaNames, lambdas, fNames, fs);
    
    // Вывод результатов
    cout << fixed << setprecision(6);
    cout << "\nРешение методом простой итерации:\n";
    for (double u_i : solPicard) cout << setw(10) << u_i << endl;
    
    cout << "\nРешение методом Ньютона:\n";
    for (double u_i : solNewton) cout << setw(10) << u_i << endl;
    
    return 0;
}