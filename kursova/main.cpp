#define _USE_MATH_DEFINES
#include <cmath>
#include <functional>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
using namespace std;

// --- Структуры ---
struct Node {
    double x, y;
    int globalIndex;
};
struct RectElement {
    int nodeIds[9];
};

// --- Алиасы типов ---
using real = double;
using integer = int;
using Func3D = std::function<double(double, double, double)>;

// --- Генерация прямоугольной сетки ---
void GenerateRectGrid(int nx, int ny, real hx, real hy, vector<Node>& nodes, vector<RectElement>& elements) {
    int nodeCountX = 2 * nx + 1;
    int nodeCountY = 2 * ny + 1;
    nodes.clear();
    elements.clear();
    for (int j = 0; j < nodeCountY; ++j) {
        for (int i = 0; i < nodeCountX; ++i) {
            Node node;
            node.x = i * (hx / (2.0 * nx));
            node.y = j * (hy / (2.0 * ny));
            node.globalIndex = j * nodeCountX + i;
            nodes.push_back(node);
        }
    }
    for (int ey = 0; ey < ny; ++ey) {
        for (int ex = 0; ex < nx; ++ex) {
            int i0 = 2 * ex;
            int j0 = 2 * ey;
            int base = j0 * nodeCountX + i0;
            RectElement elem;
            elem.nodeIds[0] = base;
            elem.nodeIds[1] = base + 1;
            elem.nodeIds[2] = base + 2;
            elem.nodeIds[3] = base + nodeCountX;
            elem.nodeIds[4] = base + nodeCountX + 1;
            elem.nodeIds[5] = base + nodeCountX + 2;
            elem.nodeIds[6] = base + 2 * nodeCountX;
            elem.nodeIds[7] = base + 2 * nodeCountX + 1;
            elem.nodeIds[8] = base + 2 * nodeCountX + 2;
            elements.push_back(elem);
        }
    }
}

// --- Одномерные лагранжевые полиномы для Q2 ---
inline double L0(double x) { return 0.5 * x * (x - 1.0); } // x=-1
inline double L1(double x) { return (1.0 - x) * (1.0 + x); } // x=0
inline double L2(double x) { return 0.5 * x * (x + 1.0); } // x=+1
inline double dL0(double x) { return x - 0.5; }
inline double dL1(double x) { return -2.0 * x; }
inline double dL2(double x) { return x + 0.5; }

// --- Биквадратичные базисные функции ---
void BiquadraticBasis(double xi, double eta, double N[9]) {
    N[0] = L0(xi) * L0(eta);
    N[1] = L1(xi) * L0(eta);
    N[2] = L2(xi) * L0(eta);
    N[3] = L0(xi) * L1(eta);
    N[4] = L1(xi) * L1(eta);
    N[5] = L2(xi) * L1(eta);
    N[6] = L0(xi) * L2(eta);
    N[7] = L1(xi) * L2(eta);
    N[8] = L2(xi) * L2(eta);
}

void BiquadraticBasis_dXi(double xi, double eta, double dN[9]) {
    dN[0] = dL0(xi) * L0(eta);
    dN[1] = dL1(xi) * L0(eta);
    dN[2] = dL2(xi) * L0(eta);
    dN[3] = dL0(xi) * L1(eta);
    dN[4] = dL1(xi) * L1(eta);
    dN[5] = dL2(xi) * L1(eta);
    dN[6] = dL0(xi) * L2(eta);
    dN[7] = dL1(xi) * L2(eta);
    dN[8] = dL2(xi) * L2(eta);
}

void BiquadraticBasis_dEta(double xi, double eta, double dN[9]) {
    dN[0] = L0(xi) * dL0(eta);
    dN[1] = L1(xi) * dL0(eta);
    dN[2] = L2(xi) * dL0(eta);
    dN[3] = L0(xi) * dL1(eta);
    dN[4] = L1(xi) * dL1(eta);
    dN[5] = L2(xi) * dL1(eta);
    dN[6] = L0(xi) * dL2(eta);
    dN[7] = L1(xi) * dL2(eta);
    dN[8] = L2(xi) * dL2(eta);
}

// --- Квадратурные точки Гаусса ---
const double gauss_pts[3] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
const double gauss_wts[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};

// --- Сборка локальных матриц ---
void AssembleLocalMatrices(const Node nodes[9], double lambda[9], double Ke[9][9], double Me[9][9]) {
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            Ke[i][j] = Me[i][j] = 0.0;
    for (int gp = 0; gp < 3; ++gp) {
        double xi = gauss_pts[gp];
        double wx = gauss_wts[gp];
        for (int gq = 0; gq < 3; ++gq) {
            double eta = gauss_pts[gq];
            double wy = gauss_wts[gq];
            double N[9], dN_dxi[9], dN_deta[9];
            BiquadraticBasis(xi, eta, N);
            BiquadraticBasis_dXi(xi, eta, dN_dxi);
            BiquadraticBasis_dEta(xi, eta, dN_deta);
            double dx_dxi = 0, dx_deta = 0, dy_dxi = 0, dy_deta = 0;
            for (int k = 0; k < 9; ++k) {
                dx_dxi += nodes[k].x * dN_dxi[k];
                dx_deta += nodes[k].x * dN_deta[k];
                dy_dxi += nodes[k].y * dN_dxi[k];
                dy_deta += nodes[k].y * dN_deta[k];
            }
            double J = dx_dxi * dy_deta - dx_deta * dy_dxi;
            if (fabs(J) < 1e-12) throw runtime_error("Singular Jacobian");
            double dN_dx[9], dN_dy[9];
            for (int k = 0; k < 9; ++k) {
                dN_dx[k] = (dN_dxi[k] * dy_deta - dN_deta[k] * dy_dxi) / J;
                dN_dy[k] = (-dN_dxi[k] * dx_deta + dN_deta[k] * dx_dxi) / J;
            }
            double lmb = 0.0;
            for (int k = 0; k < 9; ++k) lmb += lambda[k] * N[k];
            for (int i = 0; i < 9; ++i) {
                for (int j = 0; j < 9; ++j) {
                    Ke[i][j] += lmb * (dN_dx[i] * dN_dx[j] + dN_dy[i] * dN_dy[j]) * wx * wy * fabs(J);
                    Me[i][j] += N[i] * N[j] * wx * wy * fabs(J);
                }
            }
        }
    }
}

// --- Сборка глобальных матриц (плотный формат) ---
void AssembleGlobalSystem(const vector<Node>& nodes, const vector<RectElement>& elements,
                         const vector<double>& lambda, vector<vector<double>>& K, vector<vector<double>>& M) {
    int N = nodes.size();
    K.assign(N, vector<double>(N, 0.0));
    M.assign(N, vector<double>(N, 0.0));
    for (const auto& elem : elements) {
        Node elemNodes[9];
        double elemLambda[9];
        for (int i = 0; i < 9; ++i) {
            elemNodes[i] = nodes[elem.nodeIds[i]];
            elemLambda[i] = lambda[elem.nodeIds[i]];
        }
        double Ke[9][9], Me[9][9];
        AssembleLocalMatrices(elemNodes, elemLambda, Ke, Me);
        for (int i = 0; i < 9; ++i) {
            int gi = elem.nodeIds[i];
            for (int j = 0; j < 9; ++j) {
                int gj = elem.nodeIds[j];
                K[gi][gj] += Ke[i][j];
                M[gi][gj] += Me[i][j];
            }
        }
    }
}

// --- Поиск граничных узлов ---
vector<int> FindBoundaryNodes(const vector<Node>& nodes, int nx, int ny) {
    int nodeCountX = 2 * nx + 1;
    vector<int> boundary;
    for (int j = 0; j <= 2 * ny; ++j) {
        boundary.push_back(j * nodeCountX); // левая граница
        boundary.push_back(j * nodeCountX + 2 * nx); // правая граница
    }
    for (int i = 1; i < 2 * nx; ++i) {
        boundary.push_back(i); // нижняя граница
        boundary.push_back(2 * ny * nodeCountX + i); // верхняя граница
    }
    sort(boundary.begin(), boundary.end());
    boundary.erase(unique(boundary.begin(), boundary.end()), boundary.end());
    return boundary;
}

// --- Применение граничных условий Дирихле ---
void ApplyDirichletBC(vector<vector<double>>& A, vector<double>& F,
                      const vector<int>& boundaryNodes, const vector<double>& values) {
    for (size_t idx = 0; idx < boundaryNodes.size(); ++idx) {
        int node = boundaryNodes[idx];
        double val = values[idx];
        for (int j = 0; j < A.size(); ++j) {
            A[node][j] = 0.0;
            A[j][node] = 0.0;
        }
        A[node][node] = 1.0;
        F[node] = val;
    }
}

// --- Сборка правой части ---
void AssembleGlobalRHS(const vector<Node>& nodes, const vector<RectElement>& elements,
                       const Func3D& f, double t, vector<double>& F) {
    int N = nodes.size();
    F.assign(N, 0.0);
    for (const auto& elem : elements) {
        Node elemNodes[9];
        for (int i = 0; i < 9; ++i) elemNodes[i] = nodes[elem.nodeIds[i]];
        double Fe[9] = {0.0};
        for (int gp = 0; gp < 3; ++gp) {
            double xi = gauss_pts[gp];
            double wx = gauss_wts[gp];
            for (int gq = 0; gq < 3; ++gq) {
                double eta = gauss_pts[gq];
                double wy = gauss_wts[gq];
                double N[9];
                BiquadraticBasis(xi, eta, N);
                double x = 0, y = 0;
                for (int k = 0; k < 9; ++k) {
                    x += elemNodes[k].x * N[k];
                    y += elemNodes[k].y * N[k];
                }
                double dx_dxi = 0, dx_deta = 0, dy_dxi = 0, dy_deta = 0;
                double dN_dxi_arr[9], dN_deta_arr[9];
                BiquadraticBasis_dXi(xi, eta, dN_dxi_arr);
                BiquadraticBasis_dEta(xi, eta, dN_deta_arr);
                for (int k = 0; k < 9; ++k) {
                    dx_dxi += elemNodes[k].x * dN_dxi_arr[k];
                    dx_deta += elemNodes[k].x * dN_deta_arr[k];
                    dy_dxi += elemNodes[k].y * dN_dxi_arr[k];
                    dy_deta += elemNodes[k].y * dN_deta_arr[k];
                }
                double J = dx_dxi * dy_deta - dx_deta * dy_dxi;
                double fval = f(x, y, t);
                for (int i = 0; i < 9; ++i)
                    Fe[i] += N[i] * fval * wx * wy * fabs(J);
            }
        }
        for (int i = 0; i < 9; ++i)
            F[elem.nodeIds[i]] += Fe[i];
    }
}

// --- Итерационный решатель: Метод сопряжённых градиентов (CG) ---
vector<double> SolveCG(const vector<vector<double>>& A, const vector<double>& b, double tol=1e-10, int maxit=10000) {
    int N = A.size();
    vector<double> x(N, 0.0), r = b, p = b, Ap(N);
    double rsold = 0.0;
    for (int i = 0; i < N; ++i) rsold += r[i] * r[i];
    for (int it = 0; it < maxit; ++it) {
        // Ap = A*p
        for (int i = 0; i < N; ++i) {
            Ap[i] = 0.0;
            for (int j = 0; j < N; ++j) Ap[i] += A[i][j] * p[j];
        }
        double alpha = rsold;
        double pAp = 0.0;
        for (int i = 0; i < N; ++i) pAp += p[i] * Ap[i];
        if (fabs(pAp) < 1e-20) break;
        alpha /= pAp;
        for (int i = 0; i < N; ++i) x[i] += alpha * p[i];
        for (int i = 0; i < N; ++i) r[i] -= alpha * Ap[i];
        double rsnew = 0.0;
        for (int i = 0; i < N; ++i) rsnew += r[i] * r[i];
        if (sqrt(rsnew) < tol) break;
        for (int i = 0; i < N; ++i) p[i] = r[i] + (rsnew / rsold) * p[i];
        rsold = rsnew;
    }
    return x;
}
// --- Решатель СЛАУ: выбирает CG для симметричных матриц, иначе Гаусс ---
vector<double> SolveLinearSystem(const vector<vector<double>>& A, const vector<double>& b) {
    int N = A.size();
    // Проверка на симметричность
    bool isSym = true;
    for (int i = 0; i < N && isSym; ++i)
        for (int j = 0; j < N; ++j)
            if (fabs(A[i][j] - A[j][i]) > 1e-12) { isSym = false; break; }
    if (isSym) {
        // Используем МСГ
        return SolveCG(A, b);
    }
    vector<vector<double>> A_copy = A;
    vector<double> b_copy = b;
    vector<double> x(N);
    for (int k = 0; k < N; ++k) {
        int max_row = k;
        double max_val = fabs(A_copy[k][k]);
        for (int i = k + 1; i < N; ++i) {
            if (fabs(A_copy[i][k]) > max_val) {
                max_val = fabs(A_copy[i][k]);
                max_row = i;
            }
        }
        if (max_val < 1e-12) throw runtime_error("Matrix is singular");
        swap(A_copy[k], A_copy[max_row]);
        swap(b_copy[k], b_copy[max_row]);
        for (int i = k + 1; i < N; ++i) {
            double factor = A_copy[i][k] / A_copy[k][k];
            for (int j = k; j < N; ++j) {
                A_copy[i][j] -= factor * A_copy[k][j];
            }
            b_copy[i] -= factor * b_copy[k];
        }
    }
    for (int i = N - 1; i >= 0; --i) {
        x[i] = b_copy[i];
        for (int j = i + 1; j < N; ++j) {
            x[i] -= A_copy[i][j] * x[j];
        }
        x[i] /= A_copy[i][i];
    }
    return x;
}

// --- Трехслойная схема ---
void SolveThreeLayerScheme(const vector<vector<double>>& K, const vector<vector<double>>& M,
                           const vector<double>& F, const vector<double>& u_prev,
                           const vector<double>& u_curr, vector<double>& u_next, double dt,
                           const vector<int>& boundary, const vector<double>& bc_values) {
    int N = K.size();
    vector<double> rhs(N, 0.0);
    // Обнуляем вклады u_curr и u_prev для граничных узлов
    vector<double> u_curr_adj = u_curr;
    vector<double> u_prev_adj = u_prev;
    for (int node : boundary) {
        u_curr_adj[node] = 0.0;
        u_prev_adj[node] = 0.0;
    }
    // Формируем правую часть
    for (int i = 0; i < N; ++i) {
        rhs[i] = F[i];
        for (int j = 0; j < N; ++j) {
            rhs[i] += M[i][j] * (2.0 * u_curr_adj[j] - u_prev_adj[j]) / (dt * dt);
        }
    }
    vector<vector<double>> A(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = M[i][j] / (dt * dt) + K[i][j];
        }
    }
    ApplyDirichletBC(A, rhs, boundary, bc_values);
    u_next = SolveLinearSystem(A, rhs);
    for (double val : u_next) {
        if (!isfinite(val)) throw runtime_error("NaN or inf in solution");
    }
}

// --- Основной расчетный цикл ---
void RunThreeLayerFEM_Full() {
    int nx = 10, ny = 10;
    double hx = 1.0, hy = 1.0;
    vector<Node> nodes;
    vector<RectElement> elements;
    GenerateRectGrid(nx, ny, hx, hy, nodes, elements);
    int N = nodes.size();
    int NE = elements.size();
    ofstream f_nodes("nodes.txt");
    for (const auto& node : nodes) {
        f_nodes << node.globalIndex << " " << node.x << " " << node.y << "\n";
    }
    f_nodes.close();
    ofstream f_elem("elements.txt");
    for (int i = 0; i < NE; ++i) {
        f_elem << i;
        for (int j = 0; j < 9; ++j) f_elem << " " << elements[i].nodeIds[j];
        f_elem << "\n";
    }
    f_elem.close();
    vector<double> lambda(N, 1.0);
    vector<vector<double>> K(N, vector<double>(N, 0.0));
    vector<vector<double>> M(N, vector<double>(N, 0.0));
    AssembleGlobalSystem(nodes, elements, lambda, K, M);
    Func3D f = [](double x, double y, double t) {
        double pi2 = M_PI * M_PI;
        return sin(M_PI * x) * sin(M_PI * y) * (-cos(t) + 2 * pi2 * cos(t));
    };
    Func3D initial = [](double x, double y, double t) {
        return sin(M_PI * x) * sin(M_PI * y) * cos(t);
    };
    vector<int> boundary = FindBoundaryNodes(nodes, nx, ny);
    vector<double> F(N);
    double T = 1.0;
    int Nt = 10;
    double dt = T / Nt;
    vector<vector<double>> U(Nt + 1, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) {
        U[0][i] = initial(nodes[i].x, nodes[i].y, 0.0);
        U[1][i] = initial(nodes[i].x, nodes[i].y, dt);
    }
    for (int n = 1; n < Nt; ++n) {
        double t = (n + 1) * dt;
        AssembleGlobalRHS(nodes, elements, f, t, F);
        vector<double> bc_values(boundary.size(), 0.0); // Нулевые граничные условия
        SolveThreeLayerScheme(K, M, F, U[n-1], U[n], U[n+1], dt, boundary, bc_values);
    }
    ofstream fout("output.txt");
    for (int i = 0; i < N; ++i) {
        for (int n = 0; n <= Nt; ++n) {
            fout << setw(16) << U[n][i];
        }
        fout << "\n";
    }
    fout.close();
}

// --- Тест на порядок аппроксимации по времени (полиномиальные функции) ---
void TestOrderTimePoly() {
    Func3D u_exact = [](double x, double y, double t) {
        return sin(M_PI * x) * sin(M_PI * y) * pow(t, 3);
    };
    Func3D u_t = [](double x, double y, double t) {
        return 3 * sin(M_PI * x) * sin(M_PI * y) * pow(t, 2);
    };
    Func3D u_tt = [](double x, double y, double t) {
        return 6 * sin(M_PI * x) * sin(M_PI * y) * t;
    };
    Func3D f_rhs = [](double x, double y, double t) {
        double pi2 = M_PI * M_PI;
        return sin(M_PI * x) * sin(M_PI * y) * (6 * t + 2 * pi2 * pow(t, 3));
    };
    cout << "[Test 3.1] Таблица порядка аппроксимации по времени (полиномиальные функции):\n";
    cout << "Nt      rel_error       order\n";
    vector<int> Nt_values = {8, 16, 32, 64, 128, 256};
    vector<double> errors;
    for (int Nt : Nt_values) {
        int nx = 10, ny = 10; // для чистоты порядка по времени
        double hx = 1.0, hy = 1.0;
        vector<Node> nodes;
        vector<RectElement> elements;
        GenerateRectGrid(nx, ny, hx, hy, nodes, elements);
        int N = nodes.size();
        vector<double> lambda(N, 1.0);
        vector<vector<double>> K(N, vector<double>(N, 0.0));
        vector<vector<double>> M(N, vector<double>(N, 0.0));
        AssembleGlobalSystem(nodes, elements, lambda, K, M);
        // Масслимпинг ОТКЛЮЧЁН для теста порядка по времени
        // for (int i = 0; i < N; ++i)
        //     for (int j = 0; j < N; ++j)
        //         if (i != j) M[i][j] = 0.0;
        vector<int> boundary = FindBoundaryNodes(nodes, nx, ny);
        double T = 1.0;
        double dt = T / Nt;
        vector<double> u_prev(N), u_curr(N), u_next(N);
        for (int i = 0; i < N; ++i) {
            u_prev[i] = u_exact(nodes[i].x, nodes[i].y, 0.0);
            u_curr[i] = u_prev[i] + dt * u_t(nodes[i].x, nodes[i].y, 0.0) + 0.5 * dt * dt * u_tt(nodes[i].x, nodes[i].y, 0.0);
        }
        // Жёстко задаём граничные условия на начальных слоях
        for (int b : boundary) {
            u_prev[b] = 0.0;
            u_curr[b] = 0.0;
        }
        vector<vector<double>> U(Nt + 1, vector<double>(N));
        U[0] = u_prev;
        U[1] = u_curr;
        for (int n = 1; n < Nt; ++n) {
            double t = (n + 1) * dt;
            vector<double> F(N);
            AssembleGlobalRHS(nodes, elements, f_rhs, t, F); // F_{n+1}
            vector<double> bc_values(boundary.size(), 0.0); // Нулевые граничные условия
            SolveThreeLayerScheme(K, M, F, U[n-1], U[n], U[n+1], dt, boundary, bc_values);
        }
        double err = 0, norm = 0;
        for (int i = 0; i < N; ++i) {
            double uex = u_exact(nodes[i].x, nodes[i].y, T);
            err += pow(U[Nt][i] - uex, 2);
            norm += pow(uex, 2);
        }
        err = sqrt(err / N);
        norm = sqrt(norm / N);
        double rel_err = (norm > 1e-10) ? err / norm : err;
        errors.push_back(rel_err);
        double order = 0.0;
        if (errors.size() > 1) {
            order = log(errors[errors.size()-2] / errors[errors.size()-1]) / log(2.0);
        }
        cout << setw(8) << Nt << setw(16) << rel_err << setw(16) << order << endl;
    }
}

// --- Тест на порядок сходимости по времени (неполиномиальные функции) ---
void TestOrderTimeNonPoly() {
    Func3D u_exact = [](double x, double y, double t) {
        return exp(-t) * sin(M_PI * x) * sin(M_PI * y);
    };
    Func3D u_t = [](double x, double y, double t) {
        return -exp(-t) * sin(M_PI * x) * sin(M_PI * y);
    };
    Func3D u_tt = [](double x, double y, double t) {
        return exp(-t) * sin(M_PI * x) * sin(M_PI * y);
    };
    Func3D f_rhs = [](double x, double y, double t) {
        double pi2 = M_PI * M_PI;
        return exp(-t) * sin(M_PI * x) * sin(M_PI * y) * (1 + 2 * pi2);
    };
    cout << "\n[Test 3.2] Таблица порядка сходимости по времени (неполиномиальные функции):\n";
    cout << "Nt      rel_error       order\n";
    vector<int> Nt_values = {8, 16, 32, 64, 128, 256};
    vector<double> errors;
    for (int Nt : Nt_values) {
        int nx = 10, ny = 10; // для чистоты порядка по времени
        double hx = 1.0, hy = 1.0;
        vector<Node> nodes;
        vector<RectElement> elements;
        GenerateRectGrid(nx, ny, hx, hy, nodes, elements);
        int N = nodes.size();
        vector<double> lambda(N, 1.0);
        vector<vector<double>> K(N, vector<double>(N, 0.0));
        vector<vector<double>> M(N, vector<double>(N, 0.0));
        AssembleGlobalSystem(nodes, elements, lambda, K, M);
        // Масслимпинг ОТКЛЮЧЁН для теста порядка по времени
        // for (int i = 0; i < N; ++i)
        //     for (int j = 0; j < N; ++j)
        //         if (i != j) M[i][j] = 0.0;
        vector<int> boundary = FindBoundaryNodes(nodes, nx, ny);
        double T = 0.2; // для устойчивости
        double dt = T / Nt;
        vector<double> u_prev(N), u_curr(N), u_next(N);
        for (int i = 0; i < N; ++i) {
            u_prev[i] = u_exact(nodes[i].x, nodes[i].y, 0.0);
            u_curr[i] = u_prev[i] + dt * u_t(nodes[i].x, nodes[i].y, 0.0) + 0.5 * dt * dt * u_tt(nodes[i].x, nodes[i].y, 0.0);
        }
        for (int b : boundary) {
            u_prev[b] = 0.0;
            u_curr[b] = 0.0;
        }
        vector<vector<double>> U(Nt + 1, vector<double>(N));
        U[0] = u_prev;
        U[1] = u_curr;
        for (int n = 1; n < Nt; ++n) {
            double t = (n + 1) * dt;
            vector<double> F(N);
            AssembleGlobalRHS(nodes, elements, f_rhs, t, F); // F_{n+1}
            vector<double> bc_values(boundary.size(), 0.0); // Нулевые граничные условия
            SolveThreeLayerScheme(K, M, F, U[n-1], U[n], U[n+1], dt, boundary, bc_values);
        }
        double err = 0, norm = 0;
        for (int i = 0; i < N; ++i) {
            double uex = u_exact(nodes[i].x, nodes[i].y, T);
            err += pow(U[Nt][i] - uex, 2);
            norm += pow(uex, 2);
        }
        err = sqrt(err / N);
        norm = sqrt(norm / N);
        double rel_err = (norm > 1e-10) ? err / norm : err;
        errors.push_back(rel_err);
        double order = 0.0;
        if (errors.size() > 1) {
            order = log(errors[errors.size()-2] / errors[errors.size()-1]) / log(2.0);
        }
        cout << setw(8) << Nt << setw(16) << rel_err << setw(16) << order << endl;
    }
}

// --- Диагностический тест 1: чисто временная задача (K=0) ---
void TestTimeOnly() {
    cout << "\n[DiagTest 1] Чисто временная задача (K=0):\n";
    // u_tt = f(t), u(0)=1, u_t(0)=0, f(t)=0, u(t)=1
    int N = 1;
    vector<vector<double>> K(N, vector<double>(N, 0.0));
    vector<vector<double>> M(N, vector<double>(N, 1.0));
    vector<double> u_prev(N, 1.0), u_curr(N, 1.0), u_next(N);
    vector<int> boundary = {0};
    vector<double> bc_values = {1.0};
    double T = 1.0;
    vector<int> Nt_values = {8, 16, 32, 64, 128, 256};
    vector<double> errors;
    for (int Nt : Nt_values) {
        double dt = T / Nt;
        u_prev[0] = 1.0;
        u_curr[0] = 1.0;
        for (int n = 1; n < Nt; ++n) {
            vector<double> F(N, 0.0);
            SolveThreeLayerScheme(K, M, F, u_prev, u_curr, u_next, dt, boundary, bc_values);
            u_prev = u_curr;
            u_curr = u_next;
        }
        double err = fabs(u_curr[0] - 1.0);
        errors.push_back(err);
        double order = 0.0;
        if (errors.size() > 1)
            order = log(errors[errors.size()-2]/errors[errors.size()-1])/log(2.0);
        cout << "Nt=" << setw(4) << Nt << "  error=" << scientific << err << "  order=" << order << endl;
    }
}

// --- Диагностический тест 2: постоянное решение (u=const, f=0) ---
void TestConstantSolution() {
    cout << "\n[DiagTest 2] Постоянное решение (u=const, f=0):\n";
    int nx = 2, ny = 2;
    double hx = 1.0, hy = 1.0;
    vector<Node> nodes;
    vector<RectElement> elements;
    GenerateRectGrid(nx, ny, hx, hy, nodes, elements);
    int N = nodes.size();
    vector<double> lambda(N, 1.0);
    vector<vector<double>> K(N, vector<double>(N, 0.0));
    vector<vector<double>> M(N, vector<double>(N, 0.0));
    AssembleGlobalSystem(nodes, elements, lambda, K, M);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            if (i != j) M[i][j] = 0.0;
    vector<int> boundary = FindBoundaryNodes(nodes, nx, ny);
    double T = 1.0;
    vector<int> Nt_values = {8, 16, 32, 64, 128, 256};
    vector<double> errors;
    for (int Nt : Nt_values) {
        double dt = T / Nt;
        vector<double> u_prev(N, 2.0), u_curr(N, 2.0), u_next(N);
        for (int b : boundary) { u_prev[b]=2.0; u_curr[b]=2.0; }
        for (int n = 1; n < Nt; ++n) {
            vector<double> F(N, 0.0);
            vector<double> bc_values(boundary.size(), 2.0);
            SolveThreeLayerScheme(K, M, F, u_prev, u_curr, u_next, dt, boundary, bc_values);
            u_prev = u_curr;
            u_curr = u_next;
        }
        double err = 0;
        for (int i = 0; i < N; ++i) err += pow(u_curr[i]-2.0,2);
        err = sqrt(err/N);
        errors.push_back(err);
        double order = 0.0;
        if (errors.size() > 1)
            order = log(errors[errors.size()-2]/errors[errors.size()-1])/log(2.0);
        cout << "Nt=" << setw(4) << Nt << "  error=" << scientific << err << "  order=" << order << endl;
    }
}

// --- Диагностический тест 3: билинейный базис (Q1) ---
void AssembleLocalMatricesQ1(const Node nodes[4], double lambda[4], double Ke[4][4], double Me[4][4]) {
    // Q1: 4 узла, стандартные формулы
    const double xi_pts[2] = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};
    const double wts[2] = {1.0, 1.0};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            Ke[i][j] = Me[i][j] = 0.0;
    for (int gp = 0; gp < 2; ++gp) {
        double xi = xi_pts[gp];
        double wx = wts[gp];
        for (int gq = 0; gq < 2; ++gq) {
            double eta = xi_pts[gq];
            double wy = wts[gq];
            double N[4] = {0.25*(1-xi)*(1-eta), 0.25*(1+xi)*(1-eta), 0.25*(1+xi)*(1+eta), 0.25*(1-xi)*(1+eta)};
            double dN_dxi[4] = {-0.25*(1-eta), 0.25*(1-eta), 0.25*(1+eta), -0.25*(1+eta)};
            double dN_deta[4] = {-0.25*(1-xi), -0.25*(1+xi), 0.25*(1+xi), 0.25*(1-xi)};
            double dx_dxi=0, dx_deta=0, dy_dxi=0, dy_deta=0;
            for (int k=0;k<4;++k) {
                dx_dxi += nodes[k].x * dN_dxi[k];
                dx_deta += nodes[k].x * dN_deta[k];
                dy_dxi += nodes[k].y * dN_dxi[k];
                dy_deta += nodes[k].y * dN_deta[k];
            }
            double J = dx_dxi*dy_deta - dx_deta*dy_dxi;
            double dN_dx[4], dN_dy[4];
            for (int k=0;k<4;++k) {
                dN_dx[k] = (dN_dxi[k]*dy_deta - dN_deta[k]*dy_dxi)/J;
                dN_dy[k] = (-dN_dxi[k]*dx_deta + dN_deta[k]*dx_dxi)/J;
            }
            double lmb = 0.0;
            for (int k=0;k<4;++k) lmb += lambda[k]*N[k];
            for (int i=0;i<4;++i)
                for (int j=0;j<4;++j) {
                    Ke[i][j] += lmb*(dN_dx[i]*dN_dx[j]+dN_dy[i]*dN_dy[j])*wx*wy*fabs(J);
                    Me[i][j] += N[i]*N[j]*wx*wy*fabs(J);
                }
        }
    }
}
void GenerateRectGridQ1(int nx, int ny, double hx, double hy, vector<Node>& nodes, vector<RectElement>& elements) {
    int nodeCountX = nx+1, nodeCountY = ny+1;
    nodes.clear(); elements.clear();
    for (int j=0;j<nodeCountY;++j)
        for (int i=0;i<nodeCountX;++i) {
            Node node; node.x = i*hx/nx; node.y = j*hy/ny; node.globalIndex = j*nodeCountX+i; nodes.push_back(node);
        }
    for (int ey=0;ey<ny;++ey)
        for (int ex=0;ex<nx;++ex) {
            int i0=ex, j0=ey, base=j0*nodeCountX+i0;
            RectElement elem;
            elem.nodeIds[0]=base; elem.nodeIds[1]=base+1; elem.nodeIds[2]=base+nodeCountX+1; elem.nodeIds[3]=base+nodeCountX;
            elements.push_back(elem);
        }
}
void AssembleGlobalSystemQ1(const vector<Node>& nodes, const vector<RectElement>& elements, const vector<double>& lambda, vector<vector<double>>& K, vector<vector<double>>& M) {
    int N=nodes.size(); K.assign(N,vector<double>(N,0.0)); M.assign(N,vector<double>(N,0.0));
    for (const auto& elem:elements) {
        Node elemNodes[4]; double elemLambda[4];
        for (int i=0;i<4;++i) { elemNodes[i]=nodes[elem.nodeIds[i]]; elemLambda[i]=lambda[elem.nodeIds[i]]; }
        double Ke[4][4], Me[4][4];
        AssembleLocalMatricesQ1(elemNodes, elemLambda, Ke, Me);
        for (int i=0;i<4;++i) for (int j=0;j<4;++j) {
            K[elem.nodeIds[i]][elem.nodeIds[j]] += Ke[i][j];
            M[elem.nodeIds[i]][elem.nodeIds[j]] += Me[i][j];
        }
    }
}
vector<int> FindBoundaryNodesQ1(const vector<Node>& nodes, int nx, int ny) {
    int nodeCountX=nx+1; vector<int> boundary;
    for (int j=0;j<=ny;++j) { boundary.push_back(j*nodeCountX); boundary.push_back(j*nodeCountX+nx); }
    for (int i=1;i<nx;++i) { boundary.push_back(i); boundary.push_back(ny*nodeCountX+i); }
    sort(boundary.begin(),boundary.end()); boundary.erase(unique(boundary.begin(),boundary.end()),boundary.end());
    return boundary;
}
void TestQ1Basis() {
    cout << "\n[DiagTest 3] Билинейный базис (Q1):\n";
    auto u_exact = [](double x,double y,double t){return sin(M_PI*x)*sin(M_PI*y)*pow(t,3);};
    auto u_t = [](double x,double y,double t){return 3*sin(M_PI*x)*sin(M_PI*y)*pow(t,2);};
    auto u_tt = [](double x,double y,double t){return 6*sin(M_PI*x)*sin(M_PI*y)*t;};
    auto f_rhs = [](double x,double y,double t){double pi2=M_PI*M_PI;return sin(M_PI*x)*sin(M_PI*y)*(6*t+2*pi2*pow(t,3));};
    vector<int> Nt_values = {8,16,32,64,128,256};
    vector<double> errors;
    for (int Nt:Nt_values) {
        int nx=10, ny=10; double hx=1.0, hy=1.0;
        vector<Node> nodes; vector<RectElement> elements;
        GenerateRectGridQ1(nx,ny,hx,hy,nodes,elements);
        int N=nodes.size(); vector<double> lambda(N,1.0);
        vector<vector<double>> K, M;
        AssembleGlobalSystemQ1(nodes,elements,lambda,K,M);
        for (int i=0;i<N;++i) for (int j=0;j<N;++j) if (i!=j) M[i][j]=0.0;
        vector<int> boundary=FindBoundaryNodesQ1(nodes,nx,ny);
        double T=1.0, dt=T/Nt;
        vector<double> u_prev(N), u_curr(N), u_next(N);
        for (int i=0;i<N;++i) {
            u_prev[i]=u_exact(nodes[i].x,nodes[i].y,0.0);
            u_curr[i]=u_prev[i]+dt*u_t(nodes[i].x,nodes[i].y,0.0)+0.5*dt*dt*u_tt(nodes[i].x,nodes[i].y,0.0);
        }
        for (int b:boundary) {u_prev[b]=0.0;u_curr[b]=0.0;}
        vector<vector<double>> U(Nt+1,vector<double>(N)); U[0]=u_prev; U[1]=u_curr;
        for (int n=1;n<Nt;++n) {
            double t=(n+1)*dt;
            vector<double> F(N);
            for (int i=0;i<N;++i) F[i]=0.0;
            for (const auto& elem:elements) {
                Node elemNodes[4]; for (int i=0;i<4;++i) elemNodes[i]=nodes[elem.nodeIds[i]];
                double Fe[4]={0.0};
                for (int gp=0;gp<2;++gp) {
                    double xi=-1.0/sqrt(3.0)+(gp)*2.0/sqrt(3.0);
                    for (int gq=0;gq<2;++gq) {
                        double eta=-1.0/sqrt(3.0)+(gq)*2.0/sqrt(3.0);
                        double Nq[4]={0.25*(1-xi)*(1-eta),0.25*(1+xi)*(1-eta),0.25*(1+xi)*(1+eta),0.25*(1-xi)*(1+eta)};
                        double x=0,y=0;
                        for (int k=0;k<4;++k) {x+=elemNodes[k].x*Nq[k];y+=elemNodes[k].y*Nq[k];}
                        double J=0.25;
                        double fval=f_rhs(x,y,t);
                        for (int i=0;i<4;++i) Fe[i]+=Nq[i]*fval*J;
                    }
                }
                for (int i=0;i<4;++i) F[elem.nodeIds[i]]+=Fe[i];
            }
            vector<double> bc_values(boundary.size(),0.0);
            SolveThreeLayerScheme(K,M,F,U[n-1],U[n],U[n+1],dt,boundary,bc_values);
        }
        double err=0,norm=0;
        for (int i=0;i<N;++i) {double uex=u_exact(nodes[i].x,nodes[i].y,T);err+=pow(U[Nt][i]-uex,2);norm+=pow(uex,2);}
        err=sqrt(err/N); norm=sqrt(norm/N);
        double rel_err=(norm>1e-10)?err/norm:err;
        errors.push_back(rel_err);
        double order=0.0;
        if (errors.size()>1) order=log(errors[errors.size()-2]/errors[errors.size()-1])/log(2.0);
        cout<<"Nt="<<setw(4)<<Nt<<"  error="<<scientific<<rel_err<<"  order="<<order<<endl;
    }
}

// --- Диагностический тест 4: уменьшение пространственного шага ---
void TestSpaceStep() {
    cout << "\n[DiagTest 4] Влияние пространственного шага (Nx,Ny) на ошибку при фиксированном Nt=256:\n";
    vector<int> nxy = {2, 4, 8};
    int Nt = 256;
    double T = 1.0;
    auto u_exact = [](double x,double y,double t){return sin(M_PI*x)*sin(M_PI*y)*pow(t,3);};
    auto u_t = [](double x,double y,double t){return 3*sin(M_PI*x)*sin(M_PI*y)*pow(t,2);};
    auto u_tt = [](double x,double y,double t){return 6*sin(M_PI*x)*sin(M_PI*y)*t;};
    auto f_rhs = [](double x,double y,double t){double pi2=M_PI*M_PI;return sin(M_PI*x)*sin(M_PI*y)*(6*t+2*pi2*pow(t,3));};
    for (int n : nxy) {
        int nx=n, ny=n; double hx=1.0, hy=1.0;
        vector<Node> nodes; vector<RectElement> elements;
        GenerateRectGrid(nx,ny,hx,hy,nodes,elements);
        int N=nodes.size(); vector<double> lambda(N,1.0);
        vector<vector<double>> K(N,vector<double>(N,0.0)), M(N,vector<double>(N,0.0));
        AssembleGlobalSystem(nodes,elements,lambda,K,M);
        for (int i=0;i<N;++i) for (int j=0;j<N;++j) if (i!=j) M[i][j]=0.0;
        vector<int> boundary=FindBoundaryNodes(nodes,nx,ny);
        double dt=T/Nt;
        vector<double> u_prev(N), u_curr(N), u_next(N);
        for (int i=0;i<N;++i) {
            u_prev[i]=u_exact(nodes[i].x,nodes[i].y,0.0);
            u_curr[i]=u_prev[i]+dt*u_t(nodes[i].x,nodes[i].y,0.0)+0.5*dt*dt*u_tt(nodes[i].x,nodes[i].y,0.0);
        }
        for (int b:boundary) {u_prev[b]=0.0;u_curr[b]=0.0;}
        vector<vector<double>> U(Nt+1,vector<double>(N)); U[0]=u_prev; U[1]=u_curr;
        for (int n=1;n<Nt;++n) {
            double t=(n+1)*dt;
            vector<double> F(N); AssembleGlobalRHS(nodes,elements,f_rhs,t,F);
            vector<double> bc_values(boundary.size(),0.0);
            SolveThreeLayerScheme(K,M,F,U[n-1],U[n],U[n+1],dt,boundary,bc_values);
        }
        double err=0,norm=0;
        for (int i=0;i<N;++i) {double uex=u_exact(nodes[i].x,nodes[i].y,T);err+=pow(U[Nt][i]-uex,2);norm+=pow(uex,2);}
        err=sqrt(err/N); norm=sqrt(norm/N);
        double rel_err=(norm>1e-10)?err/norm:err;
        cout<<"Nx=Ny="<<setw(3)<<n<<"  rel_error="<<scientific<<rel_err<<endl;
    }
}

// --- Диагностический тест 5: сравнение с/без масслимитации ---
void TestMassLumping() {
    cout << "\n[DiagTest 5] Масслимитация: сравнение с/без диагонализации матрицы масс:\n";
    int nx=10, ny=10, Nt=128; double hx=1.0, hy=1.0, T=1.0;
    auto u_exact = [](double x,double y,double t){return sin(M_PI*x)*sin(M_PI*y)*pow(t,3);};
    auto u_t = [](double x,double y,double t){return 3*sin(M_PI*x)*sin(M_PI*y)*pow(t,2);};
    auto u_tt = [](double x,double y,double t){return 6*sin(M_PI*x)*sin(M_PI*y)*t;};
    auto f_rhs = [](double x,double y,double t){double pi2=M_PI*M_PI;return sin(M_PI*x)*sin(M_PI*y)*(6*t+2*pi2*pow(t,3));};
    vector<Node> nodes; vector<RectElement> elements;
    GenerateRectGrid(nx,ny,hx,hy,nodes,elements);
    int N=nodes.size(); vector<double> lambda(N,1.0);
    vector<vector<double>> K(N,vector<double>(N,0.0)), M(N,vector<double>(N,0.0)), Mfull;
    AssembleGlobalSystem(nodes,elements,lambda,K,M);
    Mfull = M;
    vector<int> boundary=FindBoundaryNodes(nodes,nx,ny);
    double dt=T/Nt;
    // Без масслимитации
    vector<double> u_prev(N), u_curr(N), u_next(N);
    for (int i=0;i<N;++i) {
        u_prev[i]=u_exact(nodes[i].x,nodes[i].y,0.0);
        u_curr[i]=u_prev[i]+dt*u_t(nodes[i].x,nodes[i].y,0.0)+0.5*dt*dt*u_tt(nodes[i].x,nodes[i].y,0.0);
    }
    for (int b:boundary) {u_prev[b]=0.0;u_curr[b]=0.0;}
    vector<vector<double>> U(Nt+1,vector<double>(N)); U[0]=u_prev; U[1]=u_curr;
    for (int n=1;n<Nt;++n) {
        double t=(n+1)*dt;
        vector<double> F(N); AssembleGlobalRHS(nodes,elements,f_rhs,t,F);
        vector<double> bc_values(boundary.size(),0.0);
        SolveThreeLayerScheme(K,Mfull,F,U[n-1],U[n],U[n+1],dt,boundary,bc_values);
    }
    double err_full=0,norm=0;
    for (int i=0;i<N;++i) {double uex=u_exact(nodes[i].x,nodes[i].y,T);err_full+=pow(U[Nt][i]-uex,2);norm+=pow(uex,2);}
    err_full=sqrt(err_full/N); norm=sqrt(norm/N);
    double rel_err_full=(norm>1e-10)?err_full/norm:err_full;
    // С масслимитацией
    for (int i=0;i<N;++i) for (int j=0;j<N;++j) if (i!=j) M[i][j]=0.0;
    u_prev.clear(); u_curr.clear(); u_next.clear(); U.clear();
    u_prev.resize(N); u_curr.resize(N); U.resize(Nt+1,vector<double>(N));
    for (int i=0;i<N;++i) {
        u_prev[i]=u_exact(nodes[i].x,nodes[i].y,0.0);
        u_curr[i]=u_prev[i]+dt*u_t(nodes[i].x,nodes[i].y,0.0)+0.5*dt*dt*u_tt(nodes[i].x,nodes[i].y,0.0);
    }
    for (int b:boundary) {u_prev[b]=0.0;u_curr[b]=0.0;}
    U[0]=u_prev; U[1]=u_curr;
    for (int n=1;n<Nt;++n) {
        double t=(n+1)*dt;
        vector<double> F(N); AssembleGlobalRHS(nodes,elements,f_rhs,t,F);
        vector<double> bc_values(boundary.size(),0.0);
        SolveThreeLayerScheme(K,M,F,U[n-1],U[n],U[n+1],dt,boundary,bc_values);
    }
    double err_lump=0;
    for (int i=0;i<N;++i) {double uex=u_exact(nodes[i].x,nodes[i].y,T);err_lump+=pow(U[Nt][i]-uex,2);}
    err_lump=sqrt(err_lump/N);
    double rel_err_lump=(norm>1e-10)?err_lump/norm:err_lump;
    cout<<"Без масслимитации: rel_error="<<scientific<<rel_err_full<<endl;
    cout<<"С масслимитацией:  rel_error="<<scientific<<rel_err_lump<<endl;
}

// --- Диагностический тест 6: аналитика на одном элементе (nx=1,ny=1) ---
void TestAnalyticElement() {
    cout << "\n[DiagTest 6] Аналитика на одном элементе (nx=1,ny=1):\n";
    int nx=1, ny=1, Nt=128; double hx=1.0, hy=1.0, T=1.0;
    auto u_exact = [](double x,double y,double t){return sin(M_PI*x)*sin(M_PI*y)*pow(t,3);};
    auto u_t = [](double x,double y,double t){return 3*sin(M_PI*x)*sin(M_PI*y)*pow(t,2);};
    auto u_tt = [](double x,double y,double t){return 6*sin(M_PI*x)*sin(M_PI*y)*t;};
    auto f_rhs = [](double x,double y,double t){double pi2=M_PI*M_PI;return sin(M_PI*x)*sin(M_PI*y)*(6*t+2*pi2*pow(t,3));};
    vector<Node> nodes; vector<RectElement> elements;
    GenerateRectGrid(nx,ny,hx,hy,nodes,elements);
    int N=nodes.size(); vector<double> lambda(N,1.0);
    vector<vector<double>> K(N,vector<double>(N,0.0)), M(N,vector<double>(N,0.0));
    AssembleGlobalSystem(nodes,elements,lambda,K,M);
    for (int i=0;i<N;++i) for (int j=0;j<N;++j) if (i!=j) M[i][j]=0.0;
    vector<int> boundary=FindBoundaryNodes(nodes,nx,ny);
    double dt=T/Nt;
    vector<double> u_prev(N), u_curr(N), u_next(N);
    for (int i=0;i<N;++i) {
        u_prev[i]=u_exact(nodes[i].x,nodes[i].y,0.0);
        u_curr[i]=u_prev[i]+dt*u_t(nodes[i].x,nodes[i].y,0.0)+0.5*dt*dt*u_tt(nodes[i].x,nodes[i].y,0.0);
    }
    for (int b:boundary) {u_prev[b]=0.0;u_curr[b]=0.0;}
    vector<vector<double>> U(Nt+1,vector<double>(N)); U[0]=u_prev; U[1]=u_curr;
    for (int n=1;n<Nt;++n) {
        double t=(n+1)*dt;
        vector<double> F(N); AssembleGlobalRHS(nodes,elements,f_rhs,t,F);
        vector<double> bc_values(boundary.size(),0.0);
        SolveThreeLayerScheme(K,M,F,U[n-1],U[n],U[n+1],dt,boundary,bc_values);
    }
    double err=0,norm=0;
    for (int i=0;i<N;++i) {double uex=u_exact(nodes[i].x,nodes[i].y,T);err+=pow(U[Nt][i]-uex,2);norm+=pow(uex,2);}
    err=sqrt(err/N); norm=sqrt(norm/N);
    double rel_err=(norm>1e-10)?err/norm:err;
    cout<<"rel_error="<<scientific<<rel_err<<endl;
}

// --- Диагностический тест: Проверка K*1 для Q2 ---
void TestStiffnessMatrixConstantQ2() {
    cout << "\n[DiagTest K*1 Q2] Проверка K*1 для Q2 на сетке 2x2:" << endl;
    int nx = 2, ny = 2;
    double hx = 1.0, hy = 1.0;
    vector<Node> nodes;
    vector<RectElement> elements;
    GenerateRectGrid(nx, ny, hx, hy, nodes, elements);
    int N = nodes.size();
    vector<double> lambda(N, 1.0);
    vector<vector<double>> K(N, vector<double>(N, 0.0));
    vector<vector<double>> M(N, vector<double>(N, 0.0));
    AssembleGlobalSystem(nodes, elements, lambda, K, M);
    vector<double> ones(N, 1.0);
    vector<double> K1(N, 0.0);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            K1[i] += K[i][j] * ones[j];
    double max_abs = 0.0;
    for (int i = 0; i < N; ++i) max_abs = max(max_abs, fabs(K1[i]));
    cout << "Максимум |K*1| = " << scientific << max_abs << endl;
    cout << "Вектор K*1: ";
    for (int i = 0; i < N; ++i) cout << K1[i] << " ";
    cout << endl;
}


// --- Вызов всех диагностических тестов в main ---
int main() {
    try {
        cout << "Запуск теста 3.1 (полиномиальные функции)..." << endl;
        TestOrderTimePoly();
        cout << "\nЗапуск теста 3.2 (неполиномиальные функции)..." << endl;
        TestOrderTimeNonPoly();
        cout << "\nЗапуск основного цикла..." << endl;
        RunThreeLayerFEM_Full();
        // Диагностические тесты:
        TestTimeOnly();
        TestConstantSolution();
        TestQ1Basis();
        TestSpaceStep();
        TestMassLumping();
        TestAnalyticElement();
        TestStiffnessMatrixConstantQ2();
        cout << "Программа завершена успешно." << endl;
    } catch (const exception& e) {
        cout << "Ошибка: " << e.what() << endl;
        return 1;
    }
    return 0;
}