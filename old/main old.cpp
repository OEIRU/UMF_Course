#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdexcept>
using namespace std;

class Solver {
public:
    static bool Gauss_Seidel(const vector<double>& Diag,
                             const vector<double>& Lower,
                             const vector<double>& Upper,
                             const vector<double>& Right,
                             vector<double>& U,
                             int max_iter = 10000,
                             double tol = 1e-6) {
        const int n = Diag.size();
        double error;
        for (int iter = 0; iter < max_iter; ++iter) {
            error = 0.0;
            for (int i = 0; i < n; ++i) {
                double sum = Right[i];
                if (i > 0) sum -= Lower[i-1] * U[i-1];
                if (i < n-1) sum -= Upper[i] * U[i+1];
                
                double new_val = sum / Diag[i];
                error += abs(new_val - U[i]);
                U[i] = new_val;
            }
            // Логирование ошибки на каждой итерации
            cout << "Iteration " << iter << ": Error = " << error << endl;
            if (error < tol) {
                cout << "Converged in " << iter << " iterations" << endl;
                return true;
            }
        }
        cerr << "No convergence" << endl;
        return false;
    }
};

class MKR {
    vector<double> Diag, Lower, Upper, Right, U;
    vector<bool> is_boundary;
    int nx, ny;
    double hx, hy;
    double lambda, alpha, beta;
    int index(int i, int j) const { return i*ny + j; }

    void mark_P_shape() {
        const int cut_x = nx/4;
        const int cut_y = ny/4;
        
        for (int i = cut_x; i < 3*cut_x; ++i) {
            for (int j = cut_y; j < ny - cut_y; ++j) {
                if (j > cut_y && j < ny - cut_y - 1) continue;
                is_boundary[index(i, j)] = true;
            }
        }
    }

public:
    MKR(const string& param_file) {
        ifstream in(param_file);
        if (!in.is_open()) {
            throw runtime_error("Cannot open parameter file");
        }
        in >> nx >> ny >> hx >> hy >> lambda >> alpha >> beta;
        const int n = nx * ny;
        Diag.resize(n, 0.0);
        Lower.resize(n-1, 0.0);
        Upper.resize(n-1, 0.0);
        Right.resize(n, 0.0);
        U.resize(n, 0.0);
        is_boundary.resize(n, false);

        // Разметка П-образной области
        mark_P_shape();

        // Построение матрицы системы
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                const int idx = index(i, j);
                
                if (is_boundary[idx]) {
                    // Третье краевое условие
                    Diag[idx] = alpha;
                    Right[idx] = beta;
                } else {
                    // Внутренние узлы
                    Diag[idx] = -2.0/(hx*hx) - 2.0/(hy*hy);
                    if (i > 0) Lower[idx-ny] = 1.0/(hx*hx);
                    if (i < nx-1) Upper[idx] = 1.0/(hx*hx);
                    if (j > 0) Lower[idx-1] += 1.0/(hy*hy);
                    if (j < ny-1) Upper[idx] += 1.0/(hy*hy);
                }
            }
        }

        // Логирование матрицы после её построения
        cout << "Matrix constructed:" << endl;
        //cout << "Diag: ";
        //for (double val : Diag) cout << val << " ";
        //cout << endl;
        //cout << "Lower: ";
        //for (double val : Lower) cout << val << " ";
        //cout << endl;
        //cout << "Upper: ";
        //for (double val : Upper) cout << val << " ";
        //cout << endl;
        //cout << "Right: ";
        //for (double val : Right) cout << val << " ";
        //cout << endl;
    }

    void solve() {
        if (!Solver::Gauss_Seidel(Diag, Lower, Upper, Right, U)) {
            throw runtime_error("Solver failed");
        }
    }

    void save_results(const string& filename) {
        ofstream out(filename);
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                out << U[index(i, j)] << " ";
            }
            out << "\n";
        }
    }

    size_t get_memory_usage() const {
        return (Diag.capacity() + Lower.capacity() + Upper.capacity() +
                Right.capacity() + U.capacity()) * sizeof(double) +
               is_boundary.capacity()/8;
    }
};

int main() {
    try {
        MKR mkr("params.txt");
        mkr.solve();
        mkr.save_results("result.txt");
        
        cout << "Memory used: " 
             << mkr.get_memory_usage()/1e6 << " MB" << endl;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}