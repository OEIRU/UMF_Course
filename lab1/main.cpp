#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// --- Глобальные параметры ---
const int nx = 21, ny = 21;
const double xmin = 0.0, xmax = 1.0;
const double ymin = 0.0, ymax = 1.0;
const double hx = (xmax - xmin) / (nx - 1);
const double hy = (ymax - ymin) / (ny - 1);
const double lambda_ = 1.0;
const double gamma_ = 1.0;
const double eps = 1e-8;
const int max_iter = 10000;
const double alpha = 1.0; // параметр для условий 3-го рода

struct TestCase {
  string name;
  double (*u_exact)(double, double);
  double (*f)(double, double);
  int (*getTypeOfB)(int);
  double (*getValueOfB)(int, double, double);
};

bool is_in_pi_domain(int i, int j) {
  int i1 = nx / 3, i2 = 2 * nx / 3;
  int j1 = ny / 3, j2 = 2 * ny / 3;
  return !(i > i1 && i < i2 && j > j1 && j < j2);
}

int idx(int i, int j) { return j * nx + i; }

void assemble_slae(const TestCase &test, vector<double> &di, vector<double> &dl, vector<double> &du,
                   vector<double> &df, vector<bool> &mask) {
  int N = nx * ny;
  di.assign(N, 0.0);
  dl.assign(N, 0.0);
  du.assign(N, 0.0);
  df.assign(N, 0.0);
  mask.assign(N, false);

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      int k = idx(i, j);
      if (!is_in_pi_domain(i, j)) continue;
      mask[k] = true;

      double x = xmin + i * hx;
      double y = ymin + j * hy;

      bool boundary = false;
      for (auto [di_, dj_] : { pair{-1, 0}, {1, 0}, {0, -1}, {0, 1} }) {
        int ni = i + di_, nj = j + dj_;
        if (ni < 0 || ni >= nx || nj < 0 || nj >= ny || !is_in_pi_domain(ni, nj)) {
          boundary = true;
          break;
        }
      }

      if (boundary) {
        int side = 1; // default boundary id
        int type = test.getTypeOfB(side);
        if (type == 1) {
          di[k] = 1.0;
          df[k] = test.getValueOfB(side, x, y);
        } else if (type == 3) {
          di[k] = 1.0 + alpha;
          df[k] = alpha * test.getValueOfB(side, x, y);
        }
        continue;
      }

      double hx2 = hx * hx, hy2 = hy * hy;
      double a = lambda_ / hx2;
      double b = lambda_ / hy2;
      double c = 2 * (a + b) + gamma_;

      di[k] = c;
      df[k] = test.f(x, y);

      if (i > 0 && is_in_pi_domain(i - 1, j)) dl[k - 1] = -a;
      if (i < nx - 1 && is_in_pi_domain(i + 1, j)) du[k] = -a;
      if (j > 0 && is_in_pi_domain(i, j - 1)) df[k] += -b * 0.0; // предварительно 0
      if (j < ny - 1 && is_in_pi_domain(i, j + 1)) df[k] += -b * 0.0; // предварительно 0
    }
  }
}

double solve_gauss_seidel(const vector<double> &di, const vector<double> &dl,
                          const vector<double> &du, const vector<double> &df,
                          vector<double> &x, const vector<bool> &mask) {
  int N = di.size();
  x.assign(N, 0.0);
  double r = 1.0, res0 = 0.0;
  for (int i = 0; i < N; ++i) {
    if (mask[i]) res0 += df[i] * df[i];
  }
  for (int iter = 0; iter < max_iter && r > eps; ++iter) {
    double res = 0.0;
    for (int i = 0; i < N; ++i) {
      if (!mask[i]) continue;
      double sum = di[i] * x[i];
      if (i > 0) sum += dl[i - 1] * x[i - 1];
      if (i < N - 1) sum += du[i] * x[i + 1];
      double delta = df[i] - sum;
      x[i] += delta / di[i];
      res += delta * delta;
    }
    r = sqrt(res / res0);
    if (iter % 100 == 0) cout << "iter=" << iter << ", resid=" << r << endl;
  }
  return r;
}

void run_test(const TestCase &test) {
  int N = nx * ny;
  vector<double> di, dl, du, df, x;
  vector<bool> mask;
  assemble_slae(test, di, dl, du, df, mask);
  double final_r = solve_gauss_seidel(di, dl, du, df, x, mask);

  ofstream out("result_" + test.name + ".txt");
  out << "Test: " << test.name << ", resid: " << final_r << "\n";
  out << setw(8) << "x" << setw(8) << "y"
      << setw(20) << "u_exact"
      << setw(20) << "u_numeric"
      << setw(20) << "abs_error" << endl;

  for (int j = ny - 1; j >= 0; --j) {
    for (int i = 0; i < nx; ++i) {
      int k = idx(i, j);
      if (!mask[k]) continue;
      double xx = xmin + i * hx;
      double yy = ymin + j * hy;
      double u_ex = test.u_exact(xx, yy);
      out << setw(8) << xx << setw(8) << yy
          << setw(20) << setprecision(12) << u_ex
          << setw(20) << setprecision(12) << x[k]
          << setw(20) << setprecision(12) << fabs(u_ex - x[k]) << endl;
    }
  }
  out.close();
  cout << "Тест завершён: " << test.name << endl;
}

int main() {
  auto getType = [](int b) -> int { return b == 1 ? 1 : 3; };
  auto getVal = [](int b, double x, double y) -> double {
    return x + y; // пока единый интерфейс, переопределяется для каждого теста
  };

  vector<TestCase> tests = {
    {"u = x + y",
     [](double x, double y) { return x + y; },
     [](double x, double y) { return x + y; }, getType, getVal},
    {"u = x^2 + y^2",
     [](double x, double y) { return x * x + y * y; },
     [](double x, double y) { return x * x + y * y - 4.0; }, getType, getVal},
    {"u = x^3 + y^3",
     [](double x, double y) { return pow(x, 3) + pow(y, 3); },
     [](double x, double y) { return pow(x, 3) - pow(y, 3) - 6 * x - 6 * y; }, getType, getVal},
    {"u = x^4 + y^4",
     [](double x, double y) { return pow(x, 4) + pow(y, 4); },
     [](double x, double y) { return pow(x, 4) + pow(y, 4) - 12 * x * x - 12 * y * y; }, getType, getVal},
    {"u = cos(2x + 2y)",
     [](double x, double y) { return cos(2 * x + 2 * y); },
     [](double x, double y) { return 9 * cos(2 * x + 2 * y); }, getType, getVal}
  };

  for (const auto &test : tests) {
    run_test(test);
  }
  return 0;
}
