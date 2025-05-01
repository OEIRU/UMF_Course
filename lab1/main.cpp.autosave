#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional>
#include <string>
#include <tuple>
#include <fstream>
#include <cstdlib>
using namespace std;

// Типы узлов сетки
enum class NodeType { INNER, EDGE, FICTITIOUS };

// Типы краевых условий
enum class BoundaryCondition { FIRST, THIRD };

// Структура для представления узла сетки
struct Node { 
    double x, y;  // Координаты узла
    NodeType type;  // Тип узла
};

// Базовый класс геометрии
class Figure {
public:
    // Проверяет, находится ли точка внутри области
    virtual bool isInside(double x, double y, double eps) const = 0;
    
    // Возвращает тип краевого условия в заданной точке
    virtual BoundaryCondition getBoundaryCondition(double x, double y, double eps) const = 0;
    
    // Возвращает значение краевого условия
    virtual double getBoundaryValue(double x, double y) const = 0;
    
    // Возвращает коэффициенты для уравнения
    virtual double getLambda(double x, double y) const = 0;
    virtual double getBeta(double x, double y) const = 0;
    virtual double getGamma(double x, double y) const = 0;
};

// Класс для П-образной области
class PiShapedFigure : public Figure {
public:
    int gType = 1;  // Тип краевого условия
    double g0 = 0.0;  // Константное значение для краевого условия
    double lambdaVal = 1.0, betaVal = 1.0, gammaVal = 0.0;  // Коэффициенты

    // Проверяет, находится ли точка (x,y) внутри области
    bool isInside(double x, double y, double eps) const override {
        // Проверяем три части П-образной области:
        // левую вертикальную часть, правую вертикальную часть и верхнюю горизонтальную часть
        bool left  = (x >= -eps && x <= 0.5 + eps && y >= -eps && y <= 2 + eps);
        bool right = (x >= 1.5 - eps && x <= 2 + eps && y >= -eps && y <= 2 + eps);
        bool top   = (x >= -eps && x <= 2 + eps && y >= 2 - eps && y <= 2.5 + eps);
        return left || right || top;
    }

    // Возвращает тип краевого условия в заданной точке
    BoundaryCondition getBoundaryCondition(double x, double y, double eps) const override {
        // Если точка находится на верхней границе (y=2.5), то краевое условие третьего рода,
        // иначе первое
        return (fabs(y - 2.5) < eps) ? BoundaryCondition::THIRD : BoundaryCondition::FIRST;
    }

    // Возвращает значение краевого условия
    double getBoundaryValue(double x, double) const override {
        // Если тип 2, возвращаем константу g0, иначе sin(πx)
        return (gType == 2) ? g0 : sin(M_PI * x);
    }

    // Возвращает коэффициенты для уравнения
    double getLambda(double, double) const override { return lambdaVal; }
    double getBeta(double, double)  const override { return betaVal; }
    double getGamma(double, double) const override { return gammaVal; }

    // Установка параметров краевых условий
    void setDirichletType(int t) { gType = t; }
    void setDirichletValue(double v) { g0 = v; }
    void setRobinParams(double l, double b, double g) { 
        lambdaVal = l; betaVal = b; gammaVal = g; 
    }
};

// Обработчик краевых условий
class BoundaryConditionHandler {
public:
    // Применяет краевое условие первого рода (Дирихле)
    void applyFirstCondition(int, double& di, double& b, double ug) {
        di = 1.0;  // Диагональный элемент матрицы
        b = ug;    // Правая часть
    }

    // Применяет краевое условие третьего рода (Робина)
    void applyThirdCondition(int, double& di, double& b, double h,
                             double lambda, double beta, double gamma) {
        // Обновляем диагональный элемент и правую часть для условия Робина
        di += (lambda / h + beta);
        b += gamma;
    }
};

// Класс для создания сетки
class Grid {
public:
    vector<Node> nodes;  // Список узлов
    vector<double> hx, hy;  // Шаги по осям x и y
    Figure* fig;  // Указатель на объект фигуры

    // Конструктор для равномерной сетки
    Grid(double x0, double x1, double y0, double y1, int nx, int ny, Figure* f)
        : fig(f) {
        double dx = (x1 - x0) / (nx - 1);  // Шаг по x
        double dy = (y1 - y0) / (ny - 1);  // Шаг по y
        hx.assign(nx - 1, dx);  // Инициализируем шаги по x
        hy.assign(ny - 1, dy);  // Инициализируем шаги по y
        build(x0, y0);  // Строим сетку
    }

    // Конструктор для неравномерной сетки
    Grid(double x0, double y0,
         const vector<double>& hx_, const vector<double>& hy_, Figure* f)
        : hx(hx_), hy(hy_), fig(f) {
        build(x0, y0);  // Строим сетку
    }

private:
    // Метод построения сетки
    void build(double x0, double y0) {
        int nx = hx.size() + 1;  // Количество узлов по x
        int ny = hy.size() + 1;  // Количество узлов по y
        
        // Векторы для координат узлов
        vector<double> xs(nx), ys(ny);
        
        // Вычисляем координаты по x
        xs[0] = x0;
        for (int i = 1; i < nx; ++i) xs[i] = xs[i - 1] + hx[i - 1];
        
        // Вычисляем координаты по y
        ys[0] = y0;
        for (int j = 1; j < ny; ++j) ys[j] = ys[j - 1] + hy[j - 1];
        
        // Создаем узлы
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                double x = xs[i], y = ys[j];
                Node node{x, y, NodeType::FICTITIOUS};  // По умолчанию фиктивный
                
                // Если точка внутри фигуры, то внутренний узел
                if (fig->isInside(x, y, 1e-10)) {
                    node.type = NodeType::INNER;
                    
                    // Если узел на границе, меняем тип на EDGE
                    if (isBoundary(x, y)) node.type = NodeType::EDGE;
                }
                nodes.push_back(node);
            }
        }
    }

    // Проверяет, является ли узел граничным
    bool isBoundary(double x, double y) const {
        // Проверяем, находится ли точка за пределами области при небольших отклонениях
        return !fig->isInside(x - 1e-10, y, 1e-10) ||
               !fig->isInside(x + 1e-10, y, 1e-10) ||
               !fig->isInside(x, y - 1e-10, 1e-10) ||
               !fig->isInside(x, y + 1e-10, 1e-10);
    }

    // Возвращает индекс узла по координатам (i,j)
    int getIndex(int i, int j, int nx) const { 
        return j * nx + i; 
    }
    
    friend class EllipticSolver;
};

// Абстрактный класс итерационного решателя
class IterativeSolver {
protected:
    // Основные элементы матрицы системы
    vector<double> di, au1, au2, al1, al2, b;
public:
    // Метод решения системы
    virtual vector<double> solve(int maxIter, double tol) = 0;
    
    // Устанавливает систему уравнений
    void setSystem(const vector<double>& d, const vector<double>& u1, const vector<double>& u2,
                   const vector<double>& l1, const vector<double>& l2, const vector<double>& rhs) {
        di = d; au1 = u1; au2 = u2;
        al1 = l1; al2 = l2; b = rhs;
    }
};

// Метод Якоби для решения системы
class JacobiMethod : public IterativeSolver {
public:
    // Реализация метода Якоби
    vector<double> solve(int maxIter, double tol) override {
        vector<double> x(di.size(), 0.0), xk(di.size());
        
        for (int it = 0; it < maxIter; ++it) {
            // Выполняем итерацию Якоби
            for (size_t i = 0; i < di.size(); ++i) {
                double sum = 0;
                if (i > 0) sum += al1[i-1] * x[i-1];  // Левый сосед
                if (i < di.size()-1) sum += au1[i] * x[i+1];  // Правый сосед
                xk[i] = (b[i] - sum) / di[i];  // Новое значение
            }
            
            // Вычисляем ошибку
            double err = 0;
            for (size_t i = 0; i < di.size(); ++i) 
                err = max(err, fabs(xk[i] - x[i]));
                
            x = xk;  // Обновляем решение
            
            // Проверяем сходимость
            if (err < tol) break;
        }
        return x;
    }
};

// Класс для решения эллиптической задачи
class EllipticSolver {
    Grid grid;  // Объект сетки
    IterativeSolver* solver;  // Указатель на решатель
    BoundaryConditionHandler bc;  // Обработчик краевых условий
    
    // Массив записей матрицы для визуализации
    vector<tuple<int,int,double>> matEntries;

    // Метод для сборки матрицы системы
    void assemble(vector<double>& di, vector<double>& au1, vector<double>& al1, vector<double>& b) {
        int nx = grid.hx.size() + 1;  // Количество узлов по x
        
        // Инициализируем вектора нулями
        di.assign(grid.nodes.size(), 0.0);
        au1.assign(grid.nodes.size()-1, 0.0);
        al1.assign(grid.nodes.size()-1, 0.0);
        b.assign(grid.nodes.size(), 0.0);
        matEntries.clear();
        
        // Проходим по всем узлам
        for (int j = 0; j < grid.hy.size() + 1; ++j) {
            for (int i = 0; i < nx; ++i) {
                int idx = grid.getIndex(i, j, nx);
                
                // Если узел внутренний, формируем уравнение
                if (grid.nodes[idx].type == NodeType::INNER) {
                    double hxm = grid.hx[max(0,i-1)], hxp = grid.hx[i];  // Шаги по x
                    double hym = grid.hy[max(0,j-1)], hyp = grid.hy[j];  // Шаги по y
                    
                    // Центральный коэффициент
                    di[idx] = -2/(hxm*hxp) - 2/(hym*hyp);
                    
                    // Левый сосед
                    if (i > 0) {
                        al1[idx-1] = 2/(hxm*(hxm+hxp));
                        matEntries.emplace_back(idx, idx-1, al1[idx-1]);
                    }
                    
                    // Правый сосед
                    if (i < nx-1) {
                        au1[idx] = 2/(hxp*(hxm+hxp));
                        matEntries.emplace_back(idx, idx+1, au1[idx]);
                    }
                    
                    matEntries.emplace_back(idx, idx, di[idx]);
                } 
                // Если узел на границе, применяем краевые условия
                else if (grid.nodes[idx].type == NodeType::EDGE) {
                    double x = grid.nodes[idx].x, y = grid.nodes[idx].y;
                    auto bcType = grid.fig->getBoundaryCondition(x, y, 1e-10);
                    
                    if (bcType == BoundaryCondition::FIRST) {
                        // Краевое условие первого рода
                        bc.applyFirstCondition(idx, di[idx], b[idx], grid.fig->getBoundaryValue(x,y));
                    } else {
                        // Краевое условие третьего рода
                        bc.applyThirdCondition(idx, di[idx], b[idx], grid.hx[i],
                                              grid.fig->getLambda(x,y), grid.fig->getBeta(x,y), grid.fig->getGamma(x,y));
                    }
                    matEntries.emplace_back(idx, idx, di[idx]);
                }
            }
        }
        
        // Передаем систему решателю
        solver->setSystem(di, au1, vector<double>(), al1, vector<double>(), b);
    }

public:
    // Конструкторы
    EllipticSolver(double x0, double x1, double y0, double y1, int nx, int ny, Figure* f, IterativeSolver* s)
        : grid(x0, x1, y0, y1, nx, ny, f), solver(s) {}
        
    EllipticSolver(double x0, double y0, const vector<double>& hx, const vector<double>& hy, Figure* f, IterativeSolver* s)
        : grid(x0, y0, hx, hy, f), solver(s) {}

    // Метод решения и вывода результатов
    vector<double> solveAndPrint(int maxIter, double tol, function<double(double,double)> exact) {
        vector<double> di, au1, al1, b;
        assemble(di, au1, al1, b);
        auto sol = solver->solve(maxIter, tol);
        
        // Подготавливаем вывод
        cout << left << setw(12) << "Точка" << setw(15) << "Точное" << setw(15) << "Численное" 
             << setw(20) << "Погрешность" << "" << endl;
             
        ofstream sdf("solution_data.txt");  // Файл для данных решения
        ofstream csv("results.csv");  // CSV-файл для Excel
        
        // Заголовок CSV
        csv << "x,y,Exact,Numeric,Error" << endl;
        
        // Выводим результаты
        for (size_t idx = 0; idx < grid.nodes.size(); ++idx) {
            if (grid.nodes[idx].type != NodeType::FICTITIOUS) {
                double xv = grid.nodes[idx].x, yv = grid.nodes[idx].y;
                double ev = exact(xv, yv);  // Точное значение
                double err = fabs(sol[idx] - ev);  // Ошибка
                
                // Вывод в консоль
                cout << left << setw(12) << "(" + to_string(xv) + "," + to_string(yv) + ")"
                     << setw(15) << setprecision(6) << ev
                     << setw(15) << setprecision(6) << sol[idx]
                     << setw(20) << scientific << setprecision(2) << err << " " << endl;
                
                // Запись в файлы
                sdf << xv << " " << yv << " " << sol[idx] << "\n";
                csv << xv << "," << yv << "," << ev << "," << sol[idx] << "," << err << endl;
            }
        }
        
        sdf.close();
        csv.close();
        cout << "Результаты также сохранены в файле results.csv для импорта в Excel." << endl;
        return sol;
    }

    // Метод для визуализации матрицы
    void dumpMatrix(const string& prefix) {
        string dataFile = prefix + ".dat";
        ofstream df(dataFile);
        
        // Сохраняем элементы матрицы
        for (auto& t : matEntries) {
            int i,j; double v;
            tie(i,j,v) = t;
            df << i << " " << j << " " << v << "\n";
        }
        df.close();
        
        string script = prefix + ".gnuplot";
        ofstream sf(script);
        
        // Скрипт для gnuplot
        sf << "set title 'Sparsity pattern'\n"
           << "set xlabel 'Column'\n"
           << "set ylabel 'Row'\n";
        sf << "plot '" << dataFile << "' using 2:1 with points pt 7 ps 0.5 title ''\n";
        sf << "pause -1\n";

        sf.close();
        system(("gnuplot " + script).c_str());
    }
};

// Метод визуализации решения
int plot_output() {
    // Проверяем существование файла
    ifstream f("solution_data.txt");
    if (!f.good()) {
        cerr << "Ошибка: Файл solution_data.txt не найден. Сначала выполните расчет (пункт 5)." << endl;
        return 1;
    }

    ofstream script("plot_script.gp");
    script << "set pm3d map\n"
           << "splot 'solution_data.txt' using 1:2:3 with points pt 7 ps 0.5\n"
           << "pause -1\n";
    script.close();

    // Универсальный вызов gnuplot (работает и на Windows, и на Linux)
    system("gnuplot plot_script.gp");

    remove("plot_script.gp");
    return 0;
}

// Метод визуализации области
int plot_domain(PiShapedFigure& fig, int nx, int ny) {
    ofstream df("domain_data.txt");
    double x0 = 0, x1 = 2, y0 = 0, y1 = 2.5;
    double dx = (x1 - x0) / (nx - 1), dy = (y1 - y0) / (ny - 1);
    
    // Генерируем точки для визуализации
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double x = x0 + i * dx;
            double y = y0 + j * dy;
            df << x << " " << y << " " << (fig.isInside(x, y, 1e-10) ? 1 : 0) << "\n";
        }
        df << "\n";
    }
    df.close();
    
    ofstream sf("plot_domain.gp");
    sf << "set pm3d map\n";
    sf << "splot 'domain_data.txt' using 1:2:3 with pm3d\n";
    sf << "pause -1\n";
    sf.close();
    
    system("gnuplot plot_domain.gp");
    //system("/usr/bin/gnuplot plot_domain.gp");
    remove("plot_domain.gp");
    remove("domain_data.txt");
    return 0;
}

// Основная функция
int main() {
    PiShapedFigure figure;  // Объект П-образной области
    JacobiMethod jacobi;    // Решатель методом Якоби
    
    // Параметры по умолчанию
    int nx = 10, ny = 10, maxIter = 1000, gridType = 1;
    double tol = 1e-6, rx = 1.0, ry = 1.0;
    
    // Главный цикл обработки команд
    while (true) {
        // Меню команд
        cout << "1. Размер сетки\n"
             << "2. Параметры решателя\n"
             << "3. Тип сетки\n"
             << "4. Параметры BC\n"
             << "5. Решить и вывести\n"
             << "6. Визуализировать матрицу\n"
             << "7. Визуализировать решение\n"
             << "8. Визуализировать область\n"
             << "9. Выход\n";
             
        int c; 
        cin >> c;
        
        // Обработка команд
        if (c == 1) {  // Размер сетки
            cout << "nx ny: ";
            cin >> nx >> ny;
        } else if (c == 2) {  // Параметры решателя
            cout << "maxIter tol: ";
            cin >> maxIter >> tol;
        } else if (c == 3) {  // Тип сетки
            cout << "1-uniform,2-ratio: ";
            cin >> gridType;
            if (gridType == 2) {
                cout << "rx ry: ";
                cin >> rx >> ry;
            }
        } else if (c == 4) {  // Параметры краевых условий
            int dt;
            cout << "Dirichlet 1-sin(pi x),2-const: ";
            cin >> dt;
            figure.setDirichletType(dt);
            
            if (dt == 2) {
                double v;
                cout << "g= ";
                cin >> v;
                figure.setDirichletValue(v);
            }
            
            double l, b, g;
            cout << "Robin lambda beta gamma: ";
            cin >> l >> b >> g;
            figure.setRobinParams(l, b, g);
        } else if (c == 5 || c == 6) {  // Решение задачи
            EllipticSolver* sol;
            
            // Создаем решатель в зависимости от типа сетки
            if (gridType == 1) {
                sol = new EllipticSolver(0, 2, 0, 2.5, nx, ny, &figure, &jacobi);
            } else {
                int sx = nx - 1, sy = ny - 1;
                vector<double> hx(sx), hy(sy);
                double Lx = 2, Ly = 2.5;
                
                // Вычисляем шаги для неравномерной сетки
                double dx0 = (Lx * (rx - 1)) / (pow(rx, sx) - 1);
                double dy0 = (Ly * (ry - 1)) / (pow(ry, sy) - 1);
                
                for (int i = 0; i < sx; ++i) 
                    hx[i] = (fabs(rx - 1) < 1e-12) ? Lx/sx : dx0 * pow(rx, i);
                    
                for (int j = 0; j < sy; ++j) 
                    hy[j] = (fabs(ry - 1) < 1e-12) ? Ly/sy : dy0 * pow(ry, j);
                    
                sol = new EllipticSolver(0, 0, hx, hy, &figure, &jacobi);
            }
            
            // Решаем и выводим результаты
            sol->solveAndPrint(maxIter, tol, [](double x, double y){ 
                return sin(M_PI * x) * sinh(M_PI * y); 
            });
            
            // Если нужно визуализировать матрицу
            if (c == 6) {
                sol->dumpMatrix("matrix");
            }
            
            delete sol;
        } else if (c == 7) {  // Визуализация решения
            plot_output();
        } else if (c == 8) {  // Визуализация области
            plot_domain(figure, nx, ny);
        } else if (c == 9) {  // Выход
            break;
        }
    }
    
    return 0;
}
