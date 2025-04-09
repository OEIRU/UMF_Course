#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <functional>

const double PI = 3.14159265358979;

struct area_element_t
{
    int func_number;
    int x1, y1, x2, y2;
    area_element_t() {}
    area_element_t(int _func_number, int _x1, int _y1, int _x2, int _y2)
    {
        func_number = _func_number;
        x1 = _x1;
        y1 = _y1;
        x2 = _x2;
        y2 = _y2;
    }
};

struct area_t
{
    std::vector<double> x_lines, y_lines;
    std::vector<int> Ix, Iy; // Индексы значений area_t.xy_lines в grid_t.xy
    std::vector<area_element_t> elems;
};

struct parameters_t
{
    double lambda, gamma;
};

struct grid_t
{
    std::vector<double> x, y;
};

struct SLAE
{
    std::vector<double> di1, di2, diMain, di4, di5, b;
    std::vector<int> ig;
    std::vector<std::function<double(double, double)>> F{
        [](double x, double y)
        { return 2 * (x + y); },
        [](double x, double y)
        { return 0; },
        [](double x, double y)
        { return 0; }};
};

struct functionsBC
{
    std::vector<std::function<double(double, double)>> firstBC{
        [](double x, double y)
        { return 1 + y; }, // 1
        [](double x, double y)
        { return x + 4; }, // 2
        [](double x, double y)
        { return 2 + y; }, // 3
        [](double x, double y)
        { return x + 3; }, // 4
        [](double x, double y)
        { return 3 + y; }, // 5
        [](double x, double y)
        { return x + 4; }, // 6
        [](double x, double y)
        { return 4 + y; }, // 7
        [](double x, double y)
        { return x + 2; }, // 8
        [](double x, double y)
        { return 2 + y; }, // 9
        [](double x, double y)
        { return x + 1; }, // 10
        [](double x, double y)
        { return 2 + y; }, // 11
        [](double x, double y)
        { return x + 1; } // 12
    };

    std::vector<std::function<double(double, double)>> secondBC{
        [](double x, double y)
        { return -cos(x + y); },
        [](double x, double y)
        { return 0; }};
};

void read_area(std::string filename, area_t &eA);
void read_params(std::string filename, parameters_t &par);
void read_boundary_conds(std::string filename, std::vector<std::vector<int>> &bC);
void read_grid(std::string filename, area_t &eA, grid_t &g);
int returnNumberOfEstimatedSubArea(area_t &area, int p, int s, int &l);
bool isFictiousNode(area_t &area, int p, int s);
void clearFictiousNodes(SLAE &slae, area_t &area, int nX, int nY);
void createPortrait(SLAE &slae, int nX, int nY);
void addElemInGlobalMatrix(SLAE &slae, int i, int j, double elem);
void calcGlobalMatrixAndVector(SLAE &slae, grid_t &g, area_t &area, parameters_t &par);
void calcBoundaryConditions(area_t &area, grid_t &g, SLAE &slae, parameters_t &par, std::vector<std::vector<int>> &bC, functionsBC &func);
void stringMatrixInNull(SLAE &slae, int i);
void methodOfIterations(SLAE &slae, std::vector<double> &x, std::vector<double> &xk, double eps, int maxIter, double w);
double normVector(std::vector<double> &a);
double calcDiscrepancy(SLAE &slae, std::vector<double> &x);
void output(area_t &area, std::vector<double> &u, grid_t &g); // !
double uReal(double x, double y);
