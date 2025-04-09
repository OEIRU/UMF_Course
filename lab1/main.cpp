#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <functional>
#include "function.h"

int main()
{
    const double eps = 1e-15, w = 1.2;
    const int maxIter = 1000;

    area_t area;
    parameters_t params;
    std::vector<std::vector<int>> boundary_conds;
    grid_t grid{};

    SLAE slae{};
    functionsBC func{};
    std::vector<double> u{};

    read_area("in/area.txt", area);
    read_params("in/params.txt", params);
    read_boundary_conds("in/boundary_conds.txt", boundary_conds);
    read_grid("in/grid.txt", area, grid);

    createPortrait(slae, grid.x.size(), grid.y.size());
    calcGlobalMatrixAndVector(slae, grid, area, params);
    calcBoundaryConditions(area, grid, slae, params, boundary_conds, func);
    clearFictiousNodes(slae, area, grid.x.size(), grid.y.size());
    methodOfIterations(slae, u, u, eps, maxIter, w);
    output(area, u, grid);

    return 0;
}
