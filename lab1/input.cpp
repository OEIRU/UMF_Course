#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <functional>
#include "function.h"

void read_area(std::string filename, area_t &area)
{
    std::ifstream in(filename);
    int count = 0;

    in >> count;
    area.x_lines.resize(count);
    for (int i = 0; i < count; i++)
        in >> area.x_lines[i];

    in >> count;
    area.y_lines.resize(count);
    for (int i = 0; i < count; i++)
        in >> area.y_lines[i];

    in >> count;
    area.elems.resize(count);
    for (int l = 0; l < count; l++)
    {
        in >> area.elems[l].func_number >> area.elems[l].x1 >> area.elems[l].x2 >> area.elems[l].y1 >> area.elems[l].y2;
        area.elems[l].x1--;
        area.elems[l].x2--;
        area.elems[l].y1--;
        area.elems[l].y2--;
    }

    in.close();
}

void read_params(std::string filename, parameters_t &params)
{
    std::ifstream in(filename);
    in >> params.lambda >> params.gamma;
    in.close();
}

void read_boundary_conds(std::string filename, std::vector<std::vector<int>> &boundary_conds)
{
    std::ifstream in(filename);
    int count = 0;

    in >> count;
    boundary_conds.resize(count);
    for (int i = 0; i < count; i++)
    {
        boundary_conds[i].resize(6);
        in >> boundary_conds[i][0] >> boundary_conds[i][1] >> boundary_conds[i][2] >> boundary_conds[i][3] >> boundary_conds[i][4] >> boundary_conds[i][5];
        boundary_conds[i][2]--;
        boundary_conds[i][3]--;
        boundary_conds[i][4]--;
        boundary_conds[i][5]--;
    }

    in.close();
}

void read_grid(std::string filename, area_t &area, grid_t &grid)
{
    const int nX = area.x_lines.size() - 1, nY = area.y_lines.size() - 1;
    area.Ix.resize(nX + 1);
    area.Iy.resize(nY + 1);
    std::ifstream fin(filename);

    int nXk = 0;
    grid.x.resize(1, area.x_lines[0]);
    for (int i = 0, j = 1; i < nX; i++, j++)
    {
        int countInterval = 0;
        double q = 0, step = 0;
        fin >> countInterval >> q;
        nXk += countInterval;
        grid.x.resize(nXk + 1);

        if (q != 1)
        {
            double sumProgression = (pow(q, countInterval) - 1.0) / (q - 1.0);
            step = (area.x_lines[i + 1] - area.x_lines[i]) / sumProgression;
            int jk = 1;
            for (j; j < nXk; j++, jk++)
                grid.x[j] = area.x_lines[i] + step * (pow(q, jk) - 1.0) / (q - 1.0);
        }
        else
        {
            step = (area.x_lines[i + 1] - area.x_lines[i]) / countInterval;
            int jk = 1;
            for (j; j < nXk; j++, jk++)
                grid.x[j] = area.x_lines[i] + step * jk;
        }
        grid.x[j] = area.x_lines[i + 1];
        area.Ix[i + 1] = j;
    }

    int nYk = 0;
    grid.y.resize(1, area.y_lines[0]);
    for (int i = 0, j = 1; i < nY; i++, j++)
    {
        int countInterval = 0;
        double q = 0, step = 0;
        fin >> countInterval >> q;
        nYk += countInterval;
        grid.y.resize(nYk + 1);

        if (q != 1)
        {
            double sumProgression = (pow(q, countInterval) - 1.0) / (q - 1.0);
            step = (area.y_lines[i + 1] - area.y_lines[i]) / sumProgression;
            int jk = 1;
            for (j; j < nYk; j++, jk++)
                grid.y[j] = area.y_lines[i] + step * (pow(q, jk) - 1.0) / (q - 1.0);
        }
        else
        {
            step = (area.y_lines[i + 1] - area.y_lines[i]) / countInterval;
            int jk = 1;
            for (j; j < nYk; j++, jk++)
                grid.y[j] = area.y_lines[i] + step * jk;
        }
        grid.y[j] = area.y_lines[i + 1];
        area.Iy[i + 1] = j;
    }
    fin.close();
}
