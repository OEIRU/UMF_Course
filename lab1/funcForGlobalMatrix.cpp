#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <functional>
#include "function.h"

void createPortrait(SLAE &slae, int nX, int nY)
{
    auto &di1 = slae.di1, &di2 = slae.di2, &diMain = slae.diMain, &di4 = slae.di4, &di5 = slae.di5, &b = slae.b;
    auto &ig = slae.ig;
    const int n = nX * nY;

    ig = {-(nX), -1, 1, nX};
    di1.resize(n - nX);
    di2.resize(n - 1);
    diMain.resize(n);
    di4.resize(n - 1);
    di5.resize(n - nX);
    b.resize(n);
}

void addElemInGlobalMatrix(SLAE &slae, int i, int j, double elem)
{
    auto &di1 = slae.di1, &di2 = slae.di2, &diMain = slae.diMain, &di4 = slae.di4, &di5 = slae.di5;
    auto &ig = slae.ig;

    if (i == j)
        diMain[i] += elem;
    else
    {
        int temp = j - i;
        if (temp == ig[0])
            di1[j] += elem;
        if (temp == ig[1])
            di2[j] += elem;
        if (temp == ig[2])
            di4[i] += elem;
        if (temp == ig[3])
            di5[i] += elem;
    }
}

void calcGlobalMatrixAndVector(SLAE &slae, grid_t &g, area_t &area, parameters_t &par)
{
    auto &di1 = slae.di1, &di2 = slae.di2, &diMain = slae.diMain, &di4 = slae.di4, &di5 = slae.di5;
    auto &ig = slae.ig;
    auto &x = g.x, &y = g.y, &b = slae.b;
    auto &lambda = par.lambda, &gamma = par.gamma;
    const int nX = x.size(), nY = y.size();
    int l = 0;
    auto &F = slae.F;

    for (int s = 1; s < nY - 1; s++)
    {
        for (int p = 1; p < nX - 1; p++)
        {
            int numArea = returnNumberOfEstimatedSubArea(area, p, s, l);

            if (numArea != -1)
            {
                std::vector<int> globalNumbers = {nX * (s - 1) + p, nX * s + p - 1, nX * s + p, nX * s + p + 1, nX * (s + 1) + p};
                double hxLeft = x[p] - x[p - 1], hxRight = x[p + 1] - x[p],
                       hyLeft = y[s] - y[s - 1], hyRight = y[s + 1] - y[s];

                double kxLeft = lambda * 2.0 / (hxLeft * (hxLeft + hxRight)), kx = lambda * 2.0 / (hxLeft * hxRight),
                       kxRight = lambda * 2.0 / (hxRight * (hxLeft + hxRight)), kyLeft = lambda * 2.0 / (hyLeft * (hyLeft + hyRight)),
                       ky = lambda * 2.0 / (hyLeft * hyRight), kyRight = lambda * 2.0 / (hyRight * (hyLeft + hyRight));
                addElemInGlobalMatrix(slae, globalNumbers[2], globalNumbers[0], -kyLeft);
                addElemInGlobalMatrix(slae, globalNumbers[2], globalNumbers[1], -kxLeft);
                addElemInGlobalMatrix(slae, globalNumbers[2], globalNumbers[2], kx + ky + gamma);
                addElemInGlobalMatrix(slae, globalNumbers[2], globalNumbers[3], -kxRight);
                addElemInGlobalMatrix(slae, globalNumbers[2], globalNumbers[4], -kyRight);
                b[globalNumbers[2]] += F[numArea - 1](x[p], y[s]);
            }
        }
    }
}

void clearFictiousNodes(SLAE &slae, area_t &area, int nX, int nY)
{
    auto &diMain = slae.diMain;

    for (int s = 0; s < nY; s++)
    {
        for (int p = 0; p < nX; p++)
        {
            int globalNum = s * nX + p;
            if (isFictiousNode(area, p, s))
                diMain[globalNum] = 1;
        }
    }
}

void calcBoundaryConditions(area_t &area, grid_t &g, SLAE &slae, parameters_t &par, std::vector<std::vector<int>> &bC, functionsBC &func)
{
    auto &ig = slae.ig;
    auto &x = g.x, &y = g.y, &di1 = slae.di1, &di2 = slae.di2,
         &diMain = slae.diMain, &di4 = slae.di4, &di5 = slae.di5, &b = slae.b;
    auto &lambda = par.lambda;
    auto &firstBC = func.firstBC, &secondBC = func.secondBC;
    const int countCond = bC.size(), nX = x.size(), nY = y.size(), n = nX * nY;

    for (int k = 0; k < countCond; k++)
    {
        const int typeCond = bC[k][0];
        if (typeCond == 2)
        {
            int numFunc = bC[k][1] - 1, p0 = bC[k][2], p1 = bC[k][3], s0 = bC[k][4], s1 = bC[k][5];
            if (p0 == p1)
            {
                int p = area.Ix[p0], s = area.Iy[s0];
                s1 = area.Iy[s1];
                for (s; s < s1; s++)
                {
                    int globalNum1 = s * nX + p, globalNum2 = s * nX + p + 1;
                    stringMatrixInNull(slae, globalNum1);
                    double temp = lambda / (x[p + 1] - x[p]);
                    addElemInGlobalMatrix(slae, globalNum1, globalNum1, temp);
                    addElemInGlobalMatrix(slae, globalNum1, globalNum2, -temp);
                    b[globalNum1] = secondBC[numFunc](x[p], y[s]);
                }
            }

            if (s0 == s1)
            {
                int p = area.Ix[p0], s = area.Iy[s0];
                p1 = area.Ix[p1];
                for (p; p < p1; p++)
                {
                    int globalNum1 = s * nX + p, globalNum2 = (s + 1) * nX + p;
                    stringMatrixInNull(slae, globalNum1);
                    double temp = lambda / (y[s + 1] - y[s]);
                    addElemInGlobalMatrix(slae, globalNum1, globalNum1, temp);
                    addElemInGlobalMatrix(slae, globalNum1, globalNum2, -temp);
                    b[globalNum1] = secondBC[numFunc](x[p], y[s]);
                }
            }
        }
    }

    for (int k = 0; k < countCond; k++)
    {
        int typeCond = bC[k][0];
        if (typeCond == 1)
        {
            int numFunc = bC[k][1] - 1, p0 = bC[k][2], p1 = bC[k][3], s0 = bC[k][4], s1 = bC[k][5];
            if (p0 == p1)
            {
                int p = area.Ix[p0], s = area.Iy[s0];
                s1 = area.Iy[s1];
                for (s; s <= s1; s++)
                {
                    int globalNum = s * nX + p;
                    if (globalNum + ig[0] >= 0)
                        di1[globalNum + ig[0]] = 0;
                    if (globalNum + ig[1] >= 0)
                        di2[globalNum + ig[1]] = 0;
                    diMain[globalNum] = 1;
                    if (globalNum + ig[2] < n)
                        di4[globalNum] = 0;
                    if (globalNum + ig[3] < n)
                        di5[globalNum] = 0;
                    b[globalNum] = firstBC[numFunc](x[p], y[s]);
                }
            }

            if (s0 == s1)
            {
                int p = area.Ix[p0], s = area.Iy[s0];
                p1 = area.Ix[p1];
                for (p; p <= p1; p++)
                {
                    int globalNum = s * nX + p;
                    if (globalNum + ig[0] >= 0)
                        di1[globalNum + ig[0]] = 0;
                    if (globalNum + ig[1] >= 0)
                        di2[globalNum + ig[1]] = 0;
                    diMain[globalNum] = 1;
                    if (globalNum + ig[2] < n)
                        di4[globalNum] = 0;
                    if (globalNum + ig[3] < n)
                        di5[globalNum] = 0;
                    b[globalNum] = firstBC[numFunc](x[p], y[s]);
                }
            }
        }
    }
}

void stringMatrixInNull(SLAE &slae, int i)
{
    auto &di1 = slae.di1, &di2 = slae.di2, &diMain = slae.diMain, &di4 = slae.di4, &di5 = slae.di5, &b = slae.b;
    auto &ig = slae.ig;
    const int n = diMain.size();

    if (i + ig[0] >= 0)
        di1[i + ig[0]] = 0;
    if (i + ig[1] >= 0)
        di2[i + ig[1]] = 0;
    diMain[i] = 0;
    if (i + ig[2] < n)
        di4[i] = 0;
    if (i + ig[3] < n)
        di5[i] = 0;
    b[i] = 0;
}
