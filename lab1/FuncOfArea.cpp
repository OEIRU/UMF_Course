#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <functional>
#include "function.h"

int returnNumberOfEstimatedSubArea(area_t &area, int p, int s, int &l)
{
    const int L = area.elems.size(), L1 = l;

    for (l; l < L; l++)
    {
        int mx0 = area.Ix[area.elems[l].x1], mx1 = area.Ix[area.elems[l].x2],
            my0 = area.Iy[area.elems[l].y1], my1 = area.Iy[area.elems[l].y2], m = area.elems[l].func_number;
        if (mx0 <= p && p <= mx1 && mx0 <= (p + 1) && (p + 1) <= mx1 &&
            my0 <= s && s <= my1 && my0 <= (s + 1) && (s + 1) <= my1)
            return m;
    }

    for (l = 0; l < L1; l++)
    {
        int mx0 = area.Ix[area.elems[l].x1], mx1 = area.Ix[area.elems[l].x2],
            my0 = area.Iy[area.elems[l].y1], my1 = area.Iy[area.elems[l].y2], m = area.elems[l].func_number;
        if (mx0 <= p && p <= mx1 && mx0 <= (p + 1) && (p + 1) <= mx1 &&
            my0 <= s && s <= my1 && my0 <= (s + 1) && (s + 1) <= my1)
            return m;
    }
    return -1;
}

bool isFictiousNode(area_t &area, int p, int s)
{
    const int L = area.elems.size();

    for (int l = 0; l < L; l++)
    {
        int mx0 = area.Ix[area.elems[l].x1], mx1 = area.Ix[area.elems[l].x2],
            my0 = area.Iy[area.elems[l].y1], my1 = area.Iy[area.elems[l].y2], m = area.elems[l].func_number;
        if (mx0 <= p && p <= mx1 &&
            my0 <= s && s <= my1)
            return false;
    }
    return true;
}
