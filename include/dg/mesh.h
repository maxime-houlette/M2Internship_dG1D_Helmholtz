#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <cassert>


struct Mesh
{
    size_t N;
    double L;
    double h;
    std::vector<double> positions;
    Mesh(size_t nb_elements, double length);
};
