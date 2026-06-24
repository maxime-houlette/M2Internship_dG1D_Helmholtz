#include "mesh.h"

Mesh::Mesh(size_t nb_elements, double length)
{
    N = nb_elements;
    L = length;
    h = L / N;
    positions.reserve(nb_elements + 1);

    for (size_t i = 0; i < nb_elements + 1; i++)
        positions.push_back(i * h);
}