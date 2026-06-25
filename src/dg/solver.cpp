#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include "solver.h"
#include "mass.h"
#include "quadrature.h"
using namespace std;

typedef Eigen::Triplet<complex<double>> TripletComplex;

Solver::Solver(const Mesh &input_mesh, size_t degree, double input_kappa, double input_tau)
    : mesh(input_mesh), ref_element(degree), kappa(input_kappa)
{
    /* ---------- MEMORY ALLOCATION ---------- */
    size_t n_local_dofs = ref_element.degree + 1;
    size_t n_elements = mesh.N; 
    size_t total_dofs = n_elements * n_local_dofs;

    A.resize(2 * total_dofs, 2 * total_dofs);
    u.resize(2 * total_dofs);
    b.resize(2 * total_dofs);

    vector<TripletComplex> tripletList;
    size_t nnz_estime = n_elements * n_local_dofs * n_local_dofs * 3;
    tripletList.reserve(nnz_estime);

    /* ---------- SYSTEM CREATION ---------- */
    assemble_system();
}

void Solver::assemble_system()
{
}

size_t Solver::local_to_global(size_t element, Variable var, size_t local_dof) const
{
    size_t n_local_dofs = ref_element.degree + 1; 
    size_t block_size = 2 * n_local_dofs;
    size_t offset = block_size * element;

    if (var == Q)
        offset += n_local_dofs;
    
    return offset + local_dof;
}

void Solver::solve()
{
}
