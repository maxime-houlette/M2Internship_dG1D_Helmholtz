#pragma once

#include <complex>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "mesh.h"
#include "reference_element.h"

enum Variable
{
    Q = 0,
    U = 1
};

class Solver
{
private:
    // --- Attributs cachés à l'utilisateur ---
    Mesh mesh;
    ReferenceElement ref_element;
    double kappa;
    double tau;

    Eigen::SparseMatrix<std::complex<double>> A;
    Eigen::VectorXcd b;
    Eigen::VectorXcd u;

    // --- Méthodes internes ---
    void assemble_system();
    size_t local_to_global(size_t element, Variable var, size_t local_dof) const;

public:
    // --- Constructeur ---
    Solver(const Mesh &input_mesh, size_t degree, double input_kappa, double input_tau);

    // --- Méthodes accessibles ---
    void solve();
    const Eigen::VectorXcd &get_solution() const { return u; }
};