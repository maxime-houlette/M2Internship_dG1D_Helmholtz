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

enum Problem
{
    L2 = 0,
    Primitive = 1,
    HelmHoltz = 2
};

typedef Eigen::Triplet<std::complex<double>> TripletComplex;

static std::complex<double> I(0.0, 1.0);

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

    std::vector<TripletComplex> tripletList;

    // --- Méthodes internes ---
    void assemble_L2();
    void assemble_primitive();
    void assemble_Helmholtz();
    void assemble_rhs();
    size_t local_to_global(size_t element, Variable var, size_t local_dof) const;

public:
    // --- Constructeur ---
    Solver(const Mesh &input_mesh, size_t degree, double input_kappa, double input_tau, Problem problem);

    // --- Méthodes accessibles ---
    void solve();
    void post_process();
    void plot_u();
    void compute_L2_error();
    const Eigen::VectorXcd &get_solution() const { return u; }
};