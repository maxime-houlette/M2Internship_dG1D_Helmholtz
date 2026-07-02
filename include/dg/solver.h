#pragma once

#include <complex>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "mesh.h"
#include "reference_element.h"
#include <functional>
using RhsFunction = std::function<std::complex<double>(double)>;

enum Variable
{
    Q = 0,
    U = 1
};

enum Problem
{
    L2 = 0,
    Helmholtz = 1
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

    Problem problem;

    std::vector<TripletComplex> tripletList;

    RhsFunction f_rhs = nullptr;
    RhsFunction u_exact = nullptr;

    // --- Méthodes internes ---
    void assemble_mass(std::complex<double> alpha);
    void assemble_half_stiffness(double beta);
    void assemble_fluxes(double a_mean, double b_scal, double b_vec);
    void assemble_boundaries(double a_bc, double b_bc);
    void assemble_rhs();
    size_t local_to_global(size_t element, Variable var, size_t local_dof) const;

public:
    // --- Constructeur ---
    Solver(const Mesh &input_mesh,
           size_t degree,
           double input_kappa,
           double input_tau,
           RhsFunction input_f_rhs,
           RhsFunction input_u_exact,
           Problem problem);

    // --- Méthodes accessibles ---
    void solve();
    void post_process();
    void plot_u();
    void compute_L2_error();
    void plot_error();
    const Eigen::VectorXcd &get_solution() const { return u; }
};