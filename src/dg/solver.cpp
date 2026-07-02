#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <fstream>
#include "solver.h"
#include "mass.h"
#include "half_stiffness.h"
#include "quadrature.h"
#include <numbers>
using namespace std;

Solver::Solver(const Mesh &input_mesh, size_t degree, double input_kappa,
               double input_tau, RhsFunction input_f_rhs, RhsFunction input_u_exact, Problem input_problem)
    : mesh(input_mesh), ref_element(degree), kappa(input_kappa), tau(input_tau), problem(input_problem)
{
    /* ---------- MEMORY ALLOCATION ---------- */
    size_t n_local_dofs = ref_element.degree + 1;
    size_t n_elements = mesh.N;
    size_t total_dofs = n_elements * n_local_dofs;

    A.resize(2 * total_dofs, 2 * total_dofs);
    u.resize(2 * total_dofs);
    b = Eigen::VectorXcd::Zero(2 * total_dofs);

    size_t nnz_estime = n_elements * n_local_dofs * n_local_dofs * 3;
    tripletList.reserve(nnz_estime);

    f_rhs = input_f_rhs;
    u_exact = input_u_exact;

    /* ---------- SYSTEM CREATION ---------- */
    complex<double> alpha;
    double beta;
    double a_mean;
    double b_scal;
    double b_vec;
    double a_bc;
    double b_bc;
    if (problem == L2)
    {
        alpha = 1.;
        assemble_mass(alpha);
    }
    else
    {
        alpha = I * kappa;
        beta = -1.;
        a_mean = 1.;
        // b_scal = 0.;
        b_scal = tau / 2.;
        // b_vec = 0.;
        b_vec = 1. / (2. * tau);
        a_bc = 1.;
        // b_bc = 0.;
        b_bc = tau;
        assemble_mass(alpha);
        assemble_half_stiffness(beta);
        assemble_fluxes(a_mean, b_scal, b_vec);
        assemble_boundaries(a_bc, b_bc);
    }
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    assemble_rhs();
}

void Solver::assemble_mass(complex<double> alpha)
{
    size_t n_local_dofs = ref_element.degree + 1;
    size_t n_elements = mesh.N;

    Eigen::MatrixXcd M_loc(n_local_dofs, n_local_dofs);
    mass(ref_element, mesh.h, M_loc);

    for (size_t e = 0; e < n_elements; e++)
        for (size_t i = 0; i < n_local_dofs; i++)
            for (size_t j = 0; j < n_local_dofs; j++)
            {
                size_t I_u = local_to_global(e, U, i);
                size_t J_u = local_to_global(e, U, j);
                size_t I_q = local_to_global(e, Q, i);
                size_t J_q = local_to_global(e, Q, j);

                /* Mass u,u */
                tripletList.push_back(TripletComplex(I_u, J_u, alpha * M_loc(i, j)));

                /* Mass q,q */
                tripletList.push_back(TripletComplex(I_q, J_q, alpha * M_loc(i, j)));
            }
}

void Solver::assemble_half_stiffness(double beta)
{
    size_t n_local_dofs = ref_element.degree + 1;
    size_t n_elements = mesh.N;

    Eigen::MatrixXcd hS_loc(n_local_dofs, n_local_dofs);
    half_stiffness(ref_element, hS_loc);

    for (size_t e = 0; e < n_elements; e++)
        for (size_t i = 0; i < n_local_dofs; i++)
            for (size_t j = 0; j < n_local_dofs; j++)
            {
                /* Local to global numerotation */
                size_t I_u = local_to_global(e, U, i);
                size_t J_u = local_to_global(e, U, j);
                size_t I_q = local_to_global(e, Q, i);
                size_t J_q = local_to_global(e, Q, j);

                /* Half stiffness u,q */
                tripletList.push_back(TripletComplex(I_u, J_q, beta * hS_loc(j, i)));

                /* Half stiffness u,q */
                tripletList.push_back(TripletComplex(I_q, J_u, beta * hS_loc(j, i)));
            }
}

void Solver::assemble_fluxes(double a_mean, double b_scal, double b_vec)
{
    size_t n_local_dofs = ref_element.degree + 1;
    size_t n_elements = mesh.N;

    vector<double> phi_L(n_local_dofs);
    vector<double> phi_R(n_local_dofs);

    for (size_t i = 0; i < n_local_dofs; ++i)
    {
        phi_L[i] = ref_element.eval_basis(i, -1.0);
        phi_R[i] = ref_element.eval_basis(i, 1.0);
    }

    for (size_t e = 0; e < n_elements - 1; e++)
    {
        size_t e_L = e;
        size_t e_R = e + 1;
        for (size_t i = 0; i < n_local_dofs; i++)
            for (size_t j = 0; j < n_local_dofs; j++)
            {
                /* Local to global numerotation for the left element */
                size_t I_u_L = local_to_global(e_L, U, i);
                size_t J_u_L = local_to_global(e_L, U, j);
                size_t I_q_L = local_to_global(e_L, Q, i);
                size_t J_q_L = local_to_global(e_L, Q, j);

                /* Local to global numerotation for the right element */
                size_t I_u_R = local_to_global(e_R, U, i);
                size_t J_u_R = local_to_global(e_R, U, j);
                size_t I_q_R = local_to_global(e_R, Q, i);
                size_t J_q_R = local_to_global(e_R, Q, j);

                /* Q test function of the left element */
                tripletList.push_back(TripletComplex(I_q_L, J_u_L, a_mean * 0.5 * phi_R[j] * phi_R[i]));
                tripletList.push_back(TripletComplex(I_q_L, J_u_R, a_mean * 0.5 * phi_L[j] * phi_R[i]));
                tripletList.push_back(TripletComplex(I_q_L, J_q_L, b_vec * phi_R[j] * phi_R[i])); // 1/(2tau), sur J_q
                tripletList.push_back(TripletComplex(I_q_L, J_q_R, -b_vec * phi_L[j] * phi_R[i]));

                /* Q test function of the right element */
                tripletList.push_back(TripletComplex(I_q_R, J_u_L, -a_mean * 0.5 * phi_R[j] * phi_L[i]));
                tripletList.push_back(TripletComplex(I_q_R, J_u_R, -a_mean * 0.5 * phi_L[j] * phi_L[i]));
                tripletList.push_back(TripletComplex(I_q_R, J_q_L, -b_vec * phi_R[j] * phi_L[i]));
                tripletList.push_back(TripletComplex(I_q_R, J_q_R, b_vec * phi_L[j] * phi_L[i]));

                /* U test function of the left element */
                tripletList.push_back(TripletComplex(I_u_L, J_q_L, a_mean * 0.5 * phi_R[j] * phi_R[i]));
                tripletList.push_back(TripletComplex(I_u_L, J_q_R, a_mean * 0.5 * phi_L[j] * phi_R[i]));
                tripletList.push_back(TripletComplex(I_u_L, J_u_L, b_scal * phi_R[j] * phi_R[i])); // tau/2, sur J_u
                tripletList.push_back(TripletComplex(I_u_L, J_u_R, -b_scal * phi_L[j] * phi_R[i]));

                /* U test function of the right element */
                tripletList.push_back(TripletComplex(I_u_R, J_q_L, -a_mean * 0.5 * phi_R[j] * phi_L[i]));
                tripletList.push_back(TripletComplex(I_u_R, J_q_R, -a_mean * 0.5 * phi_L[j] * phi_L[i]));
                tripletList.push_back(TripletComplex(I_u_R, J_u_L, -b_scal * phi_R[j] * phi_L[i]));
                tripletList.push_back(TripletComplex(I_u_R, J_u_R, b_scal * phi_L[j] * phi_L[i]));
            }
    }
}

void Solver::assemble_boundaries(double a_bc, double b_bc)
{
    size_t n_local_dofs = ref_element.degree + 1;

    vector<double> phi_L(n_local_dofs);
    vector<double> phi_R(n_local_dofs);

    for (size_t i = 0; i < n_local_dofs; ++i)
    {
        phi_L[i] = ref_element.eval_basis(i, -1.0);
        phi_R[i] = ref_element.eval_basis(i, 1.0);
    }
    for (size_t i = 0; i < n_local_dofs; i++)
        for (size_t j = 0; j < n_local_dofs; j++)
        {
            /* Local to global for the boundary x=0 */
            size_t I_u_0 = local_to_global(0, U, i);
            size_t J_u_0 = local_to_global(0, U, j);
            size_t J_q_0 = local_to_global(0, Q, j);

            /* Local to global for the boundary x=L */
            size_t I_u_N = local_to_global(mesh.N - 1, U, i);
            size_t J_u_N = local_to_global(mesh.N - 1, U, j);
            size_t J_q_N = local_to_global(mesh.N - 1, Q, j);

            /* Boundary x=0 */
            tripletList.push_back(TripletComplex(I_u_0, J_q_0, -a_bc * phi_L[j] * phi_L[i]));
            tripletList.push_back(TripletComplex(I_u_0, J_u_0, b_bc * phi_L[j] * phi_L[i]));

            /* Boundary x=L */
            tripletList.push_back(TripletComplex(I_u_N, J_q_N, a_bc * phi_R[j] * phi_R[i]));
            tripletList.push_back(TripletComplex(I_u_N, J_u_N, b_bc * phi_R[j] * phi_R[i]));
        }
}

void Solver::assemble_rhs()
{
    size_t n_local_dofs = ref_element.degree + 1;
    size_t n_elements = mesh.N;
    double J = mesh.h / 2;

    for (size_t e = 0; e < n_elements; e++)
    {
        double x_c = 0.5 * (mesh.positions[e] + mesh.positions[e + 1]);
        QuadratureRule quad_rule = get_gauss_legendre(n_local_dofs);

        for (size_t q = 0; q < n_local_dofs; q++)
        {
            double xi_q = quad_rule.points[q];
            double w_q = quad_rule.weights[q];

            double x_q = x_c + J * xi_q;

            complex<double> f_val = f_rhs(x_q);

            for (size_t i = 0; i < n_local_dofs; i++)
            {
                size_t I_u = local_to_global(e, U, i);
                double val_i = ref_element.eval_basis(i, xi_q);

                b(I_u) += J * w_q * f_val * val_i;
            }
        }
    }
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
    std::cout << "LU factorization in progress..." << std::endl;
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> linear_solver;
    linear_solver.compute(A);

    if (linear_solver.info() != Eigen::Success)
    {
        std::cerr << "FATAL ERROR: LU factorization has failed."
                  << "Check the implementation of the fluxes." << std::endl;
        return;
    }

    std::cout << "Resolution of the linear system..." << std::endl;
    this->u = linear_solver.solve(b);

    if (linear_solver.info() != Eigen::Success)
    {
        std::cerr << "FATAL ERROR: Linear resolution has failed." << std::endl;
        return;
    }

    std::cout << "Success!" << std::endl;
}

void Solver::post_process()
{
    std::cout << "Exportation of the solution in solution.csv..." << std::endl;
    std::ofstream out("solution_u.csv");
    out << "x,u_real,u_imag,rhs_real,rhs_imag\n";

    // Pour chaque élément, on va échantillonner la solution (ex: 5 points par élément)
    size_t n_plot_points = 5;
    for (size_t e = 0; e < mesh.N; e++)
    {
        double x_L = mesh.positions[e];
        double x_R = mesh.positions[e + 1];

        for (size_t p = 0; p <= n_plot_points; p++)
        {
            double xi = -1.0 + 2.0 * (double)p / n_plot_points;
            double x_phys = x_L + (xi + 1.0) * (x_R - x_L) / 2.0;

            std::complex<double> u_val(0.0, 0.0);
            for (size_t i = 0; i < ref_element.degree + 1; i++)
            {
                size_t I_u = local_to_global(e, U, i);
                u_val += u(I_u) * ref_element.eval_basis(i, xi);
            }

            std::complex<double> rhs_val = f_rhs(x_phys);

            out << x_phys << ","
                << u_val.real() << ","
                << u_val.imag() << ","
                << rhs_val.real() << ","
                << rhs_val.imag() << "\n";
        }
    }
    out.close();
    std::cout << "Success!" << std::endl;
}

void Solver::plot_u()
{
    // --- NOUVEAU : Lancement de Gnuplot ---
    std::cout << "GNUPlot launching..." << std::endl;
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

    if (gnuplotPipe)
    {
        fprintf(gnuplotPipe, "set title 'Solution de l\\'équation de Helmholtz (dG 1D)'\n");
        fprintf(gnuplotPipe, "set xlabel 'Position x'\n");
        fprintf(gnuplotPipe, "set ylabel 'Amplitude'\n");
        fprintf(gnuplotPipe, "set grid\n");
        fprintf(gnuplotPipe, "set datafile separator ','\n"); // Important pour lire votre CSV

        fprintf(gnuplotPipe, "plot 'solution_u.csv' every ::1 using 1:2 with linespoints pt 7 ps 0.5 lc rgb 'blue' title 'Re(u)', ");
        fprintf(gnuplotPipe, "'' every ::1 using 1:3 with linespoints pt 7 ps 0.5 lc rgb 'orange' title 'Im(u)'\n");

        pclose(gnuplotPipe);
    }
    else
        std::cerr << "FATAL ERROR: unable to plot." << std::endl;
}

void Solver::compute_L2_error()
{
    std::cout << "L2 error computation..." << std::endl;

    double error_L2_sq = 0.0;

    size_t n_quad_points = ref_element.degree + 2;
    QuadratureRule quad_rule = get_gauss_legendre(n_quad_points);

    for (size_t e = 0; e < mesh.N; e++)
    {
        double x_L = mesh.positions[e];
        double J = mesh.h / 2.0;

        for (size_t q = 0; q < n_quad_points; q++)
        {
            double xi = quad_rule.points[q];
            double w_q = quad_rule.weights[q];

            double x_phys = x_L + (xi + 1.0) * J;

            std::complex<double> u_num(0.0, 0.0);
            for (size_t i = 0; i < ref_element.degree + 1; i++)
            {
                size_t I_u = local_to_global(e, U, i);
                u_num += u(I_u) * ref_element.eval_basis(i, xi);
            }

            if (!u_exact)
                std::cerr << "FATAL ERROR: no exact solution has been provided." << std::endl;
            std::complex<double> u_ex = u_exact(x_phys);

            std::complex<double> diff = u_num - u_ex;
            double diff_sq = std::norm(diff);

            error_L2_sq += J * w_q * diff_sq;
        }
    }
    double error_L2 = std::sqrt(error_L2_sq);

    std::string filename;
    if (problem == L2)
        filename = "error_proj_L2.csv";
    else
        filename = "error_helmholtz.csv";
    std::ifstream check_file(filename);
    bool file_exists = check_file.good();
    check_file.close();

    std::ofstream out(filename, std::ios::app);
    if (!file_exists)
    {
        out << "N,degree,kappa/Pi,error_L2\n";
    }

    out << mesh.N << ","
        << ref_element.degree << ","
        << kappa / M_PI << ","
        << std::scientific << error_L2 << "\n";

    out.close();
}

void Solver::plot_error()
{
    // --- NOUVEAU : Lancement de Gnuplot ---
    std::cout << "GNUPlot launching..." << std::endl;
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

    if (gnuplotPipe)
    {
        fprintf(gnuplotPipe, "set title 'Erreur de l equation de Helmholtz (dG 1D)'\n");
        fprintf(gnuplotPipe, "set xlabel 'Position x'\n");
        fprintf(gnuplotPipe, "set ylabel 'Error'\n");
        fprintf(gnuplotPipe, "set grid\n");
        fprintf(gnuplotPipe, "set datafile separator ','\n"); // Important pour lire votre CSV

        fprintf(gnuplotPipe, "plot 'error_helmholtz.csv' every ::1 using 1:4 with linespoints pt 7 ps 0.5 lc rgb 'blue' title 'Re(u)', ");

        pclose(gnuplotPipe);
    }
    else
        std::cerr << "FATAL ERROR: unable to plot." << std::endl;
}
