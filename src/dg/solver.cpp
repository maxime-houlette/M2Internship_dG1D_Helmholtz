#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <fstream>
#include "solver.h"
#include "mass.h"
#include "half_stiffness.h"
#include "quadrature.h"
using namespace std;

Solver::Solver(const Mesh &input_mesh, size_t degree, double input_kappa,
               double input_tau, Problem problem)
    : mesh(input_mesh), ref_element(degree), kappa(input_kappa), tau(input_tau)
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

    /* ---------- SYSTEM CREATION ---------- */
    if (problem == L2)
        assemble_L2();
    else if (problem == Primitive)
        assemble_primitive();
    else
        assemble_Helmholtz();
    assemble_rhs();
}

void Solver::assemble_L2()
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
                tripletList.push_back(TripletComplex(I_u, J_u, M_loc(i, j)));

                /* Mass q,q */
                tripletList.push_back(TripletComplex(I_q, J_q, M_loc(i, j)));
            }
    A.setFromTriplets(tripletList.begin(), tripletList.end());
}

void Solver::assemble_primitive()
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

    Eigen::MatrixXcd hS_loc(n_local_dofs, n_local_dofs);
    half_stiffness(ref_element, hS_loc);

    for (size_t e = 0; e < n_elements; e++)
        for (size_t i = 0; i < n_local_dofs; i++)
            for (size_t j = 0; j < n_local_dofs; j++)
            {
                size_t I_u = local_to_global(e, U, i);
                size_t J_u = local_to_global(e, U, j);
                size_t I_q = local_to_global(e, Q, i);
                size_t J_q = local_to_global(e, Q, j);

                /* Half stiffness u,q */
                tripletList.push_back(TripletComplex(I_u, J_q, -hS_loc(j, i)));

                /* Half stiffness u,q */
                tripletList.push_back(TripletComplex(I_q, J_u, -hS_loc(j, i)));

                /* TOD: Add fluxes */
            }
    A.setFromTriplets(tripletList.begin(), tripletList.end());
}

void Solver::assemble_Helmholtz()
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

    Eigen::MatrixXcd M_loc(n_local_dofs, n_local_dofs);
    mass(ref_element, mesh.h, M_loc);

    Eigen::MatrixXcd hS_loc(n_local_dofs, n_local_dofs);
    half_stiffness(ref_element, hS_loc);

    for (size_t e = 0; e < n_elements; e++)
        for (size_t i = 0; i < n_local_dofs; i++)
            for (size_t j = 0; j < n_local_dofs; j++)
            {
                size_t I_u = local_to_global(e, U, i);
                size_t J_u = local_to_global(e, U, j);
                size_t I_q = local_to_global(e, Q, i);
                size_t J_q = local_to_global(e, Q, j);

                /* Mass u,u */
                tripletList.push_back(TripletComplex(I_u, J_u, I * kappa * M_loc(i, j)));

                /* Mass q,q */
                tripletList.push_back(TripletComplex(I_q, J_q, I * kappa * M_loc(i, j)));

                /* Half stiffness u,q */
                tripletList.push_back(TripletComplex(I_u, J_q, -hS_loc(j, i)));

                /* Half stiffness u,q */
                tripletList.push_back(TripletComplex(I_q, J_u, -hS_loc(j, i)));

                /* TOD: Add fluxes */
            }
    A.setFromTriplets(tripletList.begin(), tripletList.end());
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

            /* TOD: sortir le choix de f de lal boucle */
            complex<double> f_val = sin(kappa * x_q);

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

    // 2. Analyse et factorisation de la matrice A
    linear_solver.compute(A);

    // 3. Vérification du succès de la factorisation (Crucial !)
    if (linear_solver.info() != Eigen::Success)
    {
        std::cerr << "FATAL ERROR: LU factorization has failed."
                  << "Check the implementation of the fluxes." << std::endl;
        return;
    }

    std::cout << "Resolution of the linear system..." << std::endl;

    // 4. Résolution de Ax = b
    this->u = linear_solver.solve(b);

    // 5. Vérification post-résolution
    if (linear_solver.info() != Eigen::Success)
    {
        std::cerr << "FATAL ERROR: Linear resolution has failed." << std::endl;
        return;
    }

    std::cout << "Success!" << std::endl;
}

void Solver::post_process()
{
    std::cout << "Exportation de la solution vers solution_u.csv..." << std::endl;
    std::ofstream out("solution_u.csv");

    // Modification : Ajout des en-têtes pour le RHS
    out << "x,u_real,u_imag,rhs_real,rhs_imag\n";

    // Pour chaque élément, on va échantillonner la solution (ex: 5 points par élément)
    size_t n_plot_points = 5;
    for (size_t e = 0; e < mesh.N; e++)
    {
        double x_L = mesh.positions[e];
        double x_R = mesh.positions[e + 1];

        for (size_t p = 0; p <= n_plot_points; p++)
        {
            // Point local entre -1 et 1
            double xi = -1.0 + 2.0 * (double)p / n_plot_points;
            // Point physique correspondant
            double x_phys = x_L + (xi + 1.0) * (x_R - x_L) / 2.0;

            // On reconstruit la valeur u(x) grâce aux fonctions de base
            std::complex<double> u_val(0.0, 0.0);
            for (size_t i = 0; i < ref_element.degree + 1; i++)
            {
                size_t I_u = local_to_global(e, U, i);
                u_val += u(I_u) * ref_element.eval_basis(i, xi);
            }

            // Évaluation de la fonction source analytique (RHS)
            // Note : Remplacez 'source_function' par la vraie fonction ou méthode
            // qui définit votre terme source f(x) dans votre code.
            std::complex<double> rhs_val = sin(kappa * x_phys);

            // Modification : Écriture des nouvelles colonnes
            out << x_phys << ","
                << u_val.real() << ","
                << u_val.imag() << ","
                << rhs_val.real() << ","
                << rhs_val.imag() << "\n";
        }
    }
    out.close();
    std::cout << "Exportation terminée." << std::endl;
}

void Solver::plot_u()
{
    // --- NOUVEAU : Lancement de Gnuplot ---
    std::cout << "Lancement de Gnuplot..." << std::endl;

    // On ouvre un 'pipe' vers Gnuplot. L'option -persistent permet
    // à la fenêtre de rester ouverte même après la fin du programme C++
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

    if (gnuplotPipe)
    {
        // On dicte les commandes à Gnuplot comme si on tapait dans son terminal
        fprintf(gnuplotPipe, "set title 'Solution de l\\'équation de Helmholtz (dG 1D)'\n");
        fprintf(gnuplotPipe, "set xlabel 'Position x'\n");
        fprintf(gnuplotPipe, "set ylabel 'Amplitude'\n");
        fprintf(gnuplotPipe, "set grid\n");
        fprintf(gnuplotPipe, "set datafile separator ','\n"); // Important pour lire votre CSV

        // La commande de tracé (similaire à vos subplots en Python)
        fprintf(gnuplotPipe, "plot 'solution_u.csv' every ::1 using 1:2 with linespoints pt 7 ps 0.5 lc rgb 'blue' title 'Re(u)', ");
        fprintf(gnuplotPipe, "'' every ::1 using 1:3 with linespoints pt 7 ps 0.5 lc rgb 'orange' title 'Im(u)'\n");

        // On ferme le flux (la commande est envoyée)
        pclose(gnuplotPipe);
    }
    else
        std::cerr << "Erreur : Impossible d'ouvrir Gnuplot. Est-il bien installé ?" << std::endl;
}

void Solver::compute_L2_error()
{
    std::cout << "Calcul de l'erreur L2..." << std::endl;

    double error_L2_sq = 0.0;

    // On prend un peu plus de points de quadrature pour intégrer
    // précisément la fonction exacte (qui n'est pas polynomiale)
    size_t n_quad_points = ref_element.degree + 2;
    QuadratureRule quad_rule = get_gauss_legendre(n_quad_points);

    for (size_t e = 0; e < mesh.N; e++)
    {
        double x_L = mesh.positions[e];
        double x_R = mesh.positions[e + 1];
        double J = mesh.h / 2.0;

        for (size_t q = 0; q < n_quad_points; q++)
        {
            double xi = quad_rule.points[q];
            double w_q = quad_rule.weights[q];

            // 1. Position physique du point de quadrature
            double x_phys = x_L + (xi + 1.0) * J;

            // 2. Reconstruire la solution numérique au point de quadrature
            std::complex<double> u_num(0.0, 0.0);
            for (size_t i = 0; i < ref_element.degree + 1; i++)
            {
                size_t I_u = local_to_global(e, U, i);
                u_num += u(I_u) * ref_element.eval_basis(i, xi);
            }

            // 3. Évaluer la solution analytique exacte
            std::complex<double> u_ex = sin(kappa * x_phys);

            // 4. Calcul de l'écart au carré
            std::complex<double> diff = u_num - u_ex;

            // L'astuce C++ : std::norm(z) renvoie (Re^2 + Im^2),
            // c'est-à-dire le module au carré |z|^2 !
            double diff_sq = std::norm(diff);

            // 5. Accumuler l'intégrale (dx = J * w_q)
            error_L2_sq += J * w_q * diff_sq;
        }
    }

    // On prend la racine carrée de la somme totale
    double error_L2 = std::sqrt(error_L2_sq);

    std::string filename = "error_proj_L2.csv";
    std::ifstream check_file(filename);
    bool file_exists = check_file.good();
    check_file.close();

    std::ofstream out(filename, std::ios::app);
    if (!file_exists)
    {
        out << "N,degree,kappa,error_L2\n";
    }

    out << mesh.N << ","
        << ref_element.degree << ","
        << kappa << ","
        << std::scientific << error_L2 << "\n";

    out.close();
}
