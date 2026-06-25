#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>

struct QuadratureRule
{
    std::vector<double> points;
    std::vector<double> weights;
};

inline QuadratureRule get_gauss_legendre(int num_points)
{
    QuadratureRule rule;
    rule.points.resize(num_points);
    rule.weights.resize(num_points);

    // Cas trivial pour 1 seul point
    if (num_points == 1)
    {
        rule.points[0] = 0.0;
        rule.weights[0] = 2.0;
        return rule;
    }

    // Construction de la matrice tridiagonale (Matrice de Jacobi)
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(num_points, num_points);
    for (int i = 0; i < num_points - 1; ++i)
    {
        double beta = (i + 1.0) / std::sqrt(4.0 * (i + 1.0) * (i + 1.0) - 1.0);
        J(i, i + 1) = beta;
        J(i + 1, i) = beta;
    }

    // Calcul des valeurs propres (points) et vecteurs propres (pour les poids)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(J);

    // Extraction des résultats
    for (int i = 0; i < num_points; ++i)
    {
        rule.points[i] = solver.eigenvalues()(i);
        // Le poids est 2 * (première composante du vecteur propre)^2
        rule.weights[i] = 2.0 * std::pow(solver.eigenvectors()(0, i), 2);
    }

    return rule;
}