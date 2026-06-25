#pragma once
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "quadrature.h"
#include "reference_element.h"

inline void half_stiffness(const ReferenceElement &ref_el, double h, Eigen::MatrixXcd &hS_loc)
{
    size_t n_dofs = ref_el.degree + 1;
    hS_loc.setZero();

    QuadratureRule quad_rule = get_gauss_legendre(n_dofs);

    for (size_t q = 0; q < n_dofs; q++)
    {
        double point_q = quad_rule.points[q];
        double weight_q = quad_rule.weights[q];
        for (size_t i = 0; i < n_dofs; i++)
            for (size_t j = 0; j < n_dofs; j++)
            {
                double val_i = ref_el.eval_basis(i, point_q);
                double val_deriv_j = ref_el.eval_basis_deriv(j, point_q);
                hS_loc(i, j) += weight_q * val_i * val_deriv_j;
            }
        }
}