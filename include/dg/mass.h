#pragma once
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "quadrature.h"
#include "reference_element.h"

inline void mass(const ReferenceElement &ref_el, double h, Eigen::MatrixXcd &M_loc)
{
    size_t n_dofs = ref_el.degree + 1;
    M_loc.setZero();

    QuadratureRule quad_rule = get_gauss_legendre(n_dofs);

    for (size_t i = 0; i < n_dofs; i++)
        for (size_t j = 0; j < n_dofs; j++)
            for (size_t q = 0; q < n_dofs; q++)
            {
                double val_i = ref_el.eval_basis(i, quad_rule.points[q]);
                double val_j = ref_el.eval_basis(j, quad_rule.points[q]);
                M_loc(i, j) += h / 2 * quad_rule.weights[q] * val_i * val_j;
            }
}