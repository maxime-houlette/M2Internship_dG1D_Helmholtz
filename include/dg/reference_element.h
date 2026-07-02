#pragma once
#include <cstddef>
#include <cmath>

struct ReferenceElement
{
    size_t degree; // Le degré p des polynômes

    ReferenceElement(std::size_t p);
    double eval_basis(std::size_t k, double xi) const;
    double eval_basis_deriv(std::size_t k, double xi) const;
};