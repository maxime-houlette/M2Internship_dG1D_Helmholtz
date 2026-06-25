#pragma once
#include <cstddef>
#include <cmath>

struct ReferenceElement
{
    size_t degree; // Le degré p des polynômes

    ReferenceElement(std::size_t p);

    // Évalue le k-ième polynôme de la base au point local xi (entre -1 et 1)
    double eval_basis(std::size_t k, double xi) const;

    // Évalue la dérivée du k-ième polynôme de la base au point local xi
    double eval_basis_deriv(std::size_t k, double xi) const;
};