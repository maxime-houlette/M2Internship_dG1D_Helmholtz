#include "reference_element.h"
using namespace std;

ReferenceElement::ReferenceElement(size_t p) : degree(p) {}

double ReferenceElement::eval_basis(std::size_t k, double xi) const
{
    if (k == 0)
        return 1.0;
    if (k == 1)
        return xi;
    return std::pow(xi, k);
}

double ReferenceElement::eval_basis_deriv(std::size_t k, double xi) const
{
    if (k == 0)
        return 0.0;
    if (k == 1)
        return 1.0;
    return k * std::pow(xi, k - 1);
}