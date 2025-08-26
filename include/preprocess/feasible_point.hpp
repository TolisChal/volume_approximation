// VolEsti (volume computation and sampling library)

// Copyright (c) 2025 Vissarion Fisikopoulos
// Copyright (c) 2025 Apostolos Chalkis
// Copyright (c) 2025 Elias Tsigaridas

// Licensed under GNU LGPL.3, see LICENCE file


#ifndef FEASIBLE_POINT_HPP
#define FEASIBLE_POINT_HPP

#include <tuple>

#include "preprocess/max_inscribed_ball.hpp"

// Using MT as to deal with both dense and sparse matrices
template <typename MT, typename VT>
VT compute_feasible_point(MT const& A, VT const& b)
{
    typedef typename VT::value_type NT;
    VT x;
    bool feasibility_only = true, converged;
    unsigned max_iters = 10000;
    const NT tol = 1e-08;
    // Compute a feasible point
    std::tie(x, std::ignore, converged) = max_inscribed_ball(A, b, max_iters, tol, feasibility_only);
    if (!converged || ((A * x).array() > b.array()).any())
    {
        std::runtime_error("The computation of a feasible point failed.");
    }
    return x;
}

#endif
