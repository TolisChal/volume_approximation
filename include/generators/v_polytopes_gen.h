// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef V_POLYTOPES_GEN_H
#define V_POLYTOPES_GEN_H

#include <exception>
#include "samplers.h"


template <class MT>
void removeRow(MT &matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);

    matrix.conservativeResize(numRows,numCols);
}

template <class Polytope, class RNGType>
Polytope random_vpoly(unsigned int dim, unsigned int k) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;
    typedef typename Polytope::PolytopePoint Point;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);

    Point p;
    typename std::vector<NT>::iterator pit;
    MT V(k, dim);
    unsigned int j;

    for (unsigned int i = 0; i < k; ++i) {
        p = get_direction<RNGType, Point, NT>(dim);
        V.row(i) = p.getCoefficients();
    }

    Polytope VP;
    VT b = VT::Ones(k);
    VP.init(dim, V, b);

    return VP;

}


template <class Polytope, class RNGType>
Polytope random_vpoly_incube(unsigned int d, unsigned int k) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::NT    NT;
    typedef typename Polytope::PolytopePoint Point;

    REAL *conv_mem;
    int *colno_mem;

    conv_mem = (REAL *) malloc(k * sizeof(*conv_mem));
    colno_mem = (int *) malloc(k * sizeof(*colno_mem));

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::random::uniform_real_distribution<> urdist1(-1, 1);

    Point p(d);
    typename std::vector<NT>::iterator pit;
    MT V(k, d);
    unsigned int j, count_row,it=0;
    std::vector<int> indices;
    Polytope VP;
    VT b = VT::Ones(k);

    for (unsigned int i = 0; i < k; ++i) {
        for (int j = 0; j < d; ++j) {
            V(i, j) = urdist1(rng);
        }
    }
    if(k==d+1){
        VP.init(d, V, b);
        return VP;
    }

    MT V2(k,d);
    V2 = V;
    indices.clear();
    while(it<20) {
        V.resize(V2.rows(), d);
        V = V2;
        for (int i = 0; i < indices.size(); ++i) {
            V.conservativeResize(V.rows()+1, d);
            for (int j = 0; j < d; ++j) {
                V(V.rows()-1, j) = urdist1(rng);
            }
        }
        indices.clear();
        V2.resize(k, d);
        V2 = V;

        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < d; ++j) {
                p.set_coord(j, V(i, j));
            }
            removeRow(V2, i);
            if (memLP_Vpoly(V2, p, conv_mem, colno_mem)){
                indices.push_back(i);
            }
            V2.resize(k, d);
            V2 = V;
        }
        if (indices.size()==0) {
            VP.init(d, V, b);
            return VP;
        }
        V2.resize(k - indices.size(), d);
        count_row =0;
        for (int i = 0; i < k; ++i) {
            if(std::find(indices.begin(), indices.end(), i) != indices.end()) {
                continue;
            } else {
                for (int j = 0; j < d; ++j) V2(count_row, j) = V(i,j);
                count_row++;
            }
        }
        it++;
    }

    VP.init(d, V2, VT::Ones(V2.rows()));
    free(colno_mem);
    free(conv_mem);

    return VP;

}


template <class Polytope, class RNGType>
Polytope dual_knapsack(int d) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef typename Polytope::PolytopePoint Point;
    typedef typename Point::FT NT;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    RNGType rng(seed);
    boost::random::uniform_int_distribution<> uidist(-100, 100);

    Polytope P, Ptemp;
    MT A;
    VT b;

    A.resize(d, 2 * d);
    b.resize(2 * d);

    A << MT::Identity(d,d), -MT::Identity(d,d);
    b.setConstant(NT(1));

    MT A2(2 * d, d);
    A2 = A.transpose();

    Ptemp.init(d, A2, b);
    Point p(d);

    while (true) {

        for (int i = 0; i < d; ++i) {
            p.set_coord(i, uidist(rng));
        }
        //p = get_point_in_Dsphere<RNGType, Point >(d, NT(5));
        if (Ptemp.is_in(p) == 0) {
            Ptemp.free_them_all();
            break;
        }

    }

    A2.conservativeResize(A2.rows()+1, A2.cols());
    A2.row(2 * d) = p.getCoefficients();

    P.init(d, A2, b);
    return P;

}


template <typename VT>
bool update_everest_indices(int n, int s, VT &indices) {

    for (int i = 0; i < n; ++i) {
        indices(i) = indices(i) + 1;
        if (indices(i) > s) {
            indices(i) = indices(i) % (s+1);
            if (i==n-1) return false;
        } else {
            return true;
        }
    }
    return true;
}

template <typename VT>
bool is_zero_vertex(VT &vertex) {

    for (int i = 0; i < vertex.size(); ++i) {
        if (vertex(i) != 0) return false;
    }
    return true;

}

template <class Polytope>
Polytope everest(int n, int s) {

    typedef typename Polytope::MT    MT;
    typedef typename Polytope::VT    VT;
    typedef Eigen::Matrix<int,Eigen::Dynamic,1> VTint;
    typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> MTint;
    typedef typename Polytope::PolytopePoint Point;
    typedef typename Point::FT NT;

    Polytope P;

    int d = n * s;
    NT nvertices = std::pow(NT(s)+1.0, NT(n) + 1.0) - NT(s) - NT(1);
    MTint V = MTint::Zero(int(nvertices), d);
    VT b = VT::Ones(int(nvertices));

    VTint indices(n), vertex(d), vertex_temp(d);
    indices.setConstant(0);


    int j, i, k, nrow = 0;
    Point v(d), u(d);

    do {
        vertex.setConstant(0);
        for (i = 0; i < n; ++i) {
            if (indices(i) > 0) {
                vertex(i * s + indices(i) - 1) = -NT(1);
            }
        }

        if (!is_zero_vertex(vertex)) {
            V.row(nrow) = vertex;
            nrow++;
        }
        for (j = 1; j < s + 1; ++j) {
            vertex_temp = vertex;
            for (i = 0; i < n; ++i) {
                vertex_temp(i * s + j - 1) += NT(1);
            }
            if (!is_zero_vertex(vertex_temp)) {
                V.row(nrow) = vertex_temp;
                nrow++;
            }
        }
    }
    while (update_everest_indices(n, s, indices));

    MT V2(int(nvertices), d);
    for (i = 0; i < V2.rows(); ++i) {
        for (j = 0; j < d; ++j) {
            V2(i,j) = NT(V(i,j));
        }
    }

    P.init(d, V2, b);
    return P;
}

#endif
