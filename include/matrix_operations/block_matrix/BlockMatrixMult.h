// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_BLOCKMATRIXMULT_H
#define VOLESTI_BLOCKMATRIXMULT_H

/// A wrap class to use Eigen dense matrices when solving Eigenvalue problems with ARPACK++
/// \tparam NT Numeric Type
template<class NT>
class BlockMatrixMult {
public:

    /// Eigen matrix type
    typedef SparseBlock<NT> MT;
    /// Eigen vector type
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

    /// The matrix
    MT *M;

    VT r;

    /// number of columns
    int n;
    /// number of rows
    int m;

    /// \return Number of rows
    int nrows() { return m;}

    /// \return Number of columns
    int ncols() { return n;}

    /// \return Number of rows
    int rows() { return m;}

    /// \return Number of columns
    int cols() { return n;}

    /// Required by ARPACK++ : Multiplies the matrix with vector v
    /// \param[in] v The input vector, for example double*
    /// \param[out] w The result of M*v
    void MultMv(NT const *v, NT* w) {
        // Declaring the vectors like this, we don't copy the values of v and after to w
        Eigen::Map<const VT> _v(v, m);
        Eigen::Map<VT> _w(w, m);

        (*M).multiply(_v, r);
        _w = r;

        //_w.noalias() = (*M).template selfadjointView< Eigen::Lower >() * _v;
    }

    /// Required by ARPACK++ : Multiplies the matrix with vector v
    /// \param[in] v The input vector, for example double*
    /// \param[out] w The result of M*v
    void perform_op(NT const *v, NT* w) {
        // Declaring the vectors like this, we don't copy the values of v and after to w
        Eigen::Map<const VT> _v(v, m);
        Eigen::Map<VT> _w(w, m);

        (*M).multiply(_v, r);
        _w = r;

        //_w.noalias() = (*M).template selfadjointView< Eigen::Lower >() * _v;
    }


    /// Constructs an object
    /// \param[in] M An Eigen Matrix
    BlockMatrixMult(MT *M) {
        this->M = M;
        n = M->cols();
        m = M->rows();
        r.setZero(m);
    }

};
#endif //VOLESTI_EIGENDENSEMATRIX_H
