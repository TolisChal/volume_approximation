// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_BLOCKMATRIXQUADEIG_H
#define VOLESTI_BLOCKMATRIXQUADEIG_H

/// A wrapper class for dense Eigen matrices in Spectra and ARPACK++
/// This class will be the wrapper to use the Spectra nonsymemmetric standard eigenvalue Cx = lx solver to
/// solve a generalized eigenvalue Ax = lBx.
/// In particular, this class represents the product @f[ C = B^-1 A @f]
///
/// \tparam NT Numeric Type
template<typename NT>
class BlockMatrixQuadEig {
public:
    /// Eigen matrix type
    typedef SparseBlock<NT> MT;
    /// Eigen vector type
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    /// The number of rows
    int _rows;
    /// The number of cols
    int _cols;

    /// Pointer to matrix A
    MT const *A;
    /// Pointer to matrix B
    MT const *B;
    /// Pointer to matrix C
    MT const *C;

    VT v;

    int m;

    /// The LU decomposition of B
    //lu_decomposition Blu;
    decompositioner<NT> takis;

    /// Constructs an object of this class and computes the LU decomposition of B.
    ///
    /// \param[in] A The matrix A
    /// \param[in] B The matrix B
    BlockMatrixQuadEig(MT const *A, MT const *B, MT const *C) : A(A), B(B), C(C) {
        //Blu.compute(*B);
        
        _rows = 2*(A->rows());
        _cols = 2*(B->cols());
        //_rows *= 2;
        //_cols *= 2;

        //std::cout<<"_rows = "<<_rows<<std::endl;
        //std::cout<<"_cols = "<<_cols<<std::endl;

        m = _cols / 2;
        //llt.compute(*C);
        takis.compute_dense_cholesky(*C);

        v.setZero(_rows);
    }

    ///Required by Spectra
    /// \return The number of rows
    int rows() {
        return _rows;
    }

    ///Required by Spectra
    /// \return The number of columns
    int cols() {
        return _cols;
    }

    /// Required by Spectra.
    /// Computes the product Cx = y, i.e. @f[ (B^-1 A)v = y@$]. But B = LU, so Ax = LUy.
    /// Let Ax = v, then LUy = v. Then Lw = v and finally Uy = w to get y;
    /// \param[in] x_in
    /// \param[out] y_out
    void perform_op(NT const * x_in, NT* y_out) {

        // Declaring the vectors like this, we don't copy the values of x_in to v
        // and next of y to y_out
        Eigen::Map<const VT> x(x_in, _cols);
        
        v.block(m, 0, m, 1) = x.block(m, 0, m, 1);

        (*A).multiply(x.block(0, 0, m, 1), v.block(0, 0, m, 1));
        //v.block(0, 0, m, 1).noalias() = (*A).template selfadjointView< Eigen::Lower >() * x.block(0, 0, m, 1);

        Eigen::Map<VT> y(y_out, _rows);
        //y.noalias() = Blu.solve(-v);

        y.block(0, 0, m, 1) = -v.block(m, 0, m, 1);

        (*B).multiply(y.block(0, 0, m, 1), v.block(0, 0, m, 1));
        //v.block(0, 0, m, 1).noalias() += (*B).template selfadjointView< Eigen::Lower >() * y.block(0, 0, m, 1);

        takis.solve_chol_dense_ls(v.block(0, 0, m, 1), y.block(m, 0, m, 1));
        //y.block(m, 0, m, 1).noalias() = llt.solve(v.block(0, 0, m, 1));
    }

    /// Required by arpack.
    /// Computes the product Cx = y, i.e. @f[ (B^-1 A)v = y@$]. But B = LU, so Ax = LUy.
    /// Let Ax = v, then LUy = v. Then Lw = v and finally Uy = w to get y;
    /// \param[in] x_in
    /// \param[out] y_out
    void MultMv(NT * x_in, NT* y_out) {

        // Declaring the vectors like this, we don't copy the values of x_in to v
        // and next of y to y_out
        Eigen::Map<const VT> x(const_cast<double*>(x_in), _rows);
        v.block(m, 0, m, 1) = x.block(m, 0, m, 1);

        (*A).multiply(x.block(0, 0, m, 1), v.block(0, 0, m, 1));
        //v.block(0, 0, m, 1).noalias() = (*A).template selfadjointView< Eigen::Lower >() * x.block(0, 0, m, 1);

        Eigen::Map<VT> y(y_out, _rows);
        //y.noalias() = Blu.solve(-v);

        y.block(0, 0, m, 1) = -v.block(m, 0, m, 1);

        (*B).multiply(y.block(0, 0, m, 1), v.block(0, 0, m, 1));
        //v.block(0, 0, m, 1).noalias() += (*B).template selfadjointView< Eigen::Lower >() * y.block(0, 0, m, 1);

        takis.solve_chol_dense_ls(v.block(0, 0, m, 1), y.block(m, 0, m, 1));
        //y.block(m, 0, m, 1).noalias() = llt.solve(v.block(0, 0, m, 1));
       
    }
};
#endif //VOLESTI_DENSEPRODUCTMATRIX_H
