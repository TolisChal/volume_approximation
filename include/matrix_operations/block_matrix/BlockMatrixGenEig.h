// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_BLOCKMATRIXGENEIG_H
#define VOLESTI_BLOCKMATRIXGENEIG_H


/// A wrapper class for dense Eigen matrices in Spectra and ARPACK++
/// This class will be the wrapper to use the Spectra nonsymemmetric standard eigenvalue Cx = lx solver to
/// solve a generalized eigenvalue Ax = lBx.
/// In particular, this class represents the product @f[ C = B^-1 A @f]
///
/// \tparam NT Numeric Type
template<typename NT>
class BlockMatrixGenEig {
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
    MT *A;
    /// Pointer to matrix B
    MT *B;

    VT v, r;

    /// The LU decomposition of B
    decompositioner<MT> takis;

    /// Constructs an object of this class and computes the LU decomposition of B.
    ///
    /// \param[in] A The matrix A
    /// \param[in] B The matrix B
    BlockMatrixGenEig(MT *_A, MT *_B) : A(_A), B(_B) {
        //Blu.compute(*B);
        takis.compute_dense_cholesky((*B)*(-1.0));
        _rows = A->rows();
        _cols = B->cols();
        v.setZero(_rows);
        r.setZero(_rows);
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
        Eigen::Map<VT> const x(const_cast<double*>(x_in), _rows);
        
        (*A).multiply(x, v);

        Eigen::Map<VT> y(y_out, _rows);
        
        takis.solve_chol_dense_ls(v, r);
        y = -r;
        //y = Blu.solve(v);
    }

    /// Required by arpack.
    /// Computes the product Cx = y, i.e. @f[ (B^-1 A)v = y@$]. But B = LU, so Ax = LUy.
    /// Let Ax = v, then LUy = v. Then Lw = v and finally Uy = w to get y;
    /// \param[in] x_in
    /// \param[out] y_out
    void MultMv(NT * x_in, NT* y_out) {

        // Declaring the vectors like this, we don't copy the values of x_in to v
        // and next of y to y_out
        Eigen::Map<VT> const x(const_cast<double*>(x_in), _rows);
        (*A).multiply(x, v);

        Eigen::Map<VT> y(y_out, _rows);
        takis.solve_chol_dense_ls(v, r);
        y = -r;
        //y = Blu.solve(v);
    }
};
#endif //VOLESTI_SYMGENPRODMATRIX_H

