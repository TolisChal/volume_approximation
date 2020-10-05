// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_DENSEPRODUCTMATRIX_H
#define VOLESTI_DENSEPRODUCTMATRIX_H

#define PARTIAL_LU_DECOMPOSITION
#define CHOLESKY__NO_PIVOTING_DECOMPOSITION

/// A wrapper class for dense Eigen matrices in Spectra and ARPACK++
/// This class will be the wrapper to use the Spectra nonsymemmetric standard eigenvalue Cx = lx solver to
/// solve a generalized eigenvalue Ax = lBx.
/// In particular, this class represents the product @f[ C = B^-1 A @f]
///
/// \tparam NT Numeric Type
template<typename NT>
class DenseProductMatrix {
public:
    /// Eigen matrix type
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
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

    VT v;

    int m;

    /// The decomposition we will use
    /// If PARTIAL_LU_DECOMPOSITION is defined, use the Eigen partial LU decomposition,
    /// otherwise full. The partial is faster but assumes that the matrix has full rank.
#if defined(PARTIAL_LU_DECOMPOSITION)
    typedef Eigen::PartialPivLU<MT> lu_decomposition;    
#else
    typedef Eigen::FullPivLU<MT> lu_decomposition;
#endif

#if defined(CHOLESKY__NO_PIVOTING_DECOMPOSITION)
    //typedef Eigen::PartialPivLU<MT> cholesky_decomposition;
    typedef Eigen::LDLT<MT, Eigen::Upper> cholesky_decomposition;
#else
    typedef Eigen::LDLT<MT> cholesky_decomposition;
#endif

    /// The LU decomposition of B
    //lu_decomposition Blu;
    cholesky_decomposition llt;

    /// Constructs an object of this class and computes the LU decomposition of B.
    ///
    /// \param[in] A The matrix A
    /// \param[in] B The matrix B
    DenseProductMatrix(MT const *A, MT const *B) : A(A), B(B) {
        //Blu.compute(*B);
        
        _rows = A->rows();
        _cols = B->cols();
        m = _cols / 2;
        llt.compute((*B).block(0, m, m, m));

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
        //Eigen::Map<VT> const x(const_cast<double*>(x_in), _rows);
        Eigen::Map<const VT> x(x_in, _cols);
        //VT const v = *A * x;
        int r = _rows / 2;
        v.block(r, 0, r, 1).noalias() = (*A).block(r, r, r, r).template selfadjointView< Eigen::Upper >() * x.block(r, 0, r, 1);
        v.block(0, 0, r, 1).noalias() = (*A).block(0, 0, r, r).template selfadjointView< Eigen::Upper >() * x.block(0, 0, r, 1);
        v = -v;

        Eigen::Map<VT> y(y_out, _rows);
        //y.noalias() = Blu.solve(-v);

        y.block(0, 0, r, 1).noalias() = llt.solve(v.block(r, 0, r, 1));
        v.block(0, 0, r, 1).noalias() -= (*B).block(0, 0, r, r) * y.block(0, 0, r, 1);
        y.block(r, 0, r, 1).noalias() = llt.solve(v.block(0, 0, r, 1));
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
        //Eigen::Map<const VT> x(x_in, _cols);
        //VT const v2 = *A * x;
        int r = _rows / 2;
    
        v.block(r, 0, r, 1).noalias() = (*A).block(r, r, r, r).template selfadjointView< Eigen::Upper >() * x.block(r, 0, r, 1);
        v.block(0, 0, r, 1).noalias() = (*A).block(0, 0, r, r).template selfadjointView< Eigen::Upper >() * x.block(0, 0, r, 1);
        v = -v;

        Eigen::Map<VT> y(y_out, _rows);
        //y.noalias() = Blu.solve(v);

        //VT yy(_rows);
        //std::cout<<"\n B = "<<(*B)<<"\n"<<std::endl;
        //std::cout<<"B.block(0, m, m, m) = "<<(*B).block(0, r, r, r)<<"\n"<<std::endl;
        //std::cout<<"v = "<<v.transpose()<<std::endl;
        y.block(0, 0, r, 1).noalias() = llt.solve(v.block(r, 0, r, 1));
        v.block(0, 0, r, 1).noalias() -= (*B).block(0, 0, r, r) * y.block(0, 0, r, 1);
        y.block(r, 0, r, 1).noalias() = llt.solve(v.block(0, 0, r, 1));
        

        //std::cout<<"y = "<<y.transpose()<<std::endl;
        //std::cout<<"yy = "<<yy.transpose()<<"\n"<<std::endl;

        //y=yy;
    }
};
#endif //VOLESTI_DENSEPRODUCTMATRIX_H
