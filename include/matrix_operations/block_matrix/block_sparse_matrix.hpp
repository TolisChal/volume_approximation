// a new data structure

#ifndef VOLESTI_BLOCKSPARSEMATRIX_H
#define VOLESTI_BLOCKSPARSEMATRIX_H





template <typename NT>
struct SparseBlock {
public:

    /// The type for Sparse Eigen Matrix
    typedef Eigen::SparseMatrix<NT> MT;
    /// The type for Eigen vector
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    /// The type for Dense Eigen Matrix
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> DMT;

    std::vector<MT> blocks;
    std::vector< std::pair<int, int> > block_limits;
    int num_of_blocks, dimension;

    SparseBlock() {}

    SparseBlock(int const& _m1, int const& _m2) {}

    SparseBlock(std::vector<MT> &_blocks) : blocks(_blocks)
    {
        typename std::vector<MT>::iterator iter = _blocks.begin();
        int r_start = 0, r_end;
        dimension = 0;
        for (; iter != _blocks.end(); iter++) 
        {
            dimension += (*iter).rows();
            r_end = (*iter).rows() - 1 + r_start;
            block_limits.push_back(std::pair<int, int>(r_start, r_end));
            r_start = r_end + 1;
            num_of_blocks++;
        }
    }

    void multiply(VT const& vec, VT &res)
    {
        typename std::vector<MT>::iterator iter_mat = blocks.begin();
        typename std::vector< std::pair<int, int> >::iterator iter_limits = block_limits.begin();
        int length;

        for (;  iter_mat != blocks.end(); iter_mat++, iter_limits++) 
        {
            length = (*iter_limits).second - (*iter_limits).first + 1;
            res.block((*iter_limits).first, 0, length, 1).noalias() = 
                             (*iter_mat).template selfadjointView< Eigen::Lower >() * vec.block((*iter_limits).first, 0, length, 1);
        }
    }

    int get_num_of_blocks() const
    {
        return num_of_blocks;
    }

    int get_dimension() const {
        return dimension;
    }

    int cols() const {
        return dimension;
    }

    int rows() const {
        return dimension;
    }

    std::vector< std::pair<int, int> > get_block_limits() const 
    {
        return block_limits;
    }

    std::vector<MT> get_all_blocks() const
    {
        return blocks;
    }

    MT get_block(unsigned int const& i) const
    {
        if (i > (blocks.size()-1)) {
            std::cout<<"too large index"<<std::endl;
            exit(-1);
        }
        return blocks[i]; 
    }

    void set_all_blocks(std::vector<MT> const& _blocks) {
        blocks = _blocks;
    }

    void set_block_limits(std::vector< std::pair<int, int> > const& _block_limits) {
        block_limits = _block_limits;
    }

    void set_num_of_blocks(int const& _num_of_blocks) {
        num_of_blocks = _num_of_blocks;
    }

    void set_dimension (int const& _d) {
        dimension = _d;
    }

    void operator= (SparseBlock const& A)
    {
        set_all_blocks(A.get_all_blocks());
        set_block_limits( A.get_block_limits());
        set_num_of_blocks( A.get_num_of_blocks());
        set_dimension( A.get_dimension());
        //this->blocks = A.get_all_blocks();
        //this->block_limits = A.get_block_limits();
        //this->num_of_blocks = A.get_num_of_blocks();
    }

    void operator*= (NT const& k)
    {
        typename std::vector<MT>::iterator iter_mat = blocks.begin();
        for (;  iter_mat != blocks.end(); iter_mat++) 
        {
            (*iter_mat) *= k;
        }
    }

    SparseBlock operator* (NT const& k)
    {
        SparseBlock B = (*this);
        B *= k;
        return B;
    }

    VT operator* (VT const& x)
    {
        VT res(x.rows());
        multiply(x, res);
        return res;
    }

    void operator+= (SparseBlock const& A)
    {
        typename std::vector<MT>::iterator iter_mat = blocks.begin();
        unsigned int counter = 0;

        for (;  iter_mat != blocks.end(); iter_mat++, counter++) 
        {
            (*iter_mat) += A.get_block(counter);
        }
    }

    SparseBlock operator+ (SparseBlock const& A)
    {
        SparseBlock B = A;
        B += (*this);
        return B;
    }

};

//template<typename NT>
//SparseBlock<NT> operator* (NT const& k, SparseBlock<NT> const A)
//{
//    SparseBlock<NT> B = A * k;
//    return B;/
//}


template <typename MT>
struct decompositioner {
    
};

template <typename NT>
struct decompositioner< SparseBlock<NT> > {
public:

    /// The type of sparse block matrix
    typedef SparseBlock<NT> SBM;
    /// The type for Sparse Eigen Matrix
    typedef Eigen::SparseMatrix<NT> SMT;
    /// The type for Eigen vector
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    /// The type for Dense Eigen Matrix
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> DMT;

    typedef Eigen::FullPivLU<DMT> dense_full_lu_decomposition;
    typedef Eigen::PartialPivLU<DMT> dense_partial_lu_decomposition; 
    typedef Eigen::LDLT<DMT, Eigen::Lower> dense_chol_decomposition;

    typedef Eigen::SparseLU<SMT, Eigen::COLAMDOrdering<int> > sparse_lu_decomposition;
    typedef Eigen::SimplicialLDLT<SMT, Eigen::Lower> sparse_chol_decomposition;
    typedef Eigen::ConjugateGradient<SMT, Eigen::Lower> sparse_symmetric_iterative_solver;

    std::vector<sparse_chol_decomposition> sparse_chols;
    std::vector<dense_chol_decomposition> dense_chols;

    std::vector< std::pair<int, int> > block_limits;

    void compute_sparse_cholesky (SBM const& A) 
    {
        block_limits = A.get_block_limits();
        int n_mat = A.get_num_of_blocks();
        for (int i=0; i < n_mat; i++) {
            sparse_chol_decomposition dec;
            dec.analyzePattern(A.get_block(i));
            dec.factorize(A.get_block(i));
            sparse_chols.push_back(dec);
        }
    }

    void compute_dense_cholesky (SBM const& A) 
    {
        block_limits = A.get_block_limits();
        int n_mat = A.get_num_of_blocks();
        for (int i=0; i < n_mat; i++) {
            dense_chol_decomposition dec;
            dec.compute(DMT(A.get_block(i)));
            dense_chols.push_back(dec);
        }
    }

    void solve_chol_dense_ls(VT const& b, VT &res)
    {
        typename std::vector< std::pair<int, int> >::iterator iter_limits = block_limits.begin();
        typename std::vector< dense_chol_decomposition >::iterator iter_dec = dense_chols.begin();
        int length;

        for (; iter_limits != block_limits.end(); iter_limits++, iter_dec++) 
        {
            length = (*iter_limits).second - (*iter_limits).first + 1;
            res.block((*iter_limits).first, 0, length, 1).noalias() = 
                             (*iter_dec).solve(b.block((*iter_limits).first, 0, length, 1));
        }
    }

    void solve_chol_sparse_ls(VT const& b, VT &res)
    {
        typename std::vector< std::pair<int, int> >::iterator iter_limits = block_limits.begin();
        typename std::vector< sparse_chol_decomposition >::iterator iter_dec = sparse_chols.begin();
        int length;

        for (; iter_limits != block_limits.end(); iter_limits++, iter_dec++) 
        {
            length = (*iter_limits).second - (*iter_limits).first + 1;
            res.block((*iter_limits).first, 0, length, 1).noalias() = 
                             (*iter_dec).solve(b.block((*iter_limits).first, 0, length, 1));
        }
    }

};

#endif
