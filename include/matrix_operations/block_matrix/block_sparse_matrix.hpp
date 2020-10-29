// a new data structure

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
    typedef Eigen::PartialPivLU<MT> dense_partial_lu_decomposition; 
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
            dec.compute(MT(A.get_block(i)));
            dense_chols.push_back(dec);
        }
    }

    void solve_chol_dense_ls(VT const& b, VT &res)
    {
        typename std::vector< std::pair<int, int> >::iterator iter_limits = block_limits.begin();
        typename std::vector< std::pair<int, int> >::iterator iter_dec = dense_chols.begin();
        int length;

        for (; iter_limits != blocks.end(); iter_limits++, iter_dec++) 
        {
            length = (*iter_limits).second - (*iter_limits).first + 1;
            res.block((*iter_limits).first, 0, length, 1).noalias() = 
                             (*iter_dec).solve(b.block((*iter_limits).first, 0, length, 1));
        }
    }

    void solve_chol_sparse_ls(VT const& b, VT &res)
    {
        typename std::vector< std::pair<int, int> >::iterator iter_limits = block_limits.begin();
        typename std::vector< std::pair<int, int> >::iterator iter_dec = sparse_chols.begin();
        int length;

        for (; iter_limits != blocks.end(); iter_limits++, iter_dec++) 
        {
            length = (*iter_limits).second - (*iter_limits).first + 1;
            res.block((*iter_limits).first, 0, length, 1).noalias() = 
                             (*iter_dec).solve(b.block((*iter_limits).first, 0, length, 1));
        }
    }

};


template <typename NT>
struct SparseBlock {
public:

    /// The type for Sparse Eigen Matrix
    typedef Eigen::SparseMatrix<NT> MT;
    /// The type for Eigen vector
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    /// The type for Dense Eigen Matrix
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> DMT;

    decompositioner<NT> takis;

    std::vector<MT> blocks;
    std::vector< std::pair<int, int> > block_limits;
    int num_of_blocks;

    SparseBlock(std::vector<MT> const& _blocks) : blocks(_blocks)
    {
        typename std::vector<MT>::iterator iter = _blocks.begin();
        int r_start = 0, r_end;
        for (; iter != _blocks.end(); iter++) 
        {
            r_end = (*iter).rows() - 1 + r_start;
            block_limits.push_back(std::pair<int, int>(r_start, r_end));
            r_start = r_end + 1;
            num_of_blocks++;
        }
    }

    void multiply(VT const& vec, VT &res) const
    {
        typename std::vector<MT>::iterator iter_mat = blocks.begin();
        typename std::vector< std::pair<int, int> >::iterator iter_limits = block_limits.begin();
        int length;

        for (;  iter != blocks.end(); iter_mat++, iter_limits++) 
        {
            length = (*iter_limits).second - (*iter_limits).first + 1;
            res.block((*iter_limits).first, 0, length, 1).noalias() = 
                             (*iter_mat) * vec.block((*iter_limits).first, 0, length, 1);
        }
    }

    int get_num_of_blocks() const {
        return num_of_blocks;
    }

    std::vector< std::pair<int, int> > get_block_limits() const {
        return block_limits;
    }

    MT get_block(unsigned int const& i) const
    {
        if (i > (blocks.size()-1)) {
            std::cout<<"too large index"<<std::endl;
            exit(-1);
        }
        return blocks[i]; 
    }

    void operator+= (SparseBlock const& A)
    {
        typename std::vector<MT>::iterator iter_mat = blocks.begin();
        //typename std::vector< std::pair<int, int> >::iterator iter_limits = block_limits.begin();
        unsigned int counter = 0;

        for (;  iter != blocks.end(); iter_mat++, counter++) 
        {
            (*iter_mat) += A.get_block(counter);
        }
    }



};