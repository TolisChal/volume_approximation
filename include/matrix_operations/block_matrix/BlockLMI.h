// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_BLOCKLMI_H
#define VOLESTI_BLOCKLMI_H

#include "BlockEigenvaluesProblems.h"


/// This class handles a linear matrix inequality of the form \[A_0 +  \sum x_i A_i\]
/// A template specialization for dense Eigen matrices and vectors
/// @tparam NT Numeric Type
template<typename NT, typename MT, typename VT>
class BlockLMI {
    public:
    /// Eigen matrix type
    //typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    /// Eigen vector type
    //typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

    /// The matrices A_0, A_i
    std::vector<MT> matrices;
    MT sum_Ai;

    /// structure to evaluate the lmi fast
    //evaluate_blocklmi<MT> lmi_evaluator;

    /// The dimension of the vector x
    unsigned int d;

    /// The size of the matrices A_i
    unsigned int m;

    /// At each column keep the m*(m+1)/2 distinct elements of each matrix A_i, i=1,...,d
    //MT vectorMatrix;

    BlockLMI(){}

    /// Creates A LMI object
    /// \param[in] matrices The matrices A_0, A_i
    BlockLMI(std::vector<MT>& matrices) 
    {
        d = matrices.size() - 1;
        m = matrices[0].rows();
        typename std::vector<MT>::iterator it = matrices.begin();

        while (it!=matrices.end()) {
            this->matrices.push_back(*it);
            it++;
        }

        
        
        //setVectorMatrix();
        //lmi_evaluator.setVectorMatrix(m, d, matrices);
    }

    void add_matrix(MT const& A) {
        matrices.push_back(A);
        d += 1;
    }


    /// \returns The dimension of vector x
    unsigned int dimension() const {
        return d;
    }

    std::vector< std::pair<int, int> > get_limits_of_blocks(){
        return matrices[0].get_block_limits();
    }

    int get_number_of_blocks() {
        return matrices[0].get_num_of_blocks();
    }

    /// \return The matrices A0, A1, ..., Ad
    std::vector<MT> getMatrices() const {
        return matrices;
    }

    /// \returns The size of the matrices
    unsigned int sizeOfMatrices() const {
        return m;
    }

    /// Evaluate A_0 + \[A_0 + \sum x_i A_i \]
    /// \param[in] x The input vector
    /// \param[out] ret The output matrix
    void evaluate(VT const & x, MT& ret, bool complete_mat = false) 
    {
        //if (!lmi_evaluator.evaluateWithoutA0(x, ret, complete_mat)) {
            //std::cout<<"sparse"<<std::endl;
        typename std::vector<MT>::iterator it = matrices.begin();// = ;
        ret = (*it);
        it++;
        int i = 0;
        for ( ; it!=matrices.end(); it++, i++){
            ret +=  ((*it) * x.coeff(i));
        }
        
    }

    /// Compute  \[x_1*A_1 + ... + x_n A_n]
    /// \param[in] x Input vector
    /// \param[out] res Output matrix
    void evaluateWithoutA0(const VT& x, MT& res, bool complete_mat = false) {
        //if (!lmi_evaluator.evaluateWithoutA0(x, res, complete_mat)) {
            //std::cout<<"sparse"<<std::endl;
        typename std::vector<MT>::iterator it = matrices.begin();// = ;
        it++;
        res = (*it) * x.coeff(0);
        it++;
        int i = 1;
        for ( ; it!=matrices.end(); it++, i++){
            res += (*it) * x.coeff(i);
        }
    }

    /// Compute the gradient of the determinant of the LMI at p
    /// \param[in] p Input parameter
    /// \param[in] Input vector: lmi(p)*e = 0, e != 0
    /// \param[out] ret The normalized gradient of the determinant of the LMI at p
    void normalizedDeterminantGradient(const VT& p, const VT& e, VT& ret) {
        //ret.setZero(d);
        VT temp_vec(d);
        // i-th coordinate of the determinant is e^T * A_i * e
        for (int i = 0; i < d; i++) {
            // todo, use iterators
            matrices[i+1].multiply(e, temp_vec);
            ret(i) = e.transpose() * temp_vec;
        }

        ret.normalize();
    }

    /// \param i An indicator to a matrix
    /// \return Pointer to A_i
    MT* const getMatrix(const int i) {
        return &(matrices[i]);
    }

    /// Prints the matrices A0, ..., An
    void print() const {
        int i = 0;

        for (auto iter = matrices.begin(); iter != matrices.end(); iter++, i++) {
            std::cout << "A" << i << "\n";
            std::cout << *iter << "\n\n";
        }
    }

    /// check if the matrix is negative semidefinite
    /// \param matrix a matrix
    /// \return Pointer to A_i
    bool isNegativeSemidefinite(MT &matrix ) const {

        BlockEigenvaluesProblems<NT, MT, VT> EigenvaluesProblem;
        NT eival = EigenvaluesProblem.findSymEigenvalue(matrix);
        return eival <= 0;
    }

    /// evaluate LMI(pos) and check if its negative semidefinite
    /// \param pos a vector of our current position
    /// \return true is LMI(pos) is negative semidefinite
    bool isNegativeSemidefinite(VT &pos) {
        MT mat(sizeOfMatrices(), sizeOfMatrices());
        evaluate(pos, mat);
        return isNegativeSemidefinite(mat);
    }

};






#endif //VOLESTI_LMI_H

