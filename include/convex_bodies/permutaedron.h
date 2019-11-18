// VolEsti (volume computation and sampling library)

// Copyright (c) 2019 Vissarion Fisikopoulos
// Copyright (c) 2019 Apostolos Chalkis


// Licensed under GNU LGPL.3, see LICENCE file

#ifndef PERMUTAEDRON_H
#define PERMUTAEDRON_H

#include <limits>

#include <iostream>
#include "solve_lp.h"
#include "permutaedron_oracles.h"

//min and max values for the Hit and Run functions

// V-Polytope class
template <class Point, class  RNGType>
class Permutaedron{
public:
    typedef Point PolytopePoint;
    typedef typename Point::FT NT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef RNGType rngtype;

private:
    MT T;  //matrix V. Each row contains a vertex
    VT b, beq, interior_point;  // vector b that contains first column of ine file
    MT A, Aeq;
    unsigned int _d;  //dimension
    REAL *conv_comb, *row, *mem_comb;
    int *colno;
    NT maxNT = std::numeric_limits<NT>::max();

public:
    Permutaedron() {}

    // return dimension
    unsigned int dimension() {
        return _d;
    }


    // this function returns 0. The main sampler requests this function to set the length of lambdas vector
    int num_of_hyperplanes() {
        return 0;
    }


    // compute the number of facets of the cyclic polytope in dimension _d with the same number of vertices
    // this is an upper bound for the number of the facets from McMullen's Upper Bound Theorem
    unsigned int upper_bound_of_hyperplanes() {
        return 2*_d;
    }


    // return the number of vertices
    int num_of_vertices() {
        return T.rows();
    }


    // return the matrix V
    MT get_mat() {
        return T;
    }


    // return the vector b
    VT get_vec() {
        return b;
    }


    // change the matrix V
    void set_mat(MT V2) {
        T = V2;
    }


    // change the vector b
    void set_vec(VT b2) {
        b = b2;
    }


    // get a specific coeff of matrix V
    NT get_mat_coeff(unsigned int i, unsigned int j) {
        return T(i,j);
    }


    // get a specific coeff of vector b
    NT get_vec_coeff(unsigned int i) {
        return b(i);
    }


    // set a specific coeff of matrix V
    void put_mat_coeff(unsigned int i, unsigned int j, NT value) {
        T(i,j) = value;
    }


    // set a specific coeff of vector b
    void put_vec_coeff(unsigned int i, NT value) {
        b(i) = value;
    }

    void comp_diam(NT &diam) {
        return;
    }


    void init(unsigned int dim) {
        _d = dim-1;
        unsigned int k = dim*dim;
        MT TT(dim,k);
        T.resize(_d, k);
        TT.setConstant(0);
        Aeq.resize(2*dim, k);
        Aeq.setConstant(0);
        beq = VT::Ones(2*dim);
        A.resize(k,k);
        A.setConstant(0);
        b = VT::Zero(k);
        for (int i = 0; i < k; ++i) {
            A(i,i) = -1.0;
        }
        for (int i = 0; i < 2*dim; ++i) {
            if(i<dim){
                for (int j = 0; j < dim; ++j) {
                    Aeq(i, dim*i+j) = 1.0;
                    TT(i, dim*i+j) = j+1;
                }
            } else {
                for (int j = 0; j < dim; ++j) {
                    Aeq(i,j*dim + i-dim) = 1.0;
                }
            }
        }
        MT Aeq2(1,dim);
        for (int l = 0; l < dim; ++l) {
            Aeq2(0,l) = l+1;
        }
        Eigen::FullPivLU<MT> lu(Aeq2);
        MT N = lu.kernel();
        std::cout<<"N = "<<N<<"\n"<<std::endl;
        T = N.transpose()*TT;
        std::cout<<"T = "<<T<<"\n"<<std::endl;

        VT bbb(1); bbb(0) = (dim*(dim+1))/2;
        VT q = Aeq2.colPivHouseholderQr().solve(bbb);
        std::cout<<"q = "<<q<<"\n"<<std::endl;
        Point center(dim);
        for (int j = 0; j < dim; ++j) {
            center.set_coord(j, j+1);
            std::cout<<center[j]<<std::endl;
        }
        std::cout<<"\n";
        mem_comb = (REAL *) malloc((T.cols()) * sizeof(*mem_comb));
        std::cout<<"get feas point"<<std::endl;
        get_feas_point(TT, A, b, Aeq, beq, center, mem_comb);
        VT qq(k);
        for (int m = 0; m < k; ++m) {
            qq(m) = *(mem_comb+m);
        }
        std::cout<<"qq = "<<qq<<"\n"<<std::endl;
        b = b - A*qq;
        beq = beq - Aeq*qq;
        interior_point.resize(_d);
        for (int i = 0; i < k; ++i) {
            qq(i) = 1.0/(NT(dim)) - qq(i);
        }
        std::cout<<"qq = "<<qq<<"\n"<<std::endl;
        std::cout<<"A*qq = "<<A*qq<<"\n"<<std::endl;
        std::cout<<"b = "<<b<<"\n"<<std::endl;
        b = b - A*qq;
        beq = beq - Aeq*qq;
        interior_point.setConstant(0);

        //std::cout<<"A = "<<A<<"\n"<<std::endl;
        std::cout<<"T = "<<T<<"\n"<<std::endl;
        std::cout<<"b = "<<b<<"\n"<<std::endl;
        conv_comb = (REAL *) malloc((T.cols()+1) * sizeof(*conv_comb));
        colno = (int *) malloc((T.cols()+1) * sizeof(*colno));
        row = (REAL *) malloc((T.cols()+1) * sizeof(*row));
    }


    // Construct matrix V which contains the vertices row-wise
    void init(std::vector<std::vector<NT> > Pin) {
        _d = Pin[0][1] - 1;
        T.resize(Pin.size() - 1, _d);
        b.resize(Pin.size() - 1);
        for (unsigned int i = 1; i < Pin.size(); i++) {
            b(i - 1) = Pin[i][0];
            for (unsigned int j = 1; j < _d + 1; j++) {
                T(i - 1, j - 1) = Pin[i][j];
            }
        }
        conv_comb = (REAL *) malloc(Pin.size() * sizeof(*conv_comb));
        colno = (int *) malloc((T.rows()+1) * sizeof(*colno));
        row = (REAL *) malloc((T.rows()+1) * sizeof(*row));
    }


    // print polytope in input format
    void print() {
#ifdef VOLESTI_DEBUG
        std::cout << " " << V.rows() << " " << _d << " float" << std::endl;
#endif
        for (unsigned int i = 0; i < T.rows(); i++) {
            for (unsigned int j = 0; j < _d; j++) {
#ifdef VOLESTI_DEBUG
                std::cout << V(i, j) << " ";
#endif
            }
#ifdef VOLESTI_DEBUG
            std::cout<<"\n";
#endif
        }
    }



    // pick d+1 random vertices until they define a full dimensional simplex and then
    // compute the chebychev ball of that simplex
    std::pair<Point,NT> ComputeInnerBall() {

        //std::pair<Point, NT> che_up = ComputeChebychevBall<NT, Point>(A, b);
        //std::cout<<"rad = "<<che_up.second<<std::endl;
        //VT p_in_B = VT::Ones((_d+1)*(_d+1))*(1.0/(NT((_d+1))));
        //VT p_in_P = T*p_in_B;
        //std::cout<<p_in_B<<"\n"<<std::endl;
        //std::cout<<A*p_in_B<<"\n"<<std::endl;
        //std::cout<<T*p_in_B<<"\n"<<std::endl;
        Point center(_d);
        for (int j = 0; j < _d; ++j) {
            center.set_coord(j, interior_point(j));
            std::cout<<center[j]<<std::endl;
        }
        std::cout<<"\n";
        //b = b - A*p_in_B;
        std::cout<<"is in center = "<<is_in(center)<<std::endl;
        std::vector<NT> temp(_d,0);
        NT radius =  maxNT, min_plus;
        std::pair<NT,NT> min_max;


        for (unsigned int i = 0; i < _d; ++i) {
            temp.assign(_d,0);
            temp[i] = 1.0;
            Point v(_d,temp.begin(), temp.end());
            min_max = intersect_double_line_permutaedron(T, A, b, Aeq, beq, center, v, conv_comb, row, colno);
            //std::cout<<"pos2 = "<<intersect_line_permutaedron(T, A, b, center, v, conv_comb, row, colno)<<std::endl;

            std::cout<<"pos = "<<min_max.first<<"minus = "<<min_max.second;
            if (radius > min_max.first) radius = min_max.first;
            if (radius > -min_max.second) radius = -min_max.second;
        }

        //throw false;
        radius = radius / std::sqrt(NT(_d));
        return std::pair<Point, NT> (center, radius);
    }


    // check if point p belongs to the convex hull of V-Polytope P
    int is_in(Point p) {
        if(memLP_permutaedron(T, A, b, Aeq, beq, p, mem_comb)){
            return -1;
        }
        return 0;
    }


    // compute intersection point of ray starting from r and pointing to v
    // with the V-polytope
    std::pair<NT,NT> line_intersect(Point r, Point v) {

        return intersect_double_line_permutaedron(T, A, b, Aeq, beq, r, v, conv_comb, row, colno);
    }


    // compute intersection point of ray starting from r and pointing to v
    // with the V-polytope
    std::pair<NT,NT> line_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av) {

        return intersect_double_line_permutaedron(T, A, b, Aeq, beq, r, v, conv_comb, row, colno);
    }

    // compute intersection point of ray starting from r and pointing to v
    // with the V-polytope
    std::pair<NT,NT> line_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av, NT &lambda_prev) {

        return intersect_double_line_permutaedron(T, A, b, Aeq, beq, r, v, conv_comb, row, colno);
    }


    std::pair<NT, int> line_positive_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av) {
        //intersect_line_permutaedron(MT &T, MT &A, VT &b, MT &Aeq, MT &beq, Point &r, Point &v,  NT *conv_comb, NT *row, int *colno)
        //return std::pair<NT, int> (intersect_line_permutaedron(T, A, b, Aeq, beq, r, v, conv_comb, row, colno), 1);
        return std::pair<NT, int> (0.5,1);
    }


    std::pair<NT, int> line_positive_intersect(Point r, Point v, std::vector<NT> &Ar, std::vector<NT> &Av,
                                               NT &lambda_prev) {
        return line_positive_intersect(r, v, Ar, Av);
    }


    // Compute the intersection of a coordinate ray
    // with the V-polytope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          unsigned int rand_coord,
                                          std::vector<NT> &lamdas) {

        std::vector<NT> temp(_d);
        temp[rand_coord]=1.0;
        Point v(_d,temp.begin(), temp.end());
        return intersect_double_line_permutaedron(T, A, b, Aeq, beq, r, v, conv_comb, row, colno);
    }


    // Compute the intersection of a coordinate ray
    // with the V-polytope
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          Point &r_prev,
                                          unsigned int rand_coord,
                                          unsigned int rand_coord_prev,
                                          std::vector<NT> &lamdas) {

        std::vector<NT> temp(_d);
        temp[rand_coord]=1.0;
        Point v(_d,temp.begin(), temp.end());
        return intersect_double_line_permutaedron(T, A, b, Aeq, beq, r, v, conv_comb, row, colno);
    }


    // shift polytope by a point c
    void shift(VT c) {
        return;
    }


    // apply linear transformation, of square matrix T, to the V-Polytope
    void linear_transformIt(MT _T) {
        std::cout<<"_T.inverse() = "<<_T.inverse()<<std::endl;
        T = _T.inverse()*T;
        std::cout<<"\nT = "<<T<<std::endl;
    }

    Point get_mean_of_vertices() {
        return Point(_d);
    }

    NT get_max_vert_norm() {
        return 0.0;
    }

    // consider an upper bound for the number of facets of a V-polytope
    // for each facet consider a lower bound for the distance from the origin
    // useful for CG algorithm to get the first gaussian
    std::vector<NT> get_dists(NT radius) {
        std::vector <NT> res(upper_bound_of_hyperplanes(), radius);
        return res;
    }

    void normalize() {}


    // in number_of_vertices<=20*dimension use the vertices for the rounding
    // otherwise you have to sample from the V-polytope
    template <class PointList>
    bool get_points_for_rounding (PointList &randPoints) {
        return false;
    }

    MT get_T() {
        return T;
    }

    void compute_reflection(Point &v, Point &p, int facet) {

        int count = 0, k =T.cols();
        NT sum, e = 0.0000000001;
        MT Fmat(k - _d +1, k), HypMat(_d, _d);
        VT rand_point(k), beq(k-_d+1);

        for (int i = 0; i < A.rows(); ++i) {
            sum = 0.0;
            for (int j = 0; j < k; ++j) {
                sum += A(i,j)* (*(conv_comb+j));
            }
            if (std::abs(b(i) - sum)  < e*std::abs(b(i)) && std::abs(b(i) - sum)  < e*std::abs(sum)) {
                std::cout<<"a_"<<i<<"x = "<<sum<<" b("<<i<<") = "<<b(i)<<std::endl;
                Fmat.row(count) = A.row(i);
                beq(count) = b(i);
                count++;
            }

        }
        std::cout<<"rows equal to b = "<<count<<std::endl;
        std::cout<<"\n"<<std::endl;

        Eigen::FullPivLU<MT> lu(Fmat);
        MT NS = lu.kernel();

        VT x_pr = Fmat.colPivHouseholderQr().solve(beq);
        VT x0 = T * x_pr;
        MT TT = T * NS;
        HypMat.row(0) = x0.transpose();

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);
        boost::normal_distribution<> rdist(0,1);
        NT normal;
        for (int l = 0; l < _d-1; ++l) {
            normal = 0.0;
            for (unsigned int i = 0; i < k; i++) {
                rand_point(i) = rdist(rng);
                normal += rand_point(i) * rand_point(i);
            }
            normal = 1.0 / std::sqrt(normal);
            rand_point = rand_point * normal;
            HypMat.row(l+1) = (x0 + TT * rand_point).transpose();
            //std::cout<<"Hypmat = "<<HypMat<<"\n(x0 + TT * rand_point).transpose() ="<<(x0 + TT * rand_point).transpose()<<std::endl;
        }

        VT a = HypMat.colPivHouseholderQr().solve(VT::Ones(_d));
        sum = 0.0;
        for (int i = 0; i < _d; ++i) sum += a(i)*p[i];
        if(sum<0.0) a = -1.0*a;

        a = a/a.norm();

        Point s(_d, std::vector<NT>(&a[0], a.data()+a.cols()*a.rows()));
        //Point s(_d);
        //for (int i = 0; i < _d; ++i) s.set_coord(i, a(i));

        s = ((-2.0 * v.dot(s)) * s);
        v = s + v;
    }

    void free_them_all() {
        free(row);
        free(colno);
        free(conv_comb);
    }

};

#endif