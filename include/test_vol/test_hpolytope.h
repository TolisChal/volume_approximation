// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef TEST_HPOLYTOPE_H
#define TEST_HPOLYTOPE_H

#include <limits>

#include <iostream>
#include "test_solve_lp.h"

//min and max values for the Hit and Run functions


// H-polytope class
template <typename Point>
class HPolytope{
public:
    typedef Point PolytopePoint;
    typedef typename Point::FT NT;
    typedef typename std::vector<NT>::iterator viterator;
    using RowMatrixXd = Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    typedef RowMatrixXd MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MTc;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;

private:
    MT A; //matrix A
    MTc AA;
    VT b; // vector b, s.t.: Ax<=b
    unsigned int            _d; //dimension
    //NT maxNT = 1.79769e+308;
    //NT minNT = -1.79769e+308;
    NT maxNT = std::numeric_limits<NT>::max();
    NT minNT = std::numeric_limits<NT>::lowest();
    NT inner_ball_norm;
    bool ball_hit=false;
    int facet_k;

public:
    HPolytope() {}

    // constructor: cube(d)
    HPolytope(unsigned int d): _d(d) {
        A.resize(2 * d, d);
        b.resize(2 * d);
        for (unsigned int i = 0; i < d; ++i) {
            b(i) = 1;
            for (unsigned int j = 0; j < d; ++j) {
                if (i == j) {
                    A(i, j) = 1;
                } else {
                    A(i, j) = 0;
                }
            }
        }
        for (unsigned int i = 0; i < d; ++i) {
            b(i + d) = 1;
            for (unsigned int j = 0; j < d; ++j) {
                if (i == j) {
                    A(i + d, j) = -1;
                } else {
                    A(i + d, j) = 0;
                }
            }
        }
    }


    // return dimension
    unsigned int dimension() const {
        return _d;
    }


    // return the number of facets
    int num_of_hyperplanes() const {
        return A.rows();
    }


    // return the matrix A
    MT get_mat() const {
        return A;
    }


    // return the vector b
    VT get_vec() const {
        return b;
    }


    // change the matrix A
    void set_mat(const MT &A2) {
        A = A2;
    }


    // change the vector b
    void set_vec(const VT &b2) {
        b = b2;
    }


    // set a specific coeff of matrix A
    NT get_mat_coeff(const unsigned int &i, const unsigned int &j) const {
        return A(i,j);
    }


    // get a spesific coeff of vector b
    NT get_vec_coeff(const unsigned int &i) const {
        return b(i);
    }


    // get a specific coeff of matrix A
    void put_mat_coeff(const unsigned int &i, const unsigned int &j, const NT &value) {
        A(i,j) = value;
    }


    // set a spesific coeff of vector b
    void put_vec_coeff(const unsigned int &i, const NT &value) {
        b(i) = value;
    }

    bool need_ref() const {
        return false;
    }

    Point get_mean_of_vertices() const {
        return Point(_d);
    }


    NT get_max_vert_norm() const {
        return 0.0;
    }

    void hit_ball(bool hb) {
        ball_hit = hb;
    }

    void inner_normal_norm(NT in_b) {
        inner_ball_norm = in_b;
    }
    void recompute_AA() {
        AA.noalias() = A * A.transpose();
    }

    void comp_diam(NT &diam, const NT &cheb_rad) {
        if(cheb_rad < 0.0) {
            diam = 4.0 * std::sqrt(NT(_d)) * ComputeInnerBall().second;
        } else {
            diam = 4.0 * std::sqrt(NT(_d)) * cheb_rad;
        }
    }

    void init(const unsigned int dim, const MT &_A, const VT &_b) {
        _d = dim;
        A = _A;
        AA.noalias() = A * A.transpose();
        b = _b;
        //std::cout<<AA<<std::endl;
    }

    //define matrix A and vector b, s.t. Ax<=b and the dimension
    void init(const std::vector<std::vector<NT> > &Pin) {
        _d = Pin[0][1] - 1;
        A.resize(Pin.size() - 1, _d);
        b.resize(Pin.size() - 1);
        for (unsigned int i = 1; i < Pin.size(); i++) {
            b(i - 1) = Pin[i][0];
            for (unsigned int j = 1; j < _d + 1; j++) {
                A(i - 1, j - 1) = -Pin[i][j];
            }
        }
        AA.noalias() = A * A.transpose();
    }


    // print polytope in input format
    void print() {
#ifdef VOLESTI_DEBUG
        std::cout << " " << A.rows() << " " << _d + 1 << " float" << std::endl;
#endif
        for (unsigned int i = 0; i < A.rows(); i++) {
            for (unsigned int j = 0; j < _d; j++) {
                #ifdef VOLESTI_DEBUG
                std::cout << -A(i, j) << " ";
                #endif
            }
            #ifdef VOLESTI_DEBUG
            std::cout << "<= " << b(i) << std::endl;
            #endif
        }
    }

    
    //Check if Point p is in H-polytope P:= Ax<=b
    int is_in(const Point &p) const {
        //NT sum;
        int m = A.rows();
        const NT* b_data = b.data();

        for (int i = 0; i < m; i++) {
            if (*b_data - A.row(i) * p.getCoefficients() < NT(0)) return 0;
            //sum = *b_data - A.row(i) * p.getCoefficients();
            b_data++;
            //Check if corresponding hyperplane is violated

        }
        return -1;
    }


    //Compute Chebyshev ball of H-polytope P:= Ax<=b
    //Use LpSolve library
    std::pair<Point,NT> ComputeInnerBall() {

        //lpSolve lib for the linear program
        return test_ComputeChebychevBall<NT, Point>(A, b);
    }


    // compute intersection point of ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT,NT> line_intersect(Point &r, Point &v) {

        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
        VT sum_nom, sum_denom;
        //unsigned int i, j;
        unsigned int j;
        int m = num_of_hyperplanes();

        sum_nom.noalias() = b - A * r.getCoefficients();
        sum_denom.noalias() = A * v.getCoefficients();

        NT* sum_nom_data = sum_nom.data();
        NT* sum_denom_data = sum_denom.data();

        for (int i = 0; i < m; i++) {

            if (*sum_denom_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *sum_nom_data / *sum_denom_data;
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }

            sum_nom_data++;
            sum_denom_data++;
        }
        return std::pair<NT, NT>(min_plus, max_minus);
    }


    // compute intersection points of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT,NT> line_intersect(Point &r, Point &v, VT &Ar,
            VT &Av, NT &inner_vi_ak, bool pos = false) {

        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
        VT sum_nom;
        //unsigned int i, j;
        unsigned int j;
        int m = num_of_hyperplanes(), facet;
        //viterator rit, vit, Ariter = Ar.begin(), Aviter = Av.begin();

        Ar.noalias() = A * r.getCoefficients();
        sum_nom.noalias() = b - Ar;
        Av.noalias() = A * v.getCoefficients();;

        NT* Av_data = Av.data();
        NT* sum_nom_data = sum_nom.data();

        for (int i = 0; i < m; i++) {
            if (*Av_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *sum_nom_data / *Av_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    if (pos){
                        facet = i;
                        facet_k = i;
                        inner_vi_ak = *Av_data;
                    }
                }else if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }

            Av_data++;
            sum_nom_data++;
        }
        if (pos) return std::pair<NT, NT>(min_plus, facet);
        return std::pair<NT, NT>(min_plus, max_minus);
    }

    std::pair<NT,NT> line_intersect(Point &r, Point &v, VT &Ar,
            VT &Av, const NT &lambda_prev, NT &inner_vi_ak, bool pos = false) {

        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
        VT sum_nom;
        //unsigned int i, j;
        unsigned int j;
        int m = num_of_hyperplanes(), facet;

        Ar.noalias() += lambda_prev*Av;
        sum_nom.noalias() = b - Ar;
        Av.noalias() = A * v.getCoefficients();

        NT* sum_nom_data = sum_nom.data();
        NT* Av_data = Av.data();

        for (int i = 0; i < m; i++) {
            if (*Av_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *sum_nom_data / *Av_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    if (pos){
                        facet = i;
                        facet_k = i;
                        inner_vi_ak = *Av_data;
                    }
                }else if (lamda > max_minus && lamda < 0) max_minus = lamda;
            }
            Av_data++;
            sum_nom_data++;
        }
        if (pos) return std::pair<NT, NT>(min_plus, facet);
        return std::pair<NT, NT>(min_plus, max_minus);

    }

    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(Point &r, Point &v, VT &Ar, VT &Av, NT &inner_vi_ak) {
        return line_intersect(r, v, Ar, Av, inner_vi_ak, true);
    }


    // compute intersection point of a ray starting from r and pointing to v
    // with polytope discribed by A and b
    std::pair<NT, int> line_positive_intersect(Point &r, Point &v, VT &Ar, VT &Av,
                                               NT &lambda_prev, NT &inner_vi_ak, bool new_v = false) {
        if (new_v) {
            return line_intersect(r, v, Ar, Av, lambda_prev, inner_vi_ak, true);
        }

        NT lamda = 0, min_plus = NT(maxNT), sum2;
        VT sum_nom;//
        // , sum_denom, sum2;
        //unsigned int i, j;
        unsigned int j;
        int m = num_of_hyperplanes(), facet_prev = facet_k;
        NT inner_prev = inner_vi_ak;
        //viterator vit, Ariter = Ar.begin(), Aviter = Av.begin();

        //std::cout<<"[3]facet_k = "<<facet_k<<", inner_vi_ak = "<<inner_prev<<std::endl;
        Ar.noalias() += lambda_prev*Av;
        if(!ball_hit) {
            Av.noalias() += (-2.0 * inner_prev) * AA.col(facet_prev);
        } else {
            Av.noalias() += (-2.0 * inner_prev) * (Ar / inner_ball_norm);
            //sum2 = (-2.0 * inner_prev) * ((*Ariter)/inner_ball_norm);
        }
        sum_nom.noalias() = b - Ar;

        NT* sum_nom_data = sum_nom.data();
        NT* Av_data = Av.data();

        for (int i = 0; i < m; i++) {
            if (*Av_data == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *sum_nom_data / *Av_data;
                if (lamda < min_plus && lamda > 0) {
                    min_plus = lamda;
                    //facet = i;
                    facet_k = i;
                    inner_vi_ak = *Av_data;
                }
            }
            Av_data++;
            sum_nom_data++;
        }
        return std::pair<NT, NT>(min_plus, facet_k);
    }


    //First coordinate ray intersecting convex polytope
    std::pair<NT,NT> line_intersect_coord(Point &r, const unsigned int &rand_coord,
                                          VT& lamdas) {

        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);
        VT sum_denom;
        unsigned int j;
        int m = num_of_hyperplanes();

        sum_denom = A.col(rand_coord);
        lamdas = b - A * r.getCoefficients();

        NT* lamda_data = lamdas.data();
        NT* sum_denom_data = sum_denom.data();

        for (int i = 0; i < m; i++) {

            if (*sum_denom_data == NT(0)) {
                //std::cout<<"div0"<<sum_denom<<std::endl;
                ;
            } else {
                lamda = *lamda_data * (1 / *sum_denom_data);
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;

            }
            lamda_data++;
            sum_denom_data++;
        }
        return std::pair<NT, NT>(min_plus, max_minus);
    }


    //Not the first coordinate ray intersecting convex
    std::pair<NT,NT> line_intersect_coord(Point &r,
                                          const Point &r_prev,
                                          const unsigned int rand_coord,
                                          const unsigned int rand_coord_prev,
                                          VT& lamdas) {
        NT lamda = 0, min_plus = NT(maxNT), max_minus = NT(minNT);

        int m = num_of_hyperplanes();

        lamdas += A.col(rand_coord_prev)* (r_prev[rand_coord_prev] - r[rand_coord_prev]);
        NT* data = lamdas.data();

        for (int i = 0; i < m; i++) {
            NT a = A(i, rand_coord);

            if (a == NT(0)) {
                //std::cout<<"div0"<<std::endl;
                ;
            } else {
                lamda = *data / a;
                if (lamda < min_plus && lamda > 0) min_plus = lamda;
                if (lamda > max_minus && lamda < 0) max_minus = lamda;

            }
            data++;
        }
        return std::pair<NT, NT>(min_plus, max_minus);
    }


    // Apply linear transformation, of square matrix T^{-1}, in H-polytope P:= Ax<=b
    void linear_transformIt(const MT &T) {
        A = A * T;
    }


    // shift polytope by a point c
    void shift(const VT &c){
        b = b - A*c;
    }


    // return for each facet the distance from the origin
    std::vector<NT> get_dists(const NT &radius){
        unsigned int i=0;
        std::vector <NT> dists(num_of_hyperplanes(), NT(0));
        typename std::vector<NT>::iterator disit = dists.begin();
        for ( ; disit!=dists.end(); disit++, i++)
            *disit = b(i) / A.row(i).norm();

        return dists;
    }


    // no points given for the rounding, you have to sample from the polytope
    template <typename T>
    bool get_points_for_rounding (const T &randPoints) {
        return false;
    }

    MT get_T() const {
        return A;
    }

    void normalize() {

        NT row_norm;
        for (int i = 0; i < num_of_hyperplanes(); ++i) {
            row_norm = A.row(i).norm();
            A.row(i) = A.row(i) / row_norm;
            b(i) = b(i) / row_norm;
        }

    }

    //void compute_reflection(Point &v, const Point &p, const int facet) {

        //VT a = (-2.0 * inner_vi_ak) * A.row(facet_k);
        //for (int i = 0; i < _d; ++i) v.set_coord(i, v[i] + a(i));

        //VT a = A.row(facet);
        //Point s(_d, std::vector<NT>(&a[0], a.data()+a.cols()*a.rows()));
        //s = ((-2.0 * v.dot(s)) * s);
        //v = s + v;

    //}
    
    void compute_reflection(Point &v, const Point &p, const NT &inner_vi_ak, const int &facet) {
      
      Point a((-2.0 * inner_vi_ak) * A.row(facet_k));
      v += a;
      //for (int i = 0; i < _d; ++i) v.set_coord(i, v[i] + a(i));
      
      //VT a = A.row(facet);
      //Point s(_d, std::vector<NT>(&a[0], a.data()+a.cols()*a.rows()));
      //s = ((-2.0 * v.dot(s)) * s);
      //v = s + v;
      
    }

    void free_them_all() {}

};

#endif
