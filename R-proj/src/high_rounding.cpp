// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.
//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "convex_bodies/hpolytope.h"
#include "test_vol/test_mve_rounding.hpp"
#include "extractMatPoly.h"

//' Internal rcpp function for the rounding of a convex polytope
//'
//' @param P A convex polytope (H- or V-representation or zonotope).
//' @param method Optional. The method to use for rounding, a) \code{'mve'} for the method based on mimimmum volume enclosing ellipsoid of a dataset, b) \code{'mve'} for the method based on maximum volume enclosed ellipsoid, (c) \code{'svd'} for the method based on svd decomposition.
//' @param seed Optional. A fixed seed for the number generator.
//'
//' @keywords internal
//'
//' @return A numerical matrix that describes the rounded polytope, a numerical matrix of the inverse linear transofmation that is applied on the input polytope, the numerical vector the the input polytope is shifted and the determinant of the matrix of the linear transformation that is applied on the input polytope.
// [[Rcpp::export]]
Rcpp::List test_high_rounding (Rcpp::Reference P, Rcpp::Nullable<std::string> method = R_NilValue,
                     Rcpp::Nullable<double> seed = R_NilValue){

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    unsigned int n = P.field("dimension"), walkL, type = P.field("type");
    std::string mthd = std::string("lpsolve");
    if(method.isNotNull()) {
        mthd =  Rcpp::as<std::string>(method);
    }

    Hpolytope HP;

    std::pair <Point, NT> InnerBall;
    Rcpp::NumericMatrix Mat;

    MT N = MT::Identity(n,n);
    VT N_shift = VT::Zero(n);
    NT svd_prod = 1.0;
    switch (type) {
        case 1: {
            // Hpolytope
            
                // full dimensional
            HP.init(n, Rcpp::as<MT>(P.field("A")), Rcpp::as<VT>(P.field("b")));
            HP.normalize();

            if (mthd.compare(std::string("ipm")) == 0) {
                NT tol = 0.00000001;
                std::pair<VT, NT> res;// = compute_max_inner_ball(HP.get_mat(), HP.get_vec(), 150, tol);
                InnerBall.second = res.second;
                InnerBall.first = Point(res.first);
            }else{
                InnerBall = HP.ComputeInnerBall();
            }
            std::cout<<"inner_point = "<<InnerBall.first.getCoefficients().transpose()<<std::endl;
            std::cout<<"radius = "<<InnerBall.second<<std::endl;
            break;
        }
        default: {
            throw Rcpp::exception("volesti does not support rounding for this representation currently.");
        }
    }

    std::pair< std::pair<MT, VT>, NT > round_res;
    switch (type) {
        case 1: {
            round_res = test_mve_rounding<MT, VT>(HP, InnerBall);
            Mat = extractMatPoly(HP);
            break;
        }
        default: {
            throw Rcpp::exception("volesti does not support rounding for this representation currently.");
        }
    }

    return Rcpp::List::create(Rcpp::Named("Mat") = Mat, Rcpp::Named("T") = Rcpp::wrap(round_res.first.first),
                              Rcpp::Named("shift") = Rcpp::wrap(round_res.first.second),
                              Rcpp::Named("round_value") = round_res.second);
}
