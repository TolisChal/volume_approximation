// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
//#include <boost/math/distributions/students_t.hpp>
//#include <boost/math/special_functions/erf.hpp>
#include <iterator>
#include <vector>
#include <list>
#include <math.h>
#include <chrono>
#include "test_vol/test_cartesian_kernel.h"
#include "test_vol/test_vars.h"
#include "test_vol/test_hpolytope.h"
#include "test_vol/test_ball.h"
#include "test_vol/test_ballintersectconvex.h"
//#include "rounding.h"
#include "test_vol/test_samplers.h"


//' Sample points from a convex Polytope (H-polytope, V-polytope or a zonotope) or use direct methods for uniform sampling from the unit or the canonical or an arbitrary \eqn{d}-dimensional simplex and the boundary or the interior of a \eqn{d}-dimensional hypersphere
//'
//' Sample N points with uniform or multidimensional spherical gaussian -centered in an internal point- target distribution.
//' The \eqn{d}-dimensional unit simplex is the set of points \eqn{\vec{x}\in \R^d}, s.t.: \eqn{\sum_i x_i\leq 1}, \eqn{x_i\geq 0}. The \eqn{d}-dimensional canonical simplex is the set of points \eqn{\vec{x}\in \R^d}, s.t.: \eqn{\sum_i x_i = 1}, \eqn{x_i\geq 0}.
//'
//' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope.
//' @param N The number of points that the function is going to sample from the convex polytope. The default value is \eqn{100}.
//' @param distribution Optional. A string that declares the target distribution: a) \code{'uniform'} for the uniform distribution or b) \code{'gaussian'} for the multidimensional spherical distribution. The default target distribution is uniform.
//' @param random_walk Optional. A string that declares the random walk method: a) \code{'CDHR'} for Coordinate Directions Hit-and-Run, b) \code{'RDHR'} for Random Directions Hit-and-Run, c) \code{'BaW'} for Ball Walk or d) \code{'BiW'} for Billiard walk. The default walk is \code{'BiW'} for the uniform distribution or \code{'CDHR'} for the Normal distribution.
//' @param walk_length Optional. The number of the steps for the random walk. The default value is \eqn{1} for \code{'BiW'} and \eqn{\lfloor 10 + d/10\rfloor} otherwise, where \eqn{d} is the dimension that the polytope lies.
//' @param exact A boolean parameter. It should be used for the uniform sampling from the boundary or the interior of a hypersphere centered at the origin or from the unit or the canonical or an arbitrary simplex. The arbitrary simplex has to be given as a V-polytope. For the rest well known convex bodies the dimension has to be declared and the type of body as well as the radius of the hypersphere.
//' @param body A string that declares the type of the body for the exact sampling: a) \code{'unit simplex'} for the unit simplex, b) \code{'canonical simplex'} for the canonical simplex, c) \code{'hypersphere'} for the boundary of a hypersphere centered at the origin, d) \code{'ball'} for the interior of a hypersphere centered at the origin.
//' @param parameters Optional. A list for the parameters of the methods:
//' \itemize{
//' \item{\code{variance} }{ The variance of the multidimensional spherical gaussian. The default value is 1.}
//' \item{\code{dimension} }{ An integer that declares the dimension when exact sampling is enabled for a simplex or a hypersphere.}
//' \item{\code{radius} }{ The radius of the \eqn{d}-dimensional hypersphere. The default value is \eqn{1}.}
//' \item{\code{BW_rad} }{ The radius for the ball walk.}
//' \item{\code{diameter} }{ The diameter of the polytope. It is used to set the maximum length of the trajectory in billiard walk.}
//' }
//' @param InnerPoint A \eqn{d}-dimensional numerical vector that defines a point in the interior of polytope P.
//'
//' @references \cite{R.Y. Rubinstein and B. Melamed,
//' \dQuote{Modern simulation and modeling} \emph{ Wiley Series in Probability and Statistics,} 1998.}
//' @references \cite{A Smith, Noah and W Tromble, Roy,
//' \dQuote{Sampling Uniformly from the Unit Simplex,} \emph{ Center for Language and Speech Processing Johns Hopkins University,} 2004.}
//' @references \cite{Art B. Owen,
//' \dQuote{Monte Carlo theory, methods and examples,} \emph{ Art Owen,} 2009.}
//'
//' @return A \eqn{d\times N} matrix that contains, column-wise, the sampled points from the convex polytope P.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix test_sample(Rcpp::Nullable<Rcpp::Reference> P = R_NilValue,
                                  Rcpp::Nullable<unsigned int> N = R_NilValue,
                                  Rcpp::Nullable<std::string> random_walk = R_NilValue,
                                  Rcpp::Nullable<unsigned int> walk_length = R_NilValue,
                                  Rcpp::Nullable<Rcpp::List> parameters = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericVector> InnerPoint = R_NilValue) {

    typedef double NT;
    typedef Cartesian <NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef boost::mt19937 RNGType;
    typedef HPolytope <Point> Hpolytope;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> MT;

    Hpolytope HP;

    int type, dim, numpoints;
    NT radius = 1.0, delta = -1.0, diam = -1.0;
    bool set_mean_point = false, cdhr = false, rdhr = false, ball_walk = false, gaussian = false, billiard = true;
    std::list <Point> randPoints;
    std::pair <Point, NT> InnerBall;

    numpoints = (!N.isNotNull()) ? 100 : Rcpp::as<unsigned int>(N);


    type = Rcpp::as<Rcpp::Reference>(P).field("type");
    dim = Rcpp::as<Rcpp::Reference>(P).field("dimension");
    unsigned int walkL = 1;


    Point MeanPoint;
    if (InnerPoint.isNotNull()) {
        if (Rcpp::as<Rcpp::NumericVector>(InnerPoint).size() != dim) {
            Rcpp::warning("Internal Point has to lie in the same dimension as the polytope P");
        } else {
            set_mean_point = true;
            MeanPoint = Point(dim, Rcpp::as < std::vector < NT > > (InnerPoint).begin(),
                              Rcpp::as < std::vector < NT > > (InnerPoint).end());
        }
    }
    if (walk_length.isNotNull()) walkL = Rcpp::as<unsigned int>(walk_length);

    NT a = 0.5;

    bool rand_only = false,
            NN = false,
            birk = false,
            verbose = false;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // the random engine with this seed
    RNGType rng(seed);
    boost::random::uniform_real_distribution<> urdist(0, 1);
    boost::normal_distribution<> urdist1(0, 1);

    Point shift(dim);
    if (type == 1) {

        // Hpolytope
        HP.init(dim, Rcpp::as<MT>(Rcpp::as<Rcpp::Reference>(P).field("A")),
                Rcpp::as<VT>(Rcpp::as<Rcpp::Reference>(P).field("b")));

        if (!set_mean_point || ball_walk || billiard) {
            InnerBall = HP.ComputeInnerBall();
            if (!set_mean_point) MeanPoint = InnerBall.first;
        }
        if (diam < 0.0) HP.comp_diam(diam, InnerBall.second);
        HP.shift(InnerBall.first.getCoefficients());
        shift = InnerBall.first;
        InnerBall.first = Point(dim);
        HP.normalize();
        HP.recompute_AA();
    } else {
        throw Rcpp::exception("Wrong input!");
    }

    vars <NT, RNGType> var(1, dim, walkL, 1, 0.0, 0.0, 0, 0.0, 0, InnerBall.second, diam, rng, urdist, urdist1,
                            delta, verbose, rand_only, false, NN, birk, ball_walk, cdhr, rdhr, billiard);


    Point p(dim);
    VT lamdas, Av, vec;
    lamdas.setZero(HP.num_of_hyperplanes());
    Av.setZero(HP.num_of_hyperplanes());
    vec.setZero(dim);

    test_rand_point_generator(HP, p, numpoints, var.walk_steps, randPoints, lamdas, Av, vec, var);
    //sampling_only<Point>(randPoints, HP, walkL, numpoints, gaussian, a, MeanPoint, var1, var2);


    MT RetMat(dim, numpoints);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit != randPoints.end(); rpit++, jj++)
        RetMat.col(jj) = (rpit->getCoefficients() + shift.getCoefficients());
    return Rcpp::wrap(RetMat);

}

