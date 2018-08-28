// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#include "doctest.h"
#include <unistd.h>
#include "Eigen/Eigen"
#include "volume.h"
#include <typeinfo>

template <typename NT>
NT factorial(NT n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

template <typename NT, typename FilePath>
void test_volume(FilePath f, NT expected, NT tolerance=0.1)
{
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef boost::mt19937    RNGType;
    std::ifstream inp;
    std::vector<std::vector<NT> > Pin;
    inp.open(f,std::ifstream::in);
    read_pointset(inp,Pin);
    int n = Pin[0][1]-1;
    HPolytope<Point> P;
    P.init(Pin);

    // Setup the parameters
    int walk_len=10 + n/10;
    int nexp=1, n_threads=1;
    NT e=1, err=0.0000000001;
    int rnum = std::pow(e,-2) * 400 * n * std::log(n);
    RNGType rng(std::time(0));
    boost::normal_distribution<> rdist(0,1);
    boost::random::uniform_real_distribution<>(urdist);
    boost::random::uniform_real_distribution<> urdist1(-1,1);
    
    vars<NT, RNGType> var(rnum,n,walk_len,n_threads,err,e,0,0,0,0,rng,
             urdist,urdist1,-1.0,false,false,false,false,false,false,true);

    //Compute chebychev ball//
    std::pair<Point,NT> CheBall;

    // Estimate the volume
    std::cout << "--- Testing volume of " << f << std::endl;
    std::cout << "Number type: " << typeid(NT).name() << std::endl;
    NT vol = 0;
    unsigned int const num_of_exp = 10;
    for (unsigned int i=0; i<num_of_exp; i++)
    {
        CheBall = P.ComputeInnerBall();
        vol += volume(P,var,var,CheBall);
    }
    NT error = std::abs(((vol/num_of_exp)-expected))/expected;
    std::cout << "Computed volume (average) = " << vol/num_of_exp << std::endl;
    std::cout << "Expected volume = " << expected << std::endl;
    CHECK(error < tolerance);
}


template <typename NT>
void call_test_cube(){
    test_volume<NT>("../data/cube10.ine", 1024.0);
    test_volume<NT>("../data/cube20.ine", 1048576.0);
    test_volume<NT>("../data/cube30.ine", 1073742000.0, 0.2);
}

template <typename NT>
void call_test_cross(){
    test_volume<NT>("../data/cross_10.ine", 0.0002821869);
}

template <typename NT>
void call_test_birk() {
    test_volume<NT>("../data/birk3.ine", 0.125);
    test_volume<NT>("../data/birk4.ine", 0.000970018);
    test_volume<NT>("../data/birk5.ine", 0.000000225);
    test_volume<NT>("../data/birk6.ine", 0.0000000000009455459196, 0.5);
}

template <typename NT>
void call_test_prod_simplex() {
    test_volume<NT>("../data/prod_simplex_5_5.ine", std::pow(1.0 / factorial(5.0), 2));
    test_volume<NT>("../data/prod_simplex_10_10.ine", std::pow(1.0 / factorial(10.0), 2));
    test_volume<NT>("../data/prod_simplex_15_15.ine", std::pow(1.0 / factorial(15.0), 2));
    test_volume<NT>("../data/prod_simplex_20_20.ine", std::pow(1.0 / factorial(20.0), 2));
}

template <typename NT>
void call_test_simplex() {
    test_volume<NT>("../data/simplex10.ine", 1.0 / factorial(10.0));
    test_volume<NT>("../data/simplex20.ine", 1.0 / factorial(20.0));
    test_volume<NT>("../data/simplex30.ine", 1.0 / factorial(30.0));
    test_volume<NT>("../data/simplex40.ine", 1.0 / factorial(40.0));
    test_volume<NT>("../data/simplex50.ine", 1.0 / factorial(50.0));
}

template <typename NT>
void call_test_skinny_cube() {
    test_volume<NT>("../data/skinny_cube10.ine", 102400.0);
    test_volume<NT>("../data/skinny_cube20.ine", 104857600.0);
}


TEST_CASE("cube") {
    call_test_cube<double>();
    call_test_cube<float>();
    call_test_cube<long double>();
}

TEST_CASE("cross") {
    call_test_cross<double>();
    call_test_cross<float>();
    call_test_cross<long double>();
}

TEST_CASE("birk") {
    call_test_birk<double>();
    call_test_birk<float>();
    call_test_birk<long double>();
}

TEST_CASE("prod_simplex") {
    call_test_prod_simplex<double>();
    //call_test_prod_simplex<float>();
    call_test_prod_simplex<long double>();
}

TEST_CASE("simplex") {
    call_test_simplex<double>();
    //call_test_simplex<float>();
    //call_test_simplex<long double>();
}

TEST_CASE("skinny_cube") {
    call_test_skinny_cube<double>();
    //call_test_skinny_cube<float>();
    call_test_skinny_cube<long double>();
}
