// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef TEST_RANDOM_SAMPLERS_H
#define TEST_RANDOM_SAMPLERS_H

//Eigen::EigenMultivariateNormal<double> normX_solver1(0.0, 1.0);


// Pick a random direction as a normilized vector
template <typename RNGType, typename Point, typename NT, typename VT, typename parameters>
Point test_get_direction(const unsigned int dim, VT &vec, parameters &var) {

    //boost::normal_distribution<> rdist(0,1);
    boost::normal_distribution<> nrdist = var.urdist1;//(0,1);
    RNGType &rng = var.rng;
    //boost::random::uniform_real_distribution<> rdist = var.rdist;
    //std::vector<NT> Xs(dim,0);
    //NT normal = NT(0);
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //RNGType rng(seed);
    //RNGType rng2 = var.rng;
    for (unsigned int i=0; i<dim; i++) {
        vec(i) = nrdist(rng);
        //normal += Xs[i] * Xs[i];
    }
    vec *= (1.0/vec.norm());
    //normal=1.0/std::sqrt(normal);

    //for (unsigned int i=0; i<dim; i++) {
        //Xs[i] = Xs[i] * normal;
    //}
    Point p(vec);
    return p;



    //vec<< normX_solver1.samples(dim).transpose();
    //std::cout<<vec<<std::endl;
    //Point p(vec/vec.norm());
    //return p;

}


// Pick a random point from a d-sphere
template <typename RNGType, typename Point, typename NT, typename VT, typename parameters>
Point test_get_point_on_Dsphere(const unsigned int dim, const NT &radius, VT &vec, parameters &var){
    Point p = test_get_direction<RNGType, Point, NT>(dim, vec, var);
    p = (radius == 0) ? p : radius * p;
    return p;
}


// Pick a random point from a d-ball
template <typename RNGType, typename Point, typename NT, typename VT, typename parameters>
Point test_get_point_in_Dsphere(const unsigned int dim, const NT &radius, VT &vec, parameters &var){

    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist = var.urdist;
    //boost::random::uniform_real_distribution<> urdist(0,1);
    //NT U;
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //RNGType rng2(seed);
    Point p = test_get_direction<RNGType, Point, NT>(dim, vec, var);
    NT U = urdist(rng);
    U = std::pow(U, 1.0/(NT(dim)));
    p *= (radius*U);
    return p;
}

// ----- RANDOM POINT GENERATION FUNCTIONS ------------ //


template <typename Polytope, typename PointList, typename Parameters, typename Point, typename VT>
void test_rand_point_generator(Polytope &P,
                         Point &p,   // a point to start
                         const unsigned int rnum,
                         const unsigned int walk_len,
                         PointList &randPoints,
                         VT &lamdas,
                         VT &Av,
                         VT &vec,
                         const Parameters &var)  // constants for volume
{
    //typedef typename Polytope::VT VT;
    typedef typename Point::FT NT;

    //VT lamdas, Av;
    //lamdas.setZero(P.num_of_hyperplanes());
    //Av.setZero(P.num_of_hyperplanes());

    NT lambda;

    test_billiard_walk(P, p, var.diameter, lamdas, Av, lambda, var, vec, true);
    for (unsigned int j = 0; j < walk_len-1; ++j){
        test_billiard_walk(P, p, var.diameter, lamdas, Av, lambda,  var, vec);
    }
    randPoints.push_back(p);

    //rnum--;
    for (unsigned int i = 1; i <= rnum-1; ++i) {
        for (unsigned int j = 0; j < walk_len; ++j) {
            test_billiard_walk(P, p, var.diameter, lamdas, Av, lambda,  var, vec);
        }
        randPoints.push_back(p);
    }
 
}

template <typename Polytope, typename Point, typename Parameters, typename NT, typename VT>
void test_uniform_first_point(Polytope &P,
                         Point &p,   // a point to start
                         unsigned int walk_len, // number of steps for the random walk
                         VT &lamdas,
                         VT &Av,
                         VT &vec,
                         NT &lambda,
                         const Parameters &var) {

    test_billiard_walk(P, p, var.diameter, lamdas, Av, lambda, var, vec, true);
    walk_len--;

    for (unsigned int j = 0; j < walk_len; j++) test_billiard_walk(P, p, var.diameter, lamdas, Av, lambda, var, vec);

}



template <typename Polytope, typename Point, typename Parameters, typename NT, typename VT>
void test_uniform_next_point(Polytope &P,
                        Point &p,   // a point to start
                        const unsigned int walk_len, // number of steps for the random walk
                        VT &lamdas,
                        VT &Av,
                        VT &vec,
                        NT &lambda,
                        const Parameters &var) {
    for (unsigned int j = 0; j < walk_len; j++) test_billiard_walk(P, p, var.diameter, lamdas, Av, lambda, var, vec);

}

template <class ConvexBody, class Point, class Parameters, typename NT, typename VT>
void test_billiard_walk(ConvexBody &P, Point &p, NT diameter, VT &Ar, VT &Av, NT &lambda_prev,
                   Parameters &var, VT &vec, bool first = false) {

    typedef typename Parameters::RNGType RNGType;
    unsigned int n = P.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist = var.urdist;
    NT T = urdist(rng) * diameter, inner_vi_ak;
    const NT dl = 0.995;
    vec.setZero(n);
    Point v = test_get_direction<RNGType, Point, NT>(n, vec, var), p0 = p;
    int it = 0;

    if (first) {
        std::pair<NT, int> pbpair = P.line_positive_intersect(p, v, Ar, Av, inner_vi_ak);
        if (T <= pbpair.first) {
            p += (T * v);
            lambda_prev = T;
            return;
        }
        lambda_prev = dl * pbpair.first;
        p += (lambda_prev * v);
        T -= lambda_prev;
        P.compute_reflection(v, p, inner_vi_ak, pbpair.second);
        it++;
    } else {
        std::pair<NT, int> pbpair = P.line_positive_intersect(p, v, Ar, Av, lambda_prev, inner_vi_ak, true);
        if (T <= pbpair.first) {
            p += (T * v);
            lambda_prev = T;
            return;
        }

        lambda_prev = dl * pbpair.first;
        p += (lambda_prev * v);
        T -= lambda_prev;
        P.compute_reflection(v, p, inner_vi_ak, pbpair.second);
        it++;
    }

    while (it<10*n) {
        std::pair<NT, int> pbpair = P.line_positive_intersect(p, v, Ar, Av, lambda_prev, inner_vi_ak);
        if (T <= pbpair.first) {
            p += (T * v);
            lambda_prev = T;
            break;
        }

        lambda_prev = dl * pbpair.first;
        p += (lambda_prev * v);
        T -= lambda_prev;
        P.compute_reflection(v, p, inner_vi_ak, pbpair.second);

        it++;
    }

    if(it == 10*n) p = p0;
}


#endif //RANDOM_SAMPLERS_H
