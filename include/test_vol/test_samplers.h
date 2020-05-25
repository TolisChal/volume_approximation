// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef TEST_RANDOM_SAMPLERS_H
#define TEST_RANDOM_SAMPLERS_H

//Eigen::EigenMultivariateNormal<double> normX_solver1(0.0, 1.0);


// Pick a random direction as a normilized vector
template <typename RNGType, typename Point, typename VT, typename parameters>
void test_get_direction(const unsigned int dim, Point &vec, const parameters &var) {

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
        vec.set_coord(i, nrdist(rng));
        //normal += Xs[i] * Xs[i];
    }
    vec *= (1.0 / vec.length());
    //normal=1.0/std::sqrt(normal);

    //for (unsigned int i=0; i<dim; i++) {
        //Xs[i] = Xs[i] * normal;
    //}
    //Point p(vec);
    //return p;



    //vec<< normX_solver1.samples(dim).transpose();
    //std::cout<<vec<<std::endl;
    //Point p(vec/vec.norm());
    //return p;

}


// Pick a random point from a d-sphere
template <typename RNGType, typename Point, typename NT, typename VT, typename parameters>
void test_get_point_on_Dsphere(const unsigned int dim, const NT &radius, Point &vec, const parameters &var){
    test_get_direction<RNGType, Point, NT>(dim, vec, var);
    if (radius > NT(0)) vec *= radius;
    //vec = (radius == 0) ? p : radius * p;
    //return p;
}


// Pick a random point from a d-ball
template <typename RNGType, typename Point, typename NT, typename parameters>
void test_get_point_in_Dsphere(const unsigned int dim, const NT &radius, Point &vec, const parameters &var){

    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist = var.urdist;
    //boost::random::uniform_real_distribution<> urdist(0,1);
    //NT U;
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //RNGType rng2(seed);
    test_get_direction<RNGType, Point, NT>(dim, vec, var);
    NT U = urdist(rng);
    U = std::pow(U, 1.0/(NT(dim)));
    vec *= (radius*U);
    //return p;
}

// ----- RANDOM POINT GENERATION FUNCTIONS ------------ //


template <typename Polytope, typename PointList, typename Parameters, typename Point, typename VT, typename NT>
void test_rand_point_generator(Polytope &P,
                         Point &p,   // a point to start
                         const unsigned int rnum,
                         const unsigned int walk_len,
                         PointList &randPoints,
                         VT &lamdas,
                         VT &Av,
                         Point &vec,
                               const NT &diameter,
                         Parameters &var)  // constants for volume
{
    //typedef typename Polytope::VT VT;
    //typedef typename Point::FT NT;

    //VT lamdas, Av;
    //lamdas.setZero(P.num_of_hyperplanes());
    //Av.setZero(P.num_of_hyperplanes());

    NT lambda;
    var.nsteps += double(rnum * walk_len);

    test_billiard_walk(P, p, diameter, lamdas, Av, lambda, var, vec, true);
    for (unsigned int j = 0; j < walk_len-1; ++j){
        test_billiard_walk(P, p, diameter, lamdas, Av, lambda,  var, vec);

    }
    randPoints.push_back(p);

    //rnum--;
    for (unsigned int i = 1; i <= rnum-1; ++i) {
        for (unsigned int j = 0; j < walk_len; ++j) {
            test_billiard_walk(P, p, diameter, lamdas, Av, lambda,  var, vec);
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
                         Point &vec,
                         NT &lambda,
                              const NT &diameter,
                         Parameters &var) {

    test_billiard_walk(P, p, diameter, lamdas, Av, lambda, var, vec, true);
    var.nsteps += double(walk_len);
    walk_len--;

    for (unsigned int j = 0; j < walk_len; j++){
        test_billiard_walk(P, p, diameter, lamdas, Av, lambda, var, vec);
    }

}



template <typename Polytope, typename Point, typename Parameters, typename NT, typename VT>
void test_uniform_next_point(Polytope &P,
                        Point &p,   // a point to start
                        const unsigned int walk_len, // number of steps for the random walk
                        VT &lamdas,
                        VT &Av,
                        Point &vec,
                        NT &lambda,
                        const NT &diameter,
                        Parameters &var) {
    var.nsteps += double(walk_len);
    for (unsigned int j = 0; j < walk_len; j++){
        test_billiard_walk(P, p, diameter, lamdas, Av, lambda, var, vec);
    }

}

template <class ConvexBody, class Point, class Parameters, typename NT, typename VT>
void test_billiard_walk(ConvexBody &P, Point &p, const NT &diameter, VT &Ar, VT &Av, NT &lambda_prev,
                        Parameters &var, Point &v, bool first = false) {

    typedef typename Parameters::RNGType RNGType;
    unsigned int n = P.dimension();
    RNGType &rng = var.rng;
    boost::random::uniform_real_distribution<> urdist = var.urdist;
    NT T, inner_vi_ak;
    if (var.log_length) {
        T = -std::log(NT(urdist(rng))) * diameter;
    } else {
        T = NT(urdist(rng)) * diameter;
    }
    const NT dl = 0.995;
    //vec.setZero(n);
    test_get_direction<RNGType, Point, NT>(n, v, var);
    Point p0 = p;
    int it = 0;

    if (first) {
        std::pair<NT, int> pbpair = P.line_positive_intersect(p, v, Ar, Av, inner_vi_ak);
        //var.nboracles += 1.0;
        var.noracles = var.noracles + 1.0;
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
        //var.nboracles += 1.0;
        var.noracles = var.noracles + 1.0;
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
    //std::cout<<"bref = "<<var.lw<<std::endl;
    while (it<var.lw) {
        std::pair<NT, int> pbpair = P.line_positive_intersect(p, v, Ar, Av, lambda_prev, inner_vi_ak);
        //var.nboracles += 1.0;
        var.noracles = var.noracles + 1.0;
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

    if(it == var.lw){
        std::cout<<"limit reached"<<std::endl;
        p = p0;
    }
}


#endif //RANDOM_SAMPLERS_H
