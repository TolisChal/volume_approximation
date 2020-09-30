// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_ADAPTIVE_BILLIARD_WALK_HPP
#define RANDOM_WALKS_ADAPTIVE_BILLIARD_WALK_HPP

#include "convex_bodies/ball.h"
#include "convex_bodies/ballintersectconvex.h"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/vpolytope.h"
#include "convex_bodies/vpolyintersectvpoly.h"
#include "convex_bodies/zpolytope.h"
#include "convex_bodies/zonoIntersecthpoly.h"


template<typename VT, typename NT>
bool stopping_criterion_is_true(VT const& p, VT const& v, NT const& lambda, VT const& p0, VT const& v0) 
{ 
    NT t1 = (v.dot((p0 - p))) / (v.norm()*v.norm());
    NT t2 = (v0.dot((p0-p))) / (v.dot(v0));

    if (v.dot(v0) > 0.0) {
        if ((t1 < t2)) 
        {
            return false;
        } else if (t1 > lambda && t2 > lambda) {
            return false;
        } else if (t1 < 0.0) {
            return false;
        }
    } else {
        if (t1 < 0.0 || t2 < 0.0) {
            return false;
        }
    }
    return true;
}

template<typename NT, typename VT>
VT get_point_on_the_trajectory(NT _lambda_prev, std::vector<VT> const& points_x, 
                            std::vector<VT> const& points_v, std::vector<NT> const& lamdas)
{
    VT p(points_x[0].rows());
    for (int i = 0; i<points_x.size(); i++)
    {
        if (_lambda_prev < lamdas[i])
        {
            p = points_x[i] + _lambda_prev * points_v[i];
        } else {
            _lambda_prev -= lamdas[i];
        }
    }

    return p;
}


// Billiard walk which accelarates each step for uniform distribution

struct AdaptiveBilliardWalk
{
  AdaptiveBilliardWalk(double L)
            :   param(L, true)
    {}

  AdaptiveBilliardWalk()
            :   param(0, false)
    {}

    struct parameters
    {
        parameters(double L, bool set)
                :   m_L(L), set_L(set)
        {}
        double m_L;
        bool set_L;
    };

    struct update_parameters
    {
        update_parameters()
                :   facet_prev(0), hit_ball(false), inner_vi_ak(0.0), ball_inner_norm(0.0)
        {}
        int facet_prev;
        bool hit_ball;
        double inner_vi_ak;
        double ball_inner_norm;
    };

    parameters param;


    template
            <
                    typename Polytope,
                    typename RandomNumberGenerator
            >
    struct Walk
    {
        typedef typename Polytope::PointType Point;
        typedef typename Polytope::MT MT;
        typedef typename Point::FT NT;

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P, Point const& p, RandomNumberGenerator &rng)
        {
            _update_parameters = update_parameters();
            _L = compute_diameter<GenericPolytope>
                ::template compute<NT>(P);
            _AA.noalias() = P.get_mat() * P.get_mat().transpose();
            initialize(P, p, rng);
        }

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P, Point const& p, RandomNumberGenerator &rng,
             parameters const& params)
        {
            _update_parameters = update_parameters();
            _L = params.set_L ? params.m_L
                              : compute_diameter<GenericPolytope>
                                ::template compute<NT>(P);
            _AA.noalias() = P.get_mat() * P.get_mat().transpose();
            initialize(P, p, rng);
        }

        template
                <
                        typename GenericPolytope
                >
        inline void apply(GenericPolytope const& P,
                          Point &p,   // a point to start
                          unsigned int const& walk_length,
                          RandomNumberGenerator &rng)
        {
            typedef typename GenericPolytope::VT VT;
            unsigned int n = P.dimension();
            NT T = 0.0;
            const NT dl = 0.995;
            int it;
            std::vector<VT> points_x, points_v;
            std::vector<NT> lamdas;
            
            if (P.is_in(_p) == 0){
                std::cout<<"is_in = "<<P.is_in(_p)<<std::endl;
                exit(-1);
            } 

            for (auto j=0u; j<walk_length; ++j)
            {
                T = 0.0;
                //std::cout<<"j = "<<j<<std::endl;
                _v = GetDirection<Point>::apply(n, rng);
                Point p0 = _p, v0 = _v;
                points_x.clear();
                points_v.clear();
                lamdas.clear();

                it = 0;
                //std::pair<NT, int> pbpair = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev,
                //                                                      _update_parameters);
                std::pair<NT, int> pbpair = P.line_first_positive_intersect(_p, _v, _lambdas, _Av, 
                                                                            _update_parameters);
                //if (T <= pbpair.first) {
                //    _p += (T * _v);
                //   _lambda_prev = T;
                //    continue;
                //}

                _lambda_prev = dl * pbpair.first;
                T += _lambda_prev;
                lamdas.push_back(_lambda_prev);
                points_x.push_back(_p.getCoefficients());
                points_v.push_back(_v.getCoefficients());
                _p += (_lambda_prev * _v);
                //std::cout<<"_p = "<<_p.getCoefficients()<<std::endl;
                //T -= _lambda_prev;
                P.compute_reflection(_v, _p, _update_parameters);

                if (P.is_in(_p) == 0){
                    std::cout<<"is_in = "<<P.is_in(_p)<<std::endl;
                    exit(-1);
                } 
                it++;

                while (it < 100*n)
                {
                    std::pair<NT, int> pbpair
                            = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _AA,
                                                        _update_parameters);
                    _lambda_prev = dl * pbpair.first;
                    T += _lambda_prev;
                    lamdas.push_back(_lambda_prev);
                    points_x.push_back(_p.getCoefficients());
                    points_v.push_back(_v.getCoefficients());
                    if (stopping_criterion_is_true(_p.getCoefficients(), _v.getCoefficients(), 
                                                   _lambda_prev, p0.getCoefficients(), v0.getCoefficients())) {
                         NT lam = rng.sample_urdist() * T;
                        _p = Point(get_point_on_the_trajectory(lam, points_x, points_v, lamdas));
                        //std::cout<<"_p final = "<<_p.getCoefficients()<<std::endl;
                        //p=_p;
                        break;
                    } else if (it == 100*n - 1) {
                        NT lam = rng.sample_urdist() * T;
                        _p = Point(get_point_on_the_trajectory(lam, points_x, points_v, lamdas));
                        break;
                    }
                    _p += (_lambda_prev * _v);
                    //std::cout<<"_p = "<<_p.getCoefficients()<<std::endl;
                    if (P.is_in(_p) == 0){
                        std::cout<<"is_in = "<<P.is_in(_p)<<std::endl;
                        exit(-1);
                    } 
                    P.compute_reflection(_v, _p, _update_parameters);
                    it++;
                }
                std::cout<<"it = "<<it<<std::endl;
                //if (it == 100*n) _p = p0;
            }
            p = _p;
        }

        inline void update_delta(NT L)
        {
            _L = L;
        }

    private :

        template
                <
                        typename GenericPolytope
                >
        inline void initialize(GenericPolytope const& P,
                               Point const& p,
                               RandomNumberGenerator &rng)
        {
            unsigned int n = P.dimension();
            const NT dl = 0.995;
            _lambdas.setZero(P.num_of_hyperplanes());
            _Av.setZero(P.num_of_hyperplanes());
            _p = p;
            _v = GetDirection<Point>::apply(n, rng);

            NT T = -std::log(rng.sample_urdist()) * _L;
            Point p0 = _p;
            int it = 0;

            std::pair<NT, int> pbpair
                    = P.line_first_positive_intersect(_p, _v, _lambdas, _Av, _update_parameters);
            if (T <= pbpair.first) {
                _p += (T * _v);
                _lambda_prev = T;
                return;
            }
            _lambda_prev = dl * pbpair.first;
            _p += (_lambda_prev * _v);
            T -= _lambda_prev;
            P.compute_reflection(_v, _p, _update_parameters);

            while (it <= 100*n)
            {
                std::pair<NT, int> pbpair
                        = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _AA, _update_parameters);
                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;
                    break;
                } else if (it == 100*n) {
                    _lambda_prev = rng.sample_urdist() * pbpair.first;
                    _p += (_lambda_prev * _v);
                    break;
                }
                _lambda_prev = dl * pbpair.first;
                _p += (_lambda_prev * _v);
                T -= _lambda_prev;
                P.compute_reflection(_v, _p, _update_parameters);
                it++;
            }
        }

        double _L;
        Point _p;
        Point _v;
        NT _lambda_prev;
        MT _AA;
        update_parameters _update_parameters;
        typename Point::Coeff _lambdas;
        typename Point::Coeff _Av;
    };

};


#endif



