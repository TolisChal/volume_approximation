// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_GENERIC_STATIC_BILLIARD_WALK_HPP
#define RANDOM_WALKS_GENERIC_STATIC_BILLIARD_WALK_HPP

#include "sampling/sphere.hpp"


// Billiard walk which accelarates each step for uniform distribution

struct StaticBilliardWalk
{
    StaticBilliardWalk(double L)
            :   param(L, true)
    {}

    StaticBilliardWalk()
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
        typedef typename Polytope::VT VT;
        typedef typename Point::FT NT;

        template <typename GenericPolytope>
        Walk(GenericPolytope const& P, Point const& p, RandomNumberGenerator &rng)
        {
            _update_parameters = update_parameters();
            _L = compute_diameter<GenericPolytope>
                ::template compute<NT>(P);
            _AA.noalias() = P.get_mat() * P.get_mat().transpose();
            _A = P.get_mat();
            //invHes.noalias() = P.get_mat().transpose() * P.get_mat();
            //detH = invHes.determinant();
            //invHes = invHes.inverse();
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
            unsigned int n = P.dimension();
            NT T;
            const NT dl = 0.995;
            int it;

            for (auto j=0u; j<walk_length; ++j)
            {
                _lambda_prev0 = _lambda_prev;
                _lambdas_prev = _lambdas;
                _Av_prev = _Av;


                T = -std::log(rng.sample_urdist()) * _L;
                _v = GetDirection<Point>::apply(n, rng, false);
                compute_hessian(p, _A, P.get_vec(), Hes);
                invHes = Hes.inverse();
                Eigen::LLT<MT> lltOfA(invHes); // compute the Cholesky decomposition of A
                MT L = lltOfA.matrixL(); // retrieve factor L  in the decomposition
                _v = Point(L * _v.getCoefficients());
                Point v0 = _v;
                Point p0 = _p;

                it = 0;
                std::pair<NT, int> pbpair = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, 
                                                                      _update_parameters);
                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;
                    continue;
                }

                _lambda_prev = dl * pbpair.first;
                _p += (_lambda_prev * _v);
                T -= _lambda_prev;
                P.compute_reflection(_v, _p, _update_parameters);
                it++;

                while (it < 100*n)
                {
                    std::pair<NT, int> pbpair
                            = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, _AA, _update_parameters);
                    if (T <= pbpair.first) {
                        _p += (T * _v);
                        _lambda_prev = T;
                        break;
                    }
                    _lambda_prev = dl * pbpair.first;
                    _p += (_lambda_prev * _v);
                    T -= _lambda_prev;
                    P.compute_reflection(_v, _p, _update_parameters);
                    it++;
                }
                if (it == 100*n) _p = p0;

                NT VxSVx = -0.5*(v0.getCoefficients().dot(Hes * v0.getCoefficients())), 
                    logdetHx = 0.5*std::log(Hes.determinant());

                compute_hessian(_p, P.get_mat(), P.get_vec(), Hes);

                NT VySVy = -0.5 * (_v.getCoefficients().dot(Hes * _v.getCoefficients())), 
                    logdetHy = 0.5*std::log(Hes.determinant());

                NT log_prob = VySVy + logdetHy - VxSVx - logdetHx;
                NT u_prog = std::log(rng.sample_urdist());

                if (u_prog > log_prob) {
                    _p = p0;
                    _Av = _Av_prev;
                    _lambdas = _lambdas_prev;
                    _lambda_prev = _lambda_prev0;
                }
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
            _v = GetDirection<Point>::apply(n, rng, false);
            compute_hessian(p, _A, P.get_vec(), Hes);
            invHes = Hes.inverse();
            Eigen::LLT<MT> lltOfA(invHes); // compute the Cholesky decomposition of A
            MT L = lltOfA.matrixL(); // retrieve factor L  in the decomposition
            _v = Point(L * _v.getCoefficients());

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


        inline void compute_hessian(Point const& x, MT const& A, VT const& b, MT &Hes) {
            int m = A.rows(), d = A.cols();
            Hes = MT::Zero(d, d);
            NT b_Ax;

            for (int i = 0; i < m; i++)
            {
                b_Ax = b(i) - A.row(i) * x.getCoefficients();
                Hes += (A.row(i).transpose()*A.row(i)) * (1.0 / (b_Ax * b_Ax));
            }
            
        }

        double _L;
        Point _p;
        Point _v;
        NT _lambda_prev, _lambda_prev0;
        MT _AA;
        MT _A;
        MT invHes;
        MT Hes;
        NT detH;
        std::vector<MT> RankOneMatrices;
        update_parameters _update_parameters;
        typename Point::Coeff _lambdas;
        typename Point::Coeff _Av;
        typename Point::Coeff _lambdas_prev;
        typename Point::Coeff _Av_prev;
    };

};


#endif


