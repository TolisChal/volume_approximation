// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_SUPEROPTIMIZED_GENERIC_STATIC_BILLIARD_WALK_HPP
#define RANDOM_WALKS_SUPEROPTIMIZED_GENERIC_STATIC_BILLIARD_WALK_HPP

#include "sampling/sphere.hpp"


// Billiard walk which accelarates each step for uniform distribution

struct StaticSuperOptBilliardWalk
{
    StaticSuperOptBilliardWalk(double L)
            :   param(L, true)
    {}

    StaticSuperOptBilliardWalk()
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
            _Atrans = P.get_mat().transpose();
            _b = P.get_vec();
            _Ws.setZero(P.get_mat().rows());
            _temp.setZero(P.get_mat().rows());
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
            _Atrans = P.get_mat().transpose();
            _b = P.get_vec();
            _Ws.setZero(P.get_mat().rows());
            _temp.setZero(P.get_mat().rows());
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
                //compute_hessian(_p, P.get_mat(), P.get_vec(), _Hes);
                //_invHes = _Hes.inverse();
                //std::cout<<"_invHes = "<<_invHes<< "\n _optInvHes = "<<_optInvHes<<"\n"<<std::endl;
                
                //Eigen::LLT<MT> lltOfA(_invHes); // compute the Cholesky decomposition of A
                //MT L1 = lltOfA.matrixL(), L2 = _lltOfInvHes_next.matrixL();
                //std::cout<<"lltOfA(_invHes) = "<<L1<< "\n _lltOfInvHes_next = "<<L2<<"\n"<<std::endl;
                //MT L = lltOfA.matrixL(); // retrieve factor L  in the decomposition
                _v = Point(_lltOfInvHes_next.matrixL() * _v.getCoefficients());
                Point v0 = _v;
                Point p0 = _p;

                it = 0;
                std::pair<NT, int> pbpair = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev, 
                                                                      _update_parameters);
                _temp = _b - _lambdas;
                _VxSVxOpt = -0.5*(_Av.cwiseProduct(_Av).cwiseProduct(_temp.cwiseProduct(_temp).cwiseInverse())).sum();
                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;

                    //NT VxSVx = -0.5*(v0.getCoefficients().dot(_Hes * v0.getCoefficients())), 
                    //logdetHx = 0.5*std::log(_Hes.determinant());

                    update_filter_variables(_b, _Atrans, _Ws, _lambdas_prev, _Av_prev, _lambda_prev0, 
                                            _lambdas, _Av, T, _VySVyOpt, _optInvHes, _logdetHyOpt, _lltOfInvHes_next);

                    //compute_Ws(_b, _A, _Ws, _lambdas_prev, _Av_prev, _lambda_prev0, _lambdas, _Av, T);// P.get_mat(), P.get_vec(), Hes);

                    //VxSVxOpt = -0.5 * (_v.getCoefficients().dot(Hes * _v.getCoefficients())), 
                    //    logdetHxOpt = 0.5*std::log(Hes.determinant());
                    //_temp = _b - (_lambdas + T*_Av);
                    //_VySVyOpt = -0.5*(_Av.cwiseProduct(_Av).cwiseProduct(_temp.cwiseProduct(_temp).cwiseInverse())).sum();
                    //update_determinant_cholesky_inverse(_Ws, _Atrans, _optInvHes, _logdetHyOpt, _lltOfInvHes_next);
                    //_logdetHyOpt *= 0.5;
                    
                    //compute_hessian(_p, P.get_mat(), P.get_vec(), _Hes);

                    //NT VySVy = -0.5 * (_v.getCoefficients().dot(_Hes * _v.getCoefficients())), 
                     //   logdetHy = 0.5*std::log(_Hes.determinant());

                    //NT log_prob = VySVy + logdetHy - VxSVx - logdetHx;
                    NT log_prob = _VySVyOpt + _logdetHyOpt*0.5 - _VxSVxOpt - _logdetHxOpt*0.5;
                    NT u_prog = std::log(rng.sample_urdist());

                    //std::cout<<"[1]log_prob = "<<log_prob<<", log_prob2 = "<<log_prob2<<std::endl;
                    //std::cout<<"[1]VxSVx = "<<VxSVx<<", logdetHx = "<< logdetHx<< " VySVy = "<< VySVy<< " logdetHy = "<< logdetHy<<std::endl;
                    //std::cout<<"[1]VxSVxOpt = "<< _VxSVxOpt<< ", logdetHxOpt = "<< _logdetHxOpt*0.5<< " VySVyOpt = "<< _VySVyOpt<< " _logdetHyOpt = "<< _logdetHyOpt*0.5<<"\n"<<std::endl;

                    if (u_prog > log_prob) {
                        _p = p0;
                        _Av = _Av_prev;
                        _lambdas = _lambdas_prev;
                        _lambda_prev = _lambda_prev0;
                        _lltOfInvHes_next = _lltOfInvHes_prev;
                        _logdetHyOpt = _logdetHxOpt;
                        _optInvHes = _optInvHes_prev;
                    } else {
                        _logdetHxOpt = _logdetHyOpt;
                        _lltOfInvHes_prev = _lltOfInvHes_next;
                        _optInvHes_prev = _optInvHes;
                    }


                    continue;
                }

                _lambda_prev = dl * pbpair.first;
                _p += (_lambda_prev * _v);
                T -= _lambda_prev;
                P.compute_reflection(_v, _p, _update_parameters);
                it++;

                while (it < 1000*n)
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
                /*if (it == 1000*n){
                     _p = p0;
                    _Av = _Av_prev;
                    _lambdas = _lambdas_prev;
                    _lambda_prev = _lambda_prev0;
                    continue;
                }*/

                //NT VxSVx = -0.5*(v0.getCoefficients().dot(_Hes * v0.getCoefficients())), 
                //    logdetHx = 0.5*std::log(_Hes.determinant());
                
                update_filter_variables(_b, _Atrans, _Ws, _lambdas_prev, _Av_prev, _lambda_prev0, 
                                            _lambdas, _Av, T, _VySVyOpt, _optInvHes, _logdetHyOpt, _lltOfInvHes_next);

                //_temp = _b - (_lambdas + _lambda_prev*_Av);
                //_VySVyOpt = -0.5*(_Av.cwiseProduct(_Av).cwiseProduct(_temp.cwiseProduct(_temp).cwiseInverse())).sum();

                //compute_Ws(_b, _Ws, _lambdas_prev, _Av_prev, _lambda_prev0, _lambdas, _Av, _lambda_prev);// P.get_mat(), P.get_vec(), Hes);

                //VxSVxOpt = -0.5 * (_v.getCoefficients().dot(Hes * _v.getCoefficients())), 
                //    logdetHxOpt = 0.5*std::log(Hes.determinant());

                
                //update_determinant_cholesky_inverse(_Ws, _Atrans, _optInvHes, _logdetHyOpt, _lltOfInvHes_next);
                //_logdetHyOpt *= 0.5;
                    
                //compute_hessian(_p, P.get_mat(), P.get_vec(), _Hes);

                //NT VySVy = -0.5 * (_v.getCoefficients().dot(_Hes * _v.getCoefficients())), 
                //    logdetHy = 0.5*std::log(_Hes.determinant());

                //NT log_prob = VySVy + logdetHy - VxSVx - logdetHx;
                NT log_prob = _VySVyOpt + _logdetHyOpt*0.5 - _VxSVxOpt - _logdetHxOpt*0.5;
                NT u_prog = std::log(rng.sample_urdist());

                //std::cout<<"[2]log_prob = "<<log_prob<<", log_prob2 = "<<log_prob2<<std::endl;
                //std::cout<<"[2]VxSVx = "<<VxSVx<<", logdetHx = "<< logdetHx<< " VySVy = "<< VySVy<< " logdetHy = "<< logdetHy<<std::endl;
                //std::cout<<"[2]VxSVxOpt = "<< _VxSVxOpt<< ", logdetHxOpt = "<< _logdetHxOpt*0.5<< " VySVyOpt = "<< _VySVyOpt<< " _logdetHyOpt = "<< _logdetHyOpt*0.5<<"\n"<<std::endl;

                if (u_prog > log_prob) {
                    _p = p0;
                    _Av = _Av_prev;
                    _lambdas = _lambdas_prev;
                    _lambda_prev = _lambda_prev0;
                    _lltOfInvHes_next = _lltOfInvHes_prev;
                    _logdetHyOpt = _logdetHxOpt;
                    _optInvHes = _optInvHes_prev;
                } else {
                    _logdetHxOpt = _logdetHyOpt;
                    _lltOfInvHes_prev = _lltOfInvHes_next;
                    _optInvHes_prev = _optInvHes;
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
            compute_hessian(p, P.get_mat(), P.get_vec(), _optHes);
            _optInvHes = _optHes.inverse();
            _lltOfInvHes_next = Eigen::LLT<MT>(_optInvHes);
            //compute_hessian(p, P.get_mat(), P.get_vec(), _Hes);
            //_invHes = _Hes.inverse();
            //Eigen::LLT<MT> lltOfA(_invHes); // compute the Cholesky decomposition of A
            //MT L = lltOfA.matrixL(); // retrieve factor L  in the decomposition
            _v = Point(_lltOfInvHes_next.matrixL() * _v.getCoefficients());

            NT T = -std::log(rng.sample_urdist()) * _L;
            Point p0 = _p;
            int it = 0;

            std::pair<NT, int> pbpair
                    = P.line_first_positive_intersect(_p, _v, _lambdas, _Av, _update_parameters);
            if (T <= pbpair.first) {
                _p += (T * _v);
                _lambda_prev = T;

                compute_hessian(_p, P.get_mat(), P.get_vec(), _optHes);
                _optInvHes = _optHes.inverse();
                _logdetHxOpt = std::log(_optHes.determinant());
                //std::cout<<"_p = "<<_p.getCoefficients().transpose()<<", _logdetHxOpt = "<< _logdetHxOpt<< std::endl;
                _logdetHyOpt = _logdetHxOpt;
                _lltOfInvHes_prev = Eigen::LLT<MT>(_optInvHes);
                _lltOfInvHes_next = _lltOfInvHes_prev;
                _optInvHes_prev = _optInvHes;

                return;
            }
            _lambda_prev = dl * pbpair.first;
            _p += (_lambda_prev * _v);
            T -= _lambda_prev;
            P.compute_reflection(_v, _p, _update_parameters);

            while (it <= 1000*n)
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
            compute_hessian(_p, P.get_mat(), P.get_vec(), _optHes);
            _optInvHes = _optHes.inverse();
            _logdetHxOpt = std::log(_optHes.determinant());
            //std::cout<<"_p = "<<_p.getCoefficients().transpose()<<", _logdetHxOpt = "<< _logdetHxOpt<< std::endl;
            _logdetHyOpt = _logdetHxOpt;
            _lltOfInvHes_prev = Eigen::LLT<MT>(_optInvHes);
            _lltOfInvHes_next = _lltOfInvHes_prev;
            _optInvHes_prev = _optInvHes;
        }

        inline void update_filter_variables(VT const& b, MT const& Atrans, VT &Ws, VT const& lambdas1, VT const& Av1, NT const& t1, 
                                            VT const& lambdas2, VT const& Av2, NT const& t2, NT &VySVyOpt, MT &optInvHes, NT &logdetHyOpt, Eigen::LLT<MT> &lltOfInvHes_next)
        {
            NT g;
            VT v, u;
            MT M;

            _temp = b - (lambdas1 + t1*Av1);
            VT p1 = _temp.cwiseProduct(_temp);
            _temp = b - (lambdas2 + t2*Av2);
            VySVyOpt = -0.5*(Av2.cwiseProduct(Av2).cwiseProduct(_temp.cwiseProduct(_temp).cwiseInverse())).sum();
            VT p2 = _temp.cwiseProduct(_temp);

            Ws.noalias() = (p1-p2).cwiseProduct((p1.cwiseProduct(p2)).cwiseInverse());

            for (int i = 0; i < Ws.size(); i++)
            {
                v = Atrans.col(i);// * std::sqrt(std::abs(Ws(i)));
                //std::cout<<"v*v' = "<<v*v.transpose()* Ws(i)<<", wiai*ai^T = "<<Ws(i)*(Atrans.col(i) * Atrans.col(i).transpose())<<std::endl;
                //std::cout<<"v = "<<v.transpose()<< ", optInvHes*v = "<<optInvHes*v<<std::endl;
                g = (v.cwiseProduct(optInvHes*v)).sum()* Ws(i);
                //std::cout<<"g = "<<g<<", v'*H*v = "<<v.transpose()*optInvHes*v*Ws(i)<<std::endl;
                logdetHyOpt += std::log(1.0 + g);
                g = -(1.0 / (1.0 + g));
                u.noalias() = (optInvHes * v);// *(std::sqrt(std::abs(Ws(i))) * sgn(Ws(i)));
                //std::cout<<"u = "<<u.transpose()<<std::endl;
                M.noalias() = u*u.transpose()* Ws(i);
                //std::cout<<"M = \n"<<M<<std::endl;
                optInvHes += g*M;
                lltOfInvHes_next.rankUpdate(u* std::sqrt(std::abs(Ws(i))), g*sgn(Ws(i)));
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

        inline void compute_Ws(VT const& b, VT &Ws, VT const& lambdas1, VT const& Av1, NT const& t1,
                               VT const& lambdas2, VT const& Av, NT const& t) {
            
            _temp = b - (lambdas1 + t1*Av1);
            VT p1 = _temp.cwiseProduct(_temp);
            _temp = b - (lambdas2 + t*Av);
            VT p2 = _temp.cwiseProduct(_temp);

            Ws.noalias() = (p1-p2).cwiseProduct((p1.cwiseProduct(p2)).cwiseInverse());

        }

        inline void update_determinant_cholesky_inverse(VT const& Ws, MT const& Atrans, MT &optInvHes,
                                                        NT &logdetHyOpt, Eigen::LLT<MT> &lltOfInvHes_next) {
            NT g;
            VT v, u;
            MT M;
            std::cout<<"Ws.size() = "<<Ws.size()<<std::endl;
            for (int i = 0; i < Ws.size(); i++)
            {
                v = Atrans.col(i);// * std::sqrt(std::abs(Ws(i)));
                std::cout<<"v*v' = "<<v*v.transpose()* Ws(i)<<", wiai*ai^T = "<<Ws(i)*(Atrans.col(i) * Atrans.col(i).transpose())<<std::endl;
                std::cout<<"v = "<<v.transpose()<< ", optInvHes*v = "<<optInvHes*v<<std::endl;
                g = (v.cwiseProduct(optInvHes*v)).sum()* Ws(i);
                std::cout<<"g = "<<g<<", v'*H*v = "<<v.transpose()*optInvHes*v*Ws(i)<<std::endl;
                logdetHyOpt += std::log(1.0 + g);
                g = -(1.0 / (1.0 + g));
                u.noalias() = (optInvHes * v);// *(std::sqrt(std::abs(Ws(i))) * sgn(Ws(i)));
                std::cout<<"u = "<<u.transpose()<<std::endl;
                M.noalias() = u*u.transpose()* Ws(i);
                std::cout<<"M = \n"<<M<<std::endl;
                optInvHes += g*M;
                lltOfInvHes_next.rankUpdate(u* std::sqrt(std::abs(Ws(i))), g*sgn(Ws(i)));
            }
        }

        template <typename T> int sgn(T val) {
            return (T(0) < val) - (val < T(0));
        }

        double _L;
        Point _p;
        Point _v;
        NT _lambda_prev, _lambda_prev0;
        MT _AA, _Atrans;
        MT _optInvHes, _optInvHes_prev;
        MT _optHes;
        VT _b, _Ws, _temp;
        NT _VxSVxOpt, _VySVyOpt, _logdetHxOpt, _logdetHyOpt;
        Eigen::LLT<MT> _lltOfInvHes_prev, _lltOfInvHes_next;
        update_parameters _update_parameters;
        typename Point::Coeff _lambdas;
        typename Point::Coeff _Av;
        typename Point::Coeff _lambdas_prev;
        typename Point::Coeff _Av_prev;
    };

};


#endif



