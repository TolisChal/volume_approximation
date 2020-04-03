// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

#ifndef TEST_BALL_ANNEALING_H
#define TEST_BALL_ANNEALING_H

#define MAX_ITER 10
#define TOL 0.00000000001

//An implementation of Welford's algorithm for mean and variance.
template <typename NT>
std::pair<NT, NT> test_getMeanVariance(std::vector<NT>& vec) {
    NT mean = 0, M2 = 0, variance = 0, delta;
    typedef typename std::vector<NT>::iterator viterator;

    unsigned int i=0;
    viterator vecit = vec.begin();
    for( ; vecit!=vec.end(); vecit++, i++){
        delta = *vecit - mean;
        mean += delta / (i + 1);
        M2 += delta * (*vecit - mean);
        variance = M2 / (i + 1);
    }

    return std::pair<NT, NT> (mean, variance);
}


template <typename Point, typename ConvexBody, typename PointList, typename NT>
bool test_check_convergence(ConvexBody &P, PointList &randPoints, const NT &lb, const NT &ub, bool &too_few, NT &ratio,
                      const int &nu, NT alpha, const bool &precheck, const bool &lastball) {

    std::vector<NT> ratios;
    std::pair<NT,NT> mv;
    int m = randPoints.size()/nu, i = 1;
    NT T, rs, alpha_check = 0.01;
    size_t countsIn = 0;

    for(typename std::list<Point>::iterator pit=randPoints.begin(); pit!=randPoints.end(); ++pit, i++){

        if (P.is_in(*pit)==-1) countsIn++;
        //if(precheck) var.nmoracles += 1.0;
        if (i % m == 0) {
            ratios.push_back(NT(countsIn)/m);
            countsIn = 0;
            if (ratios.size()>1 && precheck) {
                boost::math::students_t dist(ratios.size() - 1);
                mv = test_getMeanVariance(ratios);
                ratio = mv.first;
                rs = std::sqrt(mv.second);
                T = rs * (boost::math::quantile(boost::math::complement(dist, alpha_check / 2.0))
                          / std::sqrt(NT(ratios.size())));
                if (ratio + T < lb) {
                    too_few = true;
                    return false;
                } else if (ratio - T > ub) return false;
            }
        }
    }

    //NT alpha = 0.25;
    if(precheck) alpha *= 0.5;
    mv = test_getMeanVariance(ratios);
    ratio = mv.first;
    rs = std::sqrt(mv.second);
    boost::math::students_t dist(nu - 1);
    T = rs*(boost::math::quantile(boost::math::complement(dist, alpha))
            / std::sqrt(NT(nu)));
    if (ratio > lb + T) {
        if (lastball) return true;
        if ((precheck && ratio < ub - T) || (!precheck && ratio < ub + T)) return true;
        return false;
    }
    too_few = true;
    return false;

}


template <typename RNGType, typename Polytope, typename ball, typename NT, typename Point, typename parameters>
bool test_get_first_ball(Polytope &P, ball &B0, NT &ratio, NT rad1, const NT &lb, const NT &ub, const NT &alpha,
                  NT &rmax, Point &vec, const parameters &var){

    //typedef typename Polytope::PolytopePoint Point;
    //typedef typename Polytope::VT VT;
    int n = P.dimension(), iter = 1;
    bool bisection_int = false, pass = false, too_few = false;
    std::list<Point> randPoints;
    //Point p(n);
    //VT vec;
    //vec.setZero(n);

    if(rmax>0.0) {
        for (int i = 0; i < 1200; ++i) {
            test_get_point_in_Dsphere<RNGType, Point>(n, rmax, vec, var);
            randPoints.push_back(vec);
        }
        pass = test_check_convergence<Point>(P, randPoints, lb, ub, too_few, ratio, 10, alpha, true, false);
        if (pass || !too_few) {
            B0 = ball(Point(n), rmax*rmax);
            return true;
        }
        bisection_int = true;
    } else {
        rmax = 2 * std::sqrt(NT(n)) * rad1;
    }
    NT radius = rad1;
    //std::cout<<"rmax = "<<rmax<<" rad1 = "<<rad1<<std::endl;

    while(!bisection_int) {

        randPoints.clear();
        too_few = false;

        for (int i = 0; i < 1200; ++i){
            test_get_point_in_Dsphere<RNGType, Point>(n, rmax, vec, var);
            randPoints.push_back(vec);
        }

        if(test_check_convergence<Point>(P, randPoints, lb, ub, too_few, ratio, 10, alpha, true, false)) {
            B0 = ball(Point(n), rmax*rmax);
            return true;
        }

        if (too_few) break;
        rad1 = rmax;
        rmax = rmax + 2*std::sqrt(NT(n))*radius;
    }

    NT rad_med, rad0=rad1, rad_m = rmax;

    while(iter <= MAX_ITER) {

        rad_med = 0.5*(rad1+rmax);
        randPoints.clear();
        too_few = false;
        //std::cout << "rmax = " << rmax << " rad1 = " << rad1 << " rad_med = " << rad_med << std::endl;

        for (int i = 0; i < 1200; ++i){
            test_get_point_in_Dsphere<RNGType, Point>(n, rad_med, vec, var);
            randPoints.push_back(vec);
        }

        if(test_check_convergence<Point>(P, randPoints, lb, ub, too_few, ratio, 10, alpha, true, false)) {
            B0 = ball(Point(n), rad_med*rad_med);
            return true;
        }

        if (too_few) {
            rmax = rad_med;
        } else {
            rad1 = rad_med;
        }

        if(rmax-rad1 < TOL) {
            rad1 = rad0;
            rmax = rad_m;
            iter++;
        }

    }
    return false;
}


template <typename Point, typename ball, typename PointList, typename NT>
bool test_get_next_zonoball(std::vector<ball> &BallSet, PointList &randPoints, NT rad_min, std::vector<NT> &ratios,
                       const NT &lb, const NT &ub, NT &alpha, const int &nu){

    int n = (*randPoints.begin()).dimension(), iter = 1;
    bool too_few;
    NT radmax = 0.0, rad, pnorm, ratio;

    for (typename PointList::iterator rpit = randPoints.begin();  rpit!=randPoints.end(); ++rpit) {
        pnorm = (*rpit).squared_length();
        if (pnorm > radmax) radmax = pnorm;
    }
    ball Biter;
    radmax=std::sqrt(radmax);
    NT rad0 = rad_min, rad_m = radmax;

    while (iter <= MAX_ITER) {
        rad = 0.5 * (rad_min + radmax);
        Biter = ball(Point(n), rad * rad);
        too_few = false;
        //std::cout << "rmax = " << radmax << " rad1 = " << rad_min << " rad_med = " << rad << std::endl;

        if (test_check_convergence<Point>(Biter, randPoints, lb, ub, too_few, ratio, nu, alpha, false, false)) {
            BallSet.push_back(Biter);
            ratios.push_back(ratio);
            return true;
        }

        if (too_few) {
            rad_min = rad;
        } else {
            radmax = rad;
        }

        if(radmax-rad_min < TOL) {
            rad_min = rad0;
            radmax = rad_m;
            iter++;
        }

    }
    return false;
}


template <typename PolyBall, typename RNGType,class ball, typename Polytope, typename Parameters, typename NT, typename Point, typename VT>
bool test_get_sequence_of_polyballs(Polytope &P, std::vector<ball> &BallSet, std::vector<NT> &ratios, const int &Ntot, const int &nu,
                               const NT &lb, const NT &ub, NT radius, NT &alpha, const Parameters &var, NT &rmax, Point &vec,
                                    Point &q, VT &lamdas, VT &Av) {

    //typedef typename Polytope::PolytopePoint Point;
    typedef typename Polytope::MT MT;
    //typedef typename Polytope::VT VT;

    bool fail;
    int n = P.dimension();
    NT ratio, ratio0, diameter = var.diameter;
    std::list<Point> randPoints;
    ball B0;
    //q.set_to_zero();
    PolyBall zb_it;

    if( !test_get_first_ball<RNGType>(P, B0, ratio, radius, lb, ub, alpha, rmax, vec, var) ) {
        return false;
    }
    //std::cout<<"first ball computed"<<std::endl;

    //VT lamdas, Av;
    //lamdas.setZero(P.num_of_hyperplanes());
    //Av.setZero(P.num_of_hyperplanes());

    q.set_to_zero();
    ratio0 = ratio;
    test_rand_point_generator(P, q, Ntot, var.walk_steps, randPoints, lamdas, Av, vec, diameter, var);
    //std::cout<<Ntot<<" points sampled from P"<<std::endl;

    if (test_check_convergence<Point>(B0, randPoints, lb, ub, fail, ratio, nu, alpha, false, true)) {
        ratios.push_back(ratio);
        BallSet.push_back(B0);
        ratios.push_back(ratio0);
        return true;
    }
    //std::cout<<"not the last ball, ratio = "<<ratio<<std::endl;
    if ( !test_get_next_zonoball<Point>(BallSet, randPoints, B0.radius(), ratios, lb, ub, alpha, nu) ){
        return false;
    }
    //std::cout<<"number of balls = "<<BallSet.size()+1<<std::endl;

    while (true) {
        zb_it = PolyBall(P, BallSet[BallSet.size()-1]);
        q.set_to_zero();
        randPoints.clear();
        zb_it.comp_diam(diameter, 0.0);
        test_rand_point_generator(zb_it, q, Ntot, var.walk_steps, randPoints, lamdas, Av, vec, diameter, var);
        //std::cout<<"N points sampled from BP"<<std::endl;

        if (test_check_convergence<Point>(B0, randPoints, lb, ub, fail, ratio, nu, alpha, false, true)) {
            ratios.push_back(ratio);
            BallSet.push_back(B0);
            ratios.push_back(ratio0);
            return true;
        }
        if ( !test_get_next_zonoball<Point>(BallSet, randPoints, B0.radius(), ratios, lb, ub, alpha, nu) ) {
            return false;
        }
        //std::cout<<"number of balls = "<<BallSet.size()+1<<std::endl;
    }
}


#endif


