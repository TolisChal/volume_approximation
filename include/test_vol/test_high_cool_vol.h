// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2019 Vissarion Fisikopoulos
// Copyright (c) 2018-2019 Apostolos Chalkis

#ifndef TEST_HIGH_COOL_VOL_H
#define TEST_HIGH_COOL_VOL_H


template <typename NT>
NT log_gamma_function(NT x) {
    if (x <= NT(100)) return std::log(tgamma(x));
    return (std::log(x - NT(1)) + log_gamma_function(x - NT(1)));
}

template <typename Polytope, typename Point, typename UParameters, typename AParameters, typename NT>
NT test_cooling_balls_high(Polytope &P, const UParameters &var, const AParameters &var_ban, std::pair<Point,NT> &InnerBall, bool only_phases = false) {

    typedef Ball <Point> ball;
    typedef BallIntersectPolytope <Polytope, ball> PolyBall;
    typedef typename UParameters::RNGType RNGType;
    typedef typename Polytope::VT VT;
    typedef std::list <Point> PointList;

    int n = var.n, win_len = var_ban.win_len, N = var_ban.N, nu = var_ban.nu;
    bool verbose = var.verbose, round = var.round, window2 = var_ban.window2;
    NT lb = var_ban.lb, ub = var_ban.ub, prob = var_ban.p, rmax = var_ban.rmax, radius = InnerBall.second,
            round_value = 1.0, e = var.error, alpha = var_ban.alpha, diam = var.diameter;

    std::vector <ball> BallSet;
    std::vector <NT> ratios;
    Point p = InnerBall.first, vec(n);
    //P.normalize();

    /*if (round) {
#ifdef VOLESTI_DEBUG
        if(verbose) std::cout<<"\nRounding.."<<std::endl;
#endif
        double tstart1 = (double) clock() / (double) CLOCKS_PER_SEC;
        //std::pair <NT, NT> res_round = rounding_min_ellipsoid(P, InnerBall, var);
        double tstop1 = (double) clock() / (double) CLOCKS_PER_SEC;
#ifdef VOLESTI_DEBUG
        if(verbose) std::cout << "Rounding time = " << tstop1 - tstart1 << std::endl;
#endif
        //round_value = res_round.first;
        std::pair <Point, NT> res = P.ComputeInnerBall();
        c = res.first;
        radius = res.second;
        P.normalize();
        if (var.bill_walk) {
            P.comp_diam(var.diameter, radius);
        }
        if (var.ball_walk){
            var.delta = 4.0 * radius / NT(n);
        }
    }*/

    // Save the radius of the Chebychev ball
    //var.che_rad = radius;
    // Move the chebychev center to the origin and apply the same shifting to the polytope
    P.shift(p.getCoefficients());
    P.normalize();
    P.recompute_AA();
    p.set_to_zero();//Zero(n);

    VT lamdas, Av;//, vec;
    Av.setZero(P.num_of_hyperplanes());
    lamdas.setZero(P.num_of_hyperplanes());

    if ( !test_get_sequence_of_polyballs<PolyBall, RNGType>(P, BallSet, ratios, N * nu, nu, lb, ub, radius, alpha, var,
                                                            rmax, vec, p, lamdas, Av) ){
        return -1.0;
    }
    //diam = var.diameter;

    //NT vol = (std::pow(M_PI, n / 2.0) * (std::pow((*(BallSet.end() - 1)).radius(), n))) / (tgamma(n / 2.0 + 1));
    NT vol = (NT(n)/2.0 * std::log(M_PI)) + NT(n)*std::log((*(BallSet.end() - 1)).radius()) - log_gamma_function(NT(n) / 2.0 + 1);
    if (verbose) std::cout<<"rad of B_m = "<<(*(BallSet.end() - 1)).radius()<<", vol of B_m = "<<vol<<" "<<tgamma(n / 2.0 + 1)<<std::endl;
    //return vol;

    int mm = BallSet.size() + 1;
    nballs = NT(mm - 1);
    if (only_phases) {
        P.free_them_all();
        return vol;
    }
    prob = std::pow(prob, 1.0 / NT(mm));
    NT er0 = e / (2.0 * std::sqrt(NT(mm))), er1 = (e * std::sqrt(4.0 * NT(mm) - 1)) / (2.0 * std::sqrt(NT(mm)));

    vol += std::log(test_esti_ratio_interval<RNGType, Point>(*(BallSet.end() - 1), P, *(ratios.end() - 1), er0, win_len, 1200, prob,
                                                    vec, p, diam, var, lamdas, Av, true, (*(BallSet.end() - 1)).radius()));
    //std::cout<<"vol = "<<vol <<std::endl;

    PolyBall Pb;
    typename std::vector<ball>::iterator balliter = BallSet.begin();
    typename std::vector<NT>::iterator ratioiter = ratios.begin();

    er1 = er1 / std::sqrt(NT(mm) - 1.0);

    if (*ratioiter != 1) {
        vol += std::log(NT(1) / test_esti_ratio_interval<RNGType, Point>(P, *balliter, *ratioiter, er1,
                                                                         win_len, N * nu, prob, vec, p, diam, var,
                                                                         lamdas, Av));
    }
    for ( ; balliter < BallSet.end() - 1; ++balliter, ++ratioiter) {
        Pb = PolyBall(P, *balliter);
        Pb.comp_diam(diam, 0.0);
        vol +=  std::log(NT(1) / test_esti_ratio_interval<RNGType, Point>(Pb, *(balliter + 1), *(ratioiter + 1), er1,
                win_len, N * nu, prob, vec, p, diam, var, lamdas, Av));
        //std::cout<<"vol = "<<vol <<std::endl;
    }

    P.free_them_all();
    //std::cout<<"vol = "<<vol <<std::endl;
    return std::exp(vol) * round_value;

}

#endif
