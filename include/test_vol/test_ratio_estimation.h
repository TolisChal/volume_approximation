// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis


#ifndef TEST_RATIO_ESTIMATION_H
#define TEST_RATIO_ESTIMATION_H

#define MAX_ITER_ESTI 10000000

template <typename NT>
bool test_check_max_error(const NT &a, const NT &b, const NT &error) {

    if((b-a)/a<error/2.0) {
        return true;
    }
    return false;

}


template <typename RNGType, typename Point, typename PolyBall1, typename PolyBall2, typename NT, typename Parameters, typename VT>
NT test_esti_ratio(PolyBall1 &Pb1, PolyBall2 &Pb2, const NT &ratio, const NT &error, const int &W,
        const int &Ntot, const NT &diameter, Point &p, Point &vec, Parameters &var, VT &lamdas, VT &Av,
        bool isball = false, NT radius = 0.0) {

    int n = var.n, min_index = W-1, max_index = W-1, index = 0, iter = 1;
    bool verbose = var.verbose;
    NT min_val = std::numeric_limits<NT>::lowest(), max_val = std::numeric_limits<NT>::max(), val, lambda;
    size_t totCount = Ntot, countIn = Ntot * ratio;
    std::vector<NT> last_W(W,0);
    //typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    //VT lamdas, Av;// vec;
    //lamdas.setZero(Pb1.num_of_hyperplanes());
    //Av.setZero(Pb1.num_of_hyperplanes());
    //Point vec(n);

    //std::list<Point> randPoints;
    typename std::vector<NT>::iterator minmaxIt;
    //typename std::list<Point>::iterator rpit;
    //Point p(n);
    p.set_to_zero();
    if(!var.ball_walk && !isball){
        test_uniform_first_point(Pb1,p,var.walk_steps,lamdas,Av,vec,lambda, diameter, var);
    }

    while(iter <= MAX_ITER_ESTI){
        iter++;

        if (isball) {
            test_get_point_in_Dsphere<RNGType>(n, radius, p, var);
        } else {
            test_uniform_next_point(Pb1, p, var.walk_steps, lamdas, Av, vec, lambda, diameter, var);
        }
        if(Pb2.is_in(p)==-1) countIn++;

        totCount++;
        val = NT(countIn) / NT(totCount);
        last_W[index] = val;

        if(val<=min_val){
            min_val = val;
            min_index = index;
        }else if(min_index==index){
            minmaxIt = std::min_element(last_W.begin(), last_W.end());
            min_val = *minmaxIt;
            min_index = std::distance(last_W.begin(), minmaxIt);
        }

        if(val>=max_val){
            max_val = val;
            max_index = index;
        }else if(max_index==index){
            minmaxIt = std::max_element(last_W.begin(), last_W.end());
            max_val = *minmaxIt;
            max_index = std::distance(last_W.begin(), minmaxIt);
        }

        if( (max_val-min_val)/max_val<=error/2.0 ){
            if (verbose) std::cout << "final rejection ratio = " << val << " | total points = " << totCount << std::endl;
            return val;
        }

        index = index%W+1;
        if(index==W) index=0;

    }
    if (verbose) std::cout << "final rejection ratio = " << val << " | total points = " << totCount << std::endl;
    return val;
}


template <typename RNGType, typename Point, typename PolyBall1, typename PolyBall2, typename NT, typename Parameters, typename VT>
NT test_esti_ratio_interval(PolyBall1 &Pb1, PolyBall2 &Pb2, const NT &ratio, const NT &error, const int &W,
        const int &Ntot, const NT &prob, Point &q, Point &p, const NT &diameter, Parameters &var, VT &lamdas, VT &Av,
        bool isball = false, NT radius = 0.0) {

    int n = var.n, index = 0, iter = 1;
    bool verbose = var.verbose;
    std::vector<NT> last_W(W);
    //typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    //VT lamdas, Av;//, vec;
    //Av.setZero(Pb1.num_of_hyperplanes());
    //lamdas.setZero(Pb1.num_of_hyperplanes());
    //vec.setZero(n);
    NT val, sum_sq=0.0, sum=0.0, lambda;
    size_t totCount = Ntot, countIn = Ntot * ratio;
    //std::cout<<"countIn = "<<countIn<<", totCount = "<<totCount<<std::endl;

    p.set_to_zero();//(n);

    if(!isball) test_uniform_first_point(Pb1, p, 1, lamdas, Av, q, lambda, diameter, var);
    for (int i = 0; i < W; ++i) {

        if (isball) {
            test_get_point_in_Dsphere<RNGType>(n, radius, p, var);
        } else {
            test_uniform_next_point(Pb1, p, var.walk_steps, lamdas, Av, q, lambda, diameter, var);
        }
        if (Pb2.is_in(p) == -1) countIn = countIn + 1;

        totCount = totCount + 1;
        val = NT(countIn) / NT(totCount);
        sum += val;
        sum_sq += val * val;
        last_W[index] = val;
        index = index % W + 1;

        if (index == W) index = 0;
    }

    boost::math::normal dist(0.0, 1.0);
    NT zp = boost::math::quantile(boost::math::complement(dist, (1.0 - prob)/2.0)), m=sum/NT(W), s;

    while(iter <= MAX_ITER_ESTI) {
        iter++;

        if (isball) {
            test_get_point_in_Dsphere<RNGType>(n, radius, p, var);
        } else {
            test_uniform_next_point(Pb1, p, var.walk_steps, lamdas, Av, q, lambda, diameter, var);
        }
        if (Pb2.is_in(p) == -1) countIn = countIn + 1;

        totCount = totCount + 1;
        val = NT(countIn) / NT(totCount);

        m = (m - last_W[index] / NT(W)) + val / NT(W);
        sum_sq = (sum_sq - last_W[index] * last_W[index]) + val * val;
        sum = (sum - last_W[index]) + val;
        s = std::sqrt((sum_sq + NT(W) * m * m - 2.0 * m * sum) / NT(W));
        last_W[index] = val;
        index = index % W + 1;

        if (index == W) index = 0;
        if (test_check_max_error(val - zp * s, val + zp * s, error)) {
            if (verbose) std::cout << "final rejection ratio = " << val << " | total points = " << totCount << std::endl;
            return val;
        }

    }
    return val;

}

#endif

