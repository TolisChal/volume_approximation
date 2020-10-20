// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_COMPUTE_INITIAL_POINT_HPP
#define VOLESTI_COMPUTE_INITIAL_POINT_HPP


#include "random_walks/boltzmann_hmc_walk.hpp"


/// A magic number!
/// when estimating the diameter of the spectrahedron,
/// sample 20 + sqrt(dimension) points to estimate it
//#define CONSTANT_1 20


template<typename NT, typename MT, typename Spectra, typename LMI, typename Point>
void get_inner_point(LMI lmi, Point &p, Point const& obj) {

    int m = lmi.sizeOfMatrices();
    typedef Eigen::Triplet<NT> T;
    std::vector<T> tripletList;
    std::cout << "hi1"<<std::endl;

    for (int j=0; j<m; j++) {
        tripletList.push_back(T(j, j, 1.0));
    }
    std::cout << "hi2"<<std::endl;
    MT A(m,m);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    std::cout << "hi2.5"<<std::endl;
    lmi.add_matrix(A);
    std::cout << "hi3"<<std::endl;
    int d = lmi.dimension();

    Spectra spectrahedro(lmi);
    SimulatedAnnealingSettings<Point> settings(0.01, 1, -1, 0.25);

    Point q(d), sol(d), objectii(d);
    objectii.set_coord(d-1, -1.0);

    Eigen::SelfAdjointEigenSolver<MT> solver;
    solver.compute(-lmi.getMatrices()[0], Eigen::EigenvaluesOnly);
    NT lambda = solver.eigenvalues().minCoeff();
    std::cout << "hi4"<<std::endl;
    lambda *= 2.0;
    q.set_coord(d-1, lambda);
    //std::cout << "is exterior = "<<spectrahedro.isExterior(q.getCoefficients())<<std::endl;
    //sol = q;

    solve_for_initial_point(spectrahedro, objectii, settings, q, sol, true);
    p = sol;

}


template <typename _Spectrahedron, typename Point, typename _Settings>
double solve_for_initial_point(_Spectrahedron & spectrahedron, Point const & objectiveFunction, _Settings const & settings,
         Point const & interiorPoint, Point& solution, bool verbose = false) {

    // fetch the data types we will use
    typedef  typename _Spectrahedron::NUMERIC_TYPE NT;
    typedef  typename _Spectrahedron::MATRIX_TYPE MT;
    typedef  typename _Spectrahedron::VECTOR_TYPE VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef BoltzmannHMCWalk::Walk<_Spectrahedron, RNGType > HMC;

    // the algorithm requires the objective function to be normalized
    // we will need to remember the norm
    VT _objectiveFunctionNormed = objectiveFunction.getCoefficients();
    NT objectiveFunctionNorm = _objectiveFunctionNormed.norm();
    _objectiveFunctionNormed.normalize();
    Point objectiveFunctionNormed = Point(_objectiveFunctionNormed);

    // Estimate the diameter of the spectrahedron
    // needed for the random walk and for the simulated annealing algorithm
    std::cout << "diameter to compute"<<std::endl;
    //PrecomputedValues<MT> hmcPrecomputedValues2;
    //hmcPrecomputedValues2.set_mat_size(spectrahedron.getLMI().sizeOfMatrices());
    //Point qq = interiorPoint;
    //std::list<Point> randPoints2;
    //hmcRandomWalk.apply(spectrahedron, qq, settings.walkLength, randPoints2, hmcPrecomputedValues2);
    NT diameter = spectrahedron.estimateDiameter(CONSTANT_1 + std::sqrt(spectrahedron.dimension()), interiorPoint);
    std::cout << "diameter = "<<diameter<<std::endl;

    /******** initialization *********/
    solution = interiorPoint;
    // the minimum till last iteration
    NT currentMin = objectiveFunction.dot(solution);
    int stepsCount = 0;
    // initial temperature must be the diameter of the body
    NT temperature = diameter;
    // after each iteration, temperature = temperature * tempDecreaseFactor
    NT tempDecreaseFactor = 1.0 - static_cast<NT>(1.0 / std::pow(spectrahedron.dimension(), settings.k));

    // initialize random walk;
    RNGType rng(spectrahedron.dimension());
    typename HMC::Settings hmc_settings = typename HMC::Settings(settings.walkLength, rng, objectiveFunctionNormed, temperature, diameter);
    HMC hmcRandomWalk = HMC(hmc_settings);

    // this data structure help us move computations between function calls
    PrecomputedValues<MT> hmcPrecomputedValues;
    hmcPrecomputedValues.set_mat_size(spectrahedron.getLMI().sizeOfMatrices());

    NT previous_min = objectiveFunction.dot(solution);
    int d = spectrahedron.dimension();

    /******** solve *********/
    // if settings.maxNumSteps is negative there is no
    // bound to the number of steps - stop
    // when desired relative error is achieved
    while (1) {

        // sample one point with current temperature
        std::list<Point> randPoints;

        // get a sample under the Boltzmann distribution
        // using the HMC random walk
        while (1) {
            //std::cout<<"point to sample"<<std::endl;
            hmcRandomWalk.apply(spectrahedron, solution, settings.walkLength, randPoints, hmcPrecomputedValues);
            //std::cout<<"point sampled"<<std::endl;

            // if the sampled point is not inside the spectrahedron (error in boundary oracle),
            // get a new one
            if (spectrahedron.isExterior(hmcPrecomputedValues.C)) {
                if (verbose) std::cout << "Sampled point outside the spectrahedron.\n";
                randPoints.clear();
                hmcPrecomputedValues.resetFlags();
            }
            else {
                //std::cout << "Sampled point inside the spectrahedron.\n";
                // update values;
                solution = randPoints.back();
                randPoints.clear();
                break;
            }
        }

        // update current value
        currentMin = objectiveFunction.dot(solution);
        ++stepsCount;

        // compute relative error
        NT relError = relativeError(previous_min, currentMin);
        previous_min = currentMin;

        if (verbose)
            std::cout << "Step: " << stepsCount << ", Temperature: " << temperature << ", Min: " << currentMin
                      << ", Relative error: " << relError << "\n";

        // check if we reached desired accuracy
        if (relError < settings.error)
            break;

        // decrease the temperature
        temperature *= tempDecreaseFactor;
        hmcRandomWalk.setTemperature(temperature);

    } /* while (stepsCount < settings.maxNumSteps || settings.maxNumSteps < 0) { */
    //interiorPoint = solution;

    // return the minimum w.r.t. the original objective function
    return currentMin*objectiveFunctionNorm;
}


#endif
