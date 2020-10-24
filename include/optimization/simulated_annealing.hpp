// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_SIMULATED_ANNEALING_HPP
#define VOLESTI_SIMULATED_ANNEALING_HPP

#include <format.hpp>

#include "optimization/sliding_window.hpp"
#include "random_walks/boltzmann_hmc_walk.hpp"
#include "random_walks/uniform_billiard_sdp.hpp"


/// Estimates the diameter of the spectrahedron. It samples points uniformly with coordinate directions
/// hit and run, and returns the maximum distance between them.
/// \tparam Point
/// \param[in] numPoints The number of points to sample for the estimation
/// \return An estimation of the diameter of the spectrahedron
template<typename MT, typename WalkType, typename NT, class RNGType, class Spectra, class Point>
NT estimateDiameterBilliard(Spectra &spectrahedron, int const numPoints, Point const & interiorPoint) 
{
    typedef typename WalkType::template Walk<Spectra, RNGType > BiW;

    RNGType rng(spectrahedron.dimension());
    typename BiW::Settings biw_settings(1, rng);
    BiW biwRandomWalk(biw_settings);

    // this data structure help us move computations between function calls
    PrecomputedValues<MT> biwPrecomputedValues;
    biwPrecomputedValues.set_mat_size(spectrahedron.getLMI().sizeOfMatrices());
    std::list<Point> randPoints;

    std::cout<<"Hi1"<<std::endl;
    biwRandomWalk.apply(spectrahedron, interiorPoint, numPoints, randPoints, biwPrecomputedValues);
    std::cout<<"Hi2"<<std::endl;
        
    // find maximum distance among points;
    NT maxDistance = 0;
    typename std::list<Point>::iterator itInner, itOuter = randPoints.begin();

    for (; itOuter!=randPoints.end() ; ++itOuter) {
        for (itInner=itOuter ; itInner!=randPoints.end() ; ++itInner) {
            NT current = itOuter->distance(*itInner);
            if (current > maxDistance) {
                maxDistance = current;
            }
        }
    }

    return maxDistance;
}


/// A magic number!
/// when estimating the diameter of the spectrahedron,
/// sample 20 + sqrt(dimension) points to estimate it
#define CONSTANT_1 20

/// Holds parameters of the algorithm
/// \tparam Point Class point
template<class Point>
struct SimulatedAnnealingSettings {
    /// The numeric type
    typedef typename Point::FT NT;

    /// Desired accuracy (relative error)
    NT error;
    /// The walk length of the HMC random walk
    int walkLength;
    /// A bound to the number of steps; if negative it is unbounded
    int maxNumSteps;

    /// Starting from an initial temperature, at each step it will decrease by a factor of
    /// \[ 1 - 1 / (dimension^k) \]. Default is 0.5
    NT k;

    SimulatedAnnealingSettings(NT const error, int const walkLength = 1, int const maxNumSteps = -1, NT const k = 0.5) : error(error),
        walkLength(walkLength), maxNumSteps(maxNumSteps), k(k) {}
};


/// Simulated Annealing algorithm for a semidefinite program
/// Minimize \[ c^T x \], s.t. LMI(x) <= 0
/// \param[in] spectrahedron A spectrahedron described by a linear matrix inequality
/// \param[in] objectiveFunction The function we minimize
/// \param[in] settings Parameters of the algorithm
/// \param[in] interiorPoint An initial feasible solution to start the algorithm
/// \param[out] solution The vector minimizing the objective function
/// \param[in] verbose True to print messages. Default is false
/// \return The best approximation to the optimal solution
template <typename _Spectrahedron, typename Point, typename _Settings>
double solve_sdp(_Spectrahedron & spectrahedron, Point const & objectiveFunction, _Settings const & settings,
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
    NT diameter = spectrahedron.estimateDiameter(CONSTANT_1 + std::sqrt(spectrahedron.dimension()), interiorPoint);

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
    typename HMC::Settings hmc_settings = typename HMC::Settings(settings.walkLength, rng, objectiveFunction, temperature, diameter);
    HMC hmcRandomWalk = HMC(hmc_settings);

    // this data structure help us move computations between function calls
    PrecomputedValues<MT> hmcPrecomputedValues;
    hmcPrecomputedValues.set_mat_size(spectrahedron.getLMI().sizeOfMatrices());

    NT previous_min = objectiveFunction.dot(solution);

    /******** solve *********/
    // if settings.maxNumSteps is negative there is no
    // bound to the number of steps - stop
    // when desired relative error is achieved
    while (stepsCount < settings.maxNumSteps || settings.maxNumSteps < 0) {

        // sample one point with current temperature
        std::list<Point> randPoints;

        // get a sample under the Boltzmann distribution
        // using the HMC random walk
        while (1) {
            hmcRandomWalk.apply(spectrahedron, solution, settings.walkLength, randPoints, hmcPrecomputedValues);

            // if the sampled point is not inside the spectrahedron (error in boundary oracle),
            // get a new one
            if (spectrahedron.isExterior(hmcPrecomputedValues.C)) {
                if (verbose) std::cout << "Sampled point outside the spectrahedron.\n";
                randPoints.clear();
                hmcPrecomputedValues.resetFlags();
            }
            else {
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

    // return the minimum w.r.t. the original objective function
    return currentMin*objectiveFunctionNorm;
}




template <typename WalkType, typename _Spectrahedron, typename Point, typename _Settings>
double solve_sdp_with_optimal(_Spectrahedron & spectrahedron, Point const & objectiveFunction, _Settings const & settings,
         Point const & interiorPoint, Point& solution, double optimal_val, bool verbose = false) {

    // fetch the data types we will use
    typedef  typename _Spectrahedron::NUMERIC_TYPE NT;
    typedef  typename _Spectrahedron::MATRIX_TYPE MT;
    typedef  typename _Spectrahedron::VECTOR_TYPE VT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    //typedef BoltzmannHMCWalk::Walk<_Spectrahedron, RNGType > HMC;
    typedef typename WalkType::template Walk<_Spectrahedron, RNGType > HMC;

    // the algorithm requires the objective function to be normalized
    // we will need to remember the norm
    VT _objectiveFunctionNormed = objectiveFunction.getCoefficients();
    NT objectiveFunctionNorm = _objectiveFunctionNormed.norm();
    _objectiveFunctionNormed.normalize();
    Point objectiveFunctionNormed = Point(_objectiveFunctionNormed);

    // Estimate the diameter of the spectrahedron
    // needed for the random walk and for the simulated annealing algorithm
    std::cout << "diameter to compute"<<std::endl;
    //NT diameter = estimateDiameterBilliard<MT, BilliardWalkSDP, NT, RNGType>(spectrahedron, CONSTANT_1 + std::sqrt(spectrahedron.dimension()), interiorPoint);
    //std::cout << "diameter = "<<diameter<<std::endl;
    //std::cout << "diaminteriorPointeter = "<<interiorPoint.getCoefficients().transpose()<<std::endl;
    //NT diameter2 = spectrahedron.estimateDiameter(CONSTANT_1 + std::sqrt(spectrahedron.dimension()), interiorPoint);
    //std::cout << "diameter2 = "<<diameter2<<std::endl;
    RNGType rng2(spectrahedron.dimension());
    //std::cout << "Hi"<<std::endl;
    VT best_point(spectrahedron.dimension());
    NT diameter = spectrahedron.estimateDiameterRDHR(CONSTANT_1 + std::sqrt(spectrahedron.dimension()), 
                                                      interiorPoint, rng2, _objectiveFunctionNormed, best_point);
    //interiorPoint = Point(best_point);
    std::cout << "final diameter = "<<diameter<<std::endl;
    std::cout<<"is p Exterior = "<<spectrahedron.isExterior(best_point)<<std::endl;
    if (spectrahedron.isExterior(best_point)==1){
        std::cout<<"is p Exterior = "<<spectrahedron.isExterior(best_point)<<std::endl;
        exit(-1);
    }

    /******** initialization *********/
    solution = Point(best_point);
    // the minimum till last iteration
    NT currentMin = objectiveFunction.dot(solution);
    int stepsCount = 0;
    // initial temperature must be the diameter of the body
    NT temperature = diameter;
    // after each iteration, temperature = temperature * tempDecreaseFactor
    NT tempDecreaseFactor = 1.0 - static_cast<NT>(1.0 / std::pow(spectrahedron.dimension(), settings.k));

    // initialize random walk;
    RNGType rng(spectrahedron.dimension());
    typename HMC::Settings hmc_settings = typename HMC::Settings(settings.walkLength, rng, 
                                                                 objectiveFunctionNormed, 
                                                                 temperature, diameter);
    HMC hmcRandomWalk = HMC(hmc_settings);
    // this data structure help us move computations between function calls
    PrecomputedValues<MT> hmcPrecomputedValues;
    hmcPrecomputedValues.set_mat_size(spectrahedron.getLMI().sizeOfMatrices());

    NT previous_min = objectiveFunction.dot(solution);
    std::cout << "Step: " << stepsCount << ", Temperature: " << temperature << ", Min: " << currentMin
                      << ", Relative error: " << std::abs(currentMin - optimal_val) / std::abs(optimal_val)
                      << ", optimal_val: " << optimal_val<< ", settings.error: " << settings.error<< "\n";
    temperature = 600.0;
    Point temp_point(spectrahedron.dimension());
    bool improved;
    /******** solve *********/
    // if settings.maxNumSteps is negative there is no
    // bound to the number of steps - stop
    // when desired relative error is achieved
    while (std::abs(currentMin - optimal_val) / std::abs(optimal_val) > settings.error) {

        // sample one point with current temperature
        std::list<Point> randPoints;

        // get a sample under the Boltzmann distribution
        // using the HMC random walk
        while (1) {
            //std::cout << "hello sampling" << "\n";
            hmcRandomWalk.apply(spectrahedron, solution, settings.walkLength, randPoints, hmcPrecomputedValues);

            // if the sampled point is not inside the spectrahedron (error in boundary oracle),
            // get a new one
            if (spectrahedron.isExterior(hmcPrecomputedValues.C)) {
                if (verbose) std::cout << "Sampled point outside the spectrahedron.\n";
                randPoints.clear();
                hmcPrecomputedValues.resetFlags();
            }
            else {
                // update values;
                temp_point = randPoints.back();
                randPoints.clear();
                if (objectiveFunction.dot(temp_point) < currentMin) {
                    solution = temp_point; 
                    improved = true;
                } else {
                    hmcPrecomputedValues.resetFlags();
                    improved = false;
                }
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
                      << ", Relative error: " << std::abs(currentMin - optimal_val) / std::abs(optimal_val)
                      << ", optimal_val: " << optimal_val<< ", settings.error: " << settings.error<< "\n";

        // check if we reached desired accuracy
        //if (relError < settings.error*0.1)
        //    break;

        // decrease the temperature
        temperature *= tempDecreaseFactor;
        hmcRandomWalk.setTemperature(temperature);

    } /* while (stepsCount < settings.maxNumSteps || settings.maxNumSteps < 0) { */

    // return the minimum w.r.t. the original objective function
    return currentMin;
}



#endif //VOLESTI_SIMULATED_ANNEALING_HPP
