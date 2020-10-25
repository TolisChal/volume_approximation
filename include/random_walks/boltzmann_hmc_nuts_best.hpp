// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_BOLTZMANN_HMC_FULL_OPT_WALK_HPP
#define VOLESTI_BOLTZMANN_HMC_FULL_OPT_WALK_HPP

#include "spectrahedron.h"
#include "generators/boost_random_number_generator.hpp"
#include "../sampling/sphere.hpp"
#include "root_finders/low_degree_polynomials.hpp"


template<typename VT, typename NT>
bool stopping_criterion_is_true(VT const& c, VT const& pi, VT const& vi, NT const& lambda,
                                VT const& p0, VT const& v0, NT const& T,  NT const& cv0, 
                                NT const& p0v0, NT const& cp0, std::vector<double> x) 
{ 
    NT a = cv0 / (2*T), b = -vi.dot(v0), c1 = p0v0 - pi.dot(v0);
    NT D = b*b - 4*a*c1, r1, r2;
    int nroots_sq = 0;
    std::vector<NT> roots;

    if (D < 0.0 && a > 0.0)  return false;
    if (D > 0.0) {
        r1 = (-b - std::sqrt(D) ) / (2 * a);
        r2 = (-b + std::sqrt(D) ) / (2 * a);
        nroots_sq = 2;
    }

    NT aa = 1.0/(2.0 * T * T), vi_norm_sq = vi.norm(), y;
    vi_norm_sq = vi_norm_sq * vi_norm_sq;
    NT bb = ( (-1.5 * c.dot(vi)) / (2*T) ) / aa, cc = ( (cp0 - c.dot(pi)) / T + vi_norm_sq ) / aa,
       dd = (pi.dot(vi) - p0.dot(vi)) / aa;
    
    int nroots = SolveP3(x, bb, cc, dd), total_roots = 0;

    if (nroots == 3) {
        
        if (nroots_sq == 2) {
            x.push_back(r1);
            x.push_back(r2);
            std::sort (x.begin(), x.end());
            total_roots = 5;
        } else {
            std::sort (x.begin(), x.end());
            total_roots = 3;
            for (int i = 0; i < 3; i++){
                if (x[i] > 0.0 && x[i] < lambda){
                    return true;
                }
            }
            y = lambda / 2.0;
            if ( (y*y*y + bb * y * y + cc * y + dd) < 0.0 ) {
                return true;
            } else {
                return false;
            }
        }
    } else if (nroots == 1) {
        if (nroots_sq == 2) {
            x[1] = r1;
            x[2] = r2;
            std::sort (x.begin(), x.end());
            total_roots = 3;
        } else {
            if (x[0] < 0.0 || x[0] > lambda) {
                y = lambda / 2.0;
                if ( (y*y*y + bb * y * y + cc * y + dd) < 0.0 ) {
                    return true;
                } else {
                    return false;
                }
            } else {
                return true;
            }
        }
    } else {
        if (nroots_sq == 2) {
            x[2] = r1;
            x.push_back(r2);
            std::sort (x.begin(), x.end());
            total_roots = 4;
        } else {
            NT g1 = x[0], g2 = x[1];
            if (g2 < g1) {
                NT temp_g = g1;
                g1 = g2;
                g2 = g1;
            }
            if (g2 < 0.0 || g1 > lambda) {
                y = lambda / 2.0;
                if ( (y*y*y + bb * y * y + cc * y + dd) < 0.0 ) {
                    return true;
                } else {
                    return false;
                }
            }
            total_roots = 2;
        }
        //std::cout<<"two roots"<<std::endl;
        return false;
    }

    if (x[0] > lambda || x[nroots - 1] < 0.0) {
        y = lambda / 2.0;
        if ( ( (y*y*y + bb * y * y + cc * y + dd) * (a * y*y + b * y + c1 ) ) < 0.0 ) {
            return true;
        } else {
            return false;
        }
    }

    
    int counter = 0;
    bool first_positive = false;
    while (counter < total_roots) {
        
        if (x[counter] > 0.0 && !first_positive) {
            first_positive = true;
            y = x[counter] / 2.0;
            if ( ( (y*y*y + bb * y * y + cc * y + dd) * (a * y*y + b * y + c1 ) ) < 0.0 ) {
                return true;
            }
        }

        if (counter != total_roots -1) {
            if (x[counter + 1] > lambda) {
                y = (x[counter] + lambda) / 2.0;
                if ( ( (y*y*y + bb * y * y + cc * y + dd) * (a * y*y + b * y + c1 ) ) < 0.0 ) {
                    return true;
                } else {
                    return false;
                }
            }
            y = (x[counter] + x[counter + 1]) / 2.0;
            if ( ( (y*y*y + bb * y * y + cc * y + dd) * (a * y*y + b * y + c1 ) ) < 0.0 ) {
                return true;
            }
        } else {
            y = (x[counter] + lambda) / 2.0;
            if ( ( (y*y*y + bb * y * y + cc * y + dd) * (a * y*y + b * y + c1 ) ) < 0.0 ) {
                return true;
            }
        }
        counter++;
    }

    return false;
}

/// The Hamiltonian Monte Carlo random walk, to sample from the Boltzmann distribution, i.e. e^(-c*x/T).
struct BoltzmannHMCFullOptWalk {
public:

    struct parameters {};
    parameters param;

    /// The implementation of the walk
    /// Currently implemented only for spectrahedra
    /// with template specialization
    ///@tparam ConvexBody a convex body
    ///@tparam RandomNumberGenerator
    template <typename ConvexBody, typename RandomNumberGenerator>
    struct Walk {
    };



    /// The implementation of the walk for spectrahedra
    template <typename NT, typename RandomNumberGenerator>
    struct Walk<Spectrahedron<NT, Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>, Eigen::Matrix<NT,Eigen::Dynamic,1> >, RandomNumberGenerator> {

        /// The type of the spectrahedron
        typedef Spectrahedron<NT, Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic>, Eigen::Matrix<NT,Eigen::Dynamic,1>> SPECTRAHEDRON;
        /// Type for internal structure of class Spectrahedron
        typedef typename SPECTRAHEDRON::PrecomputedValues PrecomputedValues;

        /// The matrix/vector types we use
        typedef typename SPECTRAHEDRON::MATRIX_TYPE MT;
        typedef typename SPECTRAHEDRON::VECTOR_TYPE VT;

        /// A struct containing the parameters for the random walk
        struct Settings {
            /// The number of points to "burn", before keeping the following as a sample
            int walk_length;
            /// For generating random numbers
            RandomNumberGenerator randomNumberGenerator;
            /// The c in the distribution
            VT c;
            /// The T in the distribution
            NT temperature;
            /// The diameter of the spectrahedron
            NT diameter;
            /// Set the number of allowed reflections at each step: #reflections < reflectionsBound * dimension
            unsigned int reflectionsBound;
            /// When determining we can move d long till we reach the boundary, we walk d*dl, for numerical stability
            NT dl;

            /// Constructs an object of Settings
            /// \param[in] walkLength The number of points to "burn", before keeping the following as a sample
            /// \param[in] rng For generating random numbers
            /// \param[in] c The c in the distribution
            /// \param[in] temperature The T in the distribution
            /// \param[in] diameter The diameter of the spectrahedron
            /// \param[in] reflectionsBound at each iteration allow reflectionsBound*dimension reflections at most
            /// \param[in] dl approach the boundary with a factor of dl, for numerical stability
            /// \return An instance of this struct
            template <typename Point>
            Settings(const int walkLength, const RandomNumberGenerator &randomNumberGenerator, const Point &c, const NT temperature, const NT diameter,
                     unsigned int reflectionsBound = 10, NT dl = 0.995) : walk_length(walkLength), randomNumberGenerator(randomNumberGenerator),
                                                                          c(c.getCoefficients()),
                                                                          temperature(temperature),
                                                                          diameter(diameter),
                                                                          reflectionsBound(reflectionsBound), dl(dl) {}

            Settings() {}
        };

        /// The parameters of the random walk
        Settings settings;
        std::vector<double> x;

        Walk() {}

        /// Constructor
        /// \param[in] settings The settings of the random walk
        Walk(Settings &settings) : settings(settings) {
            x = std::vector<double>(3);
        }

        /// Change the settings
        /// \param[in] settings The settings of the random walk
        void setSettings(Settings &settings) {
            this->settings = settings;
        }


        /// Samples random points from the spectrahedron from the Boltzmann distribution
        /// \param[in] spectrahedron A spectrahedron
        /// \param[in] interiorPoint A point in the interior of the spectrahedron
        /// \param[in] pointsNum The number of points to sample
        /// \param[out] points The list of the sampled points
        /// \tparam Point class Point with NT and VT as declared above in this class
        template <typename Point>
        void apply(SPECTRAHEDRON &spectrahedron, Point const & interiorPoint, const unsigned int pointsNum,
                    std::list<Point> &points) {
            // store intermediate results between successive calls of methods
            // of the class spectrahedron, to avoid repeating computations
            PrecomputedValues precomputedValues;
            VT p = interiorPoint.getCoefficients();

            // sample #pointsNum points
            for (unsigned int i = 1; i <= pointsNum; ++i) {
                // burn #walk_length points to get one sample
                for (unsigned int j = 0; j < settings.walk_length; ++j) {
                    getNextPoint<Point>(spectrahedron, p, precomputedValues);
                }

                // add the sample in the return list
                points.push_back(Point(p));
            }

            // the data in preComputedValues may be out of date in the next call
            precomputedValues.resetFlags();
        }

        /// Samples random points from the spectrahedron from the Boltzmann distribution
        /// \param[in] spectrahedron A spectrahedron
        /// \param[in] interiorPoint A point in the interior of the spectrahedron
        /// \param[in] pointsNum The number of points to sample
        /// \param[out] points The list of the sampled points
        /// \param[in, out] precomputedValues transfer data between sucessive calls
        /// \tparam Point class Point with NT and VT as declared above in this class
        template <typename Point>
        void apply(SPECTRAHEDRON &spectrahedron, Point const & interiorPoint, const unsigned int pointsNum,
                    std::list<Point> &points, PrecomputedValues &precomputedValues) {
            // store intermediate results between successive calls of methods
            // of the class spectrahedron, to avoid repeating computations

            VT p = interiorPoint.getCoefficients();

            // sample #pointsNum points
            for (unsigned int i = 1; i <= pointsNum; ++i) {
                // burn #walk_length points to get one sample
                for (unsigned int j = 0; j < settings.walk_length; ++j) {
                    getNextPoint<Point>(spectrahedron, p, precomputedValues);
                }

                // add the sample in the return list
                points.push_back(Point(p));
            }

            // the data in preComputedValues may be out of date in the next call
            precomputedValues.resetFlags();
            precomputedValues.computed_C = true;
        }


        /// A single step of the HMC random walk: choose a direction and walk on the trajectory for a random distance.
        /// If it hits the boundary, the trajectory is reflected. If #reflections < reflectionsBound * dimension, it returns the same point
        /// \param[in] spectrahedron A spectrahedron
        /// \param[in, out] p An interior point, and the next point in the random walk
        /// \param[in, out] precomputedValues Data for the methods of the class Spectrahedron
        /// \tparam Point
        template <typename Point>
        void getNextPoint(SPECTRAHEDRON &spectrahedron, VT &p, PrecomputedValues &precomputedValues) {

            // initialize
            RandomNumberGenerator &rng = settings.randomNumberGenerator;
            boost::random::uniform_real_distribution<> urdist(0, 1);
            const NT dl = settings.dl;
            unsigned int n = spectrahedron.dimension();
            int reflectionsNum = 0;
            int reflectionsNumBound = settings.reflectionsBound * n*2;
            VT previousPoint;
            VT p0 = p;

            // choose a distance to walk
            NT T = 0.0, norm_c = 1.0, 
               lambda_with_min_val, best_val = settings.c.dot(p);

            // The trajectory will be of the form a*t^2 + vt + p
            // where a = -c / 2*temperature
            // and at each reflection, v and p will change

            // crate vector a
            VT a = -settings.c / (2 * settings.temperature), p_best = p, p_temp(n);

            // The vector v will be a random a direction
            VT v = GetSpericalGaussian<Point>::apply(n, rng).getCoefficients();
            //VT v = GetDirection<Point>::apply(n, rng).getCoefficients();
            VT v0 = v;
            NT cv0 = settings.c.dot(v0), cp0 = settings.c.dot(p0), p0v0 = p0.dot(v0);

            // Begin a step of the random walk
            // Also, count the reflections and respect the bound
            //std::cout<<"norm-c = "<<norm_c<<std::endl;
            while (reflectionsNum < reflectionsNumBound) {

                // we are at point p and the trajectory a*t^2 + vt + p
                // find how long we can walk till we hit the boundary
                NT lambda = spectrahedron.positiveIntersection(a, v, p, precomputedValues);
                lambda_with_min_val = (settings.c.dot(v)) /((norm_c * norm_c) / settings.temperature);
                //std::cout<<"lambda = "<<lambda<<", lambda_with_min_val = "<<lambda_with_min_val<<std::endl;
                p_temp = (lambda_with_min_val * lambda_with_min_val) * a + lambda_with_min_val * v + p;
                //std::cout<<"lambda_max_obj = "<<settings.c.dot(p_temp)<<", obj_current = "<<settings.c.dot(p)<<", best = "<<best_val<<std::endl;

                // We just solved a quadratic polynomial eigenvalue system At^2 + Bt + C,
                // where A = lmi(a) - A0, B = lmi(v) - A0 and C = lmi(p)
                // and also did a linearization creating two matrices X, Y.
                // For the subsequent calls, we don't have to compute these matrices from scratch,
                // but can efficiently update them.
                // A remains the same
                // C := A*lambda^2 + B*lambda + C
                // X, Y will be updated in class Spectrahedron
                // Set the flags
                precomputedValues.computed_A = true;
                precomputedValues.computed_C = true;
                precomputedValues.computed_XY = true;

                // if we can walk the remaining distance without reaching he boundary
                if (stopping_criterion_is_true(settings.c, p, v, lambda, p0, v0, T, cv0, 
                                               p0v0, cp0, x)) {
                    // set the new point p:= (T^2)*a + T*V + p
                    //p += (T * T) * a + T * v;
                    if (lambda_with_min_val < 0.0) {
                        lambda_with_min_val = lambda * 0.95;
                        p_temp = (lambda_with_min_val * lambda_with_min_val) * a + lambda_with_min_val * v + p;
                        if (settings.c.dot(p_temp) < best_val) {
                            best_val = settings.c.dot(p_temp);
                            p_best = p_temp;
                        }
                    } else if (lambda < lambda_with_min_val){
                        lambda_with_min_val = lambda * 0.05;
                        p_temp = (lambda_with_min_val * lambda_with_min_val) * a + lambda_with_min_val * v + p;
                        if (settings.c.dot(p_temp) < best_val) {
                            best_val = settings.c.dot(p_temp);
                            p_best = p_temp;
                        }
                    } else if (lambda > lambda_with_min_val) {
                        lambda_with_min_val = lambda * 0.05;
                        p_temp = (lambda_with_min_val * lambda_with_min_val) * a + lambda_with_min_val * v + p;
                        if (settings.c.dot(p_temp) < best_val) {
                            best_val = settings.c.dot(p_temp);
                            p_best = p_temp;
                        }
                        lambda_with_min_val = lambda * 0.95;
                        p_temp = (lambda_with_min_val * lambda_with_min_val) * a + lambda_with_min_val * v + p;
                        if (settings.c.dot(p_temp) < best_val) {
                            best_val = settings.c.dot(p_temp);
                            p_best = p_temp;
                        }
                    }
                    //std::cout<<"STOP REFLECVTIONS!!\n\n";

                    // update matrix C
                    //precomputedValues.C += (T * T) * precomputedValues.A + T * precomputedValues.B;
                    break;
                }

                if (lambda_with_min_val < 0.0) {
                        lambda_with_min_val = lambda * 0.95;
                        p_temp = (lambda_with_min_val * lambda_with_min_val) * a + lambda_with_min_val * v + p;
                        if (settings.c.dot(p_temp) < best_val) {
                            best_val = settings.c.dot(p_temp);
                            p_best = p_temp;
                        }
                } else if (lambda < lambda_with_min_val){
                        lambda_with_min_val = lambda * 0.05;
                        p_temp = (lambda_with_min_val * lambda_with_min_val) * a + lambda_with_min_val * v + p;
                        if (settings.c.dot(p_temp) < best_val) {
                            best_val = settings.c.dot(p_temp);
                            p_best = p_temp;
                        }
                } else if (lambda > lambda_with_min_val) {
                        lambda_with_min_val = lambda * 0.05;
                        p_temp = (lambda_with_min_val * lambda_with_min_val) * a + lambda_with_min_val * v + p;
                        if (settings.c.dot(p_temp) < best_val) {
                            best_val = settings.c.dot(p_temp);
                            p_best = p_temp;
                        }
                        lambda_with_min_val = lambda * 0.95;
                        p_temp = (lambda_with_min_val * lambda_with_min_val) * a + lambda_with_min_val * v + p;
                        if (settings.c.dot(p_temp) < best_val) {
                            best_val = settings.c.dot(p_temp);
                            p_best = p_temp;
                        }
                    }

                // we hit the boundary and still have to walk
                // don't go all the way to the boundary, for numerical stability
                lambda *= dl;

                // save current and set new point
                previousPoint = p;
                p += (lambda * lambda) * a + lambda * v;

                // update remaining distance we must walk
                //T -= lambda;

                // update matrix C
                precomputedValues.C += (lambda * lambda) * precomputedValues.A + lambda * precomputedValues.B;

                // Set v to have the direction of the trajectory at t = lambda
                // i.e. the gradient of at^2 + vt + p, for t = lambda
                v += (lambda * 2) * a;

                // compute reflected direction
                VT reflectedTrajectory;
                spectrahedron.computeReflection(p, v, reflectedTrajectory, precomputedValues);
                v = reflectedTrajectory;
                reflectionsNum++;
            }

            //std::cout<<"set p!"<<std::endl;
            p = p_best;
            //std::cout<<"is in = "<< spectrahedron.is_in(p) <<std::endl;
            // if the #reflections exceeded the limit, don't move
            //if (reflectionsNum == reflectionsNumBound/2) std::cout<<"limit reached!"<<std::endl;
                //p = p0;
        }

        /// Sets the temperature in the distribution
        /// \param[in] temperature New value of temperature
        void setTemperature(NT temperature) {
            settings.temperature = temperature;
        }
    };

};

#endif //VOLESTI_BOLTZMANN_HMC_WALK_HPP
