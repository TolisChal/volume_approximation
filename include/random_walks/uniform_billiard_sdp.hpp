// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_UNIFORM_BILLIARD_SDP_HPP
#define VOLESTI_UNIFORM_BILLIARD_SDP_HPP

#include "spectrahedron.h"
#include "generators/boost_random_number_generator.hpp"
#include "../sampling/sphere.hpp"

template<typename VT, typename NT>
bool stopping_criterion_is_true(VT const& p, VT const& v, NT const& lambda, VT const& p0, VT const& v0) 
{ 
    NT t1 = (v.dot((p0 - p))) / (v.norm()*v.norm());
    NT t2 = (v0.dot((p0-p))) / (v.dot(v0));

    if (v.dot(v0) > 0.0) {
        if ((t1 < t2)) 
        {
            return false;
        } else if (t2 > lambda) {
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

/// The Hamiltonian Monte Carlo random walk, to sample from the Boltzmann distribution, i.e. e^(-c*x/T).
struct BilliardWalkSDP {
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
        typedef typename SPECTRAHEDRON::PrecomputedValues2 PrecomputedValues;

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
            //VT c;
            /// The T in the distribution
            //NT temperature;
            /// The diameter of the spectrahedron
            //NT diameter;
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
            //template <typename Point>
            Settings(const int walkLength, const RandomNumberGenerator &randomNumberGenerator,
                     unsigned int reflectionsBound = 10, NT dl = 0.995) : walk_length(walkLength), 
                                                                          randomNumberGenerator(randomNumberGenerator),
                                                                          reflectionsBound(reflectionsBound), dl(dl) {}

            Settings() {}
        };

        /// The parameters of the random walk
        Settings settings;

        Walk() {}

        /// Constructor
        /// \param[in] settings The settings of the random walk
        Walk(Settings &settings) : settings(settings) {}

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
            int reflectionsNumBound = settings.reflectionsBound * n;
            VT previousPoint;
            VT p0 = p, xi, vi;
            MT Ci;
            std::cout<<"Hi30"<<std::endl;
            std::cout<<"p0 = "<<p0.transpose()<<std::endl;
            std::cout<<"is p0 Exterior = "<<spectrahedron.isExterior(p0)<<std::endl;
            std::cout<<"Hi32"<<std::endl;

            // choose a distance to walk
            NT T = 0.0, ti, w2, w1, u, t_temp;

            // The trajectory will be of the form a*t^2 + vt + p
            // where a = -c / 2*temperature
            // and at each reflection, v and p will change

            // crate vector a
            //VT a = -settings.c / (2 * settings.temperature);

            // The vector v will be a random a direction
            VT v = GetDirection<Point>::apply(n, rng).getCoefficients();
            VT v0 = v;
            std::cout<<"Hi3"<<std::endl;
            precomputedValues.computed_C = false;

            // Begin a step of the random walk
            // Also, count the reflections and respect the bound
            while (reflectionsNum < reflectionsNumBound) {

                // we are at point p and the trajectory a*t^2 + vt + p
                // find how long we can walk till we hit the boundary
                std::cout<<"is p Exterior = "<<spectrahedron.isExterior(p)<<std::endl;
                std::cout<<"p = "<<p.transpose()<<std::endl;
                std::cout<<"computing oracle boundary"<<std::endl;
                NT lambda = spectrahedron.positiveIntersection(p, v, precomputedValues);
                if (lambda <=0){
                    std::cout<<"lambda = "<<lambda<<std::endl;
                    exit(-1);
                }
                std::cout<<"lambda = "<<lambda<<std::endl;
                lambda *= dl;
                VT qq = p + lambda * v;
                std::cout<<"is pi + tvi Exterior = "<<spectrahedron.isExterior(qq)<<std::endl;
                t_temp = rng.sample_urdist() * lambda;
                //w1 = T / (T + lambda);
                w2 = lambda / (T + lambda);
                std::cout<<"w2 = "<<w2<<std::endl;
                u = rng.sample_urdist();
                std::cout<<"u = "<<u<<std::endl;
                if (u <= w2) {
                    ti = t_temp;
                    xi = p;
                    vi = v;
                    Ci = precomputedValues.C;
                }
                T += lambda;

                // We just solved a quadratic polynomial eigenvalue system At^2 + Bt + C,
                // where A = lmi(a) - A0, B = lmi(v) - A0 and C = lmi(p)
                // and also did a linearization creating two matrices X, Y.
                // For the subsequent calls, we don't have to compute these matrices from scratch,
                // but can efficiently update them.
                // A remains the same
                // C := A*lambda^2 + B*lambda + C
                // X, Y will be updated in class Spectrahedron
                // Set the flags
                //precomputedValues.computed_A = true;
                precomputedValues.computed_C = true;
                //precomputedValues.computed_XY = true;

                // if we can walk the remaining distance without reaching he boundary
                if (stopping_criterion_is_true(p, v, lambda, p0, v0) && reflectionsNum > 0) {
                    std::cout<<"stopping reached"<<std::endl;
                    // set the new point p:= (T^2)*a + T*V + p
                    p = xi + ti * vi;
                    std::cout<<"is p Exterior = "<<spectrahedron.isExterior(p)<<std::endl;

                    // update matrix C
                    precomputedValues.C = Ci;
                    return;
                }

                // we hit the boundary and still have to walk
                // don't go all the way to the boundary, for numerical stability
                //lambda *= dl;

                // save current and set new point
                //previousPoint = p;
                p += lambda * v;

                // update remaining distance we must walk
                //T -= lambda;

                // update matrix C
                precomputedValues.C += lambda * precomputedValues.B;

                // Set v to have the direction of the trajectory at t = lambda
                // i.e. the gradient of at^2 + vt + p, for t = lambda
                //v += (lambda * 2) * a;

                // compute reflected direction
                VT reflectedTrajectory;
                spectrahedron.computeReflection(p, v, reflectedTrajectory, precomputedValues);
                v = reflectedTrajectory;
                std::cout<<"\n v = "<<v.transpose()<<"\n"<<std::endl;

                reflectionsNum++;
            }

            // if the #reflections exceeded the limit, don't move
            if (reflectionsNum == reflectionsNumBound){
                p = xi + ti * vi;
                precomputedValues.C = Ci;
            }
        }

        /// Sets the temperature in the distribution
        /// \param[in] temperature New value of temperature
        //void setTemperature(NT temperature) {
        //    settings.temperature = temperature;
       // }
    };



    /// The implementation of the walk for spectrahedra
    template <typename NT, typename RandomNumberGenerator>
    struct Walk<SpectrahedronBlock<NT, SparseBlock<NT>, Eigen::Matrix<NT,Eigen::Dynamic,1> >, RandomNumberGenerator> {

        /// The type of the spectrahedron
        typedef SpectrahedronBlock<NT, SparseBlock<NT>, Eigen::Matrix<NT,Eigen::Dynamic,1> > SPECTRAHEDRON;
        /// Type for internal structure of class Spectrahedron
        typedef typename SPECTRAHEDRON::PrecomputedValues2 PrecomputedValues;

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
            //VT c;
            /// The T in the distribution
            //NT temperature;
            /// The diameter of the spectrahedron
            //NT diameter;
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
            //template <typename Point>
            Settings(const int walkLength, const RandomNumberGenerator &randomNumberGenerator,
                     unsigned int reflectionsBound = 10, NT dl = 0.995) : walk_length(walkLength), 
                                                                          randomNumberGenerator(randomNumberGenerator),
                                                                          reflectionsBound(reflectionsBound), dl(dl) {}

            Settings() {}
        };

        /// The parameters of the random walk
        Settings settings;

        Walk() {}

        /// Constructor
        /// \param[in] settings The settings of the random walk
        Walk(Settings &settings) : settings(settings) {}

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
            int reflectionsNumBound = settings.reflectionsBound * n;
            VT previousPoint;
            VT p0 = p, xi, vi;
            MT Ci;
            //std::cout<<"Hi30"<<std::endl;
            //std::cout<<"p0 = "<<p0.transpose()<<std::endl;
            //std::cout<<"is p0 Exterior = "<<spectrahedron.isExterior(p0)<<std::endl;
            //std::cout<<"Hi32"<<std::endl;

            // choose a distance to walk
            NT T = 0.0, ti, w2, w1, u, t_temp;

            // The trajectory will be of the form a*t^2 + vt + p
            // where a = -c / 2*temperature
            // and at each reflection, v and p will change

            // crate vector a
            //VT a = -settings.c / (2 * settings.temperature);

            // The vector v will be a random a direction
            VT v = GetDirection<Point>::apply(n, rng).getCoefficients();
            VT v0 = v;
            //std::cout<<"Hi3"<<std::endl;
            precomputedValues.computed_C = false;

            // Begin a step of the random walk
            // Also, count the reflections and respect the bound
            while (reflectionsNum < reflectionsNumBound) {

                // we are at point p and the trajectory a*t^2 + vt + p
                // find how long we can walk till we hit the boundary
                //std::cout<<"is p Exterior = "<<spectrahedron.isExterior(p)<<std::endl;
                //std::cout<<"p = "<<p.transpose()<<std::endl;
                //std::cout<<"computing oracle boundary"<<std::endl;
                NT lambda = spectrahedron.positiveIntersection(p, v, precomputedValues);
                if (lambda <=0){
                    std::cout<<"lambda = "<<lambda<<std::endl;
                    exit(-1);
                }
                //std::cout<<"lambda = "<<lambda<<std::endl;
                lambda *= dl;
                VT qq = p + lambda * v;
                //std::cout<<"is pi + tvi Exterior = "<<spectrahedron.isExterior(qq)<<std::endl;
                t_temp = rng.sample_urdist() * lambda;
                //w1 = T / (T + lambda);
                w2 = lambda / (T + lambda);
                //std::cout<<"w2 = "<<w2<<std::endl;
                u = rng.sample_urdist();
                //std::cout<<"u = "<<u<<std::endl;
                if (u <= w2) {
                    ti = t_temp;
                    xi = p;
                    vi = v;
                    Ci = precomputedValues.C;
                }
                T += lambda;

                // We just solved a quadratic polynomial eigenvalue system At^2 + Bt + C,
                // where A = lmi(a) - A0, B = lmi(v) - A0 and C = lmi(p)
                // and also did a linearization creating two matrices X, Y.
                // For the subsequent calls, we don't have to compute these matrices from scratch,
                // but can efficiently update them.
                // A remains the same
                // C := A*lambda^2 + B*lambda + C
                // X, Y will be updated in class Spectrahedron
                // Set the flags
                //precomputedValues.computed_A = true;
                precomputedValues.computed_C = true;
                //precomputedValues.computed_XY = true;

                // if we can walk the remaining distance without reaching he boundary
                if (stopping_criterion_is_true(p, v, lambda, p0, v0) && reflectionsNum > 0) {
                    std::cout<<"stopping reached"<<std::endl;
                    // set the new point p:= (T^2)*a + T*V + p
                    p = xi + ti * vi;
                    //std::cout<<"is p Exterior = "<<spectrahedron.isExterior(p)<<std::endl;

                    // update matrix C
                    precomputedValues.C = Ci;
                    return;
                }

                // we hit the boundary and still have to walk
                // don't go all the way to the boundary, for numerical stability
                //lambda *= dl;

                // save current and set new point
                //previousPoint = p;
                p += lambda * v;

                // update remaining distance we must walk
                //T -= lambda;

                // update matrix C
                precomputedValues.C += precomputedValues.B * lambda;

                // Set v to have the direction of the trajectory at t = lambda
                // i.e. the gradient of at^2 + vt + p, for t = lambda
                //v += (lambda * 2) * a;

                // compute reflected direction
                VT reflectedTrajectory;
                spectrahedron.computeReflection(p, v, reflectedTrajectory, precomputedValues);
                v = reflectedTrajectory;
                //std::cout<<"\n v = "<<v.transpose()<<"\n"<<std::endl;

                reflectionsNum++;
            }

            // if the #reflections exceeded the limit, don't move
            if (reflectionsNum == reflectionsNumBound){
                p = xi + ti * vi;
                precomputedValues.C = Ci;
            }
        }

        /// Sets the temperature in the distribution
        /// \param[in] temperature New value of temperature
        //void setTemperature(NT temperature) {
        //    settings.temperature = temperature;
       // }
    };



    /// The implementation of the walk for spectrahedra
    template <typename NT, typename RandomNumberGenerator>
    struct Walk<Spectrahedron<NT, Eigen::SparseMatrix<NT>, Eigen::Matrix<NT,Eigen::Dynamic,1> >, RandomNumberGenerator> {

        /// The type of the spectrahedron
        typedef Spectrahedron<NT, Eigen::SparseMatrix<NT>, Eigen::Matrix<NT,Eigen::Dynamic,1>> SPECTRAHEDRON;
        /// Type for internal structure of class Spectrahedron
        typedef typename SPECTRAHEDRON::PrecomputedValues2 PrecomputedValues;

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
            //VT c;
            /// The T in the distribution
            //NT temperature;
            /// The diameter of the spectrahedron
            //NT diameter;
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
            //template <typename Point>
            Settings(const int walkLength, const RandomNumberGenerator &randomNumberGenerator,
                     unsigned int reflectionsBound = 10, NT dl = 0.95) : walk_length(walkLength), 
                                                                          randomNumberGenerator(randomNumberGenerator),
                                                                          reflectionsBound(reflectionsBound), dl(dl) {}

            Settings() {}
        };

        /// The parameters of the random walk
        Settings settings;

        Walk() {}

        /// Constructor
        /// \param[in] settings The settings of the random walk
        Walk(Settings &settings) : settings(settings) {}

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
            std::cout<<"walk_length = "<<settings.walk_length<<std::endl;
            // sample #pointsNum points
            for (unsigned int i = 1; i <= pointsNum; ++i) {
                // burn #walk_length points to get one sample
                for (unsigned int j = 0; j < settings.walk_length; ++j) {
                    std::cout<<"Hi0.5"<<std::endl;
                    getNextPoint<Point>(spectrahedron, p, precomputedValues);
                }

                VT qq = p;
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
            int reflectionsNumBound = settings.reflectionsBound * n;
            VT previousPoint;
            VT p0 = p, xi, vi, qq;
            MT Ci;
            std::cout<<"is p0 Exterior = "<<spectrahedron.isExterior(p0)<<std::endl;

            // choose a distance to walk
            NT T = 0.0, ti, w2, w1, u, t_temp;

            // The trajectory will be of the form a*t^2 + vt + p
            // where a = -c / 2*temperature
            // and at each reflection, v and p will change

            // crate vector a
            //VT a = -settings.c / (2 * settings.temperature);

            // The vector v will be a random a direction
            VT v = GetDirection<Point>::apply(n, rng).getCoefficients();
            VT v0 = v;
            std::cout<<"Hi3"<<std::endl;
            precomputedValues.computed_C = false;

            // Begin a step of the random walk
            // Also, count the reflections and respect the bound
            while (reflectionsNum < reflectionsNumBound) {

                // we are at point p and the trajectory a*t^2 + vt + p
                // find how long we can walk till we hit the boundary
                std::cout<<"is p Exterior = "<<spectrahedron.isExterior(p)<<std::endl;
                std::cout<<"p = "<<p.transpose()<<std::endl;
                std::cout<<"computing oracle boundary"<<std::endl;
                NT lambda = spectrahedron.positiveIntersection(p, v, precomputedValues);
                if (lambda <=0){
                    std::cout<<"lambda = "<<lambda<<std::endl;
                    exit(-1);
                }
                std::cout<<"lambda = "<<lambda<<std::endl;
                lambda *= dl;
                qq = p + lambda * v;
                std::cout<<"is pi + tvi Exterior = "<<spectrahedron.isExterior(qq)<<std::endl;
                t_temp = rng.sample_urdist() * lambda;
                //w1 = T / (T + lambda);
                w2 = lambda / (T + lambda);
                std::cout<<"w2 = "<<w2<<std::endl;
                u = rng.sample_urdist();
                std::cout<<"u = "<<u<<std::endl;
                if (u <= w2) {
                    ti = t_temp;
                    xi = p;
                    vi = v;
                    Ci = precomputedValues.C;
                }
                T += lambda;

                // We just solved a quadratic polynomial eigenvalue system At^2 + Bt + C,
                // where A = lmi(a) - A0, B = lmi(v) - A0 and C = lmi(p)
                // and also did a linearization creating two matrices X, Y.
                // For the subsequent calls, we don't have to compute these matrices from scratch,
                // but can efficiently update them.
                // A remains the same
                // C := A*lambda^2 + B*lambda + C
                // X, Y will be updated in class Spectrahedron
                // Set the flags
                //precomputedValues.computed_A = true;
                precomputedValues.computed_C = true;
                //precomputedValues.computed_XY = true;

                // if we can walk the remaining distance without reaching he boundary
                if (stopping_criterion_is_true(p, v, lambda, p0, v0) && reflectionsNum > 0) {
                    std::cout<<"stopping reached"<<std::endl;
                    // set the new point p:= (T^2)*a + T*V + p
                    p = xi + ti * vi;
                    std::cout<<"is p Exterior = "<<spectrahedron.isExterior(p)<<std::endl;

                    // update matrix C
                    precomputedValues.C = Ci;
                    return;
                }

                // we hit the boundary and still have to walk
                // don't go all the way to the boundary, for numerical stability
                //lambda *= dl;

                // save current and set new point
                //previousPoint = p;
                p += lambda * v;

                // update remaining distance we must walk
                //T -= lambda;

                // update matrix C
                precomputedValues.C += lambda * precomputedValues.B;

                // Set v to have the direction of the trajectory at t = lambda
                // i.e. the gradient of at^2 + vt + p, for t = lambda
                //v += (lambda * 2) * a;

                // compute reflected direction
                VT reflectedTrajectory;
                spectrahedron.computeReflection(p, v, reflectedTrajectory, precomputedValues);
                v = reflectedTrajectory;
                std::cout<<"\n v = "<<v.transpose()<<"\n"<<std::endl;

                reflectionsNum++;
            }

            // if the #reflections exceeded the limit, don't move
            if (reflectionsNum == reflectionsNumBound){
                p = xi + ti * vi;
                precomputedValues.C = Ci;
            }
        }

        /// Sets the temperature in the distribution
        /// \param[in] temperature New value of temperature
        void setTemperature(NT temperature) {
            settings.temperature = temperature;
        }
    };

};

#endif //VOLESTI_BOLTZMANN_HMC_WALK_HPP

