// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_SPECTRAHEDRONBLOCK_H
#define VOLESTI_SPECTRAHEDRONBLOCK_H

#include <format.hpp>

#include "BlockLMI.h"
#include "chrono"
#include "generators/boost_random_number_generator.hpp"
#include "../../sampling/sphere.hpp"


template <typename MT>
struct PrecomputedValues {
    
};


template <typename NT>
struct PrecomputedValues< SparseBlock<NT> > {

    typedef SparseBlock<NT> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    /// These flags indicate whether the corresponding matrices are computed
    /// if yes, we can use them and not compute them fro scratch
    bool computed_A = false;
    bool computed_C = false;
    bool computed_XY = false;

    /// The matrices the method positiveIntersection receives from its previous call
    /// if the flag first_positive_intersection is true.
    /// Matrix A is also used in coordinateIntersection
    MT A, B, C;//, X, Y;

    /// In method positive_intersect, the distance we are computing corresponds
    /// to the minimum positive eigenvalue of a quadratic eigenvalue problem.
    /// This will hold the eigenvector for that eigenvalue
    VT eigenvector;

    /// Sets all flags to false
    void resetFlags() {
        computed_XY = computed_C = computed_A = false;
    }

    void set_mat_size(int const& m) 
    {
        //matrix_initializer<MT>::initialize(A, B, C, X, Y, m);

        //A.setZero(m, m);
        //B.setZero(m, m);
        //C.setZero(m, m);

        eigenvector.setZero(m);

        //X.setZero(2*m, 2*m);
        //Y.setZero(2*m, 2*m);
    }
};



/// This class manipulates a spectrahedron, described by a LMI
/// \tparam NT Numeric Type
/// \tparam MT Matrix Type
/// \tparam VT Vector Type
template<typename NT>
class SpectrahedronBlock {
public:

    typedef SparseBlock<NT> MT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    /// The numeric/matrix/vector types we use
    typedef NT NUMERIC_TYPE;
    //typedef MT MATRIX_TYPE;
    //typedef VT VECTOR_TYPE;
    typedef PrecomputedValues<MT> PrecomputedValues2;

    /// The type of a pair of NT
    typedef std::pair<NT, NT> pairNT;

    BlockEigenvaluesProblems<NT, MT, VT> EigenvaluesProblem;

    /// The dimension of the spectrahedron
    unsigned int d;
    VT grad;

    /// The linear matrix inequality that describes the spectrahedron
    BlockLMI<NT> lmi;

    Spectrahedron() {}

    /// Creates a spectrahedron
    /// \param[in] lmi The linear matrix inequality that describes the spectrahedron
    Spectrahedron(const BlockLMI<NT>& _lmi) : lmi(_lmi) {
        d = _lmi.dimension();
        grad.setZero(d);
    }


    /// Construct the quadratic eigenvalue problem \[At^2 + Bt + C \] for positive_intersect.
    /// A = lmi(c) - A0, B = lmi(b) - A0 and C = lmi(c).
    /// \param[in] a Input vector
    /// \param[in] b Input vector
    /// \param[in] c Input vector
    /// \param[in, out] precomputedValues Holds matrices A, C
    void createMatricesForPositiveIntersection(const VT& a, const VT& b, const VT& c,
            PrecomputedValues2& precomputedValues) {
        // check if matrices A, C are ready
        // if not compute them
        if (!precomputedValues.computed_A) {
            lmi.evaluateWithoutA0(a, precomputedValues.A);
        }

        if (!precomputedValues.computed_C) {
            lmi.evaluate(c, precomputedValues.C);
        }

        // compute Matrix B
        lmi.evaluateWithoutA0(b, precomputedValues.B);
    }

    /// Construct the quadratic eigenvalue problem \[At^2 + Bt + C \] for positive_intersect.
    /// A = lmi(c) - A0, B = lmi(b) - A0 and C = lmi(c).
    /// \param[in] a Input vector
    /// \param[in] b Input vector
    /// \param[in] c Input vector
    /// \param[in, out] precomputedValues Holds matrices A, C
    void createMatricesForPositiveIntersection(const VT& p, const VT& v,
            PrecomputedValues2& precomputedValues) {
        // check if matrices A, C are ready
        // if not compute them
        if (!precomputedValues.computed_C) {
            lmi.evaluate(p, precomputedValues.C);
        }

        //if (!precomputedValues.computed_C) {
        //    lmi.evaluate(c, precomputedValues.C);
        //}

        // compute Matrix B
        lmi.evaluateWithoutA0(v, precomputedValues.B);
    }

    /// Computes the distance d we must travel on the parametrized polynomial curve \[at^2 + bt + c \],
    /// assuming we start at t=0, and we start increasing t.
    /// We construct the quadratic eigenvalue problem \[At^2 + Bt + C \],
    /// where A = lmi(c) - A0, B = lmi(b) - A0 and C = lmi(c).
    /// Then we do a linearization and transform it to the generalized problem X+lY,
    /// which we pass to an external library.
    /// \param[in] a Input vector, the coefficient of t \[t^2\]
    /// \param[in] b Input vector, the coefficient of t
    /// \param[in] c Input Vector, the constant term
    /// \param[in, out] precomputedValues Data we move between successive method calls
    /// \returns The distance d
    NT positiveIntersection(VT const & a, VT const & b, VT const & c, PrecomputedValues2& precomputedValues) 
    {
        //unsigned int matrixDim = lmi.sizeOfMatrices();

        // create matrices A, B, C
        createMatricesForPositiveIntersection(a, b, c, precomputedValues);

        // get the minimum positive eigenvalue of At^2 + Bt + C
        
        
        NT distance = EigenvaluesProblem.minPosQuadraticEigenvalue(precomputedValues.A, precomputedValues.B,
                                                                            precomputedValues.C,
                                                                            precomputedValues.eigenvector);

        return distance;
    }

    NT positiveIntersection(VT const & p, VT const & v, PrecomputedValues2& precomputedValues) 
    {
        //unsigned int matrixDim = lmi.sizeOfMatrices();

        // create matrices A, B, C
        createMatricesForPositiveIntersection(p, v, precomputedValues);
        std::cout<<"matrices created"<<std::endl;

        // get the minimum positive eigenvalue of At^2 + Bt + C
        
        std::cout<<"computing eigenvalue..."<<std::endl;
        NT distance = EigenvaluesProblem.minPosLinearEigenvalue(precomputedValues.C, precomputedValues.B,
                                                                            precomputedValues.eigenvector);
        std::cout<<"distance = "<<distance<<std::endl;

        return distance;
    }

     pairNT RandomDirectionIntersection(VT const &p, VT const &v, PrecomputedValues2& precomputedValues) {

        // prepare the generalized eigenvalue problem A+lB
        // we may not have to compute A!
        if (!precomputedValues.computed_C)
            lmi.evaluate(p, precomputedValues.C);

        lmi.evaluateWithoutA0(v, precomputedValues.B);

        return EigenvaluesProblem.symGeneralizedProblem(precomputedValues.C, precomputedValues.B);
    }


    /// Computes the distance d one must travel on the line a + tb,
    /// assuming we start at t=0 and that b has zero everywhere and 1 in its i-th coordinate.
    /// We must solve the generalized eigenvalue problem A+tB, where A = lmi(a) and B=(lmi) - A0 = A_i
    /// If the flag precomputedValues,computed_A is true, the matrix A is not computed.
    /// \param[in] a Input vector
    /// \param[in] coordinate Indicator of the i-th coordinate, 1 <= coordinate <= dimension
    /// \return The pair (positive t, negative t) for which we reach the boundary
    pairNT coordinateIntersection(VT const & a, int const coordinate, PrecomputedValues2& precomputedValues) {

        // prepare the generalized eigenvalue problem A+lB
        // we may not have to compute A!
        if (!precomputedValues.computed_A)
            lmi.evaluate(a, precomputedValues.A, true);

        
        return EigenvaluesProblem.symGeneralizedProblem(precomputedValues.A, *(lmi.getMatrix(coordinate)));
    }

    /// Computes the reflected direction at a point on the boundary of the spectrahedron.
    /// \param[in] point A point on the boundary of the spectrahedron
    /// \param[in] incomingDirection The direction of the trajectory as it hits the boundary
    /// \param[out] reflectedDirection The reflected direction
    /// \param[in] precomputedValues Must contain an eigenvalue needed to compute the reflection
    void computeReflection(VT const & point, VT const & incomingDirection, VT& reflectedDirection,
            PrecomputedValues2& precomputedValues) {

        // get the gradient of the determinant of the lmi at point
        //VT grad;
        lmi.normalizedDeterminantGradient(point, precomputedValues.eigenvector, grad);

        // compute reflected direction
        // if v is original direction and s the surface normal,
        // reflected direction = v - 2 <v,s>*s

        NT dot = 2 * incomingDirection.dot(grad);
        reflectedDirection = incomingDirection - dot * grad;
    }


    /// \return The dimension of the spectrahedron
    unsigned int dimension() const {
        return d;
    }

    /// \return The LMI describing this spectrahedron
    LMI<NT, MT, VT> getLMI() const {
        return lmi;
    }



    /// Estimates the diameter of the spectrahedron. It samples points uniformly with coordinate directions
    /// hit and run, and returns the maximum distance between them.
    /// \tparam Point
    /// \param[in] numPoints The number of points to sample for the estimation
    /// \return An estimation of the diameter of the spectrahedron
    template<class Point, class RNGType>
    NT estimateDiameterRDHR(int const numPoints, Point const & interiorPoint, RNGType &rng, VT const& c, VT &best_point) {
        //typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;

        //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        //RNGType rng(seed);
        int n = dimension();

        std::list<Point> randPoints;

        // initialize random numbers generators
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(1, d);

        PrecomputedValues2 precomputedValues;
        precomputedValues.set_mat_size(getLMI().sizeOfMatrices());
        VT p = interiorPoint.getCoefficients(), v;
        NT maxDistance = 0, best_opt = c.dot(p);
        //best_point = interiorPoint;

        // sample points with walk length set to 1
        int samplingNo=0;
        VT q;
        while (samplingNo<numPoints) {
            // uniformly select a line parallel to an axis,
            // i.e. an indicator i s.t. x_i = 1
            //int coordinate = uidist(rng);
            v = GetDirection<Point>::apply(n, rng).getCoefficients();

            // get the distances we can travel from p
            // on the line p + t* e_coordinate
            // before reaching the boundary
            //std::pair<NT, NT> distances = this->coordinateIntersection(p, coordinate, precomputedValues);
            std::pair<NT, NT> distances = this->RandomDirectionIntersection(p, v, precomputedValues);
            //std::cout << "distances.second = "<<distances.second<<std::endl;
            //std::cout << "distances.first = "<<distances.first<<"\n"<<std::endl;
            // uniformly set the new point on the segment
            // defined by the intersection points
            NT lambda = rng.sample_urdist();
            NT diff = distances.second + lambda * (distances.first - distances.second);
            if ( (distances.first - distances.second) > maxDistance ) {
                maxDistance = (distances.first - distances.second);
            }

            q = p + diff * v;
            if (this->isExterior(p)==1) {
                continue;
            }
            p=q;

            // update the precomputedValues, so we can skip
            // computations in the next call
            precomputedValues.computed_C = true;
            precomputedValues.C += diff * precomputedValues.B;
            randPoints.push_back(Point(p));
            samplingNo++;
        }

        // find maximum distance among points;
        
        typename std::list<Point>::iterator itInner, itOuter = randPoints.begin();

        for (; itOuter!=randPoints.end() ; ++itOuter) {
            NT temp_opt = c.dot((*itOuter).getCoefficients());
            if (temp_opt < best_opt) {
                best_opt = temp_opt;
                best_point = (*itOuter).getCoefficients();
            }
            for (itInner=itOuter ; itInner!=randPoints.end() ; ++itInner) {
                NT current = itOuter->distance(*itInner);
                if (current > maxDistance)
                    maxDistance = current;
            }
        }

        return maxDistance;
    }


   
    /// Estimates the diameter of the spectrahedron. It samples points uniformly with coordinate directions
    /// hit and run, and returns the maximum distance between them.
    /// \tparam Point
    /// \param[in] numPoints The number of points to sample for the estimation
    /// \return An estimation of the diameter of the spectrahedron
    template<class Point>
    NT estimateDiameter(int const numPoints, Point const & interiorPoint) {
        typedef boost::mt19937 RNGType;

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RNGType rng(seed);

        std::list<Point> randPoints;

        // initialize random numbers generators
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(1, d);

        PrecomputedValues2 precomputedValues;
        precomputedValues.set_mat_size(getLMI().sizeOfMatrices());
        VT p = interiorPoint.getCoefficients();

        // sample points with walk length set to 1
        for (int samplingNo=0 ; samplingNo<numPoints ; ++samplingNo) {
            // uniformly select a line parallel to an axis,
            // i.e. an indicator i s.t. x_i = 1
            int coordinate = uidist(rng);

            // get the distances we can travel from p
            // on the line p + t* e_coordinate
            // before reaching the boundary
            std::pair<NT, NT> distances = this->coordinateIntersection(p, coordinate, precomputedValues);

            // uniformly set the new point on the segment
            // defined by the intersection points
            NT lambda = urdist(rng);
            NT diff = distances.first + lambda * (distances.second - distances.first);

            p(coordinate - 1) = p(coordinate - 1) + diff;

            // update the precomputedValues, so we can skip
            // computations in the next call
            precomputedValues.computed_A = true;
            precomputedValues.A += diff * *(this->getLMI().getMatrix(coordinate));
            randPoints.push_back(Point(p));
        }

        // find maximum distance among points;
        NT maxDistance = 0;
        typename std::list<Point>::iterator itInner, itOuter = randPoints.begin();

        for (; itOuter!=randPoints.end() ; ++itOuter)
            for (itInner=itOuter ; itInner!=randPoints.end() ; ++itInner) {
                NT current = itOuter->distance(*itInner);
                if (current > maxDistance)
                    maxDistance = current;
            }

        return maxDistance;
    }

    /// Find out is lmi(current position) = mat is in the exterior of the spectrahedron
    /// \param mat a matrix where mat = lmi(current position)
    /// \return true if position is outside the spectrahedron
    bool isExterior(MT const & mat) {
        return !lmi.isNegativeSemidefinite(mat);
    }

    /// Find out is pos is in the exterior of the spectrahedron
    /// \param pos a vector
    /// \return true if pos is outside the spectrahedron
    bool isExterior(VT &pos) {
        return !lmi.isNegativeSemidefinite(pos);
    }
};

#endif //VOLESTI_SPECTRAHEDRON_H
