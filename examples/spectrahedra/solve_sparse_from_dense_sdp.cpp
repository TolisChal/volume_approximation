// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis

//Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.

// Licensed under GNU LGPL.3, see LICENCE file

// This examples illustrates how to solve a semidefinite program.
// It will read a semidefinite program from data/sdp_n2m3.txt, solve it and print its solution (minimal value).


//#define VOLESTI_DEBUG
#define SPARSE_PROBLEM

#include <fstream>
#include <iostream>

#include "random.hpp"
#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "SDPAFormatManager.h"
#include "optimization/simulated_annealing.hpp"
#include "random_walks/boltzmann_hmc_walk.hpp"
#include "compute_initial_point.hpp"
//#include "random_walks/boltzmann_hmc_opt_walk.hpp"
//#include "random_walks/boltzmann_hmc_full_opt_walk.hpp"


typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix <NT, Eigen::Dynamic, Eigen::Dynamic> DMT;
typedef Eigen::SparseMatrix <NT> MT;
typedef Cartesian <NT> Kernel;
typedef typename Kernel::Point Point;
typedef Spectrahedron <NT, DMT, VT> DSPECTRAHEDRON;
typedef Spectrahedron <NT, MT, VT> SPECTRAHEDRON;
typedef LMI <NT, MT, VT> lmi;
//typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
//typedef BoltzmannHMCWalk::Walk<Spectrahedron, RNGType > HMC1;
//typedef BoltzmannHMCOptWalk::Walk<Spectrahedron, RNGType > HMC2;





int main(const int argc, const char** argv) {

    SPECTRAHEDRON spectrahedron;
    Point objFunction, p;
    bool correct, verbose = false, best_on_traj = true, set_p = false;
    NT optimal_val;
    NT rel_error = 0.001, k=0.5;
    lmi Slmi;


    for(int i=1; i<argc; ++i){

        correct = false;

        if(!strcmp(argv[i],"-sdp")){
            
            std::cout<<"Reading input from file..."<<std::endl;
            std::ifstream in;
            in.open(argv[++i], std::ifstream::in);
            //SdpaFormatManager<NT> sdpaFormatManager;
            VT q;
            loadSparseSDPAFormatFile<MT, NT>(in, Slmi, q);
            std::cout<<"file ok..."<<std::endl;
            objFunction = Point(q);
            spectrahedron = SPECTRAHEDRON(Slmi);
            correct = true;
        }

        if(!strcmp(argv[i],"-inner_point")){

            std::cout<<"Reading inner point from file..."<<std::endl;

            std::ifstream inpF;
            inpF.open(argv[++i], std::ifstream::in);

            // read point
            std::string line;
            std::getline(inpF, line, '\n');

            std::vector<double> vec = readVector3(line);
            int n = vec.size();
            VT point(n);
            for (int i = 0 ; i<n; ++i)
                point(i) = vec[i];

            p = Point(point);
            set_p = true;

            //std::cout<<"p = "<<point<<std::endl;
            
            
            //std::ifstream in;
            //in.open(argv[++i], std::ifstream::in);
            //SdpaFormatManager<NT> sdpaFormatManager;
            //sdpaFormatManager.loadSDPAFormatFile(in, spectrahedron, objFunction);
            correct = true;
        }

        if(!strcmp(argv[i],"-opt_val")){
            optimal_val = atof(argv[++i]);
            correct = true;
        }

        if(!strcmp(argv[i],"-error")){
            rel_error = atof(argv[++i]);
            correct = true;
        }

        if(!strcmp(argv[i],"-no_opts")){
            best_on_traj = false;
            correct = true;
        }

        if(!strcmp(argv[i],"-k")){
            k = atof(argv[++i]);
            correct = true;
        }

        if(!strcmp(argv[i],"-verbose") || !strcmp(argv[i],"-v")){
            verbose = true;
            correct = true;
        }

        if(correct==false){
            std::cerr<<"unknown parameters \'"<<argv[i]<<
                     "\', try "<<argv[0]<<" --help"<<std::endl;
            exit(-2);
        }

    }

    // We will need an initial interior point. In this
    // spectrahedron the origin (zero point) is interior
    if (!set_p){
        get_inner_point<NT, MT, VT, SPECTRAHEDRON>(Slmi, p, objFunction);
        VT qq = p.getCoefficients();
        std::cout << "Point in spectra, isNegativeSemidefinite = "<<Slmi.isNegativeSemidefinite(qq)<<std::endl;
        //exit(-1);
    }
    //Point initialPoint(spectrahedron.getLMI().dimension());

    // First some parameters for the solver
    // desired relative error
    

    // Declare settings
    SimulatedAnnealingSettings<Point> settings(rel_error, 1, -1, k);

    // solve the program
    Point sol;
    NT min, average_runtime = 0.0;
    double tstart1, tstop1;
    std::cout << "dimension: " << spectrahedron.getLMI().dimension() <<std::endl;



    for (int i = 0; i < 5; i++) {

        //initializations
        SPECTRAHEDRON spectrahedron_temp = spectrahedron;
        SimulatedAnnealingSettings<Point> settings_temp = settings;
        Point initialPoint_temp(p.getCoefficients()), sol_temp(spectrahedron.getLMI().dimension());

        tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
        min = solve_sdp_with_optimal<BoltzmannHMCWalk>(spectrahedron_temp, objFunction, settings_temp, 
                                                           initialPoint_temp, sol_temp, optimal_val, verbose);
        tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
        average_runtime += (tstop1 - tstart1);

        //computation
        /*if (best_on_traj) {
            tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
            min = solve_sdp_with_optimal<BoltzmannHMCFullOptWalk>(spectrahedron_temp, objFunction, settings_temp, 
                                                              initialPoint_temp, sol_temp, optimal_val, verbose);
            tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
            average_runtime += (tstop1 - tstart1);
        } else {
            tstart1 = (double)clock()/(double)CLOCKS_PER_SEC;
            
            tstop1 = (double)clock()/(double)CLOCKS_PER_SEC;
            average_runtime += (tstop1 - tstart1);
        }*/

        // print solution
        std::cout << "time: " << (tstop1 - tstart1)<< ", min = "<< min << ", error = " << std::abs(optimal_val - min) / std::abs(optimal_val) <<std::endl;
        //sol_temp.print();
    }

    average_runtime = average_runtime / 5.0;
    std::cout << "average_runtime: " << average_runtime << "\n\n" << std::endl;
    std::ofstream myfile;
    myfile.open("average_runtimes.txt", std::ofstream::app);
    myfile << ", ";
    myfile << average_runtime;
    

    return 0;
}



