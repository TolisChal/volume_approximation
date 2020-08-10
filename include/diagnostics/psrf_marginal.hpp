// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

//Contributed and/or modified by Alexandros Manochis, as part of Google Summer of Code 2020 program.

//Licensed under GNU LGPL.3, see LICENCE file

/*
    This function implements a multivariate version of the Rubin & Gelman diagnostic.
    It is based on "Inference from iterative simulation using multiple sequences, 1992" by D. B. Rubin and A. Gelman

    The sample is splitted into two parts. Then the psrf of D.B. Rubin and A. Gelman is computed

    Inputs: samples, a matrix that contains sample points column-wise

    Output: The value of PSRF of D.B. Rubin and A. Gelman for each coordinate
*/

#ifndef PSRF_MARGINAL_HPP
#define PSRF_MARGINAL_HPP

template <typename NT, typename VT, typename MT>
VT perform_psrf_marginal(MT const& samples, bool fast = false)
{
    MT runs = samples.transpose();
    unsigned int N = samples.cols(), d = samples.rows();
    unsigned int N1 = N / 2;
    unsigned int N2 = N - N1;
    VT coord_samples(N), temp_col(N), results(d);
    NT mean1, mean2, mean00, sum, R, W, B;

    if (!fast)
    {
        for (int i = 0; i < d; i++)
        {
            coord_samples = runs.col(i);
            mean1 = coord_samples.block(0,0,N1,1).mean();
            mean2 = coord_samples.block(N1,0,N2,1).mean();

            sum = NT(0);
            for (int j = 0; j < N1; j++)
            {
                sum += (coord_samples(j) - mean1) * (coord_samples(j) - mean1);
            }
            W = sum / (NT(2) * (NT(N1) - NT(1)));

            sum = NT(0);
            for (int j = N1; j < N; j++)
            {
                sum += (coord_samples(j) - mean2) * (coord_samples(j) - mean2);
            }
            W += sum / (NT(2) * (NT(N2) - NT(1)));

            mean00 = coord_samples.mean();

            B = (mean1 - mean00) * (mean1 - mean00) + (mean2 - mean00) * (mean2 - mean00);
            V = ((NT(N1) - NT(1)) / NT(N1)) * W + B;
            R = V / W;
            results(i) = R;
        }
    } else 
    {
        VT sorted_samples(N), sorted_subsamples1(N1), sorted_subsamples2(N2);
        std::vector<NT> temp_col(N);
        NT alpha = 0.05;
        for (int i = 0; i < d; i++)
        {
            sorted_samples = runs.col(i);
            temp_col = std::vector<NT>(&sorted_samples[0], sorted_samples.data() + sorted_samples.cols() * 
                                       sorted_samples.rows());
            std::sort(temp_col.begin(), temp_col.end());
            sorted_samples = Eigen::Map<VT>(&temp_col[0], temp_col.size());

            int n1 = N * (alpha / NT(2)), n2 = N - N * (alpha / NT(2));

            int len_total_sequence_interval = sorted_samples(n2) - sorted_samples(n1);

            sorted_subsamples1 = runs.block(0,0,N1,1);
            temp_col.resize(N1);
            temp_col = std::vector<NT>(&sorted_subsamples1[0], sorted_subsamples1.data() + 
                                       sorted_subsamples1.cols() * sorted_subsamples1.rows());
            std::sort(temp_col.begin(), temp_col.end());
            sorted_subsamples1 = Eigen::Map<VT>(&temp_col[0], temp_col.size());

            sorted_subsamples2 = runs.block(N1,0,N2,1);
            temp_col.resize(N2);
            temp_col = std::vector<NT>(&sorted_subsamples2[0], sorted_subsamples2.data() + 
                                       sorted_subsamples2.cols() * sorted_subsamples2.rows());
            std::sort(temp_col.begin(), temp_col.end());
            sorted_subsamples2 = Eigen::Map<VT>(&temp_col[0], temp_col.size());

            n1 = N1 * (alpha / NT(2)), n2 = N1 - N1 * (alpha / NT(2));
            int len_sequence_interval1 = sorted_subsamples1(n2) - sorted_subsamples1(n1);

            n1 = N2 * (alpha / NT(2)), n2 = N2 - N2 * (alpha / NT(2));
            int len_sequence_interval2 = sorted_subsamples2(n2) - sorted_subsamples2(n1);

            R = (NT(len_total_sequence_interval)) / 
                 ((NT(len_sequence_interval1) + NT(len_sequence_interval2)) / NT(2));

            results(i) = std::abs(1.0 - R) + NT(1);
        }
    }

    return results;
}


#endif
