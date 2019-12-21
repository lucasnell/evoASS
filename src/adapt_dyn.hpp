#ifndef __SAURON_ADAPT_DYN_H
#define __SAURON_ADAPT_DYN_H


#include <RcppArmadillo.h>
#include "sim.hpp"

using namespace Rcpp;




/*
 Output info for one repetition:
 */
class OneRepInfoAD {
public:

    /*
     -----------------
     Objects that change size over time
    -----------------
    */
    // Current abundances:
    std::vector<double> N;
    // Density dependences:
    std::vector<double> A;
    // Indices to trait values for all lines (this will get updated, though):
    std::vector<uint32_t> I;
    // Current clone index:
    uint32_t clone_I;

    /*
    -----------------
    Objects that store all information (and just get added to over time)
    -----------------
    */
    // Trait values for all clones:
    std::vector<arma::rowvec> all_V;
    // All abundances from t = 1 to max_t:
    std::vector<std::vector<double>> all_N;
    // All indices from t = 1 to max_t:
    std::vector<std::vector<uint32_t>> all_I;
    // All time points sampled:
    std::vector<uint32_t> all_t;


    OneRepInfoAD() {};
    OneRepInfoAD(const std::vector<arma::rowvec>& V0,
                 const std::vector<double>& N0,
                 const uint32_t& max_clones,
                 const uint32_t& max_t,
                 const uint32_t& save_every,
                 const double& mut_sd)
        : N(N0), A(N0.size()), I(N0.size()), clone_I(0),
          all_V(V0), all_N(1), all_I(1), all_t(1),
          norm_distr(0, mut_sd) {

        N.reserve(max_clones);
        A.reserve(max_clones);
        while (clone_I < I.size()) {
            I[clone_I] = clone_I;
            clone_I++;
        }
        I.reserve(max_clones);

        all_V.reserve(V0.size() + max_clones);
        all_N[0] = N;
        uint32_t max_times = (max_t / save_every) + 2;
        all_N.reserve(max_times);
        all_I[0] = I;
        all_I.reserve(max_times);
        all_t[0] = 0;
        all_t.reserve(max_times);

    }


    void iterate(const uint32_t& t,
                 const double& f,
                 const double& a0,
                 const arma::mat& C,
                 const double& r0,
                 const double& d,
                 const double& max_t,
                 const double& min_N,
                 const double& mut_sd,
                 const double& mut_prob,
                 const uint32_t& save_every,
                 pcg64& eng) {

        // Fill in density dependences:
        A_VNI__<std::vector<double>>(A, all_V, N, I, a0, d);

        // Extinct clones (if any):
        std::vector<uint32_t> extinct;
        // Fill in abundances:
        for (uint32_t i = 0; i < A.size(); i++) {
            double r = r_V_(all_V[I[i]], f, C, r0);
            N[i] *= std::exp(r - A[i]);
            // See if it goes extinct:
            if (N[i] < min_N) {
                extinct.push_back(i);
            }
        }
        // If everything is gone, output results and stop simulations:
        if (extinct.size() == N.size()) {
            all_N.push_back(N);
            all_I.push_back(I);
            all_t.push_back(t);
            return;
        }
        // Remove extinct clones (starting at the back):
        for (uint32_t i = 0, j; i < extinct.size(); i++) {
            j = extinct.size() - i - 1;
            N.erase(N.begin() + extinct[j]);
            A.erase(A.begin() + extinct[j]);
            I.erase(I.begin() + extinct[j]);
        }

        // Seeing if I should add new clones:
        uint32_t n_clones = N.size(); // doing this bc N.size() might increase
        for (uint32_t i = 0; i < n_clones; i++) {

            double u = runif_01(eng);

            if (u < mut_prob) {

                N[i] -= (1.01 * min_N);

                N.push_back(1.01 * min_N);
                A.push_back(0);
                I.push_back(clone_I);
                clone_I++;

                all_V.push_back(all_V[I[i]]);
                arma::rowvec& new_V(all_V.back());
                for (double& v : new_V) {
                    v += norm_distr(eng);
                    if (v < 0) v *= -1; // <-- keeping traits >= 0
                }

            }

        }

        if ((t+1) % save_every == 0 || (t+1) == max_t) {
            all_N.push_back(N);
            all_I.push_back(I);
            all_t.push_back(t);
        }

        return;
    }



    // How many rows is required for this repetition?
    uint32_t n_rows() const {
        uint32_t nr = 0;
        for (const std::vector<double>& nn : all_N) nr += nn.size();
        return nr;
    }


    /*
     Fill a matrix with data from this rep.
     The matrix has the following columns:
       0  - rep
       1  - time
       2  - clone index
       3  - abundance
       ≥4 - clone trait values
     */
    void fill_matrix(arma::mat& matrix,
                     const double& rep_number,
                     uint32_t start_row) const {

        uint32_t q = all_V[0].n_elem;
        double rn = static_cast<double>(rep_number);

        for (uint32_t i = 0; i < all_N.size(); i++) {

            uint32_t nr = all_N[i].size();
            double t = static_cast<double>(all_t[i]);

            for (uint32_t j = 0; j < nr; j++) {

                matrix(start_row + j, 0) = rn;
                matrix(start_row + j, 1) = t;
                matrix(start_row + j, 2) = all_I[i][j];
                matrix(start_row + j, 3) = all_N[i][j];

                for (uint32_t k = 0; k < q; k++) {
                    matrix(start_row + j, 4 + k) = all_V[all_I[i][j]][k];
                }

            }

            start_row += nr;

        }
        return;
    }


private:

    std::normal_distribution<double> norm_distr;

};


#endif
