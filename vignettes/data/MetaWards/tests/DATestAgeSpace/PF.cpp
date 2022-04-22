// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::depends(sitmo)]]

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include <sitmo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// headers for bessel functions
#include "bessel.h"

#define ML_POSINF	R_PosInf

#define min0(x, y) (((x) <= (y)) ? (x) : (y))
#define max0(x, y) (((x) <= (y)) ? (y) : (x))

static void K_bessel_new(double *x, double *alpha, int *nb,
		     int *ize, double *bk, int *ncalc);
static void I_bessel_new(double *x, double *alpha, int *nb,
		     int *ize, double *bi, int *ncalc);
double bessel_i_ex_new(double x, double alpha, double expo, double *bi);
double bessel_k_ex_new(double x, double alpha, double expo, double *bk);
double Rf_gamma_cody(double x);


using namespace Rcpp;

// log-sum-exp function to prevent numerical overflow
double log_sum_exp(arma::vec &x, int mn = 0) {
    double maxx = max(x);
    double y = maxx + log(sum(exp(x - maxx)));
    if(mn == 1) y -= log(x.n_elem);
    return y;
}

// Poisson RNG using inverse transform method
// (to try to circumvent non thread-safe RNG in R)
int rpois_cpp (double lambda, sitmo::prng &eng) {
    if(lambda < 0.0) {
        stop("'lambda' must be > 0 in rpois\n");
    }
    double mx = sitmo::prng::max();
    double u = log(eng()) - log(mx);
    double temp = -lambda;
    double temp1 = 0.0;
    double maxt = 0.0;
    double lk = 0.0;
    int k = 0;
    while(temp < u) {
        k++;
        temp1 = k * log(lambda) - lambda;
        lk += log(k);
        temp1 -= lk;
        maxt = (temp > temp1 ? temp:temp1);
        temp = maxt + log(exp(temp - maxt) + exp(temp1 - maxt));
    }
    return k;
}

// Binomial RNG using inverse transform method
// (to try to circumvent non thread-safe RNG in R)
int rbinom_cpp (int n, double p, sitmo::prng &eng) {
    if(n < 0) {
        Rprintf("n = %d\n", n);
        stop("'n' must be >= 0 in rbinom\n");
    }
    if(p < 0.0 || p > 1.0) {
        Rprintf("p = %f\n", p);
        stop("Must have 0 <= p <= 1 in rbinom\n");
    }
    if(n == 0) return 0;
    double mx = sitmo::prng::max();
    double u = log(eng()) - log(mx);
    int k;
    double ln = 0.0, lx = 0.0, lnmx = 0.0;
    for(k = 1; k <= n; k++) {
        ln += log(k);
    }
    lnmx = ln;
    double temp = ln - lx - lnmx + n * log(1.0 - p);
    double temp1 = 0.0, maxt = 0.0;
    k = 0;
    while(temp < u) {
        k++;
        if(k > n) {
            // check for rounding errors
            if(fabs(temp) < 1e-15 && fabs(u) < 1e-15) return(n);
            stop("Error in binomial sampling: lcum = %e lu = %e k = %d n = %d\n", temp, u, k, n);
        }
        lx += log(k);
        lnmx -= log(n - (k - 1));
        temp1 = ln - lx - lnmx + k * log(p) + (n - k) * log(1.0 - p);
        maxt = (temp > temp1 ? temp:temp1);
        temp = maxt + log(exp(temp - maxt) + exp(temp1 - maxt));
    }
    return k;
}

// Multinomial RNG for n = 1 using inverse transform method
// (to try to circumvent non thread-safe RNG in R)
int rmultinom_cpp (arma::vec &p, sitmo::prng &eng) {
    if(fabs(sum(p) - 1.0) > 1e-10) {
        Rprintf("sum(p) = %e\n", sum(p));
        stop("Must have sum(p) == 1 in rmultinom\n");
    }
    if(any(p) > 1.0 || any(p) < 0.0) {
        stop("Must have all 0 <= p <= 1 in rmultinom\n");
    }
    double mx = sitmo::prng::max();
    double u = eng() / mx;
    double temp = p(0);
    int k = 0;
    while(temp < u) {
        k++;
        if(k >= p.n_elem) stop("Error in multinomial sampling\n");
        temp += p(k);
    }
    return k;
}  

// log Skellam CDF (note, no checks on inputs)
// adapted from "pskellam" source code by Patrick Brown (mistakes are mine)
double lpskellam_cpp(int x, double lambda1, double lambda2) {
    // different formulas for negative & nonnegative x (zero lambda is OK)
    double ldens = 0.0;
    if(x < 0) {
        ldens = R::pnchisq(2.0 * lambda2, -2.0 * x, 2.0 * lambda1, 1, 1);
    } else {
        ldens = R::pnchisq(2.0 * lambda1, 2.0 * (x + 1), 2.0 * lambda2, 0, 1);
    }
    if(!arma::is_finite(ldens)) {
        Rprintf("Non-finite pdensity in lpskellam = %f x = %d lambda1 = %f lambda2 = %f\n", ldens, x, lambda1, lambda2);
    }
    return ldens;
}

// log truncated Skellam density (note, no checks on inputs)
// adapted from "dskellam" source code by Patrick Brown (mistakes are mine)
double ldtskellam_cpp(int x, double lambda1, double lambda2, char *str, int LB = 0, int UB = -1, int print = 0) {
    
    // if UB < LB then returns untruncated density
    
    // check for this instance only
    if(LB > UB) {
        stop("ldt LB > UB\n");
    }
    
    if(print == 1) Rprintf("State: %s\n", str);
    
    // create array for Bessel (trying to deal with recursive
    // gc error message that may be coming from memory allocation
    // in bessel_i - hence swapping to bessel_i_ex and controlling
    // vector allocation directly here through calloc and free -
    // OK to use std::calloc here because object not going
    // back to R I think)
    double *bi = (double *) calloc (floor(abs(x)) + 1, sizeof(double));
    
    // from "skellam" source code: log(besselI(y, nu)) == y + log(besselI(y, nu, TRUE))
    double ldens = -(lambda1 + lambda2) + (x / 2.0) * (log(lambda1) - log(lambda2)) +
        log(bessel_i_ex_new(2.0 * sqrt(lambda1 * lambda2), abs(x), 2, bi)) + 2.0 * sqrt(lambda1 * lambda2);
    
    // free memory from the heap
    free(bi);
    
//    // if non-finite density, then try difference of CDFs
//    if(!arma::is_finite(ldens)) {
//        Rprintf("Trying difference of CDFs in ldtskellam\n");
//        double norm0 = 0.0, norm1 = 0.0;
//        norm0 = lpskellam_cpp(x - 1, lambda1, lambda2);
//        norm1 = lpskellam_cpp(x, lambda1, lambda2);
//        if(!arma::is_finite(norm0) && !arma::is_finite(norm1)) {
//            ldens = norm0;
//        } else {
//            ldens = norm1 + log(1.0 - exp(norm0 - norm1));
//        }
//    }
        
    // if truncated then adjust log-density
    if(UB >= LB) {
        if(x < LB || x > UB) {
            Rprintf("'x' must be in bounds in dtskellam_cpp: x = %d LB = %d UB = %d\n", x, LB, UB);
            stop("");
        }
        if(arma::is_finite(ldens)) {
            // declare variables
            int normsize = UB - LB + 1;
            if(normsize > 1) {
                double norm0 = 0.0, norm1 = 0.0;
                norm0 = lpskellam_cpp(LB - 1, lambda1, lambda2);
                norm1 = lpskellam_cpp(UB, lambda1, lambda2); 
                if(!arma::is_finite(norm0) && !arma::is_finite(norm1)) {
                    stop("Issue with truncation bounds in Skellam\n");
                }
                double lnorm = norm1 + log(1.0 - exp(norm0 - norm1));
                ldens -= lnorm;
            } else {
                ldens = 0.0;
            }
        }
    }
//    if(!arma::is_finite(ldens)) {
//        Rprintf("Non-finite density in ldtskellam = %f x = %d l1 = %f l2 = %f LB = %d UB = %d\n", ldens, x, lambda1, lambda2, LB, UB);
//    }
    return ldens;
}

// truncated Skellam sampler
int rtskellam_cpp(double lambda1, double lambda2, char *str, sitmo::prng &eng, int LB = 0, int UB = -1) {
    
    // if UB < LB then returns untruncated density
    
    // check in this setting:
    if(LB > UB) {
        stop("rdt LB > UB\n");
    }
    
    // declare variables
    int x = 0;
    double mx = sitmo::prng::max();
    
    // if bounds are the same, then return the
    // only viable value
    if(LB == UB) {
        x = LB;
    } else {
        // draw untruncated sample
        x = rpois_cpp(lambda1, eng) - rpois_cpp(lambda2, eng);
        
        // rejection sample if necessary
        if(UB > LB) {
            int k = 0;
            int ntries = 1000;
            while((x < LB || x > UB) && k < ntries) {
                x = rpois_cpp(lambda1, eng) - rpois_cpp(lambda2, eng);
                k++;
            }
            // if rejection sampling doesn't work
            // then try inverse transform sampling
            if(k == ntries) {
                // Rprintf("LB = %d UB = %d\n", LB, UB);
                double u = eng() / mx;
                k = 0;
                double xdens = exp(ldtskellam_cpp(LB + k, lambda1, lambda2, str, LB, UB, 0));
                if(!arma::is_finite(xdens)) xdens = 0.0;
                while(u > xdens && k < (UB - LB)) {
                    k++;
                    xdens += exp(ldtskellam_cpp(LB + k, lambda1, lambda2, str, LB, UB, 0));
                    if(!arma::is_finite(xdens)) xdens += 0.0;
                }
                if(k == (UB - LB + 1) && u > xdens) {
                    stop("Something wrong in truncated Skellam sampling\n");
                }
                x = LB + k;
            }
        }
    }
    return x;
}

// simulation model
void discreteStochModel(int ipart, arma::vec &pars, int tstart, int tstop, 
                        arma::imat &u1_moves, std::vector<arma::icube> &u1, arma::icube &u1_day, arma::icube &u1_night,
                        arma::imat &N_day, arma::imat &N_night, arma::mat &pinf, arma::imat &origE, arma::mat &C, sitmo::prng &eng) {
    
    // u1_moves is a matrix with columns: LADfrom, LADto
    // u1 is a 3D array with dimensions: nclasses x nages x nmoves
    //          each row of u1 must match u1_moves
    // u1_day/night are 3D arrays with dimensions: nclasses x nages x nlads
    // N_day/night are nlad x nage matrices of population counts
    // pinf is nages x nlads auxiliary matrix
    // origE is nages x nmoves auxiliary matrix
    
    // set up auxiliary matrix for counts
    arma::uword nclasses = (arma::uword) u1[ipart].n_rows;
    arma::uword nages = (arma::uword) u1[ipart].n_cols;
    int k, n;
    arma::uword i, j, l;
    
    // reconstruct day/night counts
    u1_day.zeros();
    u1_night.zeros();
    for(i = 0; i < u1_moves.n_rows; i++) {
        for(j = 0; j < nages; j++) {
            for(l = 0; l < nclasses; l++) {
                u1_day(l, j, (arma::uword) u1_moves(i, 1) - 1) += u1[ipart](l, j, i);
                u1_night(l, j, (arma::uword) u1_moves(i, 0) - 1) += u1[ipart](l, j, i);
            }
        }
    }
    
    // extract parameters
    double nu = pars(0);
    double nuA = pars(1);
    arma::vec probE(nages);
    arma::vec probEP(nages);
    arma::vec probA(nages);
    arma::vec probP(nages);
    arma::vec probI1(nages);
    arma::vec probI1H(nages);
    arma::vec probI1D(nages);
    arma::vec probI2(nages);
    arma::vec probH(nages);
    arma::vec probHD(nages);
    
    for(j = 0; j < nages; j++) {
        probE(j) = pars(j + 2);
        probEP(j) = pars(j + nages + 2);
        probA(j) = pars(j + 2 * nages + 2);
        probP(j) = pars(j + 3 * nages + 2);
        probI1(j) = pars(j + 4 * nages + 2);
        probI1H(j) = pars(j + 5 * nages + 2);
        probI1D(j) = pars(j + 6 * nages + 2);
        probI2(j) = pars(j + 7 * nages + 2);
        probH(j) = pars(j + 8 * nages + 2);
        probHD(j) = pars(j + 9 * nages + 2);
    }
    
    // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
    //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
    
    int tcurr = 0;
    
    // set up vector of number of infectives
    arma::vec uinf(nages);
    
    // set up transmission rates
    arma::mat beta(nages, 1);
    
    // set up auxiliary variables
    arma::vec mprobsE(3);
    arma::vec mprobsI1(4);
    arma::vec mprobsH(3);
    tcurr++;
    tstart++;
    
    while(tstart <= tstop) {
        
        // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
        //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
        
        // save current number of infectives for later transitions
        for(i = 0; i < u1_moves.n_rows; i++) {
            for(j = 0; j < nages; j++) {
                origE(j, i) = u1[ipart](1, j, i);
            }
        }
        
        // transmission probabilities (day), loop over LADs
        for(i = 0; i < u1_day.n_slices; i++) {
            
            // update infective counts for rate
            for(j = 0; j < nages; j++) {
                uinf(j) = (double) nuA * u1_day(2, j, i) + nu * (u1_day(4, j, i) + u1_day(5, j, i) + u1_day(7, j, i));
            }
            
            // SE
            beta = 0.7 * C * (uinf / N_day.col(i));
            for(j = 0; j < nages; j++) {
                pinf(j, i) = 1.0 - exp(-beta(j, 0));
                pinf(j, i) = (pinf(j, i) < 0.0 ? 0.0:pinf(j, i));
                pinf(j, i) = (pinf(j, i) > 1.0 ? 1.0:pinf(j, i));
            }
        }
        // transmission events (day), loop over network
        for(i = 0; i < u1_moves.n_rows; i++) {
            for(j = 0; j < nages; j++) {
                k = rbinom_cpp(u1[ipart](0, j, i), pinf(j, (arma::uword) u1_moves(i, 1) - 1), eng);
                u1[ipart](0, j, i) -= k;
                u1[ipart](1, j, i) += k;
                u1_day(0, j, (arma::uword) u1_moves(i, 1) - 1) -= k;
                u1_day(1, j, (arma::uword) u1_moves(i, 1) - 1) += k;
                u1_night(0, j, (arma::uword) u1_moves(i, 0) - 1) -= k;
                u1_night(1, j, (arma::uword) u1_moves(i, 0) - 1) += k;
            }
        }
        
        // transmission probabilities (night), loop over LADs
        for(i = 0; i < u1_night.n_slices; i++) {
            
            // update infective counts for rate
            for(j = 0; j < nages; j++) {
                uinf(j) = (double) nuA * u1_night(2, j, i) + nu * (u1_night(4, j, i) + u1_night(5, j, i) + u1_night(7, j, i));
            }
            
            // SE
            beta = 0.3 * C * (uinf / N_night.col(i));
            for(j = 0; j < nages; j++) {
                pinf(j, i) = 1.0 - exp(-beta(j, 0));
                pinf(j, i) = (pinf(j, i) < 0.0 ? 0.0:pinf(j, i));
                pinf(j, i) = (pinf(j, i) > 1.0 ? 1.0:pinf(j, i));
            }
        }
        // transmission events (night), loop over network
        for(i = 0; i < u1_moves.n_rows; i++) {
            for(j = 0; j < nages; j++) {
                k = rbinom_cpp(u1[ipart](0, j, i), pinf(j, (arma::uword) u1_moves(i, 0) - 1), eng);
                u1[ipart](0, j, i) -= k;
                u1[ipart](1, j, i) += k;
                u1_day(0, j, (arma::uword) u1_moves(i, 1) - 1) -= k;
                u1_day(1, j, (arma::uword) u1_moves(i, 1) - 1) += k;
                u1_night(0, j, (arma::uword) u1_moves(i, 0) - 1) -= k;
                u1_night(1, j, (arma::uword) u1_moves(i, 0) - 1) += k;
            }
        }
        
        // conduct transition moves, loop over age classes                
        for(j = 0; j < nages; j++) {
            
            // (These probs could be pre-calculated an stored as matrix as long
            // as rmultinom can use rows of matrices as inputs?)
            
            // transition probs out of hospital
            mprobsH(0) = probH(j) * (1.0 - probHD(j));
            mprobsH(1) = probH(j) * probHD(j);
            mprobsH(2) = 1.0 - probH(j);
            mprobsH = mprobsH / sum(mprobsH);
            
            // transition probs out of I1
            mprobsI1(0) = probI1(j) * probI1H(j);
            mprobsI1(1) = probI1(j) * (1.0 - probI1D(j) - probI1H(j));
            mprobsI1(2) = probI1(j) * probI1D(j);
            mprobsI1(3) = 1.0 - probI1(j);
            mprobsI1 = mprobsI1 / sum(mprobsI1);
          
            // transition probs out of E
            mprobsE(0) = probE(j) * (1.0 - probEP(j));
            mprobsE(1) = probE(j) * probEP(j);
            mprobsE(2) = 1.0 - probE(j);
            mprobsE = mprobsE / sum(mprobsE);
            
            // loop over network
            for(i = 0; i < u1_moves.n_rows; i++) {
                
                // H out
                n = u1[ipart](9, j, i);
                if(n > 0) {
                    for(l = 0; l < n; l++) {
                        k = rmultinom_cpp(mprobsH, eng);
                        u1[ipart](9, j, i) -= (k < 2 ? 1:0);
                        u1[ipart](10, j, i) += (k == 0 ? 1:0);
                        u1[ipart](11, j, i) += (k == 1 ? 1:0);
                        u1_day(9, j, (arma::uword) u1_moves(i, 1) - 1) -= (k < 2 ? 1:0);
                        u1_day(10, j, (arma::uword) u1_moves(i, 1) - 1) += (k == 0 ? 1:0);
                        u1_day(11, j, (arma::uword) u1_moves(i, 1) - 1) += (k == 1 ? 1:0);
                        u1_night(9, j, (arma::uword) u1_moves(i, 0) - 1) -= (k < 2 ? 1:0);
                        u1_night(10, j, (arma::uword) u1_moves(i, 0) - 1) += (k == 0 ? 1:0);
                        u1_night(11, j, (arma::uword) u1_moves(i, 0) - 1) += (k == 1 ? 1:0);
                    }
                }
                
                // I2RI
                k = rbinom_cpp(u1[ipart](7, j, i), probI2(j), eng);
                u1[ipart](7, j, i) -= k;
                u1[ipart](8, j, i) += k;
                u1_day(7, j, (arma::uword) u1_moves(i, 1) - 1) -= k;
                u1_day(8, j, (arma::uword) u1_moves(i, 1) - 1) += k;
                u1_night(7, j, (arma::uword) u1_moves(i, 0) - 1) -= k;
                u1_night(8, j, (arma::uword) u1_moves(i, 0) - 1) += k;
                
                // I1 out
                n = u1[ipart](5, j, i);
                if(n > 0) {
                    for(l = 0; l < n; l++) {
                        k = rmultinom_cpp(mprobsI1, eng);
                        u1[ipart](5, j, i) -= (k < 3 ? 1:0);
                        u1[ipart](9, j, i) += (k == 0 ? 1:0);
                        u1[ipart](7, j, i) += (k == 1 ? 1:0);
                        u1[ipart](6, j, i) += (k == 2 ? 1:0);
                        u1_day(5, j, (arma::uword) u1_moves(i, 1) - 1) -= (k < 3 ? 1:0);
                        u1_day(9, j, (arma::uword) u1_moves(i, 1) - 1) += (k == 0 ? 1:0);
                        u1_day(7, j, (arma::uword) u1_moves(i, 1) - 1) += (k == 1 ? 1:0);
                        u1_day(6, j, (arma::uword) u1_moves(i, 1) - 1) += (k == 2 ? 1:0);
                        u1_night(5, j, (arma::uword) u1_moves(i, 0) - 1) -= (k < 3 ? 1:0);
                        u1_night(9, j, (arma::uword) u1_moves(i, 0) - 1) += (k == 0 ? 1:0);
                        u1_night(7, j, (arma::uword) u1_moves(i, 0) - 1) += (k == 1 ? 1:0);
                        u1_night(6, j, (arma::uword) u1_moves(i, 0) - 1) += (k == 2 ? 1:0);
                    }
                }
                
                // PI1
                k = rbinom_cpp(u1[ipart](4, j, i), probP(j), eng);
                u1[ipart](4, j, i) -= k;
                u1[ipart](5, j, i) += k;
                u1_day(4, j, (arma::uword) u1_moves(i, 1) - 1) -= k;
                u1_day(5, j, (arma::uword) u1_moves(i, 1) - 1) += k;
                u1_night(4, j, (arma::uword) u1_moves(i, 0) - 1) -= k;
                u1_night(5, j, (arma::uword) u1_moves(i, 0) - 1) += k;
                
                // ARA
                k = rbinom_cpp(u1[ipart](2, j, i), probA(j), eng);
                u1[ipart](2, j, i) -= k;
                u1[ipart](3, j, i) += k;
                u1_day(2, j, (arma::uword) u1_moves(i, 1) - 1) -= k;
                u1_day(3, j, (arma::uword) u1_moves(i, 1) - 1) += k;
                u1_night(2, j, (arma::uword) u1_moves(i, 0) - 1) -= k;
                u1_night(3, j, (arma::uword) u1_moves(i, 0) - 1) += k;
                
                // E out
                if(origE(j, i) > 0) {
                    for(l = 0; l < origE(j, i); l++) {
                        k = rmultinom_cpp(mprobsE, eng);
                        u1[ipart](1, j, i) -= (k < 2 ? 1:0);
                        u1[ipart](2, j, i) += (k == 0 ? 1:0);
                        u1[ipart](4, j, i) += (k == 1 ? 1:0);
                        u1_day(1, j, (arma::uword) u1_moves(i, 1) - 1) -= (k < 2 ? 1:0);
                        u1_day(2, j, (arma::uword) u1_moves(i, 1) - 1) += (k == 0 ? 1:0);
                        u1_day(4, j, (arma::uword) u1_moves(i, 1) - 1) += (k == 1 ? 1:0);
                        u1_night(1, j, (arma::uword) u1_moves(i, 0) - 1) -= (k < 2 ? 1:0);
                        u1_night(2, j, (arma::uword) u1_moves(i, 0) - 1) += (k == 0 ? 1:0);
                        u1_night(4, j, (arma::uword) u1_moves(i, 0) - 1) += (k == 1 ? 1:0);
                    }
                }
            }
        }
        
        // update time 
        tcurr++;
        tstart++;
    }
    return;
}

//// [[Rcpp::export]]
//void testsitmo () {
//    // check RNGs
//    int ncores = 2;
//    arma::vec seeds(ncores);
//    for(int i = 0; i < ncores; i++) {
//        seeds(i) = R::rnorm(0.0, 100.0);
//    }
//    double mx = sitmo::prng::max();
//    Rprintf("Without adjusting seed state:\n");
//    for(int j = 0; j < 2; j++) {
//#pragma omp parallel for default(none) shared(j, ncores, seeds, mx)
//        for(int i = 0; i < ncores; i++) {
//            uint32_t coreseed = static_cast<uint32_t>(seeds((arma::uword) i));
//            sitmo::prng eng(coreseed);
//            Rprintf("j = %d i = %d runif = %f\n", j, i, eng() / mx);
//        }
//    }
//    Rprintf("With adjusting seed state:\n");
//    for(int j = 0; j < 2; j++) {
//#pragma omp parallel for default(none) shared(j, ncores, seeds, mx)
//        for(int i = 0; i < ncores; i++) {
//            uint32_t coreseed = static_cast<uint32_t>(seeds((arma::uword) i));
//            sitmo::prng eng(coreseed);
//            Rprintf("j = %d i = %d runif = %f\n", j, i, eng() / mx);
//            seeds((arma::uword) i) = eng();
//        }
//    }
//    return;
//}

//// [[Rcpp::export]]
//List testRNGs () {
//    // check RNGs
//    uint32_t coreseed = static_cast<uint32_t>(R::rnorm(0.0, 100.0));
//    sitmo::prng eng(coreseed);
//    arma::ivec resBin (10000);
//    arma::ivec resPois (10000);
//    arma::ivec resMulti (10000);
//    arma::vec probs = {0.1, 0.2, 0.05, 0.35, 0.05, 0.05, 0.1, 0.1};
//    for(arma::uword i = 0; i < 10000; i++) {
//        resBin(i) = rbinom_cpp(30, 0.89, eng);
//        resPois(i) = rpois_cpp(0.8, eng);
//        resMulti(i) = rmultinom_cpp(probs, eng);
//    }
//    return List::create(Named("bin") = resBin, _["pois"] = resPois, _["multi"] = resMulti);
//}

// [[Rcpp::export]]
List PF_cpp (arma::vec pars, arma::mat C, arma::imat data, arma::uword nclasses, arma::uword nages, arma::uword nlads, arma::imat u1_moves, 
         arma::icube u1_comb, arma::uword ndays, arma::uword npart, int MD, double a1, double a2, double b, double a_dis, 
         double b_dis, int saveAll, int PF, int ncores) {
    
    // set counters
    arma::uword i, j, l, k, t;
    int tempLB = 0;
    
    // split u1 up into different LADs
    std::vector<arma::icube> u1(npart);
    std::vector<arma::icube> u1_new(npart);
    arma::icube u1_night_full(nclasses, nages, nlads); u1_night_full.zeros();
    arma::icube u1_night_reduced(2, nages, nlads); u1_night_reduced.zeros();
    arma::imat N_day(nages, nlads); N_day.zeros();
    arma::imat N_night(nages, nlads); N_night.zeros();
    for(i = 0; i < u1_moves.n_rows; i++) {
        for(j = 0; j < nages; j++) {
            for(l = 0; l < nclasses; l++) {
                N_day(j, (arma::uword) u1_moves(i, 1) - 1) += u1_comb(l, j, i);
                N_night(j, (arma::uword) u1_moves(i, 0) - 1) += u1_comb(l, j, i);
            }
        }
    }
    for(i = 0; i < npart; i++) {
        u1[i] = u1_comb;
        u1_new[i] = u1_comb;
    }
    
    // adjust discrepancy parameters to match cohort to national
    a_dis = a_dis / ((double) u1_moves.n_rows);
    
    // set up weight vector
    arma::vec weights (npart);
    arma::ivec inds(npart);
    double wnorm = 0.0;
    
    // set up auxiliary objects    
    arma::ivec obsInc (data.n_cols); obsInc.zeros();
    
    // check which output required
    List out (npart * (ndays + 1));
    if(saveAll != 0) {
        if(saveAll == 1) {
            for(i = 0; i < npart; i++) {
                // extract just counts for DI and DH
                u1_night_reduced.zeros();
                for(l = 0; l < u1_moves.n_rows; l++) {
                    for(j = 0; j < nages; j++) {
                        u1_night_reduced(0, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](6, j, l);
                        u1_night_reduced(1, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](11, j, l);
                    }
                }
                out[i] = u1_night_reduced;
            }
        } else {
            for(i = 0; i < npart; i++) {
                u1_night_full.zeros();
                for(l = 0; l < u1_moves.n_rows; l++) {
                    for(j = 0; j < nages; j++) {
                        for(k = 0; k < nclasses; k++) {
                            u1_night_full(k, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](k, j, l);
                        }
                    }
                }
                out[i] = u1_night_full;
            }
        }
    }
    
    // initialise timer
    Timer timer;
    int timer_cnt = 0;
    double prev_time = 0.0;
    
    // sample seeds to set up thread-safe PRNGs
#ifdef _OPENMP
    omp_set_num_threads(ncores);
#endif
    arma::vec seeds(ncores);
    for(i = 0; i < ncores; i++) {
        seeds(i) = R::rnorm(0.0, 100.0);
    }
    uint32_t coreseedSerial = static_cast<uint32_t>(R::rnorm(0.0, 100.0));
    sitmo::prng engSerial(coreseedSerial);
    
    // loop over time
    double ll = 0.0;
    for(t = 0; t < ndays; t++) {
            
        // check for interrupt
        R_CheckUserInterrupt();
        
        // extract data
        if(PF == 1) obsInc = data.row(t).t();
        
        // loop over particles
#ifdef _OPENMP
#pragma omp parallel for default(none) private(j, l, tempLB) shared(seeds, npart, u1_moves, nages, nclasses, nlads, data, C, N_night, N_day, u1, u1_new, t, pars, weights, MD, a_dis, b_dis, a1, a2, b, obsInc, PF)
#endif
        for(i = 0; i < npart; i++) {
    
            // set up print string for debugging
            char str1[80];
	        std::strcpy(str1, "setup");
 
            // set up thread-safe RNG
            uint32_t coreseed = static_cast<uint32_t>(seeds(0));
#ifdef _OPENMP
            coreseed = static_cast<uint32_t>(seeds((arma::uword) omp_get_thread_num()));
#endif
            sitmo::prng eng(coreseed);
        
            // set up auxiliary objects
            arma::imat DHinc (nages, u1_moves.n_rows); DHinc.zeros();
            arma::imat DIinc (nages, u1_moves.n_rows); DIinc.zeros();
            
            arma::ivec Dtempinc (data.n_cols); Dtempinc.zeros();
            arma::vec tempdens (data.n_cols); tempdens.zeros();
            
            arma::mat pinf(nages, nlads); pinf.zeros();
            arma::imat origE(nages, u1_moves.n_rows); origE.zeros();
            arma::icube u1_day(nclasses, nages, nlads); u1_day.zeros();
            arma::icube u1_night(nclasses, nages, nlads); u1_night.zeros();
            
            // run model and return u1
            discreteStochModel((int) i, pars, t - 1, t, u1_moves, u1_new, u1_day, u1_night, N_day, N_night, pinf, origE, C, eng);
            
            // cols: c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")
            //          0,   1,   2,   3,    4,   5,    6,    7,    8,    9,   10,   11
            
            // set weights
            weights(i) = 0.0;
            
            // adjust states according to model discrepancy
            if(MD == 1) {
                // create auxiliary objects
                arma::imat RHinc (nages, u1_moves.n_rows); RHinc.zeros();
                arma::imat Hinc (nages, u1_moves.n_rows); Hinc.zeros();
                arma::imat RIinc (nages, u1_moves.n_rows); RIinc.zeros();
                arma::imat I2inc (nages, u1_moves.n_rows); I2inc.zeros();
                arma::imat I1inc (nages, u1_moves.n_rows); I1inc.zeros();
                arma::imat Pinc (nages, u1_moves.n_rows); Pinc.zeros();
                arma::imat RAinc (nages, u1_moves.n_rows); RAinc.zeros();
                arma::imat Ainc (nages, u1_moves.n_rows); Ainc.zeros();
                arma::imat Einc (nages, u1_moves.n_rows); Einc.zeros();
                
                // DH (MD on incidence)
                std::strcpy(str1, "DHinc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        DHinc(j, l) = u1_new[i](11, j, l) - u1[i](11, j, l);
                        DHinc(j, l) += rtskellam_cpp(
                            a_dis + b_dis * DHinc(j, l),
                            a_dis + b_dis * DHinc(j, l),
                            str1,
                            eng,
                            -DHinc(j, l),
                            u1[i](9, j, l) - DHinc(j, l)
                        );
                        if(DHinc(j, l) < 0) Rprintf("DHinc = %d j = %d l = %d\n", DHinc(j, l), j, l);
                    }
                }
                // update counts
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        u1_new[i](11, j, l) = u1[i](11, j, l) + DHinc(j, l);
                        if(u1_new[i](11, j, l) < 0) Rprintf("DH = %d j = %d l = %d\n", u1_new[i](11, j, l), j, l);
                    }
                }
                
                // RH given DH (MD on incidence)
                std::strcpy(str1, "RHinc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        RHinc(j, l) = u1_new[i](10, j, l) - u1[i](10, j, l);
                        RHinc(j, l) += rtskellam_cpp(
                            a_dis + b_dis * RHinc(j, l),
                            a_dis + b_dis * RHinc(j, l),
                            str1,
                            eng,
                            -RHinc(j, l),
                            u1[i](9, j, l) - DHinc(j, l) - RHinc(j, l)
                        );
                        if(RHinc(j, l) < 0) Rprintf("RHinc = %d j = %d l = %d\n", RHinc(j, l), j, l);
                    }
                }
                // update counts
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        u1_new[i](10, j, l) = u1[i](10, j, l) + RHinc(j, l);
                        if(u1_new[i](10, j, l) < 0) Rprintf("RH = %d j = %d l = %d\n", u1_new[i](10, j, l), j, l);
                    }
                }
                    
                // H given later
                std::strcpy(str1, "Hinc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        tempLB = -u1_new[i](9, j, l) + u1[i](9, j, l) - DHinc(j, l) - RHinc(j, l);
                        tempLB = (tempLB > -u1_new[i](9, j, l) ? tempLB:(-u1_new[i](9, j, l)));
                        u1_new[i](9, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](9, j, l),
                            a_dis + b_dis * u1_new[i](9, j, l),
                            str1,
                            eng,
                            tempLB,
                            u1[i](5, j, l) - u1_new[i](9, j, l) + u1[i](9, j, l) - DHinc(j, l) - RHinc(j, l)
                        );
                        Hinc(j, l) = u1_new[i](9, j, l) - u1[i](9, j, l) + DHinc(j, l) + RHinc(j, l);
                        if(Hinc(j, l) < 0) Rprintf("Hinc = %d j = %d l = %d\n", Hinc(j, l), j, l);
                        if(u1_new[i](9, j, l) < 0) Rprintf("H = %d j = %d l = %d\n", u1_new[i](9, j, l), j, l);
                    }
                }
                    
                // DI given H (MD on incidence)
                std::strcpy(str1, "DIinc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        DIinc(j, l) = u1_new[i](6, j, l) - u1[i](6, j, l);
                        DIinc(j, l) += rtskellam_cpp(
                            a_dis + b_dis * DIinc(j, l),
                            a_dis + b_dis * DIinc(j, l),
                            str1,
                            eng,
                            -DIinc(j, l),
                            u1[i](5, j, l) - Hinc(j, l) - DIinc(j, l)
                        );
                        if(DIinc(j, l) < 0) Rprintf("DIinc = %d j = %d l = %d\n", DIinc(j, l), j, l);
                    }
                }
                // update counts
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        u1_new[i](6, j, l) = u1[i](6, j, l) + DIinc(j, l);
                        if(u1_new[i](6, j, l) < 0) Rprintf("DI = %d j = %d l = %d\n", u1_new[i](6, j, l), j, l);
                    }
                }
                    
                // RI (MD on incidence)
                std::strcpy(str1, "RIinc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        RIinc(j, l) = u1_new[i](8, j, l) - u1[i](8, j, l);
                        RIinc(j, l) += rtskellam_cpp(
                            a_dis + b_dis * RIinc(j, l),
                            a_dis + b_dis * RIinc(j, l),
                            str1,
                            eng,
                            -RIinc(j, l),
                            u1[i](7, j, l) - RIinc(j, l)
                        );
                        if(RIinc(j, l) < 0) Rprintf("RIinc = %d j = %d l = %d\n", RIinc(j, l), j, l);
                    }
                }
                // update counts
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        u1_new[i](8, j, l) = u1[i](8, j, l) + RIinc(j, l);
                        if(u1_new[i](8, j, l) < 0) Rprintf("RI = %d j = %d l = %d\n", u1_new[i](8, j, l), j, l);
                    }
                }
                    
                // I2 given later
                std::strcpy(str1, "I2inc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        tempLB = -u1_new[i](7, j, l) + u1[i](7, j, l) - RIinc(j, l);
                        tempLB = (tempLB > -u1_new[i](7, j, l) ? tempLB:(-u1_new[i](7, j, l)));
                        u1_new[i](7, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](7, j, l),
                            a_dis + b_dis * u1_new[i](7, j, l),
                            str1,
                            eng,
                            tempLB,
                            u1[i](5, j, l) - Hinc(j, l) - DIinc(j, l) - u1_new[i](7, j, l) + u1[i](7, j, l) - RIinc(j, l)
                        );
                        I2inc(j, l) = u1_new[i](7, j, l) - u1[i](7, j, l) + RIinc(j, l);
                        if(I2inc(j, l) < 0) Rprintf("I2inc = %d j = %d l = %d\n", I2inc(j, l), j, l);
                        if(u1_new[i](7, j, l) < 0) Rprintf("I2 = %d j = %d l = %d\n", u1_new[i](7, j, l), j, l);
                    }
                }
                
                // I1 given later
                std::strcpy(str1, "I1inc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        tempLB = -u1_new[i](5, j, l) + u1[i](5, j, l) - I2inc(j, l) - Hinc(j, l) - DIinc(j, l);
                        tempLB = (tempLB > -u1_new[i](5, j, l) ? tempLB:(-u1_new[i](5, j, l)));
                        u1_new[i](5, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](5, j, l),
                            a_dis + b_dis * u1_new[i](5, j, l),
                            str1,
                            eng,
                            tempLB,
                            u1[i](4, j, l) - u1_new[i](5, j, l) + u1[i](5, j, l) - I2inc(j, l) - Hinc(j, l) - DIinc(j, l)
                        );
                        I1inc(j, l) = u1_new[i](5, j, l) - u1[i](5, j, l) + I2inc(j, l) + Hinc(j, l) + DIinc(j, l);
                        if(I1inc(j, l) < 0) Rprintf("I1inc = %d j = %d l = %d\n", I1inc(j, l), j, l);
                        if(u1_new[i](5, j, l) < 0) Rprintf("I1 = %d j = %d l = %d\n", u1_new[i](5, j, l), j, l);
                    }
                }
                
                // P given later
                std::strcpy(str1, "Pinc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        tempLB = -u1_new[i](4, j, l) + u1[i](4, j, l) - I1inc(j, l);
                        tempLB = (tempLB > -u1_new[i](4, j, l) ? tempLB:(-u1_new[i](4, j, l)));
                        u1_new[i](4, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](4, j, l),
                            a_dis + b_dis * u1_new[i](4, j, l),
                            str1,
                            eng,
                            tempLB,
                            u1[i](1, j, l) - u1_new[i](4, j, l) + u1[i](4, j, l) - I1inc(j, l)
                        );
                        Pinc(j, l) = u1_new[i](4, j, l) - u1[i](4, j, l) + I1inc(j, l);
                        if(Pinc(j, l) < 0) Rprintf("Pinc = %d j = %d l = %d\n", Pinc(j, l), j, l);
                        if(u1_new[i](4, j, l) < 0) Rprintf("P = %d j = %d l = %d\n", u1_new[i](4, j, l), j, l);
                    }
                }
                
                // RA (MD on incidence)
                std::strcpy(str1, "RAinc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        RAinc(j, l) = u1_new[i](3, j, l) - u1[i](3, j, l);
                        RAinc(j, l) += rtskellam_cpp(
                            a_dis + b_dis * RAinc(j, l),
                            a_dis + b_dis * RAinc(j, l),
                            str1,
                            eng,
                            -RAinc(j, l),
                            u1[i](2, j, l) - RAinc(j, l)
                        );
                        if(RAinc(j, l) < 0) Rprintf("RAinc = %d j = %d l = %d\n", RAinc(j, l), j, l);
                    }
                }
                // update counts
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        u1_new[i](3, j, l) = u1[i](3, j, l) + RAinc(j, l);
                        if(u1_new[i](3, j, l) < 0) Rprintf("RH = %d j = %d l = %d\n", u1_new[i](3, j, l), j, l);
                    }
                }
                
                // A given later
                std::strcpy(str1, "Ainc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        tempLB = -u1_new[i](2, j, l) + u1[i](2, j, l) - RAinc(j, l);
                        tempLB = (tempLB > -u1_new[i](2, j, l) ? tempLB:(-u1_new[i](2, j, l)));
                        u1_new[i](2, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](2, j, l),
                            a_dis + b_dis * u1_new[i](2, j, l),
                            str1,
                            eng,
                            tempLB,
                            u1[i](1, j, l) - Pinc(j, l) - u1_new[i](2, j, l) + u1[i](2, j, l) - RAinc(j, l)
                        );
                        Ainc(j, l) = u1_new[i](2, j, l) - u1[i](2, j, l) + RAinc(j, l);
                        if(Ainc(j, l) < 0) Rprintf("Ainc = %d j = %d l = %d\n", Ainc(j, l), j, l);
                        if(u1_new[i](2, j, l) < 0) Rprintf("A = %d j = %d l = %d\n", u1_new[i](2, j, l), j, l);
                    }
                }
                
                // E given later
                std::strcpy(str1, "Einc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        tempLB = -u1_new[i](1, j, l) + u1[i](1, j, l) - Ainc(j, l) - Pinc(j, l);
                        tempLB = (tempLB > -u1_new[i](1, j, l) ? tempLB:(-u1_new[i](1, j, l)));
                        u1_new[i](1, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](1, j, l),
                            a_dis + b_dis * u1_new[i](1, j, l),
                            str1,
                            eng,
                            tempLB,
                            u1[i](0, j, l) - u1_new[i](1, j, l) + u1[i](1, j, l) - Ainc(j, l) - Pinc(j, l)
                        );
                        Einc(j, l) = u1_new[i](1, j, l) - u1[i](1, j, l) + Ainc(j, l) + Pinc(j, l);
                        if(Einc(j, l) < 0) Rprintf("Einc = %d j = %d l = %d\n", Einc(j, l), j, l);
                        if(u1_new[i](1, j, l) < 0) Rprintf("E = %d j = %d l = %d\n", u1_new[i](1, j, l), j, l);
                    }
                }
                
                // S given later
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        u1_new[i](0, j, l) = u1[i](0, j, l) - Einc(j, l);
                        if(u1_new[i](0, j, l) < 0) Rprintf("S = %d j = %d l = %d\n", u1_new[i](0, j, l), j, l);
                    }
                }
            } else {
                if(PF == 1) {
                    // DH incidence
                    for(j = 0; j < nages; j++) {
                        for(l = 0; l < u1_moves.n_rows; l++) {
                            DHinc(j, l) = u1_new[i](11, j, l) - u1[i](11, j, l);
                            if(DHinc(j, l) < 0) Rprintf("DHinc = %d j = %d l = %d\n", DHinc(j, l), j, l);
                        }
                    }
                    
                    // DI incidence
                    for(j = 0; j < nages; j++) {
                        for(l = 0; l < u1_moves.n_rows; l++) {
                            DIinc(j, l) = u1_new[i](6, j, l) - u1[i](6, j, l);
                            if(DIinc(j, l) < 0) Rprintf("DIinc = %d j = %d l = %d\n", DIinc(j, l), j, l);
                        }
                    }
                }
            }
            
            if(PF == 1) {
                // generate data in correct format for observation error weights
                Dtempinc.zeros();
                for(l = 0; l < u1_moves.n_rows; l++) {
                    for(j = 0; j < nages; j++) {
                        Dtempinc((arma::uword) j * nlads + u1_moves(l, 0) - 1) += DIinc(j, l);
                        Dtempinc((arma::uword) nlads * nages + j * nlads + u1_moves(l, 0) - 1) += DHinc(j, l);
                    }
                }
                
                // calculate log observation error weights
                for(l = 0; l < Dtempinc.n_elem; l++) {
                    tempdens(l) = ldtskellam_cpp(
                        obsInc(l) - Dtempinc(l),
                        a1 + b * Dtempinc(l),
                        a2 + b * Dtempinc(l),
                        str1,
                        -Dtempinc(l),
                        obsInc(l),
                        0
                    );
                    if(!arma::is_finite(tempdens(l)) && tempdens(l) >= 0.0) stop("Error in OE\n");
                }
                weights(i) += sum(tempdens);
            }
            
            // advance seed
            seeds((arma::uword) omp_get_thread_num()) = eng();
        }
        
        if(PF == 1) {
            // calculate log-likelihood contribution
            ll += log_sum_exp(weights, 1);
            
            // if zero likelihood then return
            if(!arma::is_finite(ll)) {
                if(saveAll == 0) {
                    return List::create(Named("ll") = ll);
                } else {
                    return List::create(Named("ll") = ll, _["particles"] = out);
                }
            }
            
            // normalise weights
            wnorm = log_sum_exp(weights, 0);
            weights = exp(weights - wnorm);
            weights = weights / sum(weights);
            
            // resample
            for(i = 0; i < npart; i++) {
                l = (arma::uword) rmultinom_cpp(weights, engSerial);
                u1[i] = u1_new[l];
            }
            // copy in order to pass by reference
            for(i = 0; i < npart; i++) {
                u1_new[i] = u1[i];
            }
        } else {
            for(i = 0; i < npart; i++) {
                u1[i] = u1_new[i];
            }
        }
        
        // save particles if necessary
        if(saveAll != 0) {
            if(saveAll == 1) {
                for(i = 0; i < npart; i++) {
                    // extract just counts for DI and DH
                    u1_night_reduced.zeros();
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        for(j = 0; j < nages; j++) {
                            u1_night_reduced(0, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](6, j, l);
                            u1_night_reduced(1, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](11, j, l);
                        }
                    }
                    out[i + npart * (t + 1)] = u1_night_reduced;
                }
            } else {
                for(i = 0; i < npart; i++) {
                    u1_night_full.zeros();
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        for(j = 0; j < nages; j++) {
                            for(k = 0; k < nclasses; k++) {
                                u1_night_full(k, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](k, j, l);
                            }
                        }
                    }
                    out[i + npart * (t + 1)] = u1_night_full;
                }
            }
        }
        
        //calculate block run time
        timer.step("");
        NumericVector res(timer);
        
        Rprintf("t = %d / %d time = %.2f secs \n", t + 1, ndays, (res[timer_cnt] / 1e9) - prev_time);
        
        //reset timer and acceptance rate counter
        prev_time = res[timer_cnt] / 1e9;
        timer_cnt++;
    }
    if(saveAll == 0) {
        return List::create(Named("ll") = ll);
    } else {
        if(PF == 1) {
            return List::create(Named("ll") = ll, _["particles"] = out);
        } else {
            return List::create(Named("particles") = out);
        }
    }
}

/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2014 Ross Ihaka and the R Core team.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/*  DESCRIPTION --> see below */


/* From http://www.netlib.org/specfun/ribesl	Fortran translated by f2c,...
 *	------------------------------=#----	Martin Maechler, ETH Zurich
 */
 
 /* Original code taken from R-4.1.3. source code, and then
 amended by TJ McKinley, University of Exeter, UK to remove warnings
 which cause havoc when calling in parallel through OpenMP and Rcpp. Similarly
 have extracted only the necessary functions required for a specific application
 for simplicity, and converted into standalone Rcpp. Introduced errors are thus mine.*/

/* modified version of bessel_i that accepts a work array instead of
   allocating one. */
double bessel_i_ex_new(double x, double alpha, double expo, double *bi)
{
    int nb, ncalc, ize;
    double na;

#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(x) || ISNAN(alpha)) return x + alpha;
#endif
    if (x < 0) {
	Rprintf("bessel_i - out-of-range");
	return ML_NAN;
    }
    ize = (int)expo;
    na = floor(alpha);
    if (alpha < 0) {
	/* Using Abramowitz & Stegun  9.6.2 & 9.6.6
	 * this may not be quite optimal (CPU and accuracy wise) */
	return(bessel_i_ex_new(x, -alpha, expo, bi) +
	       ((alpha == na) ? 0 :
		bessel_k_ex_new(x, -alpha, expo, bi) *
		((ize == 1)? 2. : 2.*exp(-2.*x))/M_PI * sinpi(-alpha)));
    }
    nb = 1 + (int)na;/* nb-1 <= alpha < nb */
    alpha -= (double)(nb-1);
    I_bessel_new(&x, &alpha, &nb, &ize, bi, &ncalc);
    if(ncalc != nb) {/* error input */
	if(ncalc < 0)
	    Rprintf("bessel_i(%g): ncalc (=%d) != nb (=%d); alpha=%g. Arg. out of range?\n",
			     x, ncalc, nb, alpha);
	else
	    Rprintf("bessel_i(%g,nu=%g): precision lost in result\n",
			     x, alpha+(double)nb-1);
    }
    x = bi[nb-1];
    return x;
}

static void I_bessel_new(double *x, double *alpha, int *nb,
		     int *ize, double *bi, int *ncalc)
{
/* -------------------------------------------------------------------

 This routine calculates Bessel functions I_(N+ALPHA) (X)
 for non-negative argument X, and non-negative order N+ALPHA,
 with or without exponential scaling.


 Explanation of variables in the calling sequence

 X     - Non-negative argument for which
	 I's or exponentially scaled I's (I*EXP(-X))
	 are to be calculated.	If I's are to be calculated,
	 X must be less than exparg_BESS (IZE=1) or xlrg_BESS_IJ (IZE=2),
	 (see bessel.h).
 ALPHA - Fractional part of order for which
	 I's or exponentially scaled I's (I*EXP(-X)) are
	 to be calculated.  0 <= ALPHA < 1.0.
 NB    - Number of functions to be calculated, NB > 0.
	 The first function calculated is of order ALPHA, and the
	 last is of order (NB - 1 + ALPHA).
 IZE   - Type.	IZE = 1 if unscaled I's are to be calculated,
		    = 2 if exponentially scaled I's are to be calculated.
 BI    - Output vector of length NB.	If the routine
	 terminates normally (NCALC=NB), the vector BI contains the
	 functions I(ALPHA,X) through I(NB-1+ALPHA,X), or the
	 corresponding exponentially scaled functions.
 NCALC - Output variable indicating possible errors.
	 Before using the vector BI, the user should check that
	 NCALC=NB, i.e., all orders have been calculated to
	 the desired accuracy.	See error returns below.


 *******************************************************************
 *******************************************************************

 Error returns

  In case of an error,	NCALC != NB, and not all I's are
  calculated to the desired accuracy.

  NCALC < 0:  An argument is out of range. For example,
     NB <= 0, IZE is not 1 or 2, or IZE=1 and ABS(X) >= EXPARG_BESS.
     In this case, the BI-vector is not calculated, and NCALC is
     set to MIN0(NB,0)-1 so that NCALC != NB.

  NB > NCALC > 0: Not all requested function values could
     be calculated accurately.	This usually occurs because NB is
     much larger than ABS(X).  In this case, BI[N] is calculated
     to the desired accuracy for N <= NCALC, but precision
     is lost for NCALC < N <= NB.  If BI[N] does not vanish
     for N > NCALC (because it is too small to be represented),
     and BI[N]/BI[NCALC] = 10**(-K), then only the first NSIG-K
     significant figures of BI[N] can be trusted.


 Intrinsic functions required are:

     DBLE, EXP, gamma_cody, INT, MAX, MIN, REAL, SQRT


 Acknowledgement

  This program is based on a program written by David J.
  Sookne (2) that computes values of the Bessel functions J or
  I of float argument and long order.  Modifications include
  the restriction of the computation to the I Bessel function
  of non-negative float argument, the extension of the computation
  to arbitrary positive order, the inclusion of optional
  exponential scaling, and the elimination of most underflow.
  An earlier version was published in (3).

 References: "A Note on Backward Recurrence Algorithms," Olver,
	      F. W. J., and Sookne, D. J., Math. Comp. 26, 1972,
	      pp 941-947.

	     "Bessel Functions of Real Argument and Integer Order,"
	      Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp
	      125-132.

	     "ALGORITHM 597, Sequence of Modified Bessel Functions
	      of the First Kind," Cody, W. J., Trans. Math. Soft.,
	      1983, pp. 242-245.

  Latest modification: May 30, 1989

  Modified by: W. J. Cody and L. Stoltz
	       Applied Mathematics Division
	       Argonne National Laboratory
	       Argonne, IL  60439
*/

    /*-------------------------------------------------------------------
      Mathematical constants
      -------------------------------------------------------------------*/
    const static double const__ = 1.585;

    /* Local variables */
    int nend, intx, nbmx, k, l, n, nstart;
    double pold, test,	p, em, en, empal, emp2al, halfx,
	aa, bb, cc, psave, plast, tover, psavel, sum, nu, twonu;

    /*Parameter adjustments */
    --bi;
    nu = *alpha;
    twonu = nu + nu;

    /*-------------------------------------------------------------------
      Check for X, NB, OR IZE out of range.
      ------------------------------------------------------------------- */
    if (*nb > 0 && *x >= 0. &&	(0. <= nu && nu < 1.) &&
	(1 <= *ize && *ize <= 2) ) {

	*ncalc = *nb;
	if(*ize == 1 && *x > exparg_BESS) {
	    for(k=1; k <= *nb; k++)
		bi[k]=ML_POSINF; /* the limit *is* = Inf */
	    return;
	}
	if(*ize == 2 && *x > xlrg_BESS_IJ) {
	    for(k=1; k <= *nb; k++)
		bi[k]= 0.; /* The limit exp(-x) * I_nu(x) --> 0 : */
	    return;
	}
	intx = (int) (*x);/* fine, since *x <= xlrg_BESS_IJ <<< LONG_MAX */
	if (*x >= rtnsig_BESS) { /* "non-small" x ( >= 1e-4 ) */
/* -------------------------------------------------------------------
   Initialize the forward sweep, the P-sequence of Olver
   ------------------------------------------------------------------- */
	    nbmx = *nb - intx;
	    n = intx + 1;
	    en = (double) (n + n) + twonu;
	    plast = 1.;
	    p = en / *x;
	    /* ------------------------------------------------
	       Calculate general significance test
	       ------------------------------------------------ */
	    test = ensig_BESS + ensig_BESS;
	    if (intx << 1 > nsig_BESS * 5) {
		test = sqrt(test * p);
	    } else {
		test /= R_pow_di(const__, intx);
	    }
	    if (nbmx >= 3) {
		/* --------------------------------------------------
		   Calculate P-sequence until N = NB-1
		   Check for possible overflow.
		   ------------------------------------------------ */
		tover = enten_BESS / ensig_BESS;
		nstart = intx + 2;
		nend = *nb - 1;
		for (k = nstart; k <= nend; ++k) {
		    n = k;
		    en += 2.;
		    pold = plast;
		    plast = p;
		    p = en * plast / *x + pold;
		    if (p > tover) {
			/* ------------------------------------------------
			   To avoid overflow, divide P-sequence by TOVER.
			   Calculate P-sequence until ABS(P) > 1.
			   ---------------------------------------------- */
			tover = enten_BESS;
			p /= tover;
			plast /= tover;
			psave = p;
			psavel = plast;
			nstart = n + 1;
			do {
			    ++n;
			    en += 2.;
			    pold = plast;
			    plast = p;
			    p = en * plast / *x + pold;
			}
			while (p <= 1.);

			bb = en / *x;
			/* ------------------------------------------------
			   Calculate backward test, and find NCALC,
			   the highest N such that the test is passed.
			   ------------------------------------------------ */
			test = pold * plast / ensig_BESS;
			test *= .5 - .5 / (bb * bb);
			p = plast * tover;
			--n;
			en -= 2.;
			nend = min0(*nb,n);
			for (l = nstart; l <= nend; ++l) {
			    *ncalc = l;
			    pold = psavel;
			    psavel = psave;
			    psave = en * psavel / *x + pold;
			    if (psave * psavel > test) {
				goto L90;
			    }
			}
			*ncalc = nend + 1;
L90:
			--(*ncalc);
			goto L120;
		    }
		}
		n = nend;
		en = (double)(n + n) + twonu;
		/*---------------------------------------------------
		  Calculate special significance test for NBMX > 2.
		  --------------------------------------------------- */
		test = R::fmax2(test,sqrt(plast * ensig_BESS) * sqrt(p + p));
	    }
	    /* --------------------------------------------------------
	       Calculate P-sequence until significance test passed.
	       -------------------------------------------------------- */
	    do {
		++n;
		en += 2.;
		pold = plast;
		plast = p;
		p = en * plast / *x + pold;
	    } while (p < test);

L120:
/* -------------------------------------------------------------------
 Initialize the backward recursion and the normalization sum.
 ------------------------------------------------------------------- */
	    ++n;
	    en += 2.;
	    bb = 0.;
	    aa = 1. / p;
	    em = (double) n - 1.;
	    empal = em + nu;
	    emp2al = em - 1. + twonu;
	    sum = aa * empal * emp2al / em;
	    nend = n - *nb;
	    if (nend < 0) {
		/* -----------------------------------------------------
		   N < NB, so store BI[N] and set higher orders to 0..
		   ----------------------------------------------------- */
		bi[n] = aa;
		nend = -nend;
		for (l = 1; l <= nend; ++l) {
		    bi[n + l] = 0.;
		}
	    } else {
		if (nend > 0) {
		    /* -----------------------------------------------------
		       Recur backward via difference equation,
		       calculating (but not storing) BI[N], until N = NB.
		       --------------------------------------------------- */

		    for (l = 1; l <= nend; ++l) {
			--n;
			en -= 2.;
			cc = bb;
			bb = aa;
			/* for x ~= 1500,  sum would overflow to 'inf' here,
			 * and the final bi[] /= sum would give 0 wrongly;
			 * RE-normalize (aa, sum) here -- no need to undo */
			if(nend > 100 && aa > 1e200) {
			    /* multiply by  2^-900 = 1.18e-271 */
			    cc	= ldexp(cc, -900);
			    bb	= ldexp(bb, -900);
			    sum = ldexp(sum,-900);
			}
			aa = en * bb / *x + cc;
			em -= 1.;
			emp2al -= 1.;
			if (n == 1) {
			    break;
			}
			if (n == 2) {
			    emp2al = 1.;
			}
			empal -= 1.;
			sum = (sum + aa * empal) * emp2al / em;
		    }
		}
		/* ---------------------------------------------------
		   Store BI[NB]
		   --------------------------------------------------- */
		bi[n] = aa;
		if (*nb <= 1) {
		    sum = sum + sum + aa;
		    goto L230;
		}
		/* -------------------------------------------------
		   Calculate and Store BI[NB-1]
		   ------------------------------------------------- */
		--n;
		en -= 2.;
		bi[n] = en * aa / *x + bb;
		if (n == 1) {
		    goto L220;
		}
		em -= 1.;
		if (n == 2)
		    emp2al = 1.;
		else
		    emp2al -= 1.;

		empal -= 1.;
		sum = (sum + bi[n] * empal) * emp2al / em;
	    }
	    nend = n - 2;
	    if (nend > 0) {
		/* --------------------------------------------
		   Calculate via difference equation
		   and store BI[N], until N = 2.
		   ------------------------------------------ */
		for (l = 1; l <= nend; ++l) {
		    --n;
		    en -= 2.;
		    bi[n] = en * bi[n + 1] / *x + bi[n + 2];
		    em -= 1.;
		    if (n == 2)
			emp2al = 1.;
		    else
			emp2al -= 1.;
		    empal -= 1.;
		    sum = (sum + bi[n] * empal) * emp2al / em;
		}
	    }
	    /* ----------------------------------------------
	       Calculate BI[1]
	       -------------------------------------------- */
	    bi[1] = 2. * empal * bi[2] / *x + bi[3];
L220:
	    sum = sum + sum + bi[1];

L230:
	    /* ---------------------------------------------------------
	       Normalize.  Divide all BI[N] by sum.
	       --------------------------------------------------------- */
	    if (nu != 0.)
		sum *= (Rf_gamma_cody(1. + nu) * pow(*x * .5, -nu));
	    if (*ize == 1)
		sum *= exp(-(*x));
	    aa = enmten_BESS;
	    if (sum > 1.)
		aa *= sum;
	    for (n = 1; n <= *nb; ++n) {
		if (bi[n] < aa)
		    bi[n] = 0.;
		else
		    bi[n] /= sum;
	    }
	    return;
	} else { /* small x  < 1e-4 */
	    /* -----------------------------------------------------------
	       Two-term ascending series for small X.
	       -----------------------------------------------------------*/
	    aa = 1.;
	    empal = 1. + nu;
#ifdef IEEE_754
	    /* No need to check for underflow */
	    halfx = .5 * *x;
#else
	    if (*x > enmten_BESS) */
		halfx = .5 * *x;
	    else
	    	halfx = 0.;
#endif
	    if (nu != 0.)
		aa = pow(halfx, nu) / Rf_gamma_cody(empal);
	    if (*ize == 2)
		aa *= exp(-(*x));
	    bb = halfx * halfx;
	    bi[1] = aa + aa * bb / empal;
	    if (*x != 0. && bi[1] == 0.)
		*ncalc = 0;
	    if (*nb > 1) {
		if (*x == 0.) {
		    for (n = 2; n <= *nb; ++n)
			bi[n] = 0.;
		} else {
		    /* -------------------------------------------------
		       Calculate higher-order functions.
		       ------------------------------------------------- */
		    cc = halfx;
		    tover = (enmten_BESS + enmten_BESS) / *x;
		    if (bb != 0.)
			tover = enmten_BESS / bb;
		    for (n = 2; n <= *nb; ++n) {
			aa /= empal;
			empal += 1.;
			aa *= cc;
			if (aa <= tover * empal)
			    bi[n] = aa = 0.;
			else
			    bi[n] = aa + aa * bb / empal;
			if (bi[n] == 0. && *ncalc > n)
			    *ncalc = n - 1;
		    }
		}
	    }
	}
    } else { /* argument out of range */
	*ncalc = min0(*nb,0) - 1;
    }
}
/* modified version of bessel_k that accepts a work array instead of
   allocating one. */
double bessel_k_ex_new(double x, double alpha, double expo, double *bk)
{
    int nb, ncalc, ize;

#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(x) || ISNAN(alpha)) return x + alpha;
#endif
    if (x < 0) {
	Rprintf("bessel_k - out-of-range");
	return ML_NAN;
    }
    ize = (int)expo;
    if(alpha < 0)
	alpha = -alpha;
    nb = 1+ (int)floor(alpha);/* nb-1 <= |alpha| < nb */
    alpha -= (double)(nb-1);
    K_bessel_new(&x, &alpha, &nb, &ize, bk, &ncalc);
    if(ncalc != nb) {/* error input */
      if(ncalc < 0)
	Rprintf("bessel_k(%g): ncalc (=%d) != nb (=%d); alpha=%g. Arg. out of range?\n",
			 x, ncalc, nb, alpha);
      else
	Rprintf("bessel_k(%g,nu=%g): precision lost in result\n",
			 x, alpha+(double)nb-1);
    }
    x = bk[nb-1];
    return x;
}

static void K_bessel_new(double *x, double *alpha, int *nb,
		     int *ize, double *bk, int *ncalc)
{
/*-------------------------------------------------------------------

  This routine calculates modified Bessel functions
  of the third kind, K_(N+ALPHA) (X), for non-negative
  argument X, and non-negative order N+ALPHA, with or without
  exponential scaling.

  Explanation of variables in the calling sequence

 X     - Non-negative argument for which
	 K's or exponentially scaled K's (K*EXP(X))
	 are to be calculated.	If K's are to be calculated,
	 X must not be greater than XMAX_BESS_K.
 ALPHA - Fractional part of order for which
	 K's or exponentially scaled K's (K*EXP(X)) are
	 to be calculated.  0 <= ALPHA < 1.0.
 NB    - Number of functions to be calculated, NB > 0.
	 The first function calculated is of order ALPHA, and the
	 last is of order (NB - 1 + ALPHA).
 IZE   - Type.	IZE = 1 if unscaled K's are to be calculated,
		    = 2 if exponentially scaled K's are to be calculated.
 BK    - Output vector of length NB.	If the
	 routine terminates normally (NCALC=NB), the vector BK
	 contains the functions K(ALPHA,X), ... , K(NB-1+ALPHA,X),
	 or the corresponding exponentially scaled functions.
	 If (0 < NCALC < NB), BK(I) contains correct function
	 values for I <= NCALC, and contains the ratios
	 K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
 NCALC - Output variable indicating possible errors.
	 Before using the vector BK, the user should check that
	 NCALC=NB, i.e., all orders have been calculated to
	 the desired accuracy.	See error returns below.


 *******************************************************************

 Error returns

  In case of an error, NCALC != NB, and not all K's are
  calculated to the desired accuracy.

  NCALC < -1:  An argument is out of range. For example,
	NB <= 0, IZE is not 1 or 2, or IZE=1 and ABS(X) >= XMAX_BESS_K.
	In this case, the B-vector is not calculated,
	and NCALC is set to MIN0(NB,0)-2	 so that NCALC != NB.
  NCALC = -1:  Either  K(ALPHA,X) >= XINF  or
	K(ALPHA+NB-1,X)/K(ALPHA+NB-2,X) >= XINF.	 In this case,
	the B-vector is not calculated.	Note that again
	NCALC != NB.

  0 < NCALC < NB: Not all requested function values could
	be calculated accurately.  BK(I) contains correct function
	values for I <= NCALC, and contains the ratios
	K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.


 Intrinsic functions required are:

     ABS, AINT, EXP, INT, LOG, MAX, MIN, SINH, SQRT


 Acknowledgement

	This program is based on a program written by J. B. Campbell
	(2) that computes values of the Bessel functions K of float
	argument and float order.  Modifications include the addition
	of non-scaled functions, parameterization of machine
	dependencies, and the use of more accurate approximations
	for SINH and SIN.

 References: "On Temme's Algorithm for the Modified Bessel
	      Functions of the Third Kind," Campbell, J. B.,
	      TOMS 6(4), Dec. 1980, pp. 581-586.

	     "A FORTRAN IV Subroutine for the Modified Bessel
	      Functions of the Third Kind of Real Order and Real
	      Argument," Campbell, J. B., Report NRC/ERB-925,
	      National Research Council, Canada.

  Latest modification: May 30, 1989

  Modified by: W. J. Cody and L. Stoltz
	       Applied Mathematics Division
	       Argonne National Laboratory
	       Argonne, IL  60439

 -------------------------------------------------------------------
*/
    /*---------------------------------------------------------------------
     * Mathematical constants
     *	A = LOG(2) - Euler's constant
     *	D = SQRT(2/PI)
     ---------------------------------------------------------------------*/
    const static double a = .11593151565841244881;

    /*---------------------------------------------------------------------
      P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA + Euler's constant
      Coefficients converted from hex to decimal and modified
      by W. J. Cody, 2/26/82 */
    const static double p[8] = { .805629875690432845,20.4045500205365151,
	    157.705605106676174,536.671116469207504,900.382759291288778,
	    730.923886650660393,229.299301509425145,.822467033424113231 };
    const static double q[7] = { 29.4601986247850434,277.577868510221208,
	    1206.70325591027438,2762.91444159791519,3443.74050506564618,
	    2210.63190113378647,572.267338359892221 };
    /* R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2.D0*ALPHA) */
    const static double r[5] = { -.48672575865218401848,13.079485869097804016,
	    -101.96490580880537526,347.65409106507813131,
	    3.495898124521934782e-4 };
    const static double s[4] = { -25.579105509976461286,212.57260432226544008,
	    -610.69018684944109624,422.69668805777760407 };
    /* T    - Approximation for SINH(Y)/Y */
    const static double t[6] = { 1.6125990452916363814e-10,
	    2.5051878502858255354e-8,2.7557319615147964774e-6,
	    1.9841269840928373686e-4,.0083333333333334751799,
	    .16666666666666666446 };
    /*---------------------------------------------------------------------*/
    const static double estm[6] = { 52.0583,5.7607,2.7782,14.4303,185.3004, 9.3715 };
    const static double estf[7] = { 41.8341,7.1075,6.4306,42.511,1.35633,84.5096,20.};

    /* Local variables */
    int iend, i, j, k, m, ii, mplus1;
    double x2by4, twox, c, blpha, ratio, wminf;
    double d1, d2, d3, f0, f1, f2, p0, q0, t1, t2, twonu;
    double dm, ex, bk1, bk2, nu;

    ii = 0; /* -Wall */

    ex = *x;
    nu = *alpha;
    *ncalc = min0(*nb,0) - 2;
    if (*nb > 0 && (0. <= nu && nu < 1.) && (1 <= *ize && *ize <= 2)) {
	if(ex <= 0 || (*ize == 1 && ex > xmax_BESS_K)) {
	    if(ex <= 0) {
		if(ex < 0) Rprintf("K_bessel - out-of-range");
		for(i=0; i < *nb; i++)
		    bk[i] = ML_POSINF;
	    } else /* would only have underflow */
		for(i=0; i < *nb; i++)
		    bk[i] = 0.;
	    *ncalc = *nb;
	    return;
	}
	k = 0;
	if (nu < sqxmin_BESS_K) {
	    nu = 0.;
	} else if (nu > .5) {
	    k = 1;
	    nu -= 1.;
	}
	twonu = nu + nu;
	iend = *nb + k - 1;
	c = nu * nu;
	d3 = -c;
	if (ex <= 1.) {
	    /* ------------------------------------------------------------
	       Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA
			      Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA
	       ------------------------------------------------------------ */
	    d1 = 0.; d2 = p[0];
	    t1 = 1.; t2 = q[0];
	    for (i = 2; i <= 7; i += 2) {
		d1 = c * d1 + p[i - 1];
		d2 = c * d2 + p[i];
		t1 = c * t1 + q[i - 1];
		t2 = c * t2 + q[i];
	    }
	    d1 = nu * d1;
	    t1 = nu * t1;
	    f1 = log(ex);
	    f0 = a + nu * (p[7] - nu * (d1 + d2) / (t1 + t2)) - f1;
	    q0 = exp(-nu * (a - nu * (p[7] + nu * (d1-d2) / (t1-t2)) - f1));
	    f1 = nu * f0;
	    p0 = exp(f1);
	    /* -----------------------------------------------------------
	       Calculation of F0 =
	       ----------------------------------------------------------- */
	    d1 = r[4];
	    t1 = 1.;
	    for (i = 0; i < 4; ++i) {
		d1 = c * d1 + r[i];
		t1 = c * t1 + s[i];
	    }
	    /* d2 := sinh(f1)/ nu = sinh(f1)/(f1/f0)
	     *	   = f0 * sinh(f1)/f1 */
	    if (fabs(f1) <= .5) {
		f1 *= f1;
		d2 = 0.;
		for (i = 0; i < 6; ++i) {
		    d2 = f1 * d2 + t[i];
		}
		d2 = f0 + f0 * f1 * d2;
	    } else {
		d2 = sinh(f1) / nu;
	    }
	    f0 = d2 - nu * d1 / (t1 * p0);
	    if (ex <= 1e-10) {
		/* ---------------------------------------------------------
		   X <= 1.0E-10
		   Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X)
		   --------------------------------------------------------- */
		bk[0] = f0 + ex * f0;
		if (*ize == 1) {
		    bk[0] -= ex * bk[0];
		}
		ratio = p0 / f0;
		c = ex * DBL_MAX;
		if (k != 0) {
		    /* ---------------------------------------------------
		       Calculation of K(ALPHA,X)
		       and  X*K(ALPHA+1,X)/K(ALPHA,X),	ALPHA >= 1/2
		       --------------------------------------------------- */
		    *ncalc = -1;
		    if (bk[0] >= c / ratio) {
			return;
		    }
		    bk[0] = ratio * bk[0] / ex;
		    twonu += 2.;
		    ratio = twonu;
		}
		*ncalc = 1;
		if (*nb == 1)
		    return;

		/* -----------------------------------------------------
		   Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X),
		   L = 1, 2, ... , NB-1
		   ----------------------------------------------------- */
		*ncalc = -1;
		for (i = 1; i < *nb; ++i) {
		    if (ratio >= c)
			return;

		    bk[i] = ratio / ex;
		    twonu += 2.;
		    ratio = twonu;
		}
		*ncalc = 1;
		goto L420;
	    } else {
		/* ------------------------------------------------------
		   10^-10 < X <= 1.0
		   ------------------------------------------------------ */
		c = 1.;
		x2by4 = ex * ex / 4.;
		p0 = .5 * p0;
		q0 = .5 * q0;
		d1 = -1.;
		d2 = 0.;
		bk1 = 0.;
		bk2 = 0.;
		f1 = f0;
		f2 = p0;
		do {
		    d1 += 2.;
		    d2 += 1.;
		    d3 = d1 + d3;
		    c = x2by4 * c / d2;
		    f0 = (d2 * f0 + p0 + q0) / d3;
		    p0 /= d2 - nu;
		    q0 /= d2 + nu;
		    t1 = c * f0;
		    t2 = c * (p0 - d2 * f0);
		    bk1 += t1;
		    bk2 += t2;
		} while (fabs(t1 / (f1 + bk1)) > DBL_EPSILON ||
			 fabs(t2 / (f2 + bk2)) > DBL_EPSILON);
		bk1 = f1 + bk1;
		bk2 = 2. * (f2 + bk2) / ex;
		if (*ize == 2) {
		    d1 = exp(ex);
		    bk1 *= d1;
		    bk2 *= d1;
		}
		wminf = estf[0] * ex + estf[1];
	    }
	} else if (DBL_EPSILON * ex > 1.) {
	    /* -------------------------------------------------
	       X > 1./EPS
	       ------------------------------------------------- */
	    *ncalc = *nb;
	    bk1 = 1. / (M_SQRT_2dPI * sqrt(ex));
	    for (i = 0; i < *nb; ++i)
		bk[i] = bk1;
	    return;

	} else {
	    /* -------------------------------------------------------
	       X > 1.0
	       ------------------------------------------------------- */
	    twox = ex + ex;
	    blpha = 0.;
	    ratio = 0.;
	    if (ex <= 4.) {
		/* ----------------------------------------------------------
		   Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 <= X <= 4.0
		   ----------------------------------------------------------*/
		d2 = trunc(estm[0] / ex + estm[1]);
		m = (int) d2;
		d1 = d2 + d2;
		d2 -= .5;
		d2 *= d2;
		for (i = 2; i <= m; ++i) {
		    d1 -= 2.;
		    d2 -= d1;
		    ratio = (d3 + d2) / (twox + d1 - ratio);
		}
		/* -----------------------------------------------------------
		   Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward
		   recurrence and K(ALPHA,X) from the wronskian
		   -----------------------------------------------------------*/
		d2 = trunc(estm[2] * ex + estm[3]);
		m = (int) d2;
		c = fabs(nu);
		d3 = c + c;
		d1 = d3 - 1.;
		f1 = DBL_MIN;
		f0 = (2. * (c + d2) / ex + .5 * ex / (c + d2 + 1.)) * DBL_MIN;
		for (i = 3; i <= m; ++i) {
		    d2 -= 1.;
		    f2 = (d3 + d2 + d2) * f0;
		    blpha = (1. + d1 / d2) * (f2 + blpha);
		    f2 = f2 / ex + f1;
		    f1 = f0;
		    f0 = f2;
		}
		f1 = (d3 + 2.) * f0 / ex + f1;
		d1 = 0.;
		t1 = 1.;
		for (i = 1; i <= 7; ++i) {
		    d1 = c * d1 + p[i - 1];
		    t1 = c * t1 + q[i - 1];
		}
		p0 = exp(c * (a + c * (p[7] - c * d1 / t1) - log(ex))) / ex;
		f2 = (c + .5 - ratio) * f1 / ex;
		bk1 = p0 + (d3 * f0 - f2 + f0 + blpha) / (f2 + f1 + f0) * p0;
		if (*ize == 1) {
		    bk1 *= exp(-ex);
		}
		wminf = estf[2] * ex + estf[3];
	    } else {
		/* ---------------------------------------------------------
		   Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by
		   backward recurrence, for  X > 4.0
		   ----------------------------------------------------------*/
		dm = trunc(estm[4] / ex + estm[5]);
		m = (int) dm;
		d2 = dm - .5;
		d2 *= d2;
		d1 = dm + dm;
		for (i = 2; i <= m; ++i) {
		    dm -= 1.;
		    d1 -= 2.;
		    d2 -= d1;
		    ratio = (d3 + d2) / (twox + d1 - ratio);
		    blpha = (ratio + ratio * blpha) / dm;
		}
		bk1 = 1. / ((M_SQRT_2dPI + M_SQRT_2dPI * blpha) * sqrt(ex));
		if (*ize == 1)
		    bk1 *= exp(-ex);
		wminf = estf[4] * (ex - fabs(ex - estf[6])) + estf[5];
	    }
	    /* ---------------------------------------------------------
	       Calculation of K(ALPHA+1,X)
	       from K(ALPHA,X) and  K(ALPHA+1,X)/K(ALPHA,X)
	       --------------------------------------------------------- */
	    bk2 = bk1 + bk1 * (nu + .5 - ratio) / ex;
	}
	/*--------------------------------------------------------------------
	  Calculation of 'NCALC', K(ALPHA+I,X),	I  =  0, 1, ... , NCALC-1,
	  &	  K(ALPHA+I,X)/K(ALPHA+I-1,X),	I = NCALC, NCALC+1, ... , NB-1
	  -------------------------------------------------------------------*/
	*ncalc = *nb;
	bk[0] = bk1;
	if (iend == 0)
	    return;

	j = 1 - k;
	if (j >= 0)
	    bk[j] = bk2;

	if (iend == 1)
	    return;

	m = min0((int) (wminf - nu),iend);
	for (i = 2; i <= m; ++i) {
	    t1 = bk1;
	    bk1 = bk2;
	    twonu += 2.;
	    if (ex < 1.) {
		if (bk1 >= DBL_MAX / twonu * ex)
		    break;
	    } else {
		if (bk1 / ex >= DBL_MAX / twonu)
		    break;
	    }
	    bk2 = twonu / ex * bk1 + t1;
	    ii = i;
	    ++j;
	    if (j >= 0) {
		bk[j] = bk2;
	    }
	}

	m = ii;
	if (m == iend) {
	    return;
	}
	ratio = bk2 / bk1;
	mplus1 = m + 1;
	*ncalc = -1;
	for (i = mplus1; i <= iend; ++i) {
	    twonu += 2.;
	    ratio = twonu / ex + 1./ratio;
	    ++j;
	    if (j >= 1) {
		bk[j] = ratio;
	    } else {
		if (bk2 >= DBL_MAX / ratio)
		    return;

		bk2 *= ratio;
	    }
	}
	*ncalc = max0(1, mplus1 - k);
	if (*ncalc == 1)
	    bk[0] = bk2;
	if (*nb == 1)
	    return;

L420:
	for (i = *ncalc; i < *nb; ++i) { /* i == *ncalc */
#ifndef IEEE_754
	    if (bk[i-1] >= DBL_MAX / bk[i])
		return;
#endif
	    bk[i] *= bk[i-1];
	    (*ncalc)++;
	}
    }
}

/* From http://www.netlib.org/specfun/gamma	Fortran translated by f2c,...
 *	------------------------------#####	Martin Maechler, ETH Zurich
 *
 *=========== was part of	ribesl (Bessel I(.))
 *===========			~~~~~~
 */

// used in bessel_i.c and bessel_j.c, hidden if possible.

double Rf_gamma_cody(double x)
{
/* ----------------------------------------------------------------------

   This routine calculates the GAMMA function for a float argument X.
   Computation is based on an algorithm outlined in reference [1].
   The program uses rational functions that approximate the GAMMA
   function to at least 20 significant decimal digits.	Coefficients
   for the approximation over the interval (1,2) are unpublished.
   Those for the approximation for X >= 12 are from reference [2].
   The accuracy achieved depends on the arithmetic system, the
   compiler, the intrinsic functions, and proper selection of the
   machine-dependent constants.

   *******************************************************************

   Error returns

   The program returns the value XINF for singularities or
   when overflow would occur.	 The computation is believed
   to be free of underflow and overflow.

   Intrinsic functions required are:

   INT, DBLE, EXP, LOG, REAL, SIN


   References:
   [1]  "An Overview of Software Development for Special Functions",
	W. J. Cody, Lecture Notes in Mathematics, 506,
	Numerical Analysis Dundee, 1975, G. A. Watson (ed.),
	Springer Verlag, Berlin, 1976.

   [2]  Computer Approximations, Hart, Et. Al., Wiley and sons, New York, 1968.

   Latest modification: October 12, 1989

   Authors: W. J. Cody and L. Stoltz
   Applied Mathematics Division
   Argonne National Laboratory
   Argonne, IL 60439
   ----------------------------------------------------------------------*/

/* ----------------------------------------------------------------------
   Mathematical constants
   ----------------------------------------------------------------------*/
    const static double sqrtpi = .9189385332046727417803297; /* == ??? */

/* *******************************************************************

   Explanation of machine-dependent constants

   beta	- radix for the floating-point representation
   maxexp - the smallest positive power of beta that overflows
   XBIG	- the largest argument for which GAMMA(X) is representable
	in the machine, i.e., the solution to the equation
	GAMMA(XBIG) = beta**maxexp
   XINF	- the largest machine representable floating-point number;
	approximately beta**maxexp
   EPS	- the smallest positive floating-point number such that  1.0+EPS > 1.0
   XMININ - the smallest positive floating-point number such that
	1/XMININ is machine representable

   Approximate values for some important machines are:

   beta	      maxexp	     XBIG

   CRAY-1		(S.P.)	      2		8191	    966.961
   Cyber 180/855
   under NOS	(S.P.)	      2		1070	    177.803
   IEEE (IBM/XT,
   SUN, etc.)	(S.P.)	      2		 128	    35.040
   IEEE (IBM/XT,
   SUN, etc.)	(D.P.)	      2		1024	    171.624
   IBM 3033	(D.P.)	     16		  63	    57.574
   VAX D-Format	(D.P.)	      2		 127	    34.844
   VAX G-Format	(D.P.)	      2		1023	    171.489

   XINF	 EPS	    XMININ

   CRAY-1		(S.P.)	 5.45E+2465   7.11E-15	  1.84E-2466
   Cyber 180/855
   under NOS	(S.P.)	 1.26E+322    3.55E-15	  3.14E-294
   IEEE (IBM/XT,
   SUN, etc.)	(S.P.)	 3.40E+38     1.19E-7	  1.18E-38
   IEEE (IBM/XT,
   SUN, etc.)	(D.P.)	 1.79D+308    2.22D-16	  2.23D-308
   IBM 3033	(D.P.)	 7.23D+75     2.22D-16	  1.39D-76
   VAX D-Format	(D.P.)	 1.70D+38     1.39D-17	  5.88D-39
   VAX G-Format	(D.P.)	 8.98D+307    1.11D-16	  1.12D-308

   *******************************************************************

   ----------------------------------------------------------------------
   Machine dependent parameters
   ----------------------------------------------------------------------
   */


    const static double xbig = 171.624;
    /* ML_POSINF ==   const double xinf = 1.79e308;*/
    /* DBL_EPSILON = const double eps = 2.22e-16;*/
    /* DBL_MIN ==   const double xminin = 2.23e-308;*/

    /*----------------------------------------------------------------------
      Numerator and denominator coefficients for rational minimax
      approximation over (1,2).
      ----------------------------------------------------------------------*/
    const static double p[8] = {
	-1.71618513886549492533811,
	24.7656508055759199108314,-379.804256470945635097577,
	629.331155312818442661052,866.966202790413211295064,
	-31451.2729688483675254357,-36144.4134186911729807069,
	66456.1438202405440627855 };
    const static double q[8] = {
	-30.8402300119738975254353,
	315.350626979604161529144,-1015.15636749021914166146,
	-3107.77167157231109440444,22538.1184209801510330112,
	4755.84627752788110767815,-134659.959864969306392456,
	-115132.259675553483497211 };
    /*----------------------------------------------------------------------
      Coefficients for minimax approximation over (12, INF).
      ----------------------------------------------------------------------*/
    const static double c[7] = {
	-.001910444077728,8.4171387781295e-4,
	-5.952379913043012e-4,7.93650793500350248e-4,
	-.002777777777777681622553,.08333333333333333331554247,
	.0057083835261 };

    /* Local variables */
    int i, n;
    int parity;/*logical*/
    double fact, xden, xnum, y, z, yi, res, sum, ysq;

    parity = (0);
    fact = 1.;
    n = 0;
    y = x;
    if (y <= 0.) {
	/* -------------------------------------------------------------
	   Argument is negative
	   ------------------------------------------------------------- */
	y = -x;
	yi = trunc(y);
	res = y - yi;
	if (res != 0.) {
	    if (yi != trunc(yi * .5) * 2.)
		parity = (1);
	    fact = -M_PI / sinpi(res);
	    y += 1.;
	} else {
	    return(ML_POSINF);
	}
    }
    /* -----------------------------------------------------------------
       Argument is positive
       -----------------------------------------------------------------*/
    if (y < DBL_EPSILON) {
	/* --------------------------------------------------------------
	   Argument < EPS
	   -------------------------------------------------------------- */
	if (y >= DBL_MIN) {
	    res = 1. / y;
	} else {
	    return(ML_POSINF);
	}
    } else if (y < 12.) {
	yi = y;
	if (y < 1.) {
	    /* ---------------------------------------------------------
	       EPS < argument < 1
	       --------------------------------------------------------- */
	    z = y;
	    y += 1.;
	} else {
	    /* -----------------------------------------------------------
	       1 <= argument < 12, reduce argument if necessary
	       ----------------------------------------------------------- */
	    n = (int) y - 1;
	    y -= (double) n;
	    z = y - 1.;
	}
	/* ---------------------------------------------------------
	   Evaluate approximation for 1. < argument < 2.
	   ---------------------------------------------------------*/
	xnum = 0.;
	xden = 1.;
	for (i = 0; i < 8; ++i) {
	    xnum = (xnum + p[i]) * z;
	    xden = xden * z + q[i];
	}
	res = xnum / xden + 1.;
	if (yi < y) {
	    /* --------------------------------------------------------
	       Adjust result for case  0. < argument < 1.
	       -------------------------------------------------------- */
	    res /= yi;
	} else if (yi > y) {
	    /* ----------------------------------------------------------
	       Adjust result for case  2. < argument < 12.
	       ---------------------------------------------------------- */
	    for (i = 0; i < n; ++i) {
		res *= y;
		y += 1.;
	    }
	}
    } else {
	/* -------------------------------------------------------------
	   Evaluate for argument >= 12.,
	   ------------------------------------------------------------- */
	if (y <= xbig) {
	    ysq = y * y;
	    sum = c[6];
	    for (i = 0; i < 6; ++i) {
		sum = sum / ysq + c[i];
	    }
	    sum = sum / y - y + sqrtpi;
	    sum += (y - .5) * log(y);
	    res = exp(sum);
	} else {
	    return(ML_POSINF);
	}
    }
    /* ----------------------------------------------------------------------
       Final adjustments and return
       ----------------------------------------------------------------------*/
    if (parity)
	res = -res;
    if (fact != 1.)
	res = fact / res;
    return res;
}

