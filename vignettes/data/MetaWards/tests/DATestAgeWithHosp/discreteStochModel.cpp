// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix discreteStochModel(NumericVector pars, int tstart, int tstop, arma::imat u, arma::mat C) {
    
    // set up index variables
    int i = 0, j = 0, k = 0;
    
    // set up auxiliary matrix for counts
    int nclasses = u.n_rows;
    int nages = u.n_cols;
    arma::imat u1(nclasses, nages);
    for(i = 0; i < nclasses; i++) {
        for(j = 0; j < nages; j++) {
            u1(i, j) = u(i, j);
        }
    }
    u1 = u;
  
    // extract parameters
    double nu = pars[0];
    double nuA = pars[1];
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
        probE[j] = pars[j + 2];
        probEP[j] = pars[j + nages + 2];
        probA[j] = pars[j + 2 * nages + 2];
        probP[j] = pars[j + 3 * nages + 2];
        probI1[j] = pars[j + 4 * nages + 2];
        probI1H[j] = pars[j + 5 * nages + 2];
        probI1D[j] = pars[j + 6 * nages + 2];
        probI2[j] = pars[j + 7 * nages + 2];
        probH[j] = pars[j + 8 * nages + 2];
        probHD[j] = pars[j + 9 * nages + 2];
    }
    
    // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
    //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
    
    // extract population sizes
    arma::vec N (nages);
    for(i = 0; i < nages; i++) {
        N[i] = 0.0;
        for(j = 0; j < nclasses; j++) {
            N[i] += (double) u1(j, i);
        }
    }
    
    // set up output
    IntegerMatrix out(tstop - tstart + 1, nclasses * nages + 1);
    int tcurr = 0;
    out(0, 0) = tstart;
    k = 1;
    for(i = 0; i < nclasses; i++) {
        for(j = 0; j < nages; j++) {
            out(0, k) = u1(i, j);
            k++;
        }
    }
    
    // set up vector of number of infectives
    arma::vec uinf(nages);
    
    // set up transmission rates
    arma::mat beta(nages, 1);
    
    // set up auxiliary variables
    double prob = 0.0;
    IntegerVector pathE(3);
    NumericVector mprobsE(3);
    IntegerVector pathI1(4);
    NumericVector mprobsI1(4);
    IntegerVector pathH(3);
    NumericVector mprobsH(3);
    tcurr++;
    tstart++;
    
    while(tstart <= tstop) {
    
        // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
        //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
        
        // loop over age classes                
        for(j = 0; j < nages; j++) {
            
            // transition probs out of hospital
            mprobsH[0] = probH[j] * (1.0 - probHD[j]);
            mprobsH[1] = probH[j] * probHD[j];
            mprobsH[2] = 1.0 - probH[j];
            
            // H out
            rmultinom(u1(9, j), mprobsH.begin(), 3, pathH.begin());
            u1(9, j) -= (pathH[0] + pathH[1]);
            u1(10, j) += pathH[0];
            u1(11, j) += pathH[1];
            
            // I2RI
            i = R::rbinom(u1(7, j), probI2[j]);
            u1(7, j) -= i;
            u1(8, j) += i;
            
            // transition probs out of I1
            mprobsI1[0] = probI1[j] * probI1H[j];
            mprobsI1[1] = probI1[j] * (1.0 - probI1D[j] - probI1H[j]);
            mprobsI1[2] = probI1[j] * probI1D[j];
            mprobsI1[3] = 1.0 - probI1[j];
            
            // I1 out
            rmultinom(u1(5, j), mprobsI1.begin(), 4, pathI1.begin());
            u1(5, j) -= (pathI1[0] + pathI1[1] + pathI1[2]);
            u1(9, j) += pathI1[0];
            u1(7, j) += pathI1[1];
            u1(6, j) += pathI1[2];
            
            // PI1
            i = R::rbinom(u1(4, j), probP[j]);
            u1(4, j) -= i;
            u1(5, j) += i;
            
            // ARA
            i = R::rbinom(u1(2, j), probA[j]);
            u1(2, j) -= i;
            u1(3, j) += i;
            
            // transition probs out of E
            mprobsE[0] = probE[j] * (1.0 - probEP[j]);
            mprobsE[1] = probE[j] * probEP[j];
            mprobsE[2] = 1.0 - probE[j];
            
            // E out
            rmultinom(u1(1, j), mprobsE.begin(), 3, pathE.begin());
            u1(1, j) -= (pathE[0] + pathE[1]);
            u1(2, j) += pathE[0];
            u1(4, j) += pathE[1];
        }
        
        // update infective counts for rate
        for(j = 0; j < nages; j++) {
            uinf[j] = (double) nuA * u1(2, j) + nu * (u1(4, j) + u1(5, j) + u1(7, j));
        }
        
        // SE
        beta = C * uinf / N;
        for(j = 0; j < nages; j++) {
            // Rprintf("t = %d beta[%d] = %.25f\n", tstart, j, beta[j]);
            prob = 1 - exp(-beta(j, 0));
            i = R::rbinom(u1(0, j), prob);
            u1(0, j) -= i;
            u1(1, j) += i;
        }
        
        // record output
        out(tcurr, 0) = tstart;
        k = 1;
        for(i = 0; i < nclasses; i++) {
            for(j = 0; j < nages; j++) {
                out(tcurr, k) = u1(i, j);
                k++;
            }
        }
        
        // update time 
        tcurr++;
        tstart++;
    }
    return out;
}

