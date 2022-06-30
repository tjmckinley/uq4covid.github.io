// [[Rcpp::depends(RcppArmadillo, sitmo)]]

// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include <sitmo.h>
#include <iostream>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// BEGIN declaring external constants

#include "tnormConsts.h"

// END declaring external constants

const double invroot2pi = {
  1/pow(2*M_PI, 0.5)
};

double devroye(double a, sitmo::prng &eng){
  bool accept = false;
  double x = 0;
  double mx = sitmo::prng::max();
  while(!accept){
    arma::vec u(2);
    u(0) = eng() / mx;
    u(1) = eng() / mx;
    x = -(1/a)*log(u(0)) + a;
    accept = u(1)*exp(a*(a/2 - x)) <= exp(-0.5*pow(x,2));
  }
  return x;
}

double devroye_ab(double a, double b, sitmo::prng &eng){
  if(a >= b){
    stop("'a' must be less than 'b'");
  }
  double lambda = ( (b>0) ? a : b);
  double eal = exp(-a*lambda);
  double ebl = exp(-b*lambda);
  double mx = sitmo::prng::max();
  bool accept = false;
  double x = 0;
  while(!accept){
    arma::vec u(2);
    u(0) = eng() / mx;
    u(1) = eng() / mx;
    x = -(1/lambda)*log(eal - u(0)*(eal-ebl));
    accept = u(1)*exp(a*(a/2 - x)) <= exp(-0.5*pow(x,2));
  }
  return x;
}

double Direct(double a, double b, sitmo::prng &eng){
  double x = 0;
  bool accept = false;
  double mx = sitmo::prng::max();
  while(!accept){
    int i = rand() % (2*tN+2) + 1;
    if(i < (2*tN + 2) && i > 1){
      double u = eng() / mx;
      double y = yVec[i-2]*u;
      if(y <= ybar[i-2]){
        accept = true;
        x = xVec[i-2] + delta[i-2]*u;
      }
      else{
        double v = eng() / mx;
        x = xVec[i-2] + d[i-2]*v;
        accept = y <= invroot2pi * exp(-0.5*pow(x,2));
      }
    }
    else{
      arma::vec u(2);
      u(0) = eng() / mx;
      u(1) = eng() / mx;
      double lambda = xVec[2*tN];
      x = -(1/lambda)*log(u(0)) + lambda;
      accept = u(1)*exp(lambda*(lambda/2 - x)) <= exp(-0.5*pow(x,2));
      if(accept){
        if(i<2){
          x = -x;
        }
      }
    }
    if(accept){
      accept = (x >= a) && (x <= b);
    } 
  }
  return x;
}

double rtnorm_lower_one(double a, sitmo::prng &eng){
  //TN(0,1,a, Inf)
  if(a <= amin){
    return Direct(a, std::numeric_limits<double>::infinity(), eng);
  }
  else if(a >= amax){
    return devroye(a, eng);
  }
  int ia = ja[static_cast<int>(floor(a/h)) - static_cast<int>(ceil(amin/h))];
  double x = 0;
  bool accept = false;
  double mx = sitmo::prng::max();
  while(!accept){
    int i = rand() %  (2*tN+1) + 1  - ia;
    i = i+ia;
    if(i > 2*tN ){
      x = devroye(xVec[2*tN], eng);
      accept = true;
    }
    else if(i < ia+2){
      double u = eng() / mx;
      x = xVec[i-1] + d[i-1]*u;
      if(x >= a){
        double v = eng() / mx;
        double y = yVec[i-1]*v;
        accept = y <= invroot2pi * exp(-0.5*pow(x,2));
      }
    }
    else{
      double u = eng() / mx;
      double y = u*yVec[i-1];
      if(y <= ybar[i-1]){
        accept = true;
        x = xVec[i-1] + u*delta[i-1];
      }
      else{
        double v = eng() / mx;
        x = xVec[i-1] + d[i-1]*v;
        accept = y <= invroot2pi * exp(-0.5*pow(x,2));
      }
    }
  }
  return x;
}

double rtnorm_upper_one(double b, sitmo::prng &eng){
  //TN(0,1,-Inf, b)
  return (-rtnorm_lower_one(-b, eng));
}

double Direct_b(double a, double b, sitmo::prng &eng){
  bool accept = false;
  double x;
  while(!accept){
    x = rtnorm_upper_one(b, eng);
    accept = x >= a;
  }
  return x;
}

double Direct_a(double a, double b, sitmo::prng &eng){
  bool accept = false;
  double x;
  while(!accept){
    x = rtnorm_lower_one(a, eng);
    accept = x <= b;
  }
  return x;
}

double rtnorm_one(double a, double b, sitmo::prng &eng){
  //n number of samples of TN(0,1, a, b)
  if(a >= b){
    stop("'a' must be less than 'b'");
  }
  if(std::isinf(a)){
    if(std::isinf(b)){
      return Direct(a,b, eng);
    }
    else{
      return rtnorm_upper_one(b, eng);
    }
  }
  else if(std::isinf(b)){
    return rtnorm_lower_one(a, eng);
  }
  else{
    if(a <= amin){
      if(b < amax){
        return Direct_b(a, b, eng);
      }
      else{
        return Direct(a,b, eng);
      }
    }
    else if(a >= amax){
      return Direct_a(a,b, eng);
    }
    else if(b >= amax){
      return Direct_a(a, b, eng);
    }
  }
  int ia = ja[static_cast<int>(floor(a/h)) - static_cast<int>(ceil(amin/h))];
  int ib = ja[static_cast<int>(ceil(b/h)) - static_cast<int>(ceil(amin/h))];
  if(ib-ia < 5){
    return devroye_ab(a, b, eng);
  }
  bool accept=false;
  double mx = sitmo::prng::max();
  double x;
  while(!accept){
    int i = rand() %  ib + 1  - ia;
    i = i+ia;
    if((i < ia+2) || (i > ib-2)){
      double u = eng() / mx;
      x = xVec[i-1] + d[i-1]*u;
      if(x >= a && x<=b){
        double v = eng() / mx;
        double y = yVec[i-1]*v;
        accept = y <= invroot2pi * exp(-0.5*pow(x,2));
      }
    }
    else{
      double u = eng() / mx;
      double y = u*yVec[i-1];
      if(y <= ybar[i-1]){
        accept = true;
        x = xVec[i-1] + u*delta[i-1];
      }
      else{
        double v = eng() / mx;
        x = xVec[i-1] + d[i-1]*v;
        accept = y <= invroot2pi * exp(-0.5*pow(x,2));
      }
    }
  }
  return x;
}

arma::vec rtnorm(int n, double a, double b, sitmo::prng &eng){
  arma::vec x(n);
  for(int i=0; i < n; ++i){
    x(i) = rtnorm_one(a, b, eng);
  }
  return x;
}

// log truncated discrete Gaussian p.d.f.
double ldtnorm_cpp(int x, double mu, double sigma, double LB, double UB) {
    if(std::isinf(LB)) {
        stop("Lower bound of truncated Gaussian must be finite currently\n");
    }
    double temp1 = R::pnorm(x - 0.5, mu, sigma, 1, 1);
    double temp2 = R::pnorm(x + 0.5, mu, sigma, 1, 1);
    double ldens = temp2 + log(1.0 - exp(temp1 - temp2));
    if(!arma::is_finite(ldens)) {
        temp1 = R::pnorm(x - 0.5, mu, sigma, 0, 1);
        temp2 = R::pnorm(x + 0.5, mu, sigma, 0, 1);
        ldens = temp1 + log(1.0 - exp(temp2 - temp1));
    }
    if(!arma::is_finite(ldens)) stop("Something wrong in TN\n");
    if(std::isinf(UB)) {
        // normalising constant
        temp1 = R::pnorm(LB - 0.5, mu, sigma, 0, 1);
        if(!arma::is_finite(temp1)) stop("Something wrong in TN LB\n");
        ldens -= temp1;
    } else {
        // normalising constant
        temp1 = R::pnorm(LB - 0.5, mu, sigma, 1, 1);
        temp2 = R::pnorm(UB + 0.5, mu, sigma, 1, 1);
        temp1 = temp2 + log(1.0 - exp(temp1 - temp2));
        if(!arma::is_finite(temp1)) {
            temp1 = R::pnorm(LB - 0.5, mu, sigma, 0, 1);
            temp2 = R::pnorm(UB + 0.5, mu, sigma, 0, 1);
            temp2 = temp1 + log(1.0 - exp(temp2 - temp1));
            temp1 = temp2;
        }
        if(!arma::is_finite(temp1)) stop("Something wrong in TN LB UB\n");
        ldens -= temp1;
    }
    return ldens;
}

// log truncated normalising constant for discrete truncated Gaussian
double ltnormconst_cpp(double mu, double sigma, double LB, double UB) {
    if(std::isinf(LB) || std::isinf(UB)) {
        stop("Bounds of truncated Gaussian normalising constant must be finite currently\n");
    }
    if(LB > UB) {
        stop("'LB' must be > 'UB' in Gaussian normalising constant correction currently\n");
    }
    double temp1 = R::pnorm(LB - 0.5, mu, sigma, 1, 1);
    double temp2 = R::pnorm(UB + 0.5, mu, sigma, 1, 1);
    double norm = temp2 + log(1.0 - exp(temp1 - temp2));
    if(!arma::is_finite(norm)) {
        temp1 = R::pnorm(LB - 0.5, mu, sigma, 0, 1);
        temp2 = R::pnorm(UB + 0.5, mu, sigma, 0, 1);
        norm = temp1 + log(1.0 - exp(temp2 - temp1));
    }
    if(!arma::is_finite(norm)) stop("Something wrong in Gaussian normalising constant mu = %f sigma = %f LB = %f UB = %f\n", mu, sigma, LB, UB);
    return norm;
}

// truncated discrete Gaussian sampling
int rdtnorm_cpp(double mu, double sigma, double LB, double UB, sitmo::prng &eng) {
    if(std::isinf(LB)) {
        stop("Lower bound of truncated Gaussian must be finite currently\n");
    }
    int x;
    double u;
    if(std::isinf(UB)) {
        u = rtnorm_one((LB - 0.5 - mu) / sigma, UB, eng);
    } else {
        u = rtnorm_one((LB - 0.5 - mu) / sigma, (UB + 0.5 - mu) / sigma, eng);
    }
    u = u * sigma + mu;
    x = (int) round(u);
    return x;
}

// truncated Gaussian sampling
double rtnorm_cpp(double mu, double sigma, double LB, double UB, sitmo::prng &eng) {
    if(std::isinf(LB)) {
        stop("Lower bound of truncated Gaussian must be finite currently\n");
    }
    double u;
    if(std::isinf(UB)) {
        u = rtnorm_one((LB - mu) / sigma, UB, eng);
    } else {
        u = rtnorm_one((LB - mu) / sigma, (UB - mu) / sigma, eng);
    }
    u = u * sigma + mu;
    return u;
}

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
            if(fabs(temp) < 1e-10 && fabs(u) < 1e-10) return(n);
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

// function to move events from stagefrom to stageto
void moveEvent (int ipart, int MDfrom, int stagefrom, int stageto, 
    int age, int lad, arma::icube &MD, 
    arma::vec &tempinc, arma::vec &tempinc1, arma::ivec &ncohorts,
    std::vector<arma::icube> &u1_new, sitmo::prng &eng, int posonly) {
    
    // set counters
    int k;
    double tempDenom;
    int MDstagefrom = stagefrom;
    int MDstageto = stageto;
    if(MDfrom == 0) {
        MDstagefrom = stageto;
        MDstageto = stagefrom;
    }
    int add = (MD(MDstagefrom, age, lad) < 0 ? 0:1);
    
    if(posonly == 1 && add != 1) stop("Incorrect inputs for moveEvent\n");
    
    // stratified sampling across cohorts
    tempDenom = sum(tempinc);
    if(tempDenom < abs(MD(MDstagefrom, age, lad))) stop("Can't match number of events, tempDenom = %f MD = %d j = %d l = %d from = %d to = %d\n", tempDenom, MD(MDstagefrom, age, lad), age, lad, stagefrom, stageto);
    while(abs(MD(MDstagefrom, age, lad)) > 0) {
        // sample cohort
        tempinc1 = tempinc / tempDenom;
        if(fabs(sum(tempinc1) - 1.0) > 1e-10) stop("'tempinc1' doesn't sum to one\n");
        k = rmultinom_cpp(tempinc1, eng);
        
        // move event from stagefrom to stageto in kth cohort
        u1_new[ipart](stagefrom, age, k + ncohorts(lad))--;
        u1_new[ipart](stageto, age, k + ncohorts(lad))++;
        
        // update denominator and probability
        tempDenom--;
        tempinc(k)--;
        
        // update MD counters
        if(add == 0) {
            MD(MDstagefrom, age, lad)++;
            MD(MDstageto, age, lad)--;
        } else {
            MD(MDstagefrom, age, lad)--;
            if(posonly != 1) MD(MDstageto, age, lad)++;
        }
    }
    return;
}

// function to distribute incidence across cohorts
void redistribution (int ipart, int nages, int nlads, arma::icube &inc, arma::ivec &ncohorts,
    std::vector<arma::icube> &u1, std::vector<arma::icube> &u1_new, sitmo::prng &eng, int posonly) {
    
    // set counters
    int j, l, k, MD;
    double tempDenom;
    
    for(j = 0; j < nages; j++) {
        for(l = 0; l < nlads; l++) {
        
            // set up auxiliary vectors
            arma::vec tempinc (ncohorts(l + 1) - ncohorts(l)); tempinc.zeros();
            arma::vec tempinc1 (ncohorts(l + 1) - ncohorts(l)); tempinc1.zeros();
            
            // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
            //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
        
            // negative moves out of DH
            MD = inc(11, j, l);
            if(MD < 0) {
                // calculate event incidence in cohort
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1_new[ipart](11, j, k) - u1[ipart](11, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 1, 11, 9, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            // negative moves out of RH
            MD = inc(10, j, l);
            if(MD < 0) {
                // calculate event incidence in cohort
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1_new[ipart](10, j, k) - u1[ipart](10, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 1, 10, 9, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            // positive moves into DH
            MD = inc(11, j, l);
            if(MD > 0) {
                // calculate number remaining who can move
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1[ipart](9, j, k);
                    tempinc(k - ncohorts(l)) -= (u1_new[ipart](11, j, k) - u1[ipart](11, j, k));
                    tempinc(k - ncohorts(l)) -= (u1_new[ipart](10, j, k) - u1[ipart](10, j, k));
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 0, 9, 11, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            // positive moves into RH
            MD = inc(10, j, l);
            if(MD > 0) {
                // calculate number remaining who can move
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1[ipart](9, j, k);
                    tempinc(k - ncohorts(l)) -= (u1_new[ipart](11, j, k) - u1[ipart](11, j, k));
                    tempinc(k - ncohorts(l)) -= (u1_new[ipart](10, j, k) - u1[ipart](10, j, k));
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 0, 9, 10, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            
            // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
            //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
            
            // negative moves out of RI
            MD = inc(8, j, l);
            if(MD < 0) {
                // calculate event incidence in cohort
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1_new[ipart](8, j, k) - u1[ipart](8, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 1, 8, 7, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            // positive moves into RI
            MD = inc(8, j, l);
            if(MD > 0) {
                // calculate number remaining who can move
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1[ipart](7, j, k);
                    tempinc(k - ncohorts(l)) -= (u1_new[ipart](8, j, k) - u1[ipart](8, j, k));
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 0, 7, 8, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            
            // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
            //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
            
            // negative moves out of DI
            MD = inc(6, j, l);
            if(MD < 0) {
                // calculate event incidence in cohort
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1_new[ipart](6, j, k) - u1[ipart](6, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 1, 6, 5, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            // negative moves out of H
            MD = inc(9, j, l);
            if(MD < 0) {
                // calculate event incidence in cohort
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1_new[ipart](9, j, k) - u1[ipart](9, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](11, j, k) - u1[ipart](11, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](10, j, k) - u1[ipart](10, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 1, 9, 5, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            // negative moves out of I2
            MD = inc(7, j, l);
            if(MD < 0) {
                // calculate event incidence in cohort
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1_new[ipart](7, j, k) - u1[ipart](7, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](8, j, k) - u1[ipart](8, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 1, 7, 5, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            // positive moves into DI
            MD = inc(6, j, l);
            if(MD > 0) {
                // calculate number remaining who can move
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1[ipart](5, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](6, j, k) - u1[ipart](6, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](7, j, k) - u1[ipart](7, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](8, j, k) - u1[ipart](8, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](9, j, k) - u1[ipart](9, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](10, j, k) - u1[ipart](10, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](11, j, k) - u1[ipart](11, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 0, 5, 6, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            // positive moves into H
            MD = inc(9, j, l);
            if(MD > 0) {
                // calculate number remaining who can move
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1[ipart](5, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](6, j, k) - u1[ipart](6, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](7, j, k) - u1[ipart](7, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](8, j, k) - u1[ipart](8, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](9, j, k) - u1[ipart](9, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](10, j, k) - u1[ipart](10, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](11, j, k) - u1[ipart](11, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 0, 5, 9, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            // positive moves into I2
            MD = inc(7, j, l);
            if(MD > 0) {
                // calculate number remaining who can move
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1[ipart](5, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](6, j, k) - u1[ipart](6, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](7, j, k) - u1[ipart](7, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](8, j, k) - u1[ipart](8, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](9, j, k) - u1[ipart](9, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](10, j, k) - u1[ipart](10, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](11, j, k) - u1[ipart](11, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 0, 5, 7, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            
            // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
            //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
            
            // negative moves out of I1
            MD = inc(5, j, l);
            if(MD < 0) {
                // calculate event incidence in cohort
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1_new[ipart](5, j, k) - u1[ipart](5, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](6, j, k) - u1[ipart](6, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](7, j, k) - u1[ipart](7, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](8, j, k) - u1[ipart](8, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](9, j, k) - u1[ipart](9, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](11, j, k) - u1[ipart](11, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](10, j, k) - u1[ipart](10, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 1, 5, 4, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            // positive moves into I1
            MD = inc(5, j, l);
            if(MD > 0) {
                // calculate number remaining who can move
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1[ipart](4, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](5, j, k) - u1[ipart](5, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](6, j, k) - u1[ipart](6, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](7, j, k) - u1[ipart](7, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](8, j, k) - u1[ipart](8, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](9, j, k) - u1[ipart](9, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](11, j, k) - u1[ipart](11, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](10, j, k) - u1[ipart](10, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 0, 4, 5, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            
            // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
            //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
            
            // negative moves out of RA
            MD = inc(3, j, l);
            if(MD < 0) {
                // calculate event incidence in cohort
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1_new[ipart](3, j, k) - u1[ipart](3, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 1, 3, 2, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            // positive moves into RA
            MD = inc(3, j, l);
            if(MD > 0) {
                // calculate number remaining who can move
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1[ipart](2, j, k);
                    tempinc(k - ncohorts(l)) -= (u1_new[ipart](3, j, k) - u1[ipart](3, j, k));
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 0, 2, 3, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            
            // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
            //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
            
            // negative moves out of P
            MD = inc(4, j, l);
            if(MD < 0) {
                // calculate event incidence in cohort
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1_new[ipart](4, j, k) - u1[ipart](4, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](5, j, k) - u1[ipart](5, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](6, j, k) - u1[ipart](6, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](7, j, k) - u1[ipart](7, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](8, j, k) - u1[ipart](8, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](9, j, k) - u1[ipart](9, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](11, j, k) - u1[ipart](11, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](10, j, k) - u1[ipart](10, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 1, 4, 1, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            // negative moves out of A
            MD = inc(2, j, l);
            if(MD < 0) {
                // calculate event incidence in cohort
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1_new[ipart](2, j, k) - u1[ipart](2, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](3, j, k) - u1[ipart](3, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 1, 2, 1, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            // positive moves into P
            MD = inc(4, j, l);
            if(MD > 0) {
                // calculate number remaining who can move
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1[ipart](1, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](2, j, k) - u1[ipart](2, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](3, j, k) - u1[ipart](3, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](4, j, k) - u1[ipart](4, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](5, j, k) - u1[ipart](5, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](6, j, k) - u1[ipart](6, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](7, j, k) - u1[ipart](7, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](8, j, k) - u1[ipart](8, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](9, j, k) - u1[ipart](9, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](11, j, k) - u1[ipart](11, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](10, j, k) - u1[ipart](10, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 0, 1, 4, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            // positive moves into A
            MD = inc(2, j, l);
            if(MD > 0) {
                // calculate number remaining who can move
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1[ipart](1, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](2, j, k) - u1[ipart](2, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](3, j, k) - u1[ipart](3, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](4, j, k) - u1[ipart](4, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](5, j, k) - u1[ipart](5, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](6, j, k) - u1[ipart](6, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](7, j, k) - u1[ipart](7, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](8, j, k) - u1[ipart](8, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](9, j, k) - u1[ipart](9, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](11, j, k) - u1[ipart](11, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](10, j, k) - u1[ipart](10, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 0, 1, 2, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            
            // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
            //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
            
            // negative moves out of E
            MD = inc(1, j, l);
            if(MD < 0) {
                // calculate event incidence in cohort
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1_new[ipart](1, j, k) - u1[ipart](1, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](2, j, k) - u1[ipart](2, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](3, j, k) - u1[ipart](3, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](4, j, k) - u1[ipart](4, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](5, j, k) - u1[ipart](5, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](6, j, k) - u1[ipart](6, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](7, j, k) - u1[ipart](7, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](8, j, k) - u1[ipart](8, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](9, j, k) - u1[ipart](9, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](11, j, k) - u1[ipart](11, j, k);
                    tempinc(k - ncohorts(l)) += u1_new[ipart](10, j, k) - u1[ipart](10, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 1, 1, 0, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
            // positive moves into E
            MD = inc(1, j, l);
            if(MD > 0) {
                // calculate number remaining who can move
                for(k = ncohorts(l); k < ncohorts(l + 1); k++) {
                    tempinc(k - ncohorts(l)) = u1[ipart](0, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](1, j, k) - u1[ipart](1, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](2, j, k) - u1[ipart](2, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](3, j, k) - u1[ipart](3, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](4, j, k) - u1[ipart](4, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](5, j, k) - u1[ipart](5, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](6, j, k) - u1[ipart](6, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](7, j, k) - u1[ipart](7, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](8, j, k) - u1[ipart](8, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](9, j, k) - u1[ipart](9, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](11, j, k) - u1[ipart](11, j, k);
                    tempinc(k - ncohorts(l)) -= u1_new[ipart](10, j, k) - u1[ipart](10, j, k);
                    if(tempinc(k - ncohorts(l)) < 0) stop("Error in redistribution\n");
                }
                moveEvent(ipart, 0, 0, 1, j, l, inc, tempinc, tempinc1, ncohorts, u1_new, eng, posonly);
            }
        }
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
List TPF_cpp (arma::vec pars, arma::mat C, arma::imat data, arma::uword nclasses, arma::uword nages, arma::uword nlads, arma::imat u1_moves, arma::ivec ncohorts, 
         arma::icube u1_comb, arma::uword ndays, arma::uword npart, arma::mat munorm, arma::mat varnorm, double pmix, int twist, double a1, double a2, double b, double a_dis, 
         double b_dis, int saveAll, int writeExt, int returnPsi, int PF, int ncores) {
    
    // set counters
    arma::uword i, j, l, k, t = 0;
    
    // split u1 up into different LADs
    std::vector<arma::icube> u1(npart);
    std::vector<arma::icube> u1_new(npart);
    std::vector<arma::vec> twistnorm(npart);
    std::vector< std::vector<arma::vec> > tempdensx(npart);
    // set up new objects (need a better way to do this)
    std::vector<arma::vec> twistnorm1(npart);
    std::vector< std::vector<arma::vec> > tempdensx1(npart);
    for(i = 0; i < npart; i++) {
        twistnorm[i] = arma::vec(munorm.n_cols); twistnorm[i].zeros();
        tempdensx[i] = std::vector<arma::vec> (munorm.n_cols);
        tempdensx1[i] = std::vector<arma::vec> (munorm.n_cols);
    }
    arma::icube u1_night_full(nclasses + 2, nages, nlads); u1_night_full.zeros();
    arma::icube u1_night_reduced(4, nages, nlads); u1_night_reduced.zeros();
    std::vector<arma::icube> u1_night_obs(npart);
    std::vector<arma::icube> u1_night_obs1(npart);
    for(i = 0; i < npart; i++) {
        u1_night_obs[i] = arma::icube (2, nages, nlads); u1_night_obs[i].zeros();
    }
    arma::icube psi(4 * npart * ndays, nages, nlads); psi.zeros();
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
        
    // set up weight vector
    arma::vec weights (npart);
    arma::ivec inds(npart);
    double wnorm = 0.0;
    
    // set up auxiliary objects    
    arma::ivec obsInc (data.n_cols); obsInc.zeros();
    
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
    
    // check which output required
    List out (npart * (ndays + 1));
    std::ofstream file;
    char file_name[128];
    if(saveAll != 0) {
        // set up print string for debugging
        char str1[80];
        std::strcpy(str1, "save");
        double muy, sigma2y, u;
        if(saveAll == 1) {
            for(i = 0; i < npart; i++) {
                // extract just counts for DI and DH
                u1_night_reduced.zeros();
                for(l = 0; l < u1_moves.n_rows; l++) {
                    for(j = 0; j < nages; j++) {
                        u1_night_reduced(0, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](6, j, l);
                        u1_night_reduced(1, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](11, j, l);
                        // incidence
                        u1_night_reduced(2, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](6, j, l);
                        u1_night_reduced(3, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](11, j, l);
                    }
                }
                // apply observation error
                muy = a1 - a2;
                for(l = 0; l < nlads; l++) {
                    for(j = 0; j < nages; j++) {
                        sigma2y = a1 + a2 + 2.0 * b * u1_night_reduced(2, j, l);
                        u1_night_reduced(2, j, l) += rdtnorm_cpp(
                            muy, 
                            sqrt(sigma2y),
                            -u1_night_reduced(2, j, l),
                            std::numeric_limits<double>::infinity(),
                            engSerial
                        );
                        
                        sigma2y = a1 + a2 + 2.0 * b * u1_night_reduced(3, j, l);
                        u1_night_reduced(3, j, l) += rdtnorm_cpp(
                            muy, 
                            sqrt(sigma2y),
                            -u1_night_reduced(3, j, l),
                            std::numeric_limits<double>::infinity(),
                            engSerial
                        );
                        u1_night_obs[i](0, j, l) = u1_night_reduced(2, j, l);
                        u1_night_obs[i](1, j, l) = u1_night_reduced(3, j, l);
                    }
                }
                if(writeExt == 0) {
                    out[i] = u1_night_reduced;
                } else {
                    std::sprintf(file_name, "saveOut/p_%u.csv", i);
                    file.open(file_name);
                    file << "time, class, ";
                    for(j = 0; j < nages; j++) file << "age" << j + 1 << ", ";
                    file << "lad\n";
                    for(l = 0; l < nlads; l++) {
                        for(arma::uword r = 0; r < 4; r++) {
                            file << t << ", " << r << ", ";
                            for(j = 0; j < nages; j++) {
                                file << u1_night_reduced(r, j, l) << ", ";
                            }
                            file << l + 1 << "\n";
                        }
                    }
                    file.close();
                }
            }
        } else {
            for(i = 0; i < npart; i++) {
                u1_night_full.zeros();
                for(l = 0; l < u1_moves.n_rows; l++) {
                    for(j = 0; j < nages; j++) {
                        for(k = 0; k < nclasses; k++) {
                            u1_night_full(k, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](k, j, l);
                        }
                        // incidence
                        u1_night_full(nclasses, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](6, j, l);
                        u1_night_full(nclasses + 1, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](11, j, l);
                    }
                }
                // apply observation error
                muy = a1 - a2;
                for(l = 0; l < nlads; l++) {
                    for(j = 0; j < nages; j++) {
                        sigma2y = a1 + a2 + 2.0 * b * u1_night_full(nclasses, j, l);
                        u1_night_full(nclasses, j, l) += rdtnorm_cpp(
                            muy, 
                            sqrt(sigma2y),
                            -u1_night_full(nclasses, j, l),
                            std::numeric_limits<double>::infinity(),
                            engSerial
                        );
                        
                        sigma2y = a1 + a2 + 2.0 * b * u1_night_full(nclasses + 1, j, l);
                        u1_night_full(nclasses + 1, j, l) += rdtnorm_cpp(
                            muy, 
                            sqrt(sigma2y),
                            -u1_night_full(nclasses + 1, j, l),
                            std::numeric_limits<double>::infinity(),
                            engSerial
                        );
                        u1_night_obs[i](0, j, l) = u1_night_full(nclasses, j, l);
                        u1_night_obs[i](1, j, l) = u1_night_full(nclasses + 1, j, l);
                    }
                }
                if(writeExt == 0) {
                    out[i] = u1_night_full;
                } else {
                    std::sprintf(file_name, "saveOut/p_%u.csv", i);
                    file.open(file_name);
                    file << "time, class, ";
                    for(j = 0; j < nages; j++) file << "age" << j + 1 << ", ";
                    file << "lad\n";
                    for(l = 0; l < nlads; l++) {
                        for(arma::uword r = 0; r < (nclasses + 2); r++) {
                            file << t << ", " << r << ", ";
                            for(j = 0; j < nages; j++) {
                                file << u1_night_full(r, j, l) << ", ";
                            }
                            file << l + 1 << "\n";
                        }
                    }
                    file.close();
                }
            }
        }
    }
    // initialise timer
    Timer timer;
    int timer_cnt = 0;
    double prev_time = 0.0;
    
    // parameters for conditional sampling
    arma::vec condpars (pars.n_elem); condpars.zeros();
    condpars = pars;
    for(j = 0; j < nages; j++) {
        condpars(j + 4 * nages + 2) = pars(j + 4 * nages + 2) * (1.0 - pars(j + 6 * nages + 2));
        condpars(j + 4 * nages + 2) /= (1.0 - pars(j + 4 * nages + 2) * pars(j + 6 * nages + 2));
        condpars(j + 5 * nages + 2) = pars(j + 5 * nages + 2);
        condpars(j + 5 * nages + 2) /= (1.0 - pars(j + 6 * nages + 2));
        condpars(j + 6 * nages + 2) = 0.0;
        
        condpars(j + 8 * nages + 2) = pars(j + 8 * nages + 2) * (1.0 - pars(j + 9 * nages + 2));
        condpars(j + 8 * nages + 2) /= (1.0 - pars(j + 8 * nages + 2) * pars(j + 9 * nages + 2));
        condpars(j + 9 * nages + 2) = 0.0;
    }
            
    // if twisting, then calculate first normalising constants
    if(twist == 1) {
    
        // set counter for rates
        t = 0;
    
#ifdef _OPENMP
#pragma omp parallel for default(none) private(j, l) shared(seeds, npart, u1_moves, nages, nclasses, nlads, u1, pars, a_dis, b_dis, munorm, varnorm, twistnorm, tempdensx, t, pmix)
#endif
        for(i = 0; i < npart; i++) {
    
            // set up print string for debugging
            char str1[80];
	        std::strcpy(str1, "twist0");
 
            // set up thread-safe RNG
            uint32_t coreseed = static_cast<uint32_t>(seeds(0));
#ifdef _OPENMP
            coreseed = static_cast<uint32_t>(seeds((arma::uword) omp_get_thread_num()));
#endif
            sitmo::prng eng(coreseed);
        
            // set up auxiliary objects
            arma::icube u1_night(nclasses, nages, nlads); u1_night.zeros();
            
            // aggregate counts to LAD-level
            for(j = 0; j < nages; j++) {                    
                for(l = 0; l < u1_moves.n_rows; l++) {
                    u1_night(9, j, u1_moves(l, 0) - 1) += u1[i](9, j, l);
                    u1_night(5, j, u1_moves(l, 0) - 1) += u1[i](5, j, l);
                }
            }
            
            // start counter
            int w = 0;
            
            // set auxiliary variables
            double muy, sigma2y, pI1pI1D, pHpHD, weight1, weight2, weightnorm;
            
            // DI      
            for(j = 0; j < nages; j++) {
                
                // extract transition probabilities
                pI1pI1D = pars(j + 4 * nages + 2) * pars(j + 6 * nages + 2);
                for(l = 0; l < nlads; l++) {
                
                    // set up auxiliary matrix for sampling
                    tempdensx[i][w] = arma::vec (u1_night(5, j, l) + 1);

                    // loop over x values
                    for(int s = 0; s <= u1_night(5, j, l); s++) {
                        
                        // Gaussian correction
                        sigma2y = 2.0 * a_dis + 2.0 * b_dis * u1_night(5, j, l) * pI1pI1D;
                        tempdensx[i][w](s) = -0.5 * pow(s - munorm(t, w), 2.0) / (sigma2y + varnorm(t, w));
                        tempdensx[i][w](s) -= 0.5 * log(2.0 * M_PI * (sigma2y + varnorm(t, w)));
                        
                        // target normalising constant correction  
                        tempdensx[i][w](s) -= ltnormconst_cpp((double) s, sqrt(sigma2y), 0.0, u1_night(5, j, l));
                        
                        // integrate over twisting function densities
                        muy = varnorm(t, w) * s + sigma2y * munorm(t, w);
                        muy /= (varnorm(t, w) + sigma2y);
                        sigma2y = sigma2y * varnorm(t, w) / (varnorm(t, w) + sigma2y);
                        tempdensx[i][w](s) += ltnormconst_cpp(muy, sqrt(sigma2y), 0.0, u1_night(5, j, l));
                        
                        // correct for mixture
                        weight1 = log(pmix) + tempdensx[i][w](s);
                        weight2 = log(1.0 - pmix);
                        weightnorm = (weight1 > weight2 ? weight1:weight2);
                        weightnorm += log(exp(weight1 - weightnorm) + exp(weight2 - weightnorm));
                        tempdensx[i][w](s) = weightnorm;
                        
                        // simulator density
                        tempdensx[i][w](s) += R::dbinom(s, u1_night(5, j, l), pI1pI1D, 1);
                    }
                    
                    // calculate normalising constant
                    twistnorm[i](w) = log_sum_exp(tempdensx[i][w], 0);
                    
                    // increment counter
                    w++;
                }
            }
            
            // DH
            for(j = 0; j < nages; j++) {
                
                // extract transition probabilities
                pHpHD = pars(j + 8 * nages + 2) * pars(j + 9 * nages + 2);
                for(l = 0; l < nlads; l++) {
                
                    // set up auxiliary matrix for sampling
                    tempdensx[i][w] = arma::vec (u1_night(9, j, l) + 1);

                    // loop over x values
                    for(int s = 0; s <= u1_night(9, j, l); s++) {
                        
                        // Gaussian correction
                        sigma2y = 2.0 * a_dis + 2.0 * b_dis * u1_night(9, j, l) * pHpHD;
                        tempdensx[i][w](s) = -0.5 * pow(s - munorm(t, w), 2.0) / (sigma2y + varnorm(t, w));
                        tempdensx[i][w](s) -= 0.5 * log(2.0 * M_PI * (sigma2y + varnorm(t, w)));
                        
                        // target normalising constant correction  
                        tempdensx[i][w](s) -= ltnormconst_cpp((double) s, sqrt(sigma2y), 0.0, u1_night(9, j, l));
                        
                        // integrate over twisting function densities
                        muy = varnorm(t, w) * s + sigma2y * munorm(t, w);
                        muy /= (varnorm(t, w) + sigma2y);
                        sigma2y = sigma2y * varnorm(t, w) / (varnorm(t, w) + sigma2y);
                        tempdensx[i][w](s) += ltnormconst_cpp(muy, sqrt(sigma2y), 0.0, u1_night(9, j, l));
                        
                        // correct for mixture
                        weight1 = log(pmix) + tempdensx[i][w](s);
                        weight2 = log(1.0 - pmix);
                        weightnorm = (weight1 > weight2 ? weight1:weight2);
                        weightnorm += log(exp(weight1 - weightnorm) + exp(weight2 - weightnorm));
                        tempdensx[i][w](s) = weightnorm;
                        
                        // simulator density
                        tempdensx[i][w](s) += R::dbinom(s, u1_night(9, j, l), pHpHD, 1);
                    }
                    
                    // calculate normalising constant
                    twistnorm[i](w) = log_sum_exp(tempdensx[i][w], 0);
                    
                    // increment counter
                    w++;
                }
            }
        }
    }      
    
    // loop over time
    double ll = 0.0;
    for(t = 0; t < ndays; t++) {
            
        // check for interrupt
        R_CheckUserInterrupt();
        
        // extract data
        if(PF == 1) obsInc = data.row(t).t();
        
        // loop over particles
#ifdef _OPENMP
#pragma omp parallel for default(none) private(j, l, k) shared(seeds, npart, u1_moves, nages, nclasses, nlads, data, C, N_night, N_day, u1, u1_new, t, pars, weights, a_dis, b_dis, a1, a2, b, obsInc, PF, ncohorts, munorm, varnorm, twistnorm, tempdensx, condpars, twist, ndays, pmix)
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
            arma::imat DHinc (nages, nlads); DHinc.zeros();
            arma::imat DIinc (nages, nlads); DIinc.zeros();
            arma::imat DHinc1 (nages, nlads); DHinc1.zeros();
            arma::imat DIinc1 (nages, nlads); DIinc1.zeros();
            arma::imat RHinc (nages, nlads); RHinc.zeros();
            arma::imat Hinc (nages, nlads); Hinc.zeros();
            arma::imat RIinc (nages, nlads); RIinc.zeros();
            arma::imat I2inc (nages, nlads); I2inc.zeros();
            arma::imat I1inc (nages, nlads); I1inc.zeros();
            arma::imat Pinc (nages, nlads); Pinc.zeros();
            arma::imat RAinc (nages, nlads); RAinc.zeros();
            arma::imat Ainc (nages, nlads); Ainc.zeros();
            arma::imat Einc (nages, nlads); Einc.zeros();
            arma::icube tempMD (nclasses, nages, nlads); tempMD.zeros();
            
            arma::mat pinf(nages, nlads); pinf.zeros();
            arma::imat origE(nages, u1_moves.n_rows); origE.zeros();
            arma::icube u1_day(nclasses, nages, nlads); u1_day.zeros();
            arma::icube u1_night(nclasses, nages, nlads); u1_night.zeros();
            arma::icube u1_night1(nclasses, nages, nlads); u1_night1.zeros();
            
            // aggregate counts to LAD-level
            for(j = 0; j < nages; j++) {                    
                for(l = 0; l < u1_moves.n_rows; l++) {
                    for(int s = 0; s < nclasses; s++) {
                        u1_night(s, j, u1_moves(l, 0) - 1) += u1[i](s, j, l);
                    }
                }
            }
            
            // set weights
            weights(i) = 0.0;
            
            // set auxiliary variables
            double muy, sigma2y, pI1pI1D, pHpHD, weight1, weight2, weightnorm, pcondmix, u;
            
            // if twisting, then sample new observations
            if(twist == 1) {
            
                // start counter
                int w = 0;
                
                // DI           
                for(j = 0; j < nages; j++) {
                        
                    // set transition probability
                    pI1pI1D = pars(j + 4 * nages + 2) * pars(j + 6 * nages + 2);
                    for(l = 0; l < nlads; l++) {
                    
                        // sample simulator from twisted density
                        arma::vec tdensx (tempdensx[i][w].n_elem);
                        tdensx = exp(tempdensx[i][w] - twistnorm[i](w));
                        if(fabs(sum(tdensx) - 1.0) > 1e-10) {
                            Rprintf("sum(tdens) = %e\n", sum(tdensx));
                            stop("Must have sum(tdensx) == 1 in rmultinom\n");
                        }
                        tdensx = tdensx / sum(tdensx);
                        DIinc1(j, l) = rmultinom_cpp(tdensx, eng);
                        
                        if(tempdensx[i][w].n_elem != (u1_night(5, j, l) + 1)) stop("Grrr\n");
                        
                        // extract mixture probability
//                        weight1 = tempdensx[i][w](DIinc1(j, l));
//                        weight1 -= R::dbinom(DIinc1(j, l), u1_night(5, j, l), pI1pI1D, 1);
//                        weight2 = log(1.0 - pmix) - weight1;
//                        weight1 = log(1.0 - exp(weight2));
//                        pcondmix = exp(weight1);
//                        
                        // Gaussian correction
                        sigma2y = 2.0 * a_dis + 2.0 * b_dis * u1_night(5, j, l) * pI1pI1D;
                        weight1 = -0.5 * pow(DIinc1(j, l) - munorm(t, w), 2.0) / (sigma2y + varnorm(t, w));
                        weight1 -= 0.5 * log(2.0 * M_PI * (sigma2y + varnorm(t, w)));
                        
                        // target normalising constant correction  
                        weight1 -= ltnormconst_cpp((double) DIinc1(j, l), sqrt(sigma2y), 0.0, u1_night(5, j, l));
                        
                        // integrate over twisting function densities
                        muy = varnorm(t, w) * DIinc1(j, l) + sigma2y * munorm(t, w);
                        muy /= (varnorm(t, w) + sigma2y);
                        sigma2y = sigma2y * varnorm(t, w) / (varnorm(t, w) + sigma2y);
                        weight1 += ltnormconst_cpp(muy, sqrt(sigma2y), 0.0, u1_night(5, j, l));
                        pcondmix = weight1 + log(pmix);
                        
                        // adjust for mixture
                        weight1 = pcondmix;
                        weight2 = log(1.0 - pmix);
                        weightnorm = (weight1 > weight2 ? weight1:weight2);
                        weightnorm += log(exp(weight1 - weightnorm) + exp(weight2 - weightnorm));
                        pcondmix -= weightnorm;
                        pcondmix = exp(pcondmix);
                        if(pcondmix < 0.0 || pcondmix > 1.0) stop("Mixture sampling error %f\n", pcondmix);
                        
                        // sample from mixture density
                        u = eng() / sitmo::prng::max();
                        if(u < pcondmix) {
                        
                            // sample MD conditional on simulator
                            sigma2y = 2.0 * a_dis + 2.0 * b_dis * u1_night(5, j, l) * pI1pI1D;
                            muy = varnorm(t, w) * DIinc1(j, l) + sigma2y * munorm(t, w);
                            muy /= (varnorm(t, w) + sigma2y);
                            sigma2y = sigma2y * varnorm(t, w) / (varnorm(t, w) + sigma2y);
                            u = rtnorm_cpp(
                                muy, 
                                sqrt(sigma2y),
                                -0.5,
                                u1_night(5, j, l) + 0.5,
                                eng
                            );
                            DIinc(j, l) = (int) round(u);
                            if(DIinc(j, l) < 0 || DIinc(j, l) > u1_night(5, j, l)) stop("Error in DIinc DI = %d u = %f uorig = %f LB = %f UB = %f mu = %f sigma2 = %f\n", DIinc(j, l), u, (u - muy) / sqrt(sigma2y), (-0.5 - muy) / sqrt(sigma2y), (u1_night(5, j, l) + 0.5 - muy) / sqrt(sigma2y), muy, sigma2y);
                        } else {
                            
                            // sample MD conditional on simulator
                            muy = (double) DIinc1(j, l);
                            sigma2y = 2.0 * a_dis + 2.0 * b_dis * u1_night(5, j, l) * pI1pI1D;
                            u = rtnorm_cpp(
                                muy, 
                                sqrt(sigma2y),
                                -0.5,
                                u1_night(5, j, l) + 0.5,
                                eng
                            );
                            DIinc(j, l) = (int) round(u);
                            if(DIinc(j, l) < 0 || DIinc(j, l) > u1_night(5, j, l)) stop("Error in DIinc DI = %d u = %f uorig = %f LB = %f UB = %f mu = %f sigma2 = %f\n", DIinc(j, l), u, (u - muy) / sqrt(sigma2y), (-0.5 - muy) / sqrt(sigma2y), (u1_night(5, j, l) + 0.5 - muy) / sqrt(sigma2y), muy, sigma2y);
                        }
                        
                        // particle weights
                        if(t == 0) weights(i) += twistnorm[i](w);
                        
                        // observation error for given incidence
                        sigma2y = a1 + a2 + 2.0 * b * DIinc(j, l);
                        weights(i) += ldtnorm_cpp(
                            obsInc(w) - DIinc(j, l), 
                            a1 - a2, 
                            sqrt(sigma2y), 
                            -DIinc(j, l), 
                            std::numeric_limits<double>::infinity()
                        );
                        
                        // adjust for twisting densities
                        weight1 = log(pmix) + R::dnorm(u, munorm(t, w), sqrt(varnorm(t, w)), 1);
                        weight2 = log(1.0 - pmix);
                        weightnorm = (weight1 > weight2 ? weight1:weight2);
                        weightnorm += log(exp(weight1 - weightnorm) + exp(weight2 - weightnorm));
                        weights(i) -= weightnorm;
                        
                        // store incidence for redistribution
                        tempMD(6, j, l) = DIinc1(j, l);
                        
                        // increment counter
                        w++;
                    }
                }
                
                // DH       
                for(j = 0; j < nages; j++) {
                        
                    // set transition probability
                    pHpHD = pars(j + 8 * nages + 2) * pars(j + 9 * nages + 2);
                    for(l = 0; l < nlads; l++) {
                    
                        // sample simulator from twisted density
                        arma::vec tdensx (tempdensx[i][w].n_elem);
                        tdensx = exp(tempdensx[i][w] - twistnorm[i](w));
                        tdensx = tdensx / sum(tdensx);
                        DHinc1(j, l) = rmultinom_cpp(tdensx, eng);
                        
                        if(tempdensx[i][w].n_elem != (u1_night(9, j, l) + 1)) stop("Grrr DH\n");
                        
                        // extract mixture probability
//                        weight1 = tempdensx[i][w](DHinc1(j, l));
//                        weight1 -= R::dbinom(DHinc1(j, l), u1_night(9, j, l), pHpHD, 1);
//                        weight2 = log(1.0 - pmix) - weight1;
//                        weight1 = log(1.0 - exp(weight2));
//                        pcondmix = exp(weight1);

//                        
                        // Gaussian correction
                        sigma2y = 2.0 * a_dis + 2.0 * b_dis * u1_night(9, j, l) * pHpHD;
                        weight1 = -0.5 * pow(DHinc1(j, l) - munorm(t, w), 2.0) / (sigma2y + varnorm(t, w));
                        weight1 -= 0.5 * log(2.0 * M_PI * (sigma2y + varnorm(t, w)));
                        
                        // target normalising constant correction  
                        weight1 -= ltnormconst_cpp((double) DHinc1(j, l), sqrt(sigma2y), 0.0, u1_night(5, j, l));
                        
                        // integrate over twisting function densities
                        muy = varnorm(t, w) * DHinc1(j, l) + sigma2y * munorm(t, w);
                        muy /= (varnorm(t, w) + sigma2y);
                        sigma2y = sigma2y * varnorm(t, w) / (varnorm(t, w) + sigma2y);
                        weight1 += ltnormconst_cpp(muy, sqrt(sigma2y), 0.0, u1_night(9, j, l));
                        pcondmix = weight1 + log(pmix);
                        
                        // adjust for mixture
                        weight1 = pcondmix;
                        weight2 = log(1.0 - pmix);
                        weightnorm = (weight1 > weight2 ? weight1:weight2);
                        weightnorm += log(exp(weight1 - weightnorm) + exp(weight2 - weightnorm));
                        pcondmix -= weightnorm;
                        pcondmix = exp(pcondmix);
                        if(pcondmix < 0.0 || pcondmix > 1.0) stop("Mixture sampling error %f\n", pcondmix);
                        
                        // sample from mixture density
                        u = eng() / sitmo::prng::max();
                        if(u < pcondmix) {
                        
                            // sample MD conditional on simulator
                            sigma2y = 2.0 * a_dis + 2.0 * b_dis * u1_night(9, j, l) * pHpHD;
                            muy = varnorm(t, w) * DHinc1(j, l) + sigma2y * munorm(t, w);
                            muy /= (varnorm(t, w) + sigma2y);
                            sigma2y = sigma2y * varnorm(t, w) / (varnorm(t, w) + sigma2y);
                            u = rtnorm_cpp(
                                muy, 
                                sqrt(sigma2y),
                                -0.5,
                                u1_night(9, j, l) + 0.5,
                                eng
                            );
                            DHinc(j, l) = (int) round(u);
                            if(DHinc(j, l) < 0 || DHinc(j, l) > u1_night(9, j, l)) stop("Error in DHinc\n");
                        } else {
                            
                            // sample MD conditional on simulator
                            muy = (double) DHinc1(j, l);
                            sigma2y = 2.0 * a_dis + 2.0 * b_dis * u1_night(9, j, l) * pHpHD;
                            u = rtnorm_cpp(
                                muy, 
                                sqrt(sigma2y),
                                -0.5,
                                u1_night(9, j, l) + 0.5,
                                eng
                            );
                            DHinc(j, l) = (int) round(u);
                            if(DHinc(j, l) < 0 || DHinc(j, l) > u1_night(9, j, l)) stop("Error in DHinc\n");
                        }
                        
                        // particle weights
                        if(t == 0) weights(i) += twistnorm[i](w);
                        
                        // observation error for given incidence
                        sigma2y = a1 + a2 + 2.0 * b * DHinc(j, l);
                        weights(i) += ldtnorm_cpp(
                            obsInc(w) - DHinc(j, l), 
                            a1 - a2, 
                            sqrt(sigma2y), 
                            -DHinc(j, l), 
                            std::numeric_limits<double>::infinity()
                        );
                        
                        // adjust for twisting densities
                        weight1 = log(pmix) + R::dnorm(u, munorm(t, w), sqrt(varnorm(t, w)), 1);
                        weight2 = log(1.0 - pmix);
                        weightnorm = (weight1 > weight2 ? weight1:weight2);
                        weightnorm += log(exp(weight1 - weightnorm) + exp(weight2 - weightnorm));
                        weights(i) -= weightnorm;
                        
                        // store incidence for redistribution
                        tempMD(11, j, l) = DHinc1(j, l);
                        
                        // increment counter
                        w++;
                    }
                }
                
                // now redistribute simulator incidence across cohorts
                redistribution(i, nages, nlads, tempMD, ncohorts, u1, u1_new, eng, 1);
                
                // sample remaining x values from conditional simulator
                discreteStochModel((int) i, condpars, t - 1, t, u1_moves, u1_new, u1_day, u1_night1, N_day, N_night, pinf, origE, C, eng);
                
                // set model discrepancy counts for later re-distribution
                tempMD.zeros();
                for(l = 0; l < nlads; l++) {
                    for(j = 0; j < nages; j++) {
                        tempMD(6, j, l) = DIinc(j, l) - DIinc1(j, l);
                        tempMD(11, j, l) = DHinc(j, l) - DHinc1(j, l);
                    }
                }
            } else {
            
                // run model and return u1
                discreteStochModel((int) i, pars, t - 1, t, u1_moves, u1_new, u1_day, u1_night1, N_day, N_night, pinf, origE, C, eng);
            
                // cols: c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")
                //          0,   1,   2,   3,    4,   5,    6,    7,    8,    9,   10,   11
                
                // adjust states according to model discrepancy
                                
                // DH (MD on incidence)
                std::strcpy(str1, "DHinc");
                // aggregate incidence to LAD-level
                for(j = 0; j < nages; j++) {                    
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        DHinc(j, u1_moves(l, 0) - 1) += (u1_new[i](11, j, l) - u1[i](11, j, l));
                    }
                }
                for(j = 0; j < nages; j++) {
                
                    // extract transition probabilities
                    pHpHD = pars(j + 8 * nages + 2) * pars(j + 9 * nages + 2);
                    for(l = 0; l < nlads; l++) {
                    
                        // sample MD conditional on simulator
                        sigma2y = 2.0 * a_dis + 2.0 * b_dis * u1_night(9, j, l) * pHpHD;
                        muy = (double) DHinc(j, l);
                        int s = rdtnorm_cpp(
                            muy, 
                            sqrt(sigma2y),
                            0.0,
                            u1_night(9, j, l),
                            eng
                        );
                        tempMD(11, j, l) = s - DHinc(j, l);
                        DHinc(j, l) = s;
                        if(DHinc(j, l) < 0 || DHinc(j, l) > u1_night(9, j, l)) stop("DHinc error\n");
                        
                        // observation error for given incidence
                        if(PF == 1) {
                            muy = a1 - a2;
                            sigma2y = a1 + a2 + 2.0 * b * DHinc(j, l);
                            weights(i) += ldtnorm_cpp(
                                obsInc(nlads * nages + j * nlads + l) - DHinc(j, l), 
                                muy, 
                                sqrt(sigma2y), 
                                -DHinc(j, l), 
                                std::numeric_limits<double>::infinity()
                            );
                        }
                    }
                }
                
                // DI (MD on incidence)
                std::strcpy(str1, "DIinc");
                // aggregate incidence to LAD-level
                for(j = 0; j < nages; j++) {                    
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        DIinc(j, u1_moves(l, 0) - 1) += u1_new[i](6, j, l) - u1[i](6, j, l);
                    }
                }
                for(j = 0; j < nages; j++) {
                    
                    // set transition probability
                    pI1pI1D = pars(j + 4 * nages + 2) * pars(j + 6 * nages + 2);
                    for(l = 0; l < nlads; l++) {
                        
                        // sample MD conditional on simulator
                        sigma2y = 2.0 * a_dis + 2.0 * b_dis * u1_night(5, j, l) * pI1pI1D;
                        muy = (double) DIinc(j, l);
                        int s = rdtnorm_cpp(
                            muy, 
                            sqrt(sigma2y),
                            0.0,
                            u1_night(5, j, l),
                            eng
                        );
                        tempMD(6, j, l) = s - DIinc(j, l);
                        DIinc(j, l) = s;
                        if(DIinc(j, l) < 0 || DIinc(j, l) > u1_night(5, j, l)) stop("DIinc error\n");
                        
                        // observation error for given incidence
                        if(PF == 1) {
                            muy = a1 - a2;
                            sigma2y = a1 + a2 + 2.0 * b * DIinc(j, l);
                            weights(i) += ldtnorm_cpp(
                                obsInc(j * nlads + l) - DIinc(j, l), 
                                muy, 
                                sqrt(sigma2y), 
                                -DIinc(j, l), 
                                std::numeric_limits<double>::infinity()
                            );
                        }
                    }
                }
            }
            
            // calculate remaining MD terms
            
            // RH given DH (MD on incidence)
            std::strcpy(str1, "RHinc");
            // aggregate incidence to LAD-level
            for(j = 0; j < nages; j++) {                    
                for(l = 0; l < u1_moves.n_rows; l++) {
                    RHinc(j, u1_moves(l, 0) - 1) += u1_new[i](10, j, l) - u1[i](10, j, l);
                }
            }
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    sigma2y = 2.0 * a_dis + 2.0 * b_dis * RHinc(j, l);
                    tempMD(10, j, l) = rdtnorm_cpp(
                        0.0, 
                        sqrt(sigma2y),
                        -RHinc(j, l),
                        u1_night(9, j, l) - DHinc(j, l) - RHinc(j, l),
                        eng
                    );
                    RHinc(j, l) += tempMD(10, j, l);
                    if(RHinc(j, l) < 0 || RHinc(j, l) > u1_night(9, j, l) - DHinc(j, l)) stop("RHinc error\n");
                }
            }
                
            // H given later
            std::strcpy(str1, "H");
            // aggregate count to LAD-level
            for(j = 0; j < nages; j++) {
                for(l = 0; l < u1_moves.n_rows; l++) {
                    Hinc(j, u1_moves(l, 0) - 1) += u1_new[i](9, j, l);
                }
            }
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    // let Hinc be count MD
                    sigma2y = 2.0 * a_dis + 2.0 * b_dis * Hinc(j, l);
                    tempMD(9, j, l) = rdtnorm_cpp(
                        0.0, 
                        sqrt(sigma2y),
                        -Hinc(j, l) + u1_night(9, j, l) - DHinc(j, l) - RHinc(j, l),
                        u1_night(5, j, l) - DIinc(j, l) - Hinc(j, l) + u1_night(9, j, l) - DHinc(j, l) - RHinc(j, l),
                        eng
                    );
                    Hinc(j, l) += tempMD(9, j, l);
                    if(Hinc(j, l) < 0) stop("Hinc error\n");
                }
            }
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    Hinc(j, l) = Hinc(j, l) - u1_night(9, j, l) + DHinc(j, l) + RHinc(j, l);
                    if(Hinc(j, l) < 0 || Hinc(j, l) > u1_night(5, j, l) - DIinc(j, l)) stop("Hinc error1 %d\n", Hinc(j, l));
                }
            }
                
            // RI (MD on incidence)
            std::strcpy(str1, "RIinc");
            // aggregate incidence to LAD-level
            for(j = 0; j < nages; j++) {                    
                for(l = 0; l < u1_moves.n_rows; l++) {
                    RIinc(j, u1_moves(l, 0) - 1) += u1_new[i](8, j, l) - u1[i](8, j, l);
                }
            }
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    sigma2y = 2.0 * a_dis + 2.0 * b_dis * RIinc(j, l);
                    tempMD(8, j, l) = rdtnorm_cpp(
                        0.0, 
                        sqrt(sigma2y),
                        -RIinc(j, l),
                        u1_night(7, j, l) - RIinc(j, l),
                        eng
                    );
                    RIinc(j, l) += tempMD(8, j, l);
                    if(RIinc(j, l) < 0 || RIinc(j, l) > u1_night(7, j, l)) stop("RIinc error\n");
                }
            }
                
            // I2 given later
            std::strcpy(str1, "I2inc");
            // aggregate count to LAD-level
            for(j = 0; j < nages; j++) {
                for(l = 0; l < u1_moves.n_rows; l++) {
                    I2inc(j, u1_moves(l, 0) - 1) += u1_new[i](7, j, l);
                }
            }
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    sigma2y = 2.0 * a_dis + 2.0 * b_dis * I2inc(j, l);
                    tempMD(7, j, l) = rdtnorm_cpp(
                        0.0, 
                        sqrt(sigma2y),
                        -I2inc(j, l) + u1_night(7, j, l) - RIinc(j, l),
                        u1_night(5, j, l) - DIinc(j, l) - Hinc(j, l) - I2inc(j, l) + u1_night(7, j, l) - RIinc(j, l),
                        eng
                    );
                    I2inc(j, l) += tempMD(7, j, l);
                    if(I2inc(j, l) < 0) stop("I2inc error\n");
                }
            }                
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    I2inc(j, l) = I2inc(j, l) - u1_night(7, j, l) + RIinc(j, l);
                    if(I2inc(j, l) < 0 || I2inc(j, l) > u1_night(5, j, l) - DIinc(j, l) - Hinc(j, l)) stop("I2inc error\n");
                }
            }
            
            // I1 given later
            std::strcpy(str1, "I1inc");
            // aggregate count to LAD-level
            for(j = 0; j < nages; j++) {
                for(l = 0; l < u1_moves.n_rows; l++) {
                    I1inc(j, u1_moves(l, 0) - 1) += u1_new[i](5, j, l);
                }
            }
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    sigma2y = 2.0 * a_dis + 2.0 * b_dis * I1inc(j, l);
                    tempMD(5, j, l) = rdtnorm_cpp(
                        0.0, 
                        sqrt(sigma2y),
                        -I1inc(j, l) + u1_night(5, j, l) - I2inc(j, l) - DIinc(j, l) - Hinc(j, l),
                        u1_night(4, j, l) - I1inc(j, l) + u1_night(5, j, l) - I2inc(j, l) - DIinc(j, l) - Hinc(j, l),
                        eng
                    );
                    I1inc(j, l) += tempMD(5, j, l);
                    if(I1inc(j, l) < 0) stop("I1inc error\n");
                }
            }                
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    I1inc(j, l) = I1inc(j, l) - u1_night(5, j, l) + I2inc(j, l) + DIinc(j, l) + Hinc(j, l);
                    if(I1inc(j, l) < 0 || I1inc(j, l) > u1_night(4, j, l)) stop("I1inc error\n");
                }
            }
            
            // P given later
            std::strcpy(str1, "Pinc");
            // aggregate count to LAD-level
            for(j = 0; j < nages; j++) {
                for(l = 0; l < u1_moves.n_rows; l++) {
                    Pinc(j, u1_moves(l, 0) - 1) += u1_new[i](4, j, l);
                }
            }
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    sigma2y = 2.0 * a_dis + 2.0 * b_dis * Pinc(j, l);
                    tempMD(4, j, l) = rdtnorm_cpp(
                        0.0, 
                        sqrt(sigma2y),
                        -Pinc(j, l) + u1_night(4, j, l) - I1inc(j, l),
                        u1_night(1, j, l) - Pinc(j, l) + u1_night(4, j, l) - I1inc(j, l),
                        eng
                    );
                    Pinc(j, l) += tempMD(4, j, l);
                    if(Pinc(j, l) < 0) stop("Pinc error\n");
                }
            }                
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    Pinc(j, l) = Pinc(j, l) - u1_night(4, j, l) + I1inc(j, l);
                    if(Pinc(j, l) < 0 || Pinc(j, l) > u1_night(1, j, l)) stop("Pinc error\n");
                }
            }
            
            // RA (MD on incidence)
            std::strcpy(str1, "RAinc");
            // aggregate incidence to LAD-level
            for(j = 0; j < nages; j++) {                    
                for(l = 0; l < u1_moves.n_rows; l++) {
                    RAinc(j, u1_moves(l, 0) - 1) += u1_new[i](3, j, l) - u1[i](3, j, l);
                }
            }
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    sigma2y = 2.0 * a_dis + 2.0 * b_dis * RAinc(j, l);
                    tempMD(3, j, l) = rdtnorm_cpp(
                        0.0, 
                        sqrt(sigma2y),
                        -RAinc(j, l),
                        u1_night(2, j, l) - RAinc(j, l),
                        eng
                    );
                    RAinc(j, l) += tempMD(3, j, l);
                    if(RAinc(j, l) < 0 || RAinc(j, l) > u1_night(2, j, l)) stop("RAinc error\n");
                }
            }
            
            // A given later
            std::strcpy(str1, "Ainc");
            // aggregate count to LAD-level
            for(j = 0; j < nages; j++) {
                for(l = 0; l < u1_moves.n_rows; l++) {
                    Ainc(j, u1_moves(l, 0) - 1) += u1_new[i](2, j, l);
                }
            }
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    sigma2y = 2.0 * a_dis + 2.0 * b_dis * Ainc(j, l);
                    tempMD(2, j, l) = rdtnorm_cpp(
                        0.0, 
                        sqrt(sigma2y),
                        -Ainc(j, l) + u1_night(2, j, l) - RAinc(j, l),
                        u1_night(1, j, l) - Pinc(j, l) - Ainc(j, l) + u1_night(2, j, l) - RAinc(j, l),
                        eng
                    );
                    Ainc(j, l) += tempMD(2, j, l);
                }
            }                
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    Ainc(j, l) = Ainc(j, l) - u1_night(2, j, l) + RAinc(j, l);
                    if(Ainc(j, l) < 0 || Ainc(j, l) > u1_night(1, j, l) - Pinc(j, l)) stop("Ainc error\n");
                }
            }
            
            // E given later
            std::strcpy(str1, "Einc");
            // aggregate count to LAD-level
            for(j = 0; j < nages; j++) {
                for(l = 0; l < u1_moves.n_rows; l++) {
                    Einc(j, u1_moves(l, 0) - 1) += u1_new[i](1, j, l);
                }
            }
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    sigma2y = 2.0 * a_dis + 2.0 * b_dis * Einc(j, l);
                    tempMD(1, j, l) = rdtnorm_cpp(
                        0.0, 
                        sqrt(sigma2y),
                        -Einc(j, l) + u1_night(1, j, l) - Pinc(j, l) - Ainc(j, l),
                        u1_night(0, j, l) - Einc(j, l) + u1_night(1, j, l) - Pinc(j, l) - Ainc(j, l),
                        eng
                    );
                    Einc(j, l) += tempMD(1, j, l);
                    if(Einc(j, l) < 0) stop("'Einc' error\n");
                }
            }
            
            // re-distribute incidence across cohorts
            redistribution(i, nages, nlads, tempMD, ncohorts, u1, u1_new, eng, 0);
            
            // calculate new normalising constants
            if(twist == 1 && t < (ndays - 1)) {
            
                // aggregate incidence to LAD-level
                u1_night.zeros();
                for(j = 0; j < nages; j++) {                    
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        u1_night(9, j, u1_moves(l, 0) - 1) += u1_new[i](9, j, l);
                        u1_night(5, j, u1_moves(l, 0) - 1) += u1_new[i](5, j, l);
                    }
                }
                
                // start counter
                int w = 0;
                
                // set auxiliary variables
                double muy, sigma2y, pI1pI1D, pHpHD, weight1, weight2, weightnorm;
                
                // DI      
                for(j = 0; j < nages; j++) {
                    
                    // extract transition probabilities
                    pI1pI1D = pars(j + 4 * nages + 2) * pars(j + 6 * nages + 2);
                    for(l = 0; l < nlads; l++) {
                    
                        // set up auxiliary matrix for sampling
                        tempdensx[i][w] = arma::vec (u1_night(5, j, l) + 1);

                        // loop over x values
                        for(int s = 0; s <= u1_night(5, j, l); s++) {  
                            
                            // Gaussian correction
                            sigma2y = 2.0 * a_dis + 2.0 * b_dis * u1_night(5, j, l) * pI1pI1D;
                            tempdensx[i][w](s) = -0.5 * pow(s - munorm(t + 1, w), 2.0) / (sigma2y + varnorm(t + 1, w));
                            tempdensx[i][w](s) -= 0.5 * log(2.0 * M_PI * (sigma2y + varnorm(t + 1, w)));
                            
                            // target normalising constant correction  
                            tempdensx[i][w](s) -= ltnormconst_cpp((double) s, sqrt(sigma2y), 0.0, u1_night(5, j, l));
                            
                            // integrate over twisting function densities
                            muy = varnorm(t + 1, w) * s + sigma2y * munorm(t + 1, w);
                            muy /= (varnorm(t + 1, w) + sigma2y);
                            sigma2y = sigma2y * varnorm(t + 1, w) / (varnorm(t + 1, w) + sigma2y);
                            tempdensx[i][w](s) += ltnormconst_cpp(muy, sqrt(sigma2y), 0.0, u1_night(5, j, l));
                            
                            // correct for mixture
                            weight1 = log(pmix) + tempdensx[i][w](s);
                            weight2 = log(1.0 - pmix);
                            weightnorm = (weight1 > weight2 ? weight1:weight2);
                            weightnorm += log(exp(weight1 - weightnorm) + exp(weight2 - weightnorm));
                            tempdensx[i][w](s) = weightnorm;
                            
                            // simulator density
                            tempdensx[i][w](s) += R::dbinom(s, u1_night(5, j, l), pI1pI1D, 1);
                        }
                        
                        // calculate normalising constant
                        twistnorm[i](w) = log_sum_exp(tempdensx[i][w], 0);
                        
                        // update weights
                        weights(i) += twistnorm[i](w);
                        
                        // increment counter
                        w++;
                    }
                }
                
                // DH
                for(j = 0; j < nages; j++) {
                    
                    // extract transition probabilities
                    pHpHD = pars(j + 8 * nages + 2) * pars(j + 9 * nages + 2);
                    for(l = 0; l < nlads; l++) {
                    
                        // set up auxiliary matrix for sampling
                        tempdensx[i][w] = arma::vec (u1_night(9, j, l) + 1);

                        // loop over x values
                        for(int s = 0; s <= u1_night(9, j, l); s++) {   
                            
                            // Gaussian correction
                            sigma2y = 2.0 * a_dis + 2.0 * b_dis * u1_night(9, j, l) * pHpHD;
                            tempdensx[i][w](s) = -0.5 * pow(s - munorm(t + 1, w), 2.0) / (sigma2y + varnorm(t + 1, w));
                            tempdensx[i][w](s) -= 0.5 * log(2.0 * M_PI * (sigma2y + varnorm(t + 1, w)));
                            
                            // target normalising constant correction  
                            tempdensx[i][w](s) -= ltnormconst_cpp((double) s, sqrt(sigma2y), 0.0, u1_night(9, j, l));
                            
                            // integrate over twisting function densities
                            muy = varnorm(t + 1, w) * s + sigma2y * munorm(t + 1, w);
                            muy /= (varnorm(t + 1, w) + sigma2y);
                            sigma2y = sigma2y * varnorm(t + 1, w) / (varnorm(t + 1, w) + sigma2y);
                            tempdensx[i][w](s) += ltnormconst_cpp(muy, sqrt(sigma2y), 0.0, u1_night(9, j, l));
                            
                            // correct for mixture
                            weight1 = log(pmix) + tempdensx[i][w](s);
                            weight2 = log(1.0 - pmix);
                            weightnorm = (weight1 > weight2 ? weight1:weight2);
                            weightnorm += log(exp(weight1 - weightnorm) + exp(weight2 - weightnorm));
                            tempdensx[i][w](s) = weightnorm;
                            
                            // simulator density
                            tempdensx[i][w](s) += R::dbinom(s, u1_night(9, j, l), pHpHD, 1);
                        }
                        
                        // calculate normalising constant
                        twistnorm[i](w) = log_sum_exp(tempdensx[i][w], 0);
                        
                        // update weights
                        weights(i) += twistnorm[i](w);
                        
                        // increment counter
                        w++;
                    }
                }
            }
            
            // advance seed
            seeds((arma::uword) omp_get_thread_num()) = eng();
        } 
        
        // save particles if necessary
        if(saveAll != 0) {
            // set up print string for debugging
            char str1[80];
            std::strcpy(str1, "save");
            double muy, sigma2y, u;
            if(saveAll == 1) {
                for(i = 0; i < npart; i++) {
                    // extract just counts for DI and DH
                    u1_night_reduced.zeros();
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        for(j = 0; j < nages; j++) {
                            u1_night_reduced(0, j, (arma::uword) u1_moves(l, 0) - 1) += u1_new[i](6, j, l);
                            u1_night_reduced(1, j, (arma::uword) u1_moves(l, 0) - 1) += u1_new[i](11, j, l);
                            // incidence
                            u1_night_reduced(2, j, (arma::uword) u1_moves(l, 0) - 1) += (u1_new[i](6, j, l) - u1[i](6, j, l));
                            u1_night_reduced(3, j, (arma::uword) u1_moves(l, 0) - 1) += (u1_new[i](11, j, l) - u1[i](11, j, l));
                        }
                    }
                    // adjust for observation error
                    muy = a1 - a2;
                    for(l = 0; l < nlads; l++) {
                        for(j = 0; j < nages; j++) {
                            sigma2y = a1 + a2 + 2.0 * b * u1_night_reduced(2, j, l);
                            u1_night_reduced(2, j, l) += rdtnorm_cpp(
                                muy, 
                                sqrt(sigma2y),
                                -u1_night_reduced(2, j, l),
                                std::numeric_limits<double>::infinity(),
                                engSerial
                            );
                            
                            sigma2y = a1 + a2 + 2.0 * b * u1_night_reduced(3, j, l);
                            u1_night_reduced(3, j, l) += rdtnorm_cpp(
                                muy, 
                                sqrt(sigma2y),
                                -u1_night_reduced(3, j, l),
                                std::numeric_limits<double>::infinity(),
                                engSerial
                            );
                            // cumulate
                            u1_night_reduced(2, j, l) += u1_night_obs[i](0, j, l);
                            u1_night_reduced(3, j, l) += u1_night_obs[i](1, j, l);
                            u1_night_obs[i](0, j, l) = u1_night_reduced(2, j, l);
                            u1_night_obs[i](1, j, l) = u1_night_reduced(3, j, l);
                        }
                    }
                    if(writeExt == 0) {
                        out[i + npart * (t + 1)] = u1_night_reduced;
                    } else {
                        std::sprintf(file_name, "saveOut/p_%u.csv", i);
                        file.open(file_name, std::ios::app);
                        for(l = 0; l < nlads; l++) {
                            for(arma::uword r = 0; r < 4; r++) {
                                file << t + 1 << ", " << r << ", ";
                                for(j = 0; j < nages; j++) {
                                    file << u1_night_reduced(r, j, l) << ", ";
                                }
                                file << l + 1 << "\n";
                            }
                        }
                        file.close();
                    }
                }
            } else {
                for(i = 0; i < npart; i++) {
                    u1_night_full.zeros();
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        for(j = 0; j < nages; j++) {
                            for(k = 0; k < nclasses; k++) {
                                u1_night_full(k, j, (arma::uword) u1_moves(l, 0) - 1) += u1_new[i](k, j, l);
                            }
                            // incidence
                            u1_night_full(nclasses, j, (arma::uword) u1_moves(l, 0) - 1) += (u1_new[i](6, j, l) - u1[i](6, j, l));
                            u1_night_full(nclasses + 1, j, (arma::uword) u1_moves(l, 0) - 1) += (u1_new[i](11, j, l) - u1[i](11, j, l));
                        }
                    }
                    // adjust for observation error
                    muy = a1 - a2;
                    for(l = 0; l < nlads; l++) {
                        for(j = 0; j < nages; j++) {
                            sigma2y = a1 + a2 + 2.0 * b * u1_night_full(nclasses, j, l);
                            u1_night_full(nclasses, j, l) += rdtnorm_cpp(
                                muy, 
                                sqrt(sigma2y),
                                -u1_night_full(nclasses, j, l),
                                std::numeric_limits<double>::infinity(),
                                engSerial
                            );
                            
                            sigma2y = a1 + a2 + 2.0 * b * u1_night_full(nclasses + 1, j, l);
                            u1_night_full(nclasses + 1, j, l) += rdtnorm_cpp(
                                muy, 
                                sqrt(sigma2y),
                                -u1_night_full(nclasses + 1, j, l),
                                std::numeric_limits<double>::infinity(),
                                engSerial
                            );
                            // cumulate
                            u1_night_full(nclasses, j, l) += u1_night_obs[i](0, j, l);
                            u1_night_full(nclasses + 1, j, l) += u1_night_obs[i](1, j, l);
                            u1_night_obs[i](0, j, l) = u1_night_full(nclasses, j, l);
                            u1_night_obs[i](1, j, l) = u1_night_full(nclasses + 1, j, l);
                        }
                    }
                    if(writeExt == 0) {
                        out[i + npart * (t + 1)] = u1_night_full;
                    } else {
                        std::sprintf(file_name, "saveOut/p_%u.csv", i);
                        file.open(file_name, std::ios::app);
                        for(l = 0; l < nlads; l++) {
                            for(arma::uword r = 0; r < (nclasses + 2); r++) {
                                file << t + 1 << ", " << r << ", ";
                                for(j = 0; j < nages; j++) {
                                    file << u1_night_full(r, j, l) << ", ";
                                }
                                file << l + 1 << "\n";
                            }
                        }
                        file.close();
                    }
                }
            }
        }
        
        if(PF == 1) {
            // calculate log-likelihood contribution
            ll += log_sum_exp(weights, 1);
            
            // if zero likelihood then return
            if(!arma::is_finite(ll)) {
                if(saveAll == 0 || writeExt == 1) {
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
                inds(i) = (arma::uword) rmultinom_cpp(weights, engSerial);
            }
            if(returnPsi == 1) {
                // save particle summaries for twisting functions (BEFORE resampling for more variation)
                for(i = 0; i < npart; i++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        for(j = 0; j < nages; j++) {
                            psi(t * npart * 4 + i * 4, j, (arma::uword) u1_moves(l, 0) - 1) += (u1_new[i](6, j, l) - u1[i](6, j, l));
                            psi(t * npart * 4 + i * 4 + 1, j, (arma::uword) u1_moves(l, 0) - 1) += (u1_new[i](11, j, l) - u1[i](11, j, l));
                            psi(t * npart * 4 + i * 4 + 2, j, (arma::uword) u1_moves(l, 0) - 1) += u1_new[i](5, j, l);
                            psi(t * npart * 4 + i * 4 + 3, j, (arma::uword) u1_moves(l, 0) - 1) += u1_new[i](9, j, l);
                        }
                    }
                }
            }
            for(i = 0; i < npart; i++) {
                u1[i] = u1_new[inds(i)];
                u1_night_obs1[i] = u1_night_obs[inds(i)];
                if(twist == 1) {
                    twistnorm1[i] = twistnorm[inds(i)];
                    tempdensx1[i] = tempdensx[inds(i)];
                }
            }
            // copy in order to pass by reference
            for(i = 0; i < npart; i++) {
                u1_new[i] = u1[i];
                u1_night_obs[i] = u1_night_obs1[i];
                if(twist == 1) {
                    twistnorm[i] = twistnorm1[i];
                    tempdensx[i] = tempdensx1[i];
                }
            }
        } else {
            for(i = 0; i < npart; i++) {
                u1[i] = u1_new[i];
            }
        }
        
        //calculate block run time
        timer.step("");
        NumericVector res(timer);
        
        // calculate ESS
        double ESS = 0.0;
        for(i = 0; i < npart; i++) {
            ESS += pow(weights(i), 2.0);
        }
        ESS = 1.0 / ESS;
        ESS = ESS / ((double) npart);
        
        Rprintf("t = %d / %d RESS = %.2f time = %.2f secs \n", t + 1, ndays, ESS, (res[timer_cnt] / 1e9) - prev_time);
        
        //reset timer and acceptance rate counter
        prev_time = res[timer_cnt] / 1e9;
        timer_cnt++;
    }
    if(returnPsi == 1) {
        if(saveAll == 0) {
            return List::create(Named("ll") = ll, _["psi"] = psi);
        } else {
            if(writeExt == 0) {
                if(PF == 1) {
                    return List::create(Named("ll") = ll, _["particles"] = out, _["psi"] = psi);
                } else {
                    return List::create(Named("particles") = out, _["psi"] = psi);
                }
            } else {
                if(PF == 1) {
                    return List::create(Named("ll") = ll, _["psi"] = psi);
                } else {
                    return List::create(Named("psi") = psi);
                }
            }
        }
    } else {
        if(saveAll == 0) {
            return List::create(Named("ll") = ll);
        } else {
            if(writeExt == 0) {
                if(PF == 1) {
                    return List::create(Named("ll") = ll, _["particles"] = out);
                } else {
                    return List::create(Named("particles") = out);
                }
            } else {
                if(PF == 1) {
                    return List::create(Named("ll") = ll);
                } else {
                    return List::create(Named("particles") = NA_INTEGER);
                }
            }
        }
    }
}

// function to calculate rates for using in twisting function optimisation
// [[Rcpp::export]]
arma::mat TPF_rates_obs_cpp (arma::ivec data, arma::uword nages, arma::uword nlads, arma::cube psi, 
    arma::uword npart, double a1, double a2, double b, int ncores) {
    
    // set counters
    arma::uword i, j, l;
    
    // set auxiliary objects
    arma::mat twistnorm(npart * 2, data.n_elem); twistnorm.zeros();

#ifdef _OPENMP
#pragma omp parallel for default(none) private(j, l) shared(npart, nages, nlads, a1, a2, b, twistnorm, psi, data)
#endif
    for(i = 0; i < npart; i++) {
        
        // start counter
        int w = 0;
        
        // DI      
        for(j = 0; j < nages; j++) {
                
            for(l = 0; l < nlads; l++) {
                // observation error for given incidence
                double sigma2y = a1 + a2 + 2.0 * b * psi(i * 4, j, l);
                twistnorm(i * 2 + 1, w) = ldtnorm_cpp(
                    data(w) - psi(i * 4, j, l), 
                    a1 - a2, 
                    sqrt(sigma2y), 
                    -psi(i * 4, j, l), 
                    std::numeric_limits<double>::infinity()
                );
                
                // set data
                twistnorm(i * 2, w) = psi(i * 4 + 2, j, l);
                
                // increment counter
                w++;
            }
        }
        
        // DH
        for(j = 0; j < nages; j++) {
            
            for(l = 0; l < nlads; l++) {
                
                // observation error for given incidence
                double sigma2y = a1 + a2 + 2.0 * b * psi(i * 4 + 1, j, l);
                twistnorm(i * 2 + 1, w) = ldtnorm_cpp(
                    data(w) - psi(i * 4 + 1, j, l), 
                    a1 - a2, 
                    sqrt(sigma2y), 
                    -psi(i * 4 + 1, j, l), 
                    std::numeric_limits<double>::infinity()
                );
                
                // set data
                twistnorm(i * 2, w) = psi(i * 4 + 3, j, l);
                
                // increment counter
                w++;
            }
        }
    }
    return twistnorm;
}

// [[Rcpp::export]]
arma::mat TPF_rates_cpp (arma::vec pars, arma::ivec data, arma::uword nages, arma::uword nlads, arma::cube psi, 
    arma::uword npart, arma::vec munorm, arma::vec varnorm, double pmix, double a1, double a2, double b, double a_dis, double b_dis, int ncores) {
    
    // set counters
    arma::uword i, j, l;
    
    // set auxiliary objects
    arma::mat twistnorm(npart * 2, munorm.n_elem); twistnorm.zeros();

#ifdef _OPENMP
#pragma omp parallel for default(none) private(j, l) shared(npart, nages, nlads, pars, a1, a2, b, a_dis, b_dis, munorm, varnorm, twistnorm, psi, data, pmix)
#endif
    for(i = 0; i < npart; i++) {
        
        // start counter
        int w = 0;
        
        // set auxiliary variables
        double muy, sigma2y, weight1, weight2, weightnorm, pI1pI1D, pHpHD;
        
        // DI      
        for(j = 0; j < nages; j++) {
            
            // extract transition probabilities
            pI1pI1D = pars(j + 4 * nages + 2) * pars(j + 6 * nages + 2);
            for(l = 0; l < nlads; l++) {
            
                // set up auxiliary matrix for sampling
                arma::vec tempdensx (psi(i * 4 + 2, j, l) + 1); tempdensx.zeros();

                // loop over x values
                for(int s = 0; s <= psi(i * 4 + 2, j, l); s++) {
                    
                    // Gaussian correction
                    sigma2y = 2.0 * a_dis + 2.0 * b_dis * psi(i * 4 + 2, j, l) * pI1pI1D;
                    tempdensx(s) = -0.5 * pow(s - munorm(w), 2.0) / (sigma2y + varnorm(w));
                    tempdensx(s) -= 0.5 * log(2.0 * M_PI * (sigma2y + varnorm(w)));
                    
                    // target normalising constant correction  
                    tempdensx(s) -= ltnormconst_cpp((double) s, sqrt(sigma2y), 0.0, psi(i * 4 + 2, j, l));
                    
                    // integrate over twisting function densities
                    muy = varnorm(w) * s + sigma2y * munorm(w);
                    muy /= (varnorm(w) + sigma2y);
                    sigma2y = sigma2y * varnorm(w) / (varnorm(w) + sigma2y);
                    tempdensx(s) += ltnormconst_cpp(muy, sqrt(sigma2y), 0.0, psi(i * 4 + 2, j, l));
                    
                    // correct for mixture
                    weight1 = log(pmix) + tempdensx(s);
                    weight2 = log(1.0 - pmix);
                    weightnorm = (weight1 > weight2 ? weight1:weight2);
                    weightnorm += log(exp(weight1 - weightnorm) + exp(weight2 - weightnorm));
                    tempdensx(s) = weightnorm;
                    
                    // simulator density
                    tempdensx(s) += R::dbinom(s, psi(i * 4 + 2, j, l), pI1pI1D, 1);
                }
                
                // calculate normalising constant
                twistnorm(i * 2 + 1, w) = log_sum_exp(tempdensx, 0);
                
                // observation error for given incidence
                sigma2y = a1 + a2 + 2.0 * b * psi(i * 4, j, l);
                twistnorm(i * 2 + 1, w) += ldtnorm_cpp(
                    data(w) - psi(i * 4, j, l), 
                    a1 - a2, 
                    sqrt(sigma2y), 
                    -psi(i * 4, j, l), 
                    std::numeric_limits<double>::infinity()
                );
                
                // set data
                twistnorm(i * 2, w) = psi(i * 4 + 2, j, l);
                
                // increment counter
                w++;
            }
        }
        
        // DH
        for(j = 0; j < nages; j++) {
            
            // extract transition probabilities
            pHpHD = pars(j + 8 * nages + 2) * pars(j + 9 * nages + 2);
            for(l = 0; l < nlads; l++) {
            
                // set up auxiliary matrix for sampling
                arma::vec tempdensx (psi(i * 4 + 3, j, l) + 1); tempdensx.zeros();

                // loop over x values
                for(int s = 0; s <= psi(i * 4 + 3, j, l); s++) {
                    
                    // Gaussian correction
                    sigma2y = 2.0 * a_dis + 2.0 * b_dis * psi(i * 4 + 3, j, l) * pHpHD;
                    tempdensx(s) = -0.5 * pow(s - munorm(w), 2.0) / (sigma2y + varnorm(w));
                    tempdensx(s) -= 0.5 * log(2.0 * M_PI * (sigma2y + varnorm(w)));
                    
                    // target normalising constant correction  
                    tempdensx(s) -= ltnormconst_cpp((double) s, sqrt(sigma2y), 0.0, psi(i * 4 + 3, j, l));
                    
                    // integrate over twisting function densities
                    muy = varnorm(w) * s + sigma2y * munorm(w);
                    muy /= (varnorm(w) + sigma2y);
                    sigma2y = sigma2y * varnorm(w) / (varnorm(w) + sigma2y);
                    tempdensx(s) += ltnormconst_cpp(muy, sqrt(sigma2y), 0.0, psi(i * 4 + 3, j, l));
                    
                    // correct for mixture
                    weight1 = log(pmix) + tempdensx(s);
                    weight2 = log(1.0 - pmix);
                    weightnorm = (weight1 > weight2 ? weight1:weight2);
                    weightnorm += log(exp(weight1 - weightnorm) + exp(weight2 - weightnorm));
                    tempdensx(s) = weightnorm;
                    
                    // simulator density
                    tempdensx(s) += R::dbinom(s, psi(i * 4 + 3, j, l), pHpHD, 1);
                }
                
                // calculate normalising constant
                twistnorm(i * 2 + 1, w) = log_sum_exp(tempdensx, 0); 
                
                // observation error for given incidence
                sigma2y = a1 + a2 + 2.0 * b * psi(i * 4 + 1, j, l);
                twistnorm(i * 2 + 1, w) += ldtnorm_cpp(
                    data(w) - psi(i * 4 + 1, j, l), 
                    a1 - a2, 
                    sqrt(sigma2y), 
                    -psi(i * 4 + 1, j, l), 
                    std::numeric_limits<double>::infinity()
                );
                
                // set data
                twistnorm(i * 2, w) = psi(i * 4 + 3, j, l);
                
                // increment counter
                w++;
            }
        }
    }
    return twistnorm;
}

