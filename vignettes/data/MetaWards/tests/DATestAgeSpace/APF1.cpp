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
    if(!arma::is_finite(ldens)) stop("Something wrong in dTN\n");
    if(std::isinf(UB)) {
        // normalising constant
        temp1 = R::pnorm(LB - 0.5, mu, sigma, 0, 1);
        if(!arma::is_finite(temp1)) stop("Something wrong in dTN LB\n");
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
        if(!arma::is_finite(temp1)) stop("Something wrong in dTN LB UB\n");
        ldens -= temp1;
    }
    return ldens;
}

// truncated discrete Gaussian sampling
int rdtnorm_cpp(double mu, double sigma, double LB, double UB, sitmo::prng &eng) {
    if(std::isinf(LB)) {
        stop("Lower bound of truncated Gaussian must be finite currently\n");
    }
    if(LB > UB) {
        stop("'LB' can't be > 'UB' in rdtnorm_cpp\n");
    }
    if(LB == UB) {
        return (int) LB;
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

// log truncated Gaussian p.d.f.
double ltnorm_cpp(double x, double mu, double sigma, double LB, double UB) {
    if(std::isinf(LB)) {
        stop("Lower bound of truncated Gaussian must be finite currently\n");
    }
    double ldens = R::dnorm(x, mu, sigma, 1);
    if(!arma::is_finite(ldens)) stop("Something wrong in TN\n");
    if(std::isinf(UB)) {
        // normalising constant
        double temp1 = R::pnorm(LB, mu, sigma, 0, 1);
        if(!arma::is_finite(temp1)) stop("Something wrong in TN LB\n");
        ldens -= temp1;
    } else {
        // normalising constant
        double temp1 = R::pnorm(LB , mu, sigma, 1, 1);
        double temp2 = R::pnorm(UB, mu, sigma, 1, 1);
        temp1 = temp2 + log(1.0 - exp(temp1 - temp2));
        if(!arma::is_finite(temp1)) {
            temp1 = R::pnorm(LB, mu, sigma, 0, 1);
            temp2 = R::pnorm(UB, mu, sigma, 0, 1);
            temp2 = temp1 + log(1.0 - exp(temp2 - temp1));
            temp1 = temp2;
        }
        if(!arma::is_finite(temp1)) stop("Something wrong in TN LB UB\n");
        ldens -= temp1;
    }
    return ldens;
}

// expected value of truncated Gaussian p.d.f.
double exptnorm_cpp(double mu, double sigma, double LB, double UB) {
    if(std::isinf(LB)) {
        stop("Lower bound of truncated Gaussian must be finite currently\n");
    }
    if(LB == UB) {
        stop("'LB' shouldn't equal 'UB' in truncated Gaussian mean function\n");
    }
    if(std::isinf(UB)) {
        stop("Mean not yet implemented for left-truncated Gaussian\n");
    } else {
        // normalising constant
        double temp1 = R::pnorm(LB, mu, sigma, 1, 1);
        double temp2 = R::pnorm(UB, mu, sigma, 1, 1);
        temp1 = temp2 + log(1.0 - exp(temp1 - temp2));
        if(!arma::is_finite(temp1)) {
            temp1 = R::pnorm(LB, mu, sigma, 0, 1);
            temp2 = R::pnorm(UB, mu, sigma, 0, 1);
            temp2 = temp1 + log(1.0 - exp(temp2 - temp1));
            temp1 = temp2;
            if(!arma::is_finite(temp1)) stop("Something wrong in mean TN LB UB\n");
        }
        double Z = exp(temp1);
        temp1 = R::dnorm(LB, mu, sigma, 0);
        temp2 = R::dnorm(UB, mu, sigma, 0);
        temp1 = mu + (temp1 - temp2) * sigma / Z;
        if(!arma::is_finite(temp1) || temp1 < LB || temp1 > UB) {
            Rprintf("mu = %f sigma = %f LB = %f UB = %f\n", mu, sigma, LB, UB);
        }
        return temp1;
    }
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
int rbinom_cpp (int n, double p, sitmo::prng &eng, int approx = 1) {
    if(n < 0) {
        Rprintf("n = %d\n", n);
        stop("'n' must be >= 0 in rbinom\n");
    }
    if(p < 0.0 || p > 1.0) {
        Rprintf("p = %f\n", p);
        stop("Must have 0 <= p <= 1 in rbinom\n");
    }
    if(n == 0 || p <= 0.0) return 0;
    if(p >= 1.0) return(n);
    // if approximation turned on then use truncated 
    // Gaussian approximation where appropriate
    if(approx == 1) {
        if(n > 20 && (n * p) > 5 && (n * (1.0 - p)) > 5) {
            double mu = n * p;
            double sigma = sqrt(mu * (1.0 - p));
            int k = rdtnorm_cpp(mu, sigma, 0.0, (double) n, eng);
            return(k);
        }
    }
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
    while(temp < u && k < (p.n_elem - 1)) {
        k++;
        temp += p(k);
    }
    return k;
} 

// simulation model
void discreteStochModel(int ipart, int nclasses, int nages, int nlads, 
                        arma::vec &pars, int tstart, int tstop, 
                        arma::imat &u1_moves, std::vector<arma::icube> &u1, 
                        arma::mat &C, sitmo::prng &eng) {
    
    // u1_moves is a matrix with columns: LADfrom, LADto
    // u1 is a 3D array with dimensions: nclasses x nages x nmoves
    //          each row of u1 must match u1_moves
    
    // set up auxiliary matrix for counts
    int k, n;
    arma::uword i, j, l;
    
    // reconstruct day/night counts
    arma::icube u1_day(nclasses, nages, nlads); u1_day.zeros();
    arma::icube u1_night(nclasses, nages, nlads); u1_night.zeros();
    for(i = 0; i < u1_moves.n_rows; i++) {
        for(j = 0; j < nages; j++) {
            for(l = 0; l < nclasses; l++) {
                u1_day(l, j, (arma::uword) u1_moves(i, 1) - 1) += u1[ipart](l, j, i);
                u1_night(l, j, (arma::uword) u1_moves(i, 0) - 1) += u1[ipart](l, j, i);
            }
        }
    }
    
    // reconstruct population counts
    arma::imat N_day(nages, nlads); N_day.zeros();
    arma::imat N_night(nages, nlads); N_night.zeros();
    for(i = 0; i < nlads; i++) {
        for(j = 0; j < nages; j++) {
            for(l = 0; l < nclasses; l++) {
                N_day(j, i) += u1_day(l, j, i);
                N_night(j, i) += u1_night(l, j, i);
            }
        }
    }
    
    // auxiliary vectors
    arma::mat pinf(nages, nlads); pinf.zeros();
    arma::imat origE(nages, u1_moves.n_rows); origE.zeros();
    
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
        for(i = 0; i < nlads; i++) {
            
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
        for(i = 0; i < nlads; i++) {
            
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
                    }
                }
                
                // I2RI
                k = rbinom_cpp(u1[ipart](7, j, i), probI2(j), eng);
                u1[ipart](7, j, i) -= k;
                u1[ipart](8, j, i) += k;
                
                // I1 out
                n = u1[ipart](5, j, i);
                if(n > 0) {
                    for(l = 0; l < n; l++) {
                        k = rmultinom_cpp(mprobsI1, eng);
                        u1[ipart](5, j, i) -= (k < 3 ? 1:0);
                        u1[ipart](9, j, i) += (k == 0 ? 1:0);
                        u1[ipart](7, j, i) += (k == 1 ? 1:0);
                        u1[ipart](6, j, i) += (k == 2 ? 1:0);
                    }
                }
                
                // PI1
                k = rbinom_cpp(u1[ipart](4, j, i), probP(j), eng);
                u1[ipart](4, j, i) -= k;
                u1[ipart](5, j, i) += k;
                
                // ARA
                k = rbinom_cpp(u1[ipart](2, j, i), probA(j), eng);
                u1[ipart](2, j, i) -= k;
                u1[ipart](3, j, i) += k;
                
                // E out
                if(origE(j, i) > 0) {
                    for(l = 0; l < origE(j, i); l++) {
                        k = rmultinom_cpp(mprobsE, eng);
                        u1[ipart](1, j, i) -= (k < 2 ? 1:0);
                        u1[ipart](2, j, i) += (k == 0 ? 1:0);
                        u1[ipart](4, j, i) += (k == 1 ? 1:0);
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
List APF1_cpp (arma::vec pars, arma::mat C, arma::imat deathInc_lad, arma::imat deathInc_age_region,
    arma::imat hosp_nhsregion, arma::imat hospInc_age_nhsregion, arma::imat lookup, arma::imat age_lookup,
    arma::uword nclasses, arma::uword nages, arma::uword nlads, 
    arma::uword ndeathlads, arma::uword nregions, arma::uword nnhsages, arma::uword nnhsregions,  
    arma::imat u1_moves, arma::ivec ncohorts1, arma::icube u1_comb, 
    arma::icube u2_comb, arma::vec playprobs, arma::ivec ncohorts2, arma::uword ndays, 
    arma::uword npart, int niter, double a1, double a2, double b, double a_dis, double b_dis,
    double sigma2_lad, double sigma2_age_region, double sigma2_nhsregion, double sigma2_age_nhsregion,
    int saveAll, int writeExt, CharacterVector outputName, int PF, int ncores) {
    
    // set counters
    arma::uword i, j, l, k, t = 0;
    
    // split u1 up into different LADs
    std::vector<arma::icube> u1(npart);
    std::vector<arma::icube> u11(npart);
    std::vector<arma::icube> u1_new(npart);
    std::vector<arma::icube> u2(npart);
    std::vector<arma::icube> u2_new(npart);
    for(i = 0; i < npart; i++) {
        u1[i] = u1_comb;
        u1_new[i] = u1_comb;
        u2[i] = u2_comb;
        u2_new[i] = u2_comb;
    }
    
    // set up output objects
    arma::icube u_night_full(nclasses, nages, nlads); u_night_full.zeros();
    
    arma::ivec u_night_lad(ndeathlads); u_night_lad.zeros();
    arma::imat u_night_age_region(nages, nregions); u_night_age_region.zeros();
    arma::imat u_night_age_nhsregion(nnhsages, nnhsregions); u_night_age_nhsregion.zeros();
    arma::ivec u_night_nhsregion(nnhsregions); u_night_nhsregion.zeros();
    
    std::vector<arma::ivec> u_night_lad_cum(npart);
    std::vector<arma::imat> u_night_age_region_cum(npart);
    std::vector<arma::imat> u_night_age_nhsregion_cum(npart);
    std::vector<arma::ivec> u_night_vec_cum1(npart);
    std::vector<arma::imat> u_night_mat_cum1(npart);
    for(i = 0; i < npart; i++) {
        u_night_lad_cum[i] = arma::ivec(ndeathlads); u_night_lad_cum[i].zeros();
        u_night_age_region_cum[i] = arma::imat(nages, nregions); u_night_age_region_cum[i].zeros();
        u_night_age_nhsregion_cum[i] = arma::imat(nnhsages, nnhsregions); u_night_age_nhsregion_cum[i].zeros();
    }
        
    // set up weight vector
    arma::vec weights (npart);
    arma::ivec inds(npart);
    double wnorm = 0.0;
    
    // set up auxiliary objects
    arma::ivec obsInc_lad (deathInc_lad.n_cols); obsInc_lad.zeros();
    arma::ivec obsInc_age_region (deathInc_age_region.n_cols); obsInc_age_region.zeros();
    arma::ivec obsInc_age_nhsregion (hospInc_age_nhsregion.n_cols); obsInc_age_nhsregion.zeros();
    arma::ivec obs_nhsregion (hosp_nhsregion.n_cols); obs_nhsregion.zeros();
    
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
    
    // calculate number of LADs in each region
    arma::ivec nlads_region(nregions); nlads_region.zeros();
    arma::ivec nlads_nhsregion(nnhsregions); nlads_nhsregion.zeros();
    for(i = 0; i < nlads; i++) {
        if(lookup(i, 2) >= 0) {
            nlads_region(lookup(i, 2) - 1)++;
        }
        if(lookup(i, 3) >= 0) {
            nlads_nhsregion(lookup(i, 3) - 1)++;
        }
    }
    
    // check which output required
    List out (npart * (ndays + 1));
    std::ofstream file;
    char file_name[128];
    if(saveAll != 0) {
        // set up print string for debugging
        char str1[80];
        std::strcpy(str1, "save");
        double muy, sigma2y, u;
        for(i = 0; i < npart; i++) {
            // extract just counts for DI and DH
            u_night_lad.zeros();
            u_night_age_region.zeros();
            u_night_age_nhsregion.zeros();
            u_night_nhsregion.zeros();
            for(l = 0; l < u1_moves.n_rows; l++) {
                k = (arma::uword) u1_moves(l, 0) - 1;
                for(j = 0; j < nages; j++) {
                    
                    // death incidence in LADs
                    if(lookup(k, 1) >= 0) {
                        u_night_lad(lookup(k, 1) - 1) += u1[i](6, j, l);
                        u_night_lad(lookup(k, 1) - 1) += u1[i](11, j, l);
                    }
                    
                    // death incidence by age and region
                    if(lookup(k, 2) >= 0) {
                        u_night_age_region(j, lookup(k, 2) - 1) += u1[i](6, j, l);
                        u_night_age_region(j, lookup(k, 2) - 1) += u1[i](11, j, l);
                    }
                    
                    if(lookup(k, 3) >= 0) {
                        // hospital incidence by age and NHS region
                        u_night_age_nhsregion(age_lookup(j, 1) - 1, lookup(k, 3) - 1) += u1[i](9, j, l);
                        // and hospital count by NHS region
                        u_night_nhsregion(lookup(k, 3) - 1) += u1[i](9, j, l);
                    }
                }
            }
            for(l = 0; l < nlads; l++) {
                for(j = 0; j < nages; j++) {
                    
                    // death incidence in LADs
                    if(lookup(l, 1) >= 0) {
                        u_night_lad(lookup(l, 1) - 1) += u2[i](6, j, l);
                        u_night_lad(lookup(l, 1) - 1) += u2[i](11, j, l);
                    }
                    
                    // death incidence by age and region
                    if(lookup(l, 2) >= 0) {
                        u_night_age_region(j, lookup(l, 2) - 1) += u2[i](6, j, l);
                        u_night_age_region(j, lookup(l, 2) - 1) += u2[i](11, j, l);
                    }
                    
                    if(lookup(l, 3) >= 0) {
                        // hospital incidence by age and NHS region
                        u_night_age_nhsregion(age_lookup(j, 1) - 1, lookup(l, 3) - 1) += u2[i](9, j, l);
                        // and hospital count by NHS region
                        u_night_nhsregion(lookup(l, 3) - 1) += u2[i](9, j, l);
                    }
                }
            }
            
            // apply observation error
            for(l = 0; l < ndeathlads; l++) {
                // death incidence in LADs
                sigma2y = nages * (a1 + a2) + 2.0 * b * u_night_lad(l);
                sigma2y += sigma2_lad;
                muy = u_night_lad(l) + nages * (a1 - a2);
                u_night_lad(l) = rdtnorm_cpp(
                    muy, 
                    sqrt(sigma2y),
                    0,
                    std::numeric_limits<double>::infinity(),
                    engSerial
                );
            }
            for(l = 0; l < nregions; l++) {
                for(j = 0; j < nages; j++) {
                    // death incidence by age and region
                    sigma2y = nlads_region(l) * (a1 + a2) + 2.0 * b * u_night_age_region(j, l);
                    sigma2y += sigma2_age_region;
                    muy = u_night_age_region(j, l) + nlads_region(l) * (a1 - a2);
                    u_night_age_region(j, l) = rdtnorm_cpp(
                        muy, 
                        sqrt(sigma2y),
                        0,
                        std::numeric_limits<double>::infinity(),
                        engSerial
                    );
                }
            }
            for(l = 0; l < nnhsregions; l++) {
                // and hospital count by NHS region
                sigma2y = nages * nlads_nhsregion(l) * (a1 + a2) + 2.0 * b * u_night_nhsregion(l);
                sigma2y += sigma2_nhsregion;
                muy = u_night_nhsregion(l) + nages * nlads_nhsregion(l) * (a1 - a2);
                u_night_nhsregion(l) = rdtnorm_cpp(
                    muy, 
                    sqrt(sigma2y),
                    0,
                    std::numeric_limits<double>::infinity(),
                    engSerial
                );  
                for(j = 0; j < nnhsages; j++) {
                    // hospital incidence by age and NHS region
                    sigma2y = nlads_nhsregion(l) * (a1 + a2) + 2.0 * b * u_night_age_nhsregion(j, l);
                    sigma2y += sigma2_age_nhsregion;
                    muy = u_night_age_nhsregion(j, l) + nlads_nhsregion(l) * (a1 - a2);
                    u_night_age_nhsregion(j, l) = rdtnorm_cpp(
                        muy, 
                        sqrt(sigma2y),
                        0,
                        std::numeric_limits<double>::infinity(),
                        engSerial
                    );
                }
            }
            // cumulate incidence
            for(l = 0; l < ndeathlads; l++) {
                u_night_lad_cum[i](l) = u_night_lad(l);
            }
            for(l = 0; l < nregions; l++) {
                for(j = 0; j < nages; j++) {
                    u_night_age_region_cum[i](j, l) = u_night_age_region(j, l);
                }
            }
            for(l = 0; l < nnhsregions; l++) {
                for(j = 0; j < nnhsages; j++) {
                    u_night_age_nhsregion_cum[i](j, l) = u_night_age_nhsregion(j, l);
                }
            }
            if(saveAll == 2) {
                u_night_full.zeros();
                for(l = 0; l < u1_moves.n_rows; l++) {
                    for(j = 0; j < nages; j++) {
                        for(k = 0; k < nclasses; k++) {
                            u_night_full(k, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](k, j, l);
                        }
                    }
                }
                for(l = 0; l < nlads; l++) {
                    for(j = 0; j < nages; j++) {
                        for(k = 0; k < nclasses; k++) {
                            u_night_full(k, j, l) += u2[i](k, j, l);
                        }
                    }
                }
            }
            if(writeExt == 0) { 
                if(saveAll == 1) {
                    out[i] = List::create(Named("lads") = u_night_lad_cum[i], _["age_region"] = u_night_age_region_cum[i], _["nhsregion"] = u_night_nhsregion, _["age_nhsregion"] = u_night_age_nhsregion_cum[i]);
                } else {
                    out[i] = List::create(Named("full") = u_night_full, _["lads"] = u_night_lad_cum[i], _["age_region"] = u_night_age_region_cum[i], _["nhsregion"] = u_night_nhsregion, _["age_nhsregion"] = u_night_age_nhsregion_cum[i]);
                }
            } else {
                if(saveAll == 2) {
                    std::sprintf(file_name, "%s/p_%u.csv", std::string(outputName[0]).c_str(), i);
                    file.open(file_name);
                    file << "time, class, ";
                    for(j = 0; j < nages; j++) file << "age" << j + 1 << ", ";
                    file << "lad\n";
                    for(l = 0; l < nlads; l++) {
                        for(arma::uword r = 0; r < nclasses; r++) {
                            file << t << ", " << r << ", ";
                            for(j = 0; j < nages; j++) {
                                file << u_night_full(r, j, l) << ", ";
                            }
                            file << l + 1 << "\n";
                        }
                    }
                    file.close();
                }
                
                std::sprintf(file_name, "%s/p_lads_%u.csv", std::string(outputName[0]).c_str(), i);
                file.open(file_name);
                file << "time, deaths, lad\n";
                for(l = 0; l < ndeathlads; l++) {
                    file << t << ", " << u_night_lad_cum[i](l) << ", ";
                    file << l + 1 << "\n";
                }
                file.close();
                
                std::sprintf(file_name, "%s/p_age_region_%u.csv", std::string(outputName[0]).c_str(), i);
                file.open(file_name);
                file << "time, ";
                for(j = 0; j < nages; j++) file << "age" << j + 1 << ", ";
                file << "region\n";
                for(l = 0; l < nregions; l++) {
                    file << t << ", ";
                    for(j = 0; j < nages; j++) {
                        file << u_night_age_region_cum[i](j, l) << ", ";
                    }
                    file << l + 1 << "\n";
                }
                file.close();
                
                std::sprintf(file_name, "%s/p_nhsregion_%u.csv", std::string(outputName[0]).c_str(), i);
                file.open(file_name);
                file << "time, hosp, nhsregion\n";
                for(l = 0; l < nnhsregions; l++) {
                    file << t << ", ";
                    file << u_night_nhsregion(l) << ", ";
                    file << l + 1 << "\n";
                }
                file.close();
                
                std::sprintf(file_name, "%s/p_age_nhsregion_%u.csv", std::string(outputName[0]).c_str(), i);
                file.open(file_name);
                file << "time, ";
                for(j = 0; j < nnhsages; j++) file << "age" << j + 1 << ", ";
                file << "nhsregion\n";
                for(l = 0; l < nnhsregions; l++) {
                    file << t << ", ";
                    for(j = 0; j < nnhsages; j++) {
                        file << u_night_age_nhsregion_cum[i](j, l) << ", ";
                    }
                    file << l + 1 << "\n";
                }
                file.close();
            }
        }
    }
    
    // vectors for storing acceptance rates of MCMC
    arma::ivec nacc(npart); nacc.zeros();
    arma::mat MHweights_lad(npart, ndeathlads); MHweights_lad.zeros();
    
    // parameters for conditional sampling of simulator
    arma::vec condpars (pars.n_elem); condpars.zeros();
    condpars = pars;
    for(j = 0; j < nages; j++) {
        condpars(j + 4 * nages + 2) = pars(j + 4 * nages + 2) * (1.0 - pars(j + 5 * nages + 2) - pars(j + 6 * nages + 2));
        condpars(j + 4 * nages + 2) /= (1.0 - pars(j + 4 * nages + 2) * pars(j + 5 * nages + 2) - pars(j + 4 * nages + 2) * pars(j + 6 * nages + 2));
        condpars(j + 5 * nages + 2) = 0.0;
        condpars(j + 6 * nages + 2) = 0.0;
        
        condpars(j + 8 * nages + 2) = 0.0;
        condpars(j + 9 * nages + 2) = 0.0;
    }
    
    // initialise timer
    Timer timer;
    int timer_cnt = 0;
    double prev_time = 0.0;    
    
    // loop over time
    double ll = 0.0;
    for(t = 0; t < ndays; t++) {
            
        // check for interrupt
        R_CheckUserInterrupt();
        
        // extract data
        if(PF == 1) {
            obsInc_lad = deathInc_lad.row(t).t();
            obsInc_age_region = deathInc_age_region.row(t).t();
            obsInc_age_nhsregion = hospInc_age_nhsregion.row(t).t();
            obs_nhsregion = hosp_nhsregion.row(t).t();
        }
        
        // loop over particles
#ifdef _OPENMP
#pragma omp parallel for default(none) private(j, l, k) shared(seeds, npart, nages, nclasses, nlads, C, u1_moves, u1, u1_new, ncohorts1, u2, u2_new, playprobs, ncohorts2, t, pars, weights, a_dis, b_dis, a1, a2, b, obsInc_lad, obsInc_age_region, obsInc_age_nhsregion, obs_nhsregion, PF, ndays, ndeathlads, nregions, nnhsages, nnhsregions, lookup, age_lookup, sigma2_lad, sigma2_age_region, sigma2_nhsregion, sigma2_age_nhsregion, nlads_region, nlads_nhsregion, MHweights_lad)
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
            arma::mat ExpDHinc (nages, ndeathlads); ExpDHinc.zeros();
            arma::mat ExpDIinc (nages, ndeathlads); ExpDIinc.zeros();
            arma::vec yD (ndeathlads); yD.zeros();
            arma::vec muD (ndeathlads); muD.zeros();
            arma::mat muDc (nages, ndeathlads); muDc.zeros();
            arma::mat muDHc (nages, ndeathlads); muDHc.zeros();
            arma::mat muDIc (nages, ndeathlads); muDIc.zeros();
            arma::mat muDra (nages, nregions); muDra.zeros();
            arma::imat RHinc (nages, nlads); RHinc.zeros();
            arma::imat Hinc (nages, nlads); Hinc.zeros();
            arma::imat H (nages, nlads); H.zeros();
            arma::imat RIinc (nages, nlads); RIinc.zeros();
            arma::imat I2inc (nages, nlads); I2inc.zeros();
            arma::imat I1inc (nages, nlads); I1inc.zeros();
            arma::imat Pinc (nages, nlads); Pinc.zeros();
            arma::imat RAinc (nages, nlads); RAinc.zeros();
            arma::imat Ainc (nages, nlads); Ainc.zeros();
            arma::imat Einc (nages, nlads); Einc.zeros();
            arma::icube tempMD (nclasses, nages, nlads); tempMD.zeros();
            
            arma::icube u1_night(nclasses, nages, nlads); u1_night.zeros();
            
            // play movements
            for(l = 0; l < nlads; l++) {
                // loop over age and classes
                for(j = 0; j < nages; j++) {
                    for(k = 0; k < nclasses; k++) {
                        // set number in class
                        int tn = u2[i](k, j, l);
                        // distribute across connected lads
                        for(int s = ncohorts1(l); s < (ncohorts1(l) + ncohorts2(l)); s++) {
                            u1[i](k, j, s) = 0;
                        }
                        if(tn > 0) {
                            int s = ncohorts1(l);
                            int r;
                            while(s < (ncohorts1(l) + ncohorts2(l)) && tn > 0) {
                                r = rbinom_cpp(tn, playprobs(s), eng);
                                u1[i](k, j, s) = r;
                                tn -= r;
                                if(tn < 0) stop("Error in multinomial sampling of play movements\n");
                                s++;
                            }
                            if(tn != 0) stop("Non-zero tn %d %d %f\n", tn, r, playprobs(s - 1));
                        }
                    }
                }
            }
            
            // ensure player movements match between 
            // current and past movement cohorts
            u1_new[i] = u1[i];
            
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
            double muy, sigma2y, theta2y, nuy, muy1, sigma2y1;
            
            // run model and return u1
            discreteStochModel(i, nclasses, nages, nlads, pars, t - 1, t, u1_moves, u1_new, C, eng);
        
            // cols: c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")
            //          0,   1,   2,   3,    4,   5,    6,    7,    8,    9,   10,   11
            
            // adjust states according to model discrepancy
                            
            // D (MD on incidence)
            std::strcpy(str1, "Dinc");
            // aggregate incidence to LAD-level and calculate expectations
            for(j = 0; j < nages; j++) {                    
                for(l = 0; l < u1_moves.n_rows; l++) {
                    DHinc(j, u1_moves(l, 0) - 1) += (u1_new[i](11, j, l) - u1[i](11, j, l));
                    DIinc(j, u1_moves(l, 0) - 1) += (u1_new[i](6, j, l) - u1[i](6, j, l));
                }
            }
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    if(lookup(l, 1) >= 0) {
                        muy = (double) DHinc(j, l);
                        sigma2y = 2.0 * a_dis + 2.0 * b_dis * DHinc(j, l);
                        ExpDHinc(j, lookup(l, 1) - 1) = exptnorm_cpp(muy, sqrt(sigma2y), -0.5, u1_night(9, j, l) + 0.5);
                        muy = (double) DIinc(j, l);
                        sigma2y = 2.0 * a_dis + 2.0 * b_dis * DIinc(j, l);
                        ExpDIinc(j, lookup(l, 1) - 1) = exptnorm_cpp(muy, sqrt(sigma2y), -0.5, u1_night(5, j, l) + 0.5);
                    }
                }
            }
            // conditional sampling of MD terms
            for(k = 0; k < nlads; k++) {
                
                if(lookup(k, 1) >= 0) {
                    
                    // set deathlad
                    l = lookup(k, 1) - 1;
            
                    // sample y given DZ, DX
                    muy = 0.0;
                    for(j = 0; j < nages; j++) {
                        muy += ExpDHinc(j, l) + ExpDIinc(j, l);
                    }
                    theta2y = 2.0 * nages * (a1 + a2) + 2.0 * b * muy;
                    sigma2y = sigma2_lad + theta2y;
                    nuy = 2.0 * nages * (a1 - a2) + muy;
                    if(!arma::is_finite(muy) || !arma::is_finite(nuy) || !arma::is_finite(sigma2y)) {
                        Rprintf("muy = %f nuy = %f sigma2y = %f\n", muy, nuy, sigma2y);
                    }
                    yD(l) = rtnorm_one((obsInc_lad(l) - 0.5 - nuy) / sqrt(sigma2y), (obsInc_lad(l) + 0.5 - nuy) / sqrt(sigma2y), eng);
                    yD(l) = yD(l) * sqrt(sigma2y) + nuy;
                    if(yD(l) < (obsInc_lad(l) - 0.5) || yD(l) > (obsInc_lad(l) + 0.5)) {
                        stop("Error in y sampling\n");
                    }
                    
                    // calculate weights
                    if(PF == 1) {
                        MHweights_lad(i, l) = ltnorm_cpp(yD(l), nuy, sqrt(sigma2y), obsInc_lad(l) - 0.5, obsInc_lad(l) + 0.5);
                    }
                    
                    // sample mu given y
                    muy = yD(l) * theta2y + nuy * sigma2_lad;
                    muy /= (theta2y + sigma2_lad);
                    sigma2y = (sigma2_lad * theta2y) / (sigma2_lad + theta2y);
                    muD(l) = rtnorm_one(-10000000, 10000000, eng);
                    muD(l) = muD(l) * sqrt(sigma2y) + muy;
                    
                    // calculate weights
                    if(PF == 1) {
                        MHweights_lad(i, l) += R::dnorm(muD(l), muy, sqrt(sigma2y), 1);
                    }
                        
                    // now sample mean components from convolution
                    
                    // get sum of expectations
                    muy = nuy - 2.0 * nages * (a1 - a2);
                    
                    // DH means
                    int s = nages * 2;
                    double muD1 = muD(l);
                    for(j = 0; j < nages; j++) {
                        muy -= ExpDHinc(j, l);
                        sigma2y = (s - 1) * (a1 + a2) + 2.0 * b * muy;
                        nuy = (a1 + a2 + 2.0 * b * ExpDHinc(j, l)) * (muD1 - ExpDHinc(j, l) - a1 + a2);
                        nuy += sigma2y * (ExpDHinc(j, l) + a1 - a2);
                        nuy /= theta2y;
                        sigma2y *= (a1 + a2 + 2.0 * b * ExpDHinc(j, l));
                        sigma2y /= theta2y;
                        muDHc(j, l) = rtnorm_one(-10000000, 10000000, eng);
                        muDHc(j, l) = muDHc(j, l) * sqrt(sigma2y) + nuy;
                        // calculate weights
                        if(PF == 1) {
                            MHweights_lad(i, l) += R::dnorm(muDHc(j, l), nuy, sqrt(sigma2y), 1);
                        }
                        muD1 -= muDHc(j, l);
                        theta2y = (s - 1) * (a1 + a2) + 2.0 * b * muy;
                        s--;
                    }
                    for(j = 0; j < nages; j++) {
                        if(j < (nages - 1)) {
                            muy -= ExpDIinc(j, l);
                            sigma2y = (s - 1) * (a1 + a2) + 2.0 * b * muy;
                            nuy = (a1 + a2 + 2.0 * b * ExpDIinc(j, l)) * (muD1 - ExpDIinc(j, l) - a1 + a2);
                            nuy += sigma2y * (ExpDIinc(j, l) + a1 - a2);
                            nuy /= theta2y;
                            sigma2y *= (a1 + a2 + 2.0 * b * ExpDIinc(j, l));
                            sigma2y /= theta2y;
                            muDIc(j, l) = rtnorm_one(-10000000, 10000000, eng);
                            muDIc(j, l) = muDIc(j, l) * sqrt(sigma2y) + nuy;
                            // calculate weights
                            if(PF == 1) {
                                MHweights_lad(i, l) += R::dnorm(muDIc(j, l), nuy, sqrt(sigma2y), 1);
                            }
                            muD1 -= muDIc(j, l);
                            theta2y = (s - 1) * (a1 + a2) + 2.0 * b * muy;
                            s--;
                        } else {
                            muDIc(j, l) = muD1;
                        } 
                    }
                }
            }
            
            // now sample MD terms conditional on mu
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    if(lookup(l, 1) >= 0) {
                        // sample MD conditional on mu
                        arma::vec tempdensy (u1_night(9, j, l) + 1); 
                        for(int s = 0; s <= u1_night(9, j, l); s++) {
                            
                            // MD density
                            muy = (double) DHinc(j, l);
                            sigma2y = 2.0 * a_dis + 2.0 * b_dis * DHinc(j, l);
                            tempdensy(s) += ldtnorm_cpp(
                                s, 
                                muy, 
                                sqrt(sigma2y), 
                                0.0, 
                                u1_night(9, j, l)
                            );
                            
                            // sampled mean given MD
                            muy = s + a1 - a2;
                            sigma2y = a1 + a2 + 2.0 * b * s;
                            tempdensy(s) += R::dnorm(muDHc(j, lookup(l, 1) - 1), muy, sqrt(sigma2y), 1);
                        }
                        double tempnorm = log_sum_exp(tempdensy, 0);
                        tempdensy = tempdensy - tempnorm;
                        tempdensy = exp(tempdensy);
                        tempdensy = tempdensy / sum(tempdensy);
                        int r = rmultinom_cpp(tempdensy, eng);
                        tempMD(11, j, l) = r - DHinc(j, l);
                        DHinc(j, l) = r;
                        if(DHinc(j, l) < 0 || DHinc(j, l) > u1_night(9, j, l)) stop("DHinc error\n");
                        
                        // adjust weights
                        if(PF == 1) {
                            MHweights_lad(i, lookup(l, 1) - 1) += log(tempdensy(r));
                        }
                        
                        // sample MD conditional on mu
                        arma::vec tempdensy1 (u1_night(5, j, l) + 1); 
                        for(int s = 0; s <= u1_night(5, j, l); s++) {
                            
                            // MD density
                            muy = (double) DIinc(j, l);
                            sigma2y = 2.0 * a_dis + 2.0 * b_dis * DIinc(j, l);
                            tempdensy1(s) += ldtnorm_cpp(
                                s, 
                                muy, 
                                sqrt(sigma2y), 
                                0.0, 
                                u1_night(5, j, l)
                            );
                            
                            // sampled mean given MD
                            muy = s + a1 - a2;
                            sigma2y = a1 + a2 + 2.0 * b * s;
                            tempdensy1(s) += R::dnorm(muDIc(j, lookup(l, 1) - 1), muy, sqrt(sigma2y), 1);
                        }
                        tempnorm = log_sum_exp(tempdensy1, 0);
                        tempdensy1 = tempdensy1 - tempnorm;
                        tempdensy1 = exp(tempdensy1);
                        tempdensy1 = tempdensy1 / sum(tempdensy1);
                        r = rmultinom_cpp(tempdensy1, eng);
                        tempMD(6, j, l) = r - DIinc(j, l);
                        DIinc(j, l) = r;
                        if(DIinc(j, l) < 0 || DIinc(j, l) > u1_night(5, j, l)) stop("DIinc error\n");
                                                
                        // adjust weights
                        if(PF == 1) {
                            MHweights_lad(i, lookup(l, 1) - 1) += log(tempdensy1(r));
                        }
                        
                        // numerator of weights
                        if(PF == 1) {  
                            // now calculate numerator of weights
                            muy = (double) DHinc(j, l) - tempMD(11, j, l);
                            sigma2y = 2.0 * a_dis + 2.0 * b_dis * (DHinc(j, l) - tempMD(11, j, l));
                            weights(i) += ldtnorm_cpp(
                                DHinc(j, l), 
                                muy, 
                                sqrt(sigma2y),
                                0,
                                u1_night(9, j, l)
                            );
                            muy = a1 - a2 + DHinc(j, l);
                            sigma2y = a1 + a2 + 2.0 * b * DHinc(j, l);
                            weights(i) += R::dnorm(
                                muDHc(j, lookup(l, 1) - 1), 
                                muy, 
                                sqrt(sigma2y), 
                                1
                            );
                            muy = (double) DIinc(j, l) - tempMD(6, j, l);
                            sigma2y = 2.0 * a_dis + 2.0 * b_dis * (DIinc(j, l) - tempMD(6, j, l));
                            weights(i) += ldtnorm_cpp(
                                DIinc(j, l), 
                                muy, 
                                sqrt(sigma2y),
                                0,
                                u1_night(5, j, l)
                            );
                            muy = a1 - a2 + DIinc(j, l);
                            sigma2y = a1 + a2 + 2.0 * b * DIinc(j, l);
                            weights(i) += R::dnorm(
                                muDIc(j, lookup(l, 1) - 1), 
                                muy, 
                                sqrt(sigma2y), 
                                1
                            );
                        }
                    } else {
                        // sample MD conditional on simulator
                        sigma2y = 2.0 * a_dis + 2.0 * b_dis * DHinc(j, l);
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
                
                        // sample MD conditional on simulator
                        sigma2y = 2.0 * a_dis + 2.0 * b_dis * DIinc(j, l);
                        muy = (double) DIinc(j, l);
                        s = rdtnorm_cpp(
                            muy, 
                            sqrt(sigma2y),
                            0.0,
                            u1_night(5, j, l),
                            eng
                        );
                        tempMD(6, j, l) = s - DIinc(j, l);
                        DIinc(j, l) = s;
                        if(DIinc(j, l) < 0 || DIinc(j, l) > u1_night(5, j, l)) stop("DIinc error\n");
                    }
                }
            }
            
            // now sample remaining latent variables
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    if(lookup(l, 2) >= 0) {
                        muDra(j, lookup(l, 2) - 1) += muDHc(j, lookup(l, 1) - 1);
                        muDra(j, lookup(l, 2) - 1) += muDIc(j, lookup(l, 1) - 1);
                    }
                }
            }
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nregions; l++) {
                
                    // sample y given DZ, mu
                    double yDra = rtnorm_one((obsInc_age_region(j * nregions + l) - 0.5 - muDra(j, l)) / sqrt(sigma2_age_region), (obsInc_age_region(j * nregions + l) + 0.5 - muDra(j, l)) / sqrt(sigma2_age_region), eng);
                    yDra = yDra * sqrt(sigma2_age_region) + muDra(j, l);
                    if(yDra < (obsInc_age_region(j * nregions + l) - 0.5) || yDra > (obsInc_age_region(j * nregions + l) + 0.5)) {
                        stop("Error in yDra sampling\n");
                    }
                    
                    // adjust weights
                    if(PF == 1) {
                        weights(i) -= ltnorm_cpp(yDra, muDra(j, l), sqrt(sigma2_age_region), obsInc_age_region(j * nregions + l) - 0.5, obsInc_age_region(j * nregions + l) + 0.5);
                        weights(i) += R::dnorm(yDra, muDra(j, l), sqrt(sigma2_age_region), 1);
                    }
                }
            }
                    
            if(PF == 1) {
                // now adjust weights for remaining observation processes
                for(l = 0; l < ndeathlads; l++) {
                                
                    // denominator of weights
                    weights(i) -= MHweights_lad(i, l);
                    
                    // now adjust weights for remaining observation processes
                    weights(i) += R::dnorm(
                        yD(l),
                        muD(l), 
                        sqrt(sigma2_lad),
                        1
                    );
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
                    H(j, u1_moves(l, 0) - 1) += u1_new[i](9, j, l);
                }
            }
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    sigma2y = 2.0 * a_dis + 2.0 * b_dis * H(j, l);
                    tempMD(9, j, l) = rdtnorm_cpp(
                        0.0, 
                        sqrt(sigma2y),
                        -H(j, l) + u1_night(9, j, l) - DHinc(j, l) - RHinc(j, l),
                        u1_night(5, j, l) - DIinc(j, l) - H(j, l) + u1_night(9, j, l) - DHinc(j, l) - RHinc(j, l),
                        eng
                    );
                    H(j, l) += tempMD(9, j, l);
                    if(H(j, l) < 0) stop("H error\n");
                    // calculate incidence
                    Hinc(j, l) = H(j, l) - u1_night(9, j, l) + RHinc(j, l) + DHinc(j, l);
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
            redistribution(i, nages, nlads, tempMD, ncohorts1, u1, u1_new, eng, 0);
            
            // aggregate play movements back to LAD level
            u2_new[i].zeros();
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    for(k = 0; k < nclasses; k++) {
                        for(int r = ncohorts1(l); r < (ncohorts1(l) + ncohorts2(l)); r++) {
                            u2_new[i](k, j, l) += u1_new[i](k, j, r);
                        }
                    }
                }
            }
            
            if(PF == 1) {
                // set up observation errors
                arma::ivec u_night_lad1(ndeathlads); u_night_lad1.zeros();
                arma::imat u_night_age_region1(nages, nregions); u_night_age_region1.zeros();
                arma::imat u_night_age_nhsregion1(nnhsages, nnhsregions); u_night_age_nhsregion1.zeros();
                arma::ivec u_night_nhsregion1(nnhsregions); u_night_nhsregion1.zeros();
                
                for(l = 0; l < nlads; l++) {
                    for(j = 0; j < nages; j++) {
                        
                        // death incidence in LADs
                        if(lookup(l, 1) >= 0) {
                            u_night_lad1(lookup(l, 1) - 1) += DIinc(j, l);
                            u_night_lad1(lookup(l, 1) - 1) += DHinc(j, l);
                        }
                        
                        // death incidence by age and region
                        if(lookup(l, 2) >= 0) {
                            u_night_age_region1(j, lookup(l, 2) - 1) += DIinc(j, l);
                            u_night_age_region1(j, lookup(l, 2) - 1) += DHinc(j, l);
                        }
                        
                        if(lookup(l, 3) >= 0) {
                            // hospital incidence by age and NHS region
                            u_night_age_nhsregion1(age_lookup(j, 1) - 1, lookup(l, 3) - 1) += Hinc(j, l);  
                            // and hospital count by NHS region
                            u_night_nhsregion1(lookup(l, 3) - 1) += H(j, l);
                        }
                    }
                }
                
                // calculate observation error for remaining terms
                for(l = 0; l < nnhsregions; l++) {
                    // and hospital count by NHS region
                    sigma2y = nages * nlads_nhsregion(l) * (a1 + a2) + 2.0 * b * u_night_nhsregion1(l);
                    sigma2y += sigma2_nhsregion;
                    muy = u_night_nhsregion1(l) + nages * nlads_nhsregion(l) * (a1 - a2);
                    weights(i) += ldtnorm_cpp(
                        obs_nhsregion(l),
                        muy, 
                        sqrt(sigma2y),
                        0,
                        std::numeric_limits<double>::infinity()
                    );
                    for(j = 0; j < nnhsages; j++) {
                        // hospital incidence by age and NHS region
                        sigma2y = nlads_nhsregion(l) * (a1 + a2) + 2.0 * b * u_night_age_nhsregion1(j, l);
                        sigma2y += sigma2_age_nhsregion;
                        muy = u_night_age_nhsregion1(j, l) + nlads_nhsregion(l) * (a1 - a2);
                        weights(i) += ldtnorm_cpp(
                            obsInc_age_nhsregion(j * nnhsregions + l),
                            muy, 
                            sqrt(sigma2y),
                            0,
                            std::numeric_limits<double>::infinity()
                        );
                    }
                }
            }
            
            // advance seed
            seeds((arma::uword) omp_get_thread_num()) = eng();
        }
        
        if(PF == 1) {
            // calculate log-likelihood contribution
            ll += log_sum_exp(weights, 1);
            
            // if zero likelihood then return
            if(!arma::is_finite(ll)) {
                Rprintf("Finishing early due to non-finite log-likelihood estimate\n");
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
            for(i = 0; i < npart; i++) inds(i) = (arma::uword) rmultinom_cpp(weights, engSerial);
            
            for(i = 0; i < npart; i++) u11[i] = u1[inds(i)];
            for(i = 0; i < npart; i++) u1[i] = u11[i];
            
            for(i = 0; i < npart; i++) u11[i] = u1_new[inds(i)];
            for(i = 0; i < npart; i++) u1_new[i] = u11[i];
            
            for(i = 0; i < npart; i++) u11[i] = u2[inds(i)];
            for(i = 0; i < npart; i++) u2[i] = u11[i];
            
            for(i = 0; i < npart; i++) u11[i] = u2_new[inds(i)];
            for(i = 0; i < npart; i++) u2_new[i] = u11[i];
            
            for(i = 0; i < npart; i++) u_night_vec_cum1[i] = u_night_lad_cum[inds(i)];
            for(i = 0; i < npart; i++) u_night_lad_cum[i] = u_night_vec_cum1[i];
            
            for(i = 0; i < npart; i++) u_night_mat_cum1[i] = u_night_age_region_cum[inds(i)];
            for(i = 0; i < npart; i++) u_night_age_region_cum[i] = u_night_mat_cum1[i];
            
            for(i = 0; i < npart; i++) u_night_mat_cum1[i] = u_night_age_nhsregion_cum[inds(i)];
            for(i = 0; i < npart; i++) u_night_age_nhsregion_cum[i] = u_night_mat_cum1[i];
            
            // Metropolis-Hastings steps to deal with particle impoverishment
            if(niter > 0) {
#ifdef _OPENMP
#pragma omp parallel for default(none) private(j, l, k) shared(seeds, npart, nages, nclasses, nlads, C, ncohorts1, u1_moves, u1, u1_new, ncohorts2, u2, u2_new, playprobs, t, pars, a_dis, b_dis, a1, a2, b, obsInc_lad, obsInc_age_region, obsInc_age_nhsregion, obs_nhsregion, PF, ndays, ndeathlads, nregions, nnhsages, nnhsregions, lookup, age_lookup, sigma2_lad, sigma2_age_region, sigma2_nhsregion, sigma2_age_nhsregion, nlads_region, nlads_nhsregion, condpars, niter, nacc, MHweights_lad)
#endif
                for(i = 0; i < npart; i++) {
            
                    // set up print string for debugging
                    char str1[80];
	                std::strcpy(str1, "MH setup");
         
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
                    arma::imat RHinc1 (nages, nlads); RHinc1.zeros();
                    arma::imat Hinc (nages, nlads); Hinc.zeros();
                    arma::imat Hinc1 (nages, nlads); Hinc1.zeros();
                    arma::imat H (nages, nlads); H.zeros();
                    arma::imat RIinc (nages, nlads); RIinc.zeros();
                    arma::imat I2inc (nages, nlads); I2inc.zeros();
                    arma::imat I1inc (nages, nlads); I1inc.zeros();
                    arma::imat Pinc (nages, nlads); Pinc.zeros();
                    arma::imat RAinc (nages, nlads); RAinc.zeros();
                    arma::imat Ainc (nages, nlads); Ainc.zeros();
                    arma::imat Einc (nages, nlads); Einc.zeros();
                    arma::icube tempMD (nclasses, nages, nlads); tempMD.zeros();
                    arma::icube tempMD1 (4, nages, nlads); tempMD1.zeros();
                    arma::icube tempsim1 (4, nages, nlads); tempsim1.zeros();
                    
                    arma::icube u1_night(nclasses, nages, nlads); u1_night.zeros();
                    
                    // PLAY MOVEMENTS HAVE ALREADY BEEN DONE AND DON'T NEED TO BE RESAMPLED HERE
                    
                    // aggregate counts to LAD-level
                    for(j = 0; j < nages; j++) {                    
                        for(l = 0; l < u1_moves.n_rows; l++) {
                            for(int s = 0; s < nclasses; s++) {
                                u1_night(s, j, u1_moves(l, 0) - 1) += u1[i](s, j, l);
                            }
                        }
                    }
                    
                    // set auxiliary variables
                    double muy, sigma2y, pI1pI1D, pHpHD, pI1pI1H, pHpHR, acc, acccurr, accprop, u;
                
                    // cols: c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")
                    //          0,   1,   2,   3,    4,   5,    6,    7,    8,    9,   10,   11
                    
                    // adjust states according to model discrepancy
                    
                    // current likelihood
                    acccurr = 0.0;
                    
                    // set up observation errors
                    arma::ivec u_night_lad1(ndeathlads); u_night_lad1.zeros();
                    arma::imat u_night_age_region1(nages, nregions); u_night_age_region1.zeros();
                    arma::imat u_night_age_nhsregion1(nnhsages, nnhsregions); u_night_age_nhsregion1.zeros();
                    arma::ivec u_night_nhsregion1(nnhsregions); u_night_nhsregion1.zeros();
                    
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        k = (arma::uword) u1_moves(l, 0) - 1;
                        for(j = 0; j < nages; j++) {
                            
                            // death incidence in LADs
                            if(lookup(k, 1) >= 0) {
                                u_night_lad1(lookup(k, 1) - 1) += (u1_new[i](6, j, l) - u1[i](6, j, l));
                                u_night_lad1(lookup(k, 1) - 1) += (u1_new[i](11, j, l) - u1[i](11, j, l));
                            }
                            
                            // death incidence by age and region
                            if(lookup(k, 2) >= 0) {
                                u_night_age_region1(j, lookup(k, 2) - 1) += (u1_new[i](6, j, l) - u1[i](6, j, l));
                                u_night_age_region1(j, lookup(k, 2) - 1) += (u1_new[i](11, j, l) - u1[i](11, j, l));
                            }
                            
                            if(lookup(k, 3) >= 0) {
                                // hospital incidence by age and NHS region
                                u_night_age_nhsregion1(age_lookup(j, 1) - 1, lookup(k, 3) - 1) += (u1_new[i](9, j, l) - u1[i](9, j, l)) + (u1_new[i](10, j, l) - u1[i](10, j, l)) + (u1_new[i](11, j, l) - u1[i](11, j, l));
                                // and hospital count by NHS region
                                u_night_nhsregion1(lookup(k, 3) - 1) += u1[i](9, j, l);
                            }
                        }
                    }
                    
                    // calculate observation error
                    for(j = 0; l < MHweights_lad.n_cols; j++) {
                        // death incidence in LADs
                        acccurr += MHweights_lad(i, j);
                    }
                    for(l = 0; l < nregions; l++) {
                        for(j = 0; j < nages; j++) {
                            // death incidence by age and region
                            sigma2y = nlads_region(l) * (a1 + a2) + 2.0 * b * u_night_age_region1(j, l);
                            sigma2y += sigma2_age_region;
                            muy = u_night_age_region1(j, l) + nlads_region(l) * (a1 - a2);
                            acccurr += ldtnorm_cpp(
                                obsInc_age_region(j * nregions + l),
                                muy, 
                                sqrt(sigma2y),
                                0,
                                std::numeric_limits<double>::infinity()
                            );
                        }
                    }
                    for(l = 0; l < nnhsregions; l++) {
                        // and hospital count by NHS region
                        sigma2y = nages * nlads_nhsregion(l) * (a1 + a2) + 2.0 * b * u_night_nhsregion1(l);
                        sigma2y += sigma2_nhsregion;
                        muy = u_night_nhsregion1(l) + nages * nlads_nhsregion(l) * (a1 - a2);
                        acccurr += ldtnorm_cpp(
                            obs_nhsregion(l),
                            muy, 
                            sqrt(sigma2y),
                            0,
                            std::numeric_limits<double>::infinity()
                        );
                        for(j = 0; j < nnhsages; j++) {
                            // hospital incidence by age and NHS region
                            sigma2y = nlads_nhsregion(l) * (a1 + a2) + 2.0 * b * u_night_age_nhsregion1(j, l);
                            sigma2y += sigma2_age_nhsregion;
                            muy = u_night_age_nhsregion1(j, l) + nlads_nhsregion(l) * (a1 - a2);
                            acccurr += ldtnorm_cpp(
                                obsInc_age_nhsregion(j * nnhsregions + l),
                                muy, 
                                sqrt(sigma2y),
                                0,
                                std::numeric_limits<double>::infinity()
                            );
                        }
                    }
                                    
                    // run over Metropolis-Hastings steps
                    nacc(i) = 0;
                    for(int it = 0; it < niter; it++) {
                        
                        // loop over ages
                        std::strcpy(str1, "MH");
                        for(j = 0; j < nages; j++) {
                        
                            // extract transition probabilities
                            pHpHD = pars(j + 8 * nages + 2) * pars(j + 9 * nages + 2);
                            pHpHR = pars(j + 8 * nages + 2) * (1.0 - pars(j + 9 * nages + 2));
                            pI1pI1H = pars(j + 4 * nages + 2) * pars(j + 5 * nages + 2);
                            pI1pI1D = pars(j + 4 * nages + 2) * pars(j + 6 * nages + 2);
                            
                            // loop over LADs
                            for(l = 0; l < nlads; l++) {
                            
                                // sample DH from simulator
                                int r = rbinom_cpp(u1_night(9, j, l), pHpHD, eng);
                                DHinc1(j, l) = r;
                            
                                // sample MD conditional on simulator
                                sigma2y = 2.0 * a_dis + 2.0 * b_dis * DHinc1(j, l);
                                muy = (double) DHinc1(j, l);
                                int s = rdtnorm_cpp(
                                    muy, 
                                    sqrt(sigma2y),
                                    0.0,
                                    u1_night(9, j, l),
                                    eng
                                );
                                DHinc(j, l) = s;
                                
                                // sample RH from simulator given DH
                                r = rbinom_cpp(u1_night(9, j, l) - DHinc1(j, l), pHpHR / (1.0 - pHpHD), eng);
                                RHinc1(j, l) = r;
                            
                                // sample MD conditional on simulator
                                sigma2y = 2.0 * a_dis + 2.0 * b_dis * RHinc1(j, l);
                                muy = (double) RHinc1(j, l);
                                s = rdtnorm_cpp(
                                    muy, 
                                    sqrt(sigma2y),
                                    0.0,
                                    u1_night(9, j, l) - DHinc(j, l),
                                    eng
                                );
                                RHinc(j, l) = s;
                            
                                // sample DI from simulator
                                r = rbinom_cpp(u1_night(5, j, l), pI1pI1D, eng);
                                DIinc1(j, l) = r;
                            
                                // sample MD conditional on simulator
                                sigma2y = 2.0 * a_dis + 2.0 * b_dis * DIinc1(j, l);
                                muy = (double) DIinc1(j, l);
                                s = rdtnorm_cpp(
                                    muy, 
                                    sqrt(sigma2y),
                                    0.0,
                                    u1_night(5, j, l),
                                    eng
                                );
                                DIinc(j, l) = s;
                                
                                // sample H from conditional simulator
                                r = rbinom_cpp(u1_night(5, j, l) - DIinc1(j, l), pI1pI1H / (1.0 - pI1pI1D), eng);
                                Hinc1(j, l) = r;
                            
                                // sample MD conditional on simulator
                                H(j, l) = u1_night(9, j, l) + Hinc1(j, l) - DHinc1(j, l) - RHinc1(j, l);
                                sigma2y = 2.0 * a_dis + 2.0 * b_dis * H(j, l);
                                s = rdtnorm_cpp(
                                    0.0, 
                                    sqrt(sigma2y),
                                    -H(j, l) + u1_night(9, j, l) - DHinc(j, l) - RHinc(j, l),
                                    u1_night(5, j, l) - DIinc(j, l) - H(j, l) + u1_night(9, j, l) - DHinc(j, l) - RHinc(j, l),
                                    eng
                                );
                                H(j, l) += s;
                                Hinc(j, l) = H(j, l) - u1_night(9, j, l) + DHinc(j, l) + RHinc(j, l);
                            }
                        }
                    
                        // set acceptance probability
                        accprop = 0.0;
                        
                        // set up observation errors
                        u_night_lad1.zeros();
                        u_night_age_region1.zeros();
                        u_night_age_nhsregion1.zeros();
                        u_night_nhsregion1.zeros();
                        for(l = 0; l < nlads; l++) {
                            for(j = 0; j < nages; j++) {
                                // death incidence in LADs
                                if(lookup(l, 1) >= 0) {
                                    u_night_lad1(lookup(l, 1) - 1) += DIinc(j, l);
                                    u_night_lad1(lookup(l, 1) - 1) += DHinc(j, l);
                                }
                                // death incidence by age and region
                                if(lookup(l, 2) >= 0) {
                                    u_night_age_region1(j, lookup(l, 2) - 1) += DIinc(j, l);
                                    u_night_age_region1(j, lookup(l, 2) - 1) += DHinc(j, l);
                                }
                                if(lookup(l, 3) >= 0) {
                                    // hospital incidence by age and NHS region
                                    u_night_age_nhsregion1(age_lookup(j, 1) - 1, lookup(l, 3) - 1) += Hinc(j, l);
                                    // and hospital count by NHS region
                                    u_night_nhsregion1(lookup(l, 3) - 1) += H(j, l);
                                }
                            }
                        }
                    
                        // calculate observation error
                        for(l = 0; l < ndeathlads; l++) {
                            // death incidence in LADs
                            sigma2y = nages * (a1 + a2) + 2.0 * b * u_night_lad1(l);
                            sigma2y += sigma2_lad;
                            muy = u_night_lad1(l) + nages * (a1 - a2);
                            accprop += ldtnorm_cpp(
                                obsInc_lad(l),
                                muy, 
                                sqrt(sigma2y),
                                0,
                                std::numeric_limits<double>::infinity()
                            );
                        }
                        for(l = 0; l < nregions; l++) {
                            for(j = 0; j < nages; j++) {
                                // death incidence by age and region
                                sigma2y = nlads_region(l) * (a1 + a2) + 2.0 * b * u_night_age_region1(j, l);
                                sigma2y += sigma2_age_region;
                                muy = u_night_age_region1(j, l) + nlads_region(l) * (a1 - a2);
                                accprop += ldtnorm_cpp(
                                    obsInc_age_region(j * nregions + l),
                                    muy, 
                                    sqrt(sigma2y),
                                    0,
                                    std::numeric_limits<double>::infinity()
                                );
                            }
                        }
                        for(l = 0; l < nnhsregions; l++) {
                            // and hospital count by NHS region
                            sigma2y = nages * nlads_nhsregion(l) * (a1 + a2) + 2.0 * b * u_night_nhsregion1(l);
                            sigma2y += sigma2_nhsregion;
                            muy = u_night_nhsregion1(l) + nages * nlads_nhsregion(l) * (a1 - a2);
                            accprop += ldtnorm_cpp(
                                obs_nhsregion(l),
                                muy, 
                                sqrt(sigma2y),
                                0,
                                std::numeric_limits<double>::infinity()
                            );
                            for(j = 0; j < nnhsages; j++) {
                                // hospital incidence by age and NHS region
                                sigma2y = nlads_nhsregion(l) * (a1 + a2) + 2.0 * b * u_night_age_nhsregion1(j, l);
                                sigma2y += sigma2_age_nhsregion;
                                muy = u_night_nhsregion1(l) + nlads_region(l) * (a1 - a2);
                                accprop += ldtnorm_cpp(
                                    obsInc_age_nhsregion(j * nnhsregions + l),
                                    muy, 
                                    sqrt(sigma2y),
                                    0,
                                    std::numeric_limits<double>::infinity()
                                );
                            }
                        }
                                
                        // accept-reject
                        acc = accprop - acccurr;
                        u = log(eng()) - log(sitmo::prng::max());
                        if(u < acc) {
                            nacc(i)++;
                            acccurr = accprop;
                            for(j = 0; j < nages; j++) {
                                for(l = 0; l < nlads; l++) {
                                    tempsim1(0, j, l) = DIinc1(j, l);
                                    tempsim1(1, j, l) = Hinc1(j, l);
                                    tempsim1(2, j, l) = RHinc1(j, l);
                                    tempsim1(3, j, l) = DHinc1(j, l);
                                    tempMD1(0, j, l) = DIinc(j, l);
                                    tempMD1(1, j, l) = Hinc(j, l);
                                    tempMD1(2, j, l) = RHinc(j, l);
                                    tempMD1(3, j, l) = DHinc(j, l);
                                }
                            }
                        }
                    }
                    
                    // if some moves have been made then update counts
                    if(nacc(i) > 0) {
                    
                        // set new particle
                        u1_new[i] = u1[i];
                        
                        // set up simulator adjustments
                        tempMD.zeros();
                        for(l = 0; l < nlads; l++) {
                            for(j = 0; j < nages; j++) {
                                tempMD(6, j, l) = tempsim1(0, j, l);
                                tempMD(9, j, l) = tempsim1(1, j, l);
                                tempMD(10, j, l) = tempsim1(2, j, l);
                                tempMD(11, j, l) = tempsim1(3, j, l);
                            }
                        }
                
                        // now redistribute simulator incidence across cohorts
                        redistribution(i, nages, nlads, tempMD, ncohorts1, u1, u1_new, eng, 1);
                        
                        // sample remaining x values from conditional simulator
                        discreteStochModel(i, nclasses, nages, nlads, condpars, t - 1, t, u1_moves, u1_new, C, eng);
                        
                        // set model discrepancy counts for later re-distribution
                        tempMD.zeros();
                        for(l = 0; l < nlads; l++) {
                            for(j = 0; j < nages; j++) {
                                tempMD(6, j, l) = tempMD1(0, j, l) - tempsim1(0, j, l);
                                tempMD(9, j, l) = tempMD1(1, j, l) - tempMD1(2, j, l) - tempMD1(3, j, l);
                                tempMD(9, j, l) -= (tempsim1(1, j, l) - tempsim1(2, j, l) - tempsim1(3, j, l));
                                tempMD(10, j, l) = tempMD1(2, j, l) - tempsim1(2, j, l);
                                tempMD(11, j, l) = tempMD1(3, j, l) - tempsim1(3, j, l);
                                DIinc(j, l) = tempMD1(0, j, l);
                                Hinc(j, l) = tempMD1(1, j, l);
                                RHinc(j, l) = tempMD1(2, j, l);
                                DHinc(j, l) = tempMD1(3, j, l);
                            }
                        }
                        
                        // calculate remaining MD terms
                            
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
                        redistribution(i, nages, nlads, tempMD, ncohorts1, u1, u1_new, eng, 0);
                        
                        // aggregate play movements back to LAD level
                        u2_new[i].zeros();
                        for(j = 0; j < nages; j++) {
                            for(l = 0; l < nlads; l++) {
                                for(k = 0; k < nclasses; k++) {
                                    for(int r = ncohorts1(l); r < (ncohorts1(l) + ncohorts2(l)); r++) {
                                        u2_new[i](k, j, l) += u1_new[i](k, j, r);
                                    }
                                }
                            }
                        }           
                    }
                    
                    // advance seed
                    seeds((arma::uword) omp_get_thread_num()) = eng();
                }
            }
        }
        
        if(saveAll != 0) {
            // set up print string for debugging
            char str1[80];
            std::strcpy(str1, "save");
            double muy, sigma2y, u;
            for(i = 0; i < npart; i++) {
            
                // U1() BELOW ALREADY HAS THE CORRECT PLAYER MOVEMENTS SO CAN USE THE EXISTING
                // CODE WITHOUT NEEDING TO CROSS-REFERENCE U2
                
                // extract just counts for DI and DH
                u_night_lad.zeros();
                u_night_age_region.zeros();
                u_night_age_nhsregion.zeros();
                u_night_nhsregion.zeros();
                for(l = 0; l < u1_moves.n_rows; l++) {
                    k = (arma::uword) u1_moves(l, 0) - 1;
                    for(j = 0; j < nages; j++) {
                        
                        // death incidence in LADs
                        if(lookup(k, 1) >= 0) {
                            u_night_lad(lookup(k, 1) - 1) += (u1_new[i](6, j, l) - u1[i](6, j, l));
                            u_night_lad(lookup(k, 1) - 1) += (u1_new[i](11, j, l) - u1[i](11, j, l));
                        }
                        
                        // death incidence by age and region
                        if(lookup(k, 2) >= 0) {
                            u_night_age_region(j, lookup(k, 2) - 1) += (u1_new[i](6, j, l) - u1[i](6, j, l));
                            u_night_age_region(j, lookup(k, 2) - 1) += (u1_new[i](11, j, l) - u1[i](11, j, l));
                        }
                        
                        if(lookup(k, 3) >= 0) {
                            // hospital incidence by age and NHS region
                            u_night_age_nhsregion(age_lookup(j, 1) - 1, lookup(k, 3) - 1) += (u1_new[i](9, j, l) - u1[i](9, j, l)) + (u1_new[i](10, j, l) - u1[i](10, j, l)) + (u1_new[i](11, j, l) - u1[i](11, j, l));                    
                            // and hospital count by NHS region
                            u_night_nhsregion(lookup(k, 3) - 1) += u1_new[i](9, j, l);
                        }
                    }
                }
                
                nned to simulate this better
                
                // apply observation error
                for(l = 0; l < ndeathlads; l++) {
                    // death incidence in LADs
                    sigma2y = 2.0 * nages * (a1 + a2) + 2.0 * b * u_night_lad(l);
                    sigma2y += sigma2_lad;
                    muy = u_night_lad(l) + 2.0 * nages * (a1 - a2);
                    u_night_lad(l) = rdtnorm_cpp(
                        muy, 
                        sqrt(sigma2y),
                        0,
                        std::numeric_limits<double>::infinity(),
                        engSerial
                    );
                }
                for(l = 0; l < nregions; l++) {
                    for(j = 0; j < nages; j++) {
                        // death incidence by age and region
                        sigma2y = nlads_region(l) * (a1 + a2) + 2.0 * b * u_night_age_region(j, l);
                        sigma2y += sigma2_age_region;
                        muy = u_night_age_region(j, l) + nlads_region(l) * (a1 - a2);
                        u_night_age_region(j, l) = rdtnorm_cpp(
                            muy, 
                            sqrt(sigma2y),
                            0,
                            std::numeric_limits<double>::infinity(),
                            engSerial
                        );
                    }
                }
                for(l = 0; l < nnhsregions; l++) {
                    // and hospital count by NHS region
                    sigma2y = nages * nlads_nhsregion(l) * (a1 + a2) + 2.0 * b * u_night_nhsregion(l);
                    sigma2y += sigma2_nhsregion;
                    muy = u_night_nhsregion(l) + nages * nlads_nhsregion(l) * (a1 - a2);
                    u_night_nhsregion(l) = rdtnorm_cpp(
                        muy, 
                        sqrt(sigma2y),
                        0,
                        std::numeric_limits<double>::infinity(),
                        engSerial
                    );  
                    for(j = 0; j < nnhsages; j++) {
                        // hospital incidence by age and NHS region
                        sigma2y = nlads_nhsregion(l) * (a1 + a2) + 2.0 * b * u_night_age_nhsregion(j, l);
                        sigma2y += sigma2_age_nhsregion;
                        muy = u_night_age_nhsregion(j, l) + nlads_nhsregion(l) * (a1 - a2);
                        u_night_age_nhsregion(j, l) = rdtnorm_cpp(
                            muy, 
                            sqrt(sigma2y),
                            0,
                            std::numeric_limits<double>::infinity(),
                            engSerial
                        );
                    }
                }
                // cumulate incidence
                for(l = 0; l < ndeathlads; l++) {
                    u_night_lad_cum[i](l) += u_night_lad(l);
                }
                for(l = 0; l < nregions; l++) {
                    for(j = 0; j < nages; j++) {
                        u_night_age_region_cum[i](j, l) += u_night_age_region(j, l);
                    }
                }
                for(l = 0; l < nnhsregions; l++) {
                    for(j = 0; j < nnhsages; j++) {
                        u_night_age_nhsregion_cum[i](j, l) += u_night_age_nhsregion(j, l);
                    }
                }
                if(saveAll == 2) {
                    u_night_full.zeros();
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        for(j = 0; j < nages; j++) {
                            for(k = 0; k < nclasses; k++) {
                                u_night_full(k, j, (arma::uword) u1_moves(l, 0) - 1) += u1_new[i](k, j, l);
                            }
                        }
                    }
                }
                if(writeExt == 0) { 
                    if(saveAll == 1) {
                        out[i + npart * (t + 1)] = List::create(Named("lads") = u_night_lad_cum[i], _["age_region"] = u_night_age_region_cum[i], _["nhsregion"] = u_night_nhsregion, _["age_nhsregion"] = u_night_age_nhsregion_cum[i]);
                    } else {
                        out[i + npart * (t + 1)] = List::create(Named("full") = u_night_full, _["lads"] = u_night_lad_cum[i], _["age_region"] = u_night_age_region_cum[i], _["nhsregion"] = u_night_nhsregion, _["age_nhsregion"] = u_night_age_nhsregion_cum[i]);
                    }
                } else {
                    if(saveAll == 2) {
                        std::sprintf(file_name, "%s/p_%u.csv", std::string(outputName[0]).c_str(), i);
                        file.open(file_name, std::ios::app);
                        for(l = 0; l < nlads; l++) {
                            for(arma::uword r = 0; r < nclasses; r++) {
                                file << t + 1 << ", " << r << ", ";
                                for(j = 0; j < nages; j++) {
                                    file << u_night_full(r, j, l) << ", ";
                                }
                                file << l + 1 << "\n";
                            }
                        }
                        file.close();
                    }
                    
                    std::sprintf(file_name, "%s/p_lads_%u.csv", std::string(outputName[0]).c_str(), i);
                    file.open(file_name, std::ios::app);
                    for(l = 0; l < ndeathlads; l++) {
                        file << t + 1 << ", " << u_night_lad_cum[i](l) << ", ";
                        file << l + 1 << "\n";
                    }
                    file.close();
                    
                    std::sprintf(file_name, "%s/p_age_region_%u.csv", std::string(outputName[0]).c_str(), i);
                    file.open(file_name, std::ios::app);
                    for(l = 0; l < nregions; l++) {
                        file << t + 1 << ", ";
                        for(j = 0; j < nages; j++) {
                            file << u_night_age_region_cum[i](j, l) << ", ";
                        }
                        file << l + 1 << "\n";
                    }
                    file.close();
                    
                    std::sprintf(file_name, "%s/p_nhsregion_%u.csv", std::string(outputName[0]).c_str(), i);
                    file.open(file_name, std::ios::app);
                    for(l = 0; l < nnhsregions; l++) {
                        file << t + 1 << ", ";
                        file << u_night_nhsregion(l) << ", ";
                        file << l + 1 << "\n";
                    }
                    file.close();
                    
                    std::sprintf(file_name, "%s/p_age_nhsregion_%u.csv", std::string(outputName[0]).c_str(), i);
                    file.open(file_name, std::ios::app);
                    for(l = 0; l < nnhsregions; l++) {
                        file << t + 1 << ", ";
                        for(j = 0; j < nnhsages; j++) {
                            file << u_night_age_nhsregion_cum[i](j, l) << ", ";
                        }
                        file << l + 1 << "\n";
                    }
                    file.close();
                }
            }
        }
        
        // update particles for next time point
        for(i = 0; i < npart; i++) {
            u1[i] = u1_new[i];
            u2[i] = u2_new[i];
        }
        
        //calculate block run time
        timer.step("");
        NumericVector res(timer);
        
        // calculate ESS
        if(PF == 1) {
            double ESS = 0.0;
            for(i = 0; i < npart; i++) {
                ESS += pow(weights(i), 2.0);
            }
            ESS = 1.0 / ESS;
            ESS = ESS / ((double) npart);
            
            if(niter == 0) {
                Rprintf("t = %d / %d RESS = %.2f time = %.2f secs \n", t + 1, ndays, ESS, (res[timer_cnt] / 1e9) - prev_time);
            } else  {
                for(i = 0; i < npart; i++) {
                    nacc(i) = (nacc(i) > 0 ? 1:0);
                }
                double accrate = (double) sum(nacc);
                accrate /= ((double) npart);
                Rprintf("t = %d / %d RESS = %.2f nacc = %.2f time = %.2f secs \n", t + 1, ndays, ESS, accrate, (res[timer_cnt] / 1e9) - prev_time);
            }
        } else {
            Rprintf("t = %d / %d time = %.2f secs \n", t + 1, ndays, (res[timer_cnt] / 1e9) - prev_time);
        }
        
        //reset timer and acceptance rate counter
        prev_time = res[timer_cnt] / 1e9;
        timer_cnt++;
    }
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

