/******************************************************************************
 *
 * TODO: Project Title
 *
 * Copyright (C) 2012-2013 Lars Simon Zehnder. All Rights Reserved.
 * Web: -
 *
 * Author: Lars Simon Zehnder <simon.zehnder@gmail.com>
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 ******************************************************************************/

#ifndef __ALGORITHMS_H__
#define __ALGORITHMS_H__

#include <RcppArmadilloExtensions/sample.h>

template<typename T, typename U>
inline
  T sample_arma(const T &x, const int &size, 
                const bool &replace, const arma::vec &prob) 
  {
    U RcppSample = Rcpp::as<U>(Rcpp::wrap(x));    
    Rcpp::NumericVector RcppProb = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(prob));
    U RcppRet = Rcpp::RcppArmadillo::sample(RcppSample, 
                                            size, replace, RcppProb);
    T ret(RcppRet.begin(), size, false, true);
    return ret;
  }
#endif /* __ALGORITHMS_H__ */