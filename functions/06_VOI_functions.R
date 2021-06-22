# # VOI analysis functions for analysis 
# Based on:
# Jalal H, Alarid-Escudero F. A Gaussian Approximation Approach for Value of Information Analysis. Med. Decis. Making. 2017:1–16


gg_color_hue <- function(n, hue_min = 15, hue_max = 375, l = 65, c = 100) {
  hues = seq(hue_min, hue_max, length = n + 1)
  hcl(h = hues, l = l, c = c)[1:n]
}

predict.ga <- function(object, n, n0, verbose = T){
  #### Function to compute the preposterior for each of the 
  #### basis functions of the GAM model.
  #### Inputs:
  #### - object: gam object
  #### - n: scalar or vector of new sample size to compute evsi on
  #### - n0: scalar or vector of effective prior sample size
  #### - verbose: Prints the variance reduction factor for each parameter
  
  ### Name of parameters
  names.data <- colnames(object$model)
  ### Create dataframe with parameter values
  data <- data.frame(object$model[,-1])
  ## Name columns of dataframe 
  colnames(data) <- names.data[-1]
  
  ### Number of parameters
  n.params <- ncol(data)
  
  ### Sanity checks
  if(!(length(n)==1 | length(n)==n.params)){
    stop("Variable 'n' should be either a scalar or a vector 
         the same size as the number of parameters")
  }
  if(!(length(n0)==1 | length(n0)==n.params)){
    stop("Variable 'n0' should be either a scalar or a vector 
         the same size as the number of parameters")
  }
  
  ### Make n & n0 consistent with the number of parameters
  if(length(n) == 1){
    n <- rep(n, n.params)
  }
  if(length(n0) == 1){
    n0 <- rep(n0, n.params)
  }
  
  ### Compute variance reduction factor
  v.ga <- sqrt(n/(n+n0))
  if (verbose){
    print(paste("Variance reduction factor =", round(v.ga, 3)))
  }
  
  ### Number of smoothers
  n.smooth <- length(object$smooth)
  ### Number of total basis functions
  n.colX <- length(object$coefficients)
  ### Number of observations 
  n.rowX <- nrow(object$model)
  
  ### Initialize matrix for preposterior of total basis functions
  X <- matrix(NA, n.rowX, n.colX)
  X[, 1] <- 1
  
  for (k in 1:n.smooth) { # k <- 1 
    klab <- substr(object$smooth[[k]]$label, 1, 1)
    if (klab == "s"){
      Xfrag <- Predict.smooth.ga(object$smooth[[k]], data, v.ga[k])
    } else {
      Xfrag <- Predict.matrix.tensor.smooth.ga(object$smooth[[k]], data, v.ga)
    }
    X[, object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag
  }
  
  ### Coefficients of GAM model
  Beta <- coef(object)
  
  ### Compute conditional Loss
  Ltilde <- X %*% Beta
  
  return(Ltilde)
}


Predict.smooth.ga <- function (object, data, v.ga = 1) {
  #### Function to compute the preposterior for each of the 
  #### basis functions of a smooth for one parameter
  
  ### Produce basis functions for one parameter
  X <- PredictMat(object, data) # ‘mgcv’ version 1.8-17
  ## Number of observations
  n.rct <- nrow(X)
  
  ### Apply variance reduction to compute the preposterior 
  ### for each of the basis functions
  ## Vector of ones
  ones <- matrix(1, n.rct, 1)
  ## Compute phi on each of the basis function
  X.ga <- v.ga*X + (1-v.ga)*(ones %*% colMeans(X))
  
  return(X.ga)
}




Predict.matrix.tensor.smooth.ga <- function (object, 
                                             data, 
                                             v.ga = rep(1, ncol(data))){
  #### Function to compute the preposterior for each of the 
  #### basis functions for one or more parameters and calculates
  #### the tensor product if more than one parameter is selected
  #### (Heavily based on function Predict.matrix.tensor.smooth from
  #### mgcv package)
  
  m <- length(object$margin)
  X <- list()
  for (i in 1:m) { # i <- 1
    term <- object$margin[[i]]$term
    dat <- list()
    for (j in 1:length(term)) { # j <- 1
      dat[[term[j]]] <- data[[term[j]]]
    }
    X[[i]] <- if (!is.null(object$mc[i])) # before: object$mc[i]
      PredictMat(object$margin[[i]], dat, n = length(dat[[1]])) # ‘mgcv’ version 1.8-17
    else Predict.matrix(object$margin[[i]], dat)
    n.rct <- nrow(X[[i]])
  } # end for 'i'
  mxp <- length(object$XP)
  if (mxp > 0) 
    for (i in 1:mxp) if (!is.null(object$XP[[i]])) 
      X[[i]] <- X[[i]] %*% object$XP[[i]]
  
  ### Apply variance reduction to compute the preposterior 
  ### for each of the basis functions
  ## Vector of ones
  ones <- matrix(1, n.rct, 1)
  ## Initialize and fill list with preposterior of basis functions 
  ## for each parameter
  X.ga <- list()
  for (i in 1:m) { # i <- 1
    X.ga[[i]] <- v.ga[i]*X[[i]] + (1-v.ga[i])*(ones %*% colMeans(X[[i]]))
  }
  
  ### Compute tensor product
  T.ga <- tensor.prod.model.matrix(X.ga) # ‘mgcv’ version 1.8-17
  
  return(T.ga)
}


#evppi_lrmm_pop.R
#' Estimation of the Expected Value of Partial Perfect Information (EVPPI)
#' using a linear regression metamodel approach
#'
#' \code{evppi_lrmm_pop} is used to estimate the Expected Value of Partial Perfect
#' Information (EVPPI) using a linear regression metamodel approach from a
#' probabilistic sensitivity analysis (PSA) dataset.
#' @param nmb Matrix of net monetary benefits (NMB). Each column corresponds to
#' the NMB of a different strategy.
#' @param params Vector or matrix of parameters.
#' @param sel.params A vector including the column index of parameters for
#' which EVPPI should be estimated.
#' @param sel.gam Logical variable indicating if a generalized additive model,
#' GAM, (i.e., a spline model) should be fitted (Default = T). If FALSE,
#' a polynomial of degree k is fitted.
#' @param k Scalar with the order of basis functions for the spline model
#' if (sel.gam == T) or the degree of polynomial if (sel.gam != T)
#' @param verbose Logical variable indicating if estimation progress should be
#' reported.
#' @param pop population
#' @keywords Expected Value of Partial Perfect Information
#' @keywords Linear regression metamodel
#' @keywords Splines
#' @details
#' The expected value of partial pefect information (EVPPI) is the expected
#' value of perfect information from a subset of parameters of interest,
#' \eqn{\theta_I} of a cost-effectiveness analysis (CEA) of \eqn{D} different
#' strategies with parameters \eqn{\theta = \{ \theta_I, \theta_C\}}, where
#' \eqn{\theta_C} is the set of complimenatry parameters of the CEA. The
#' function \code{evppi_lrmm} computes the EVPPI of \eqn{\theta_I} from a
#' matrix of net monetary benefits \eqn{B} of the CEA. Each column of \eqn{B}
#' corresponds to the net benefit \eqn{B_d} of strategy \eqn{d}. The function
#' \code{evppi_lrmm} computes the EVPPI using a linear regression metamodel
#' approach following these steps:
#' \enumerate{
#' \item Determine the optimal strategy \eqn{d^*} from the expected net
#' benefits \eqn{\bar{B}}
#' \deqn{d^* = argmax_{d} \{\bar{B}\}}
#' \item Compute the opportunity loss for each \eqn{d} strategy, \eqn{L_d}
#' \deqn{L_d = B_d - B_{d^*}}
#' \item Estimate a linear metamodel for the opportunity loss of each \eqn{d}
#' strategy, \eqn{L_d}, by regressing them on the spline basis functions of
#' \eqn{\theta_I}, \eqn{f(\theta_I)}
#' \deqn{L_d = \beta_0 + f(\theta_I) + \epsilon,}
#' where \eqn{\epsilon} is the residual term that captures the complementary
#' parameters \eqn{\theta_C} and the difference between the original simulation
#' model and the metamodel.
#' \item Compute the EVPPI of \eqn{\theta_I} using the estimated losses for
#' each \eqn{d} strategy, \eqn{\hat{L}_d} from the linear regression metamodel
#' and applying the following equation:
#' \deqn{EVPPI_{\theta_I} = \frac{1}{K}\sum_{i=1}^{K}\max_d(\hat{L}_d)}
#' The spline model in step 3 is fitted using the `mgcv` package.
#' }
#' @return evppi A numeric vector of size one with the EVPPI of the selected
#' parameters
#' @references
#' \enumerate{
#' \item Jalal H, Alarid-Escudero F. A General Gaussian Approximation Approach
#' for Value of Information Analysis. Med Decis Making. 2018;38(2):174-188.
#' \item Strong M, Oakley JE, Brennan A. Estimating Multiparameter Partial
#' Expected Value of Perfect Information from a Probabilistic Sensitivity
#' Analysis Sample: A Nonparametric Regression Approach. Med Decis Making.
#' 2014;34(3):311–26.
#' }
#' @examples
#' ## Load mgcv package and matrixStats
#' library(mgcv)
#' library(matrixStats)
#' ## Load PSA dataset
#' data(syndX)
#' ## Net monetary benefit (NMB) matrix
#' nmb <- syndX[, 5:7]
#' ## Matrix of model parameter inputs values theta
#' theta <- syndX[, 1:4]
#' ## Optimal strategy (d*) based on the highest expected NMB
#' d.star <- which.max(colMeans(nmb))
#' d.star
#' ## Define the Loss matrix
#' loss <- nmb - nmb[, d.star]
#' ## Estimate EVPPI for parameter 1 (MeanVisitsA)
#' evppi_lrmm(nmb = nmb, params = theta, sel.params = 1, verbose = TRUE)
evppi_lrmm_pop <- function (nmb = NULL, params = NULL, sel.params = 1,
                            sel.gam = T, k = NULL,
                            verbose = F,
                            pop=1)
{
  library(mgcv, matrixStats)
  if (is.null(nmb)) {
    stop("A matrix of NMB, 'nmb', hasn't been specified")
  }
  if (is.null(params)) {
    stop("A matrix of parameters, 'params', hasn't been specified")
  }
  if (is.null(dim(nmb))) {
    stop("'nmb' must be an array with at least two strategies")
  }
  n.sel.params <- length(sel.params)
  if (is.null(dim(params))) {
    n.params <- 1
    if (sel.params > 1) {
      stop("Parameter selected is not included in the vector of parameters, 'params'")
    }
  }
  else {
    n.params <- ncol(params)
  }
  if (n.sel.params > n.params) {
    stop("Number of selected parameters exceeds the number of parameters on 'params' (the matrix or vector of parameters)")
  }
  
  ### Check for correct input k
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if(!is.null(k)){
    if(!is.wholenumber(k)){
      stop("Parameter 'k' should be either NULL or an integer")
    }
  }
  
  if (verbose) {
    print(paste("Estimating EVPPI of", n.sel.params, "parameters"))
  }
  n.sim <- nrow(nmb)
  n.strategies <- ncol(nmb)
  d.star <- which.max(colMeans(nmb))
  d.star
  Loss <- nmb - nmb[, d.star]
  lrm <- vector("list", n.strategies)
  Lhatp <- matrix(0, nrow = n.sim, ncol = n.strategies)
  if(sel.gam == T){
    if (is.null(k)){
      for (d in 1:n.strategies) {
        if(d == d.star){
          if (verbose) {
            print(paste("Strategy", d, "is d*; loss of d* = 0"))
          }
          Lhatp[, d] <- 0
        } else {
          if (verbose) {
            print(paste("Constructing metamodel for the Loss of strategy",
                        d))
          }
          if (length(sel.params) == 1) {
            if (n.params == 1) {
              lrm[[d]] <- gam(Loss[, d] ~ s(params))
            }
            else {
              lrm[[d]] <- gam(as.formula(paste("Loss[, d] ~ s(",
                                               colnames(params)[sel.params], ")")), data = params)
            }
          }
          else {
            lrm[[d]] <- gam(as.formula(paste("Loss[, d] ~ s(",
                                             paste(colnames(params[, sel.params]), collapse = ") + s("),
                                             ") + ti(", paste(colnames(params[, sel.params]),
                                                              collapse = ", "), ")")), data = params)
          }
          Lhatp[, d] <- lrm[[d]]$fitted
        }
      }
    } else{
      print(paste0("Spline with k = ", k, " basis functions selected"))
      for (d in 1:n.strategies) {
        if(d == d.star){
          if (verbose) {
            print(paste("Strategy", d, "is d*; loss of d* = 0"))
          }
          Lhatp[, d] <- 0
        } else {
          if (verbose) {
            print(paste("Constructing metamodel for the Loss of strategy",
                        d))
          }
          if (length(sel.params) == 1) {
            if (n.params == 1) {
              lrm[[d]] <- gam(Loss[, d] ~ s(params, k = k))
            }
            else {
              lrm[[d]] <- gam(as.formula(paste("Loss[, d] ~ s(",
                                               colnames(params)[sel.params], ", k = ", k,")")), data = params)
            }
          }
          else {
            lrm[[d]] <- gam(as.formula(paste("Loss[, d] ~ s(",
                                             paste(colnames(params[, sel.params]), collapse = paste0(", k = ", k,") + s(")),
                                             ",k=", k, ") + ti(", paste(colnames(params[, sel.params]),
                                                                        collapse = ", "), ", k=",3,")")), data = params)
          }
          Lhatp[, d] <- lrm[[d]]$fitted
        }
      }
    }
  } else { ## sel.gam != T
    if (is.null(k)){
      k <- 1
    }
    print(paste0("Polynomial of degree k = ", k, " selected"))
    for (d in 1:n.strategies) { # d <- 2
      if(d == d.star){
        if (verbose) {
          print(paste("Strategy", d, "is d*; loss of d* = 0"))
        }
        Lhatp[, d] <- 0
      } else {
        if (verbose) {
          print(paste("Constructing metamodel for the Loss of strategy",
                      d))
        }
        if (length(sel.params) == 1) {
          if (n.params == 1) {
            lrm[[d]] <- lm(Loss[, d] ~ polym(params, degree = k))
          }
          else {
            lrm[[d]] <- lm(as.formula(paste("Loss[, d] ~ polym(",
                                            paste(colnames(params)[sel.params], collapse = ","),
                                            ", degree = ", k, ", raw = TRUE)")), data = params)
          }
        }
        else {
          lrm[[d]] <- lm(as.formula(paste("Loss[, d] ~ polym(",
                                          paste(colnames(params[, sel.params]), collapse = ","),
                                          ", degree = ", k, ", raw = TRUE)")), data = params)
        }
        Lhatp[, d] <- lrm[[d]]$fitted
      }
    }
  }
  evppi <- mean(rowMaxs(Lhatp))*pop
  return(c(evppi = round(evppi, 1)))
}


#evppi_lrmm.R - original
########
#' Estimation of the Expected Value of Partial Perfect Information (EVPPI)
#' using a linear regression metamodel approach
#'
#' \code{evppi_lrmm} is used to estimate the Expected Value of Partial Perfect
#' Information (EVPPI) using a linear regression metamodel approach from a
#' probabilistic sensitivity analysis (PSA) dataset.
#' @param nmb Matrix of net monetary benefits (NMB). Each column corresponds to
#' the NMB of a different strategy.
#' @param params Vector or matrix of parameters.
#' @param sel.params A vector including the column index of parameters for
#' which EVPPI should be estimated.
#' @param sel.gam Logical variable indicating if a generalized additive model,
#' GAM, (i.e., a spline model) should be fitted (Default = T). If FALSE,
#' a polynomial of degree k is fitted.
#' @param k Scalar with the order of basis functions for the spline model
#' if (sel.gam == T) or the degree of polynomial if (sel.gam != T)
#' @param verbose Logical variable indicating if estimation progress should be
#' reported.
#' @keywords Expected Value of Partial Perfect Information
#' @keywords Linear regression metamodel
#' @keywords Splines
#' @details
#' The expected value of partial pefect information (EVPPI) is the expected
#' value of perfect information from a subset of parameters of interest,
#' \eqn{\theta_I} of a cost-effectiveness analysis (CEA) of \eqn{D} different
#' strategies with parameters \eqn{\theta = \{ \theta_I, \theta_C\}}, where
#' \eqn{\theta_C} is the set of complimenatry parameters of the CEA. The
#' function \code{evppi_lrmm} computes the EVPPI of \eqn{\theta_I} from a
#' matrix of net monetary benefits \eqn{B} of the CEA. Each column of \eqn{B}
#' corresponds to the net benefit \eqn{B_d} of strategy \eqn{d}. The function
#' \code{evppi_lrmm} computes the EVPPI using a linear regression metamodel
#' approach following these steps:
#' \enumerate{
#' \item Determine the optimal strategy \eqn{d^*} from the expected net
#' benefits \eqn{\bar{B}}
#' \deqn{d^* = argmax_{d} \{\bar{B}\}}
#' \item Compute the opportunity loss for each \eqn{d} strategy, \eqn{L_d}
#' \deqn{L_d = B_d - B_{d^*}}
#' \item Estimate a linear metamodel for the opportunity loss of each \eqn{d}
#' strategy, \eqn{L_d}, by regressing them on the spline basis functions of
#' \eqn{\theta_I}, \eqn{f(\theta_I)}
#' \deqn{L_d = \beta_0 + f(\theta_I) + \epsilon,}
#' where \eqn{\epsilon} is the residual term that captures the complementary
#' parameters \eqn{\theta_C} and the difference between the original simulation
#' model and the metamodel.
#' \item Compute the EVPPI of \eqn{\theta_I} using the estimated losses for
#' each \eqn{d} strategy, \eqn{\hat{L}_d} from the linear regression metamodel
#' and applying the following equation:
#' \deqn{EVPPI_{\theta_I} = \frac{1}{K}\sum_{i=1}^{K}\max_d(\hat{L}_d)}
#' The spline model in step 3 is fitted using the `mgcv` package.
#' }
#' @return evppi A numeric vector of size one with the EVPPI of the selected
#' parameters
#' @references
#' \enumerate{
#' \item Jalal H, Alarid-Escudero F. A General Gaussian Approximation Approach
#' for Value of Information Analysis. Med Decis Making. 2018;38(2):174-188.
#' \item Strong M, Oakley JE, Brennan A. Estimating Multiparameter Partial
#' Expected Value of Perfect Information from a Probabilistic Sensitivity
#' Analysis Sample: A Nonparametric Regression Approach. Med Decis Making.
#' 2014;34(3):311–26.
#' }
#' @examples
#' ## Load mgcv package and matrixStats
#' library(mgcv)
#' library(matrixStats)
#' ## Load PSA dataset
#' data(syndX)
#' ## Net monetary benefit (NMB) matrix
#' nmb <- syndX[, 5:7]
#' ## Matrix of model parameter inputs values theta
#' theta <- syndX[, 1:4]
#' ## Optimal strategy (d*) based on the highest expected NMB
#' d.star <- which.max(colMeans(nmb))
#' d.star
#' ## Define the Loss matrix
#' loss <- nmb - nmb[, d.star]
#' ## Estimate EVPPI for parameter 1 (MeanVisitsA)
#' evppi_lrmm(nmb = nmb, params = theta, sel.params = 1, verbose = TRUE)
evppi_lrmm <- function (nmb = NULL, params = NULL, sel.params = 1,
                        sel.gam = T, k = NULL,
                        verbose = F)
{
  library(mgcv, matrixStats)
  if (is.null(nmb)) {
    stop("A matrix of NMB, 'nmb', hasn't been specified")
  }
  if (is.null(params)) {
    stop("A matrix of parameters, 'params', hasn't been specified")
  }
  if (is.null(dim(nmb))) {
    stop("'nmb' must be an array with at least two strategies")
  }
  n.sel.params <- length(sel.params)
  if (is.null(dim(params))) {
    n.params <- 1
    if (sel.params > 1) {
      stop("Parameter selected is not included in the vector of parameters, 'params'")
    }
  }
  else {
    n.params <- ncol(params)
  }
  if (n.sel.params > n.params) {
    stop("Number of selected parameters exceeds the number of parameters on 'params' (the matrix or vector of parameters)")
  }
  
  ### Check for correct input k
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if(!is.null(k)){
    if(!is.wholenumber(k)){
      stop("Parameter 'k' should be either NULL or an integer")
    }
  }
  
  if (verbose) {
    print(paste("Estimating EVPPI of", n.sel.params, "parameters"))
  }
  n.sim <- nrow(nmb)
  n.strategies <- ncol(nmb)
  d.star <- which.max(colMeans(nmb))
  d.star
  Loss <- nmb - nmb[, d.star]
  lrm <- vector("list", n.strategies)
  Lhatp <- matrix(0, nrow = n.sim, ncol = n.strategies)
  if(sel.gam == T){
    if (is.null(k)){
      for (d in 1:n.strategies) {
        if(d == d.star){
          if (verbose) {
            print(paste("Strategy", d, "is d*; loss of d* = 0"))
          }
          Lhatp[, d] <- 0
        } else {
          if (verbose) {
            print(paste("Constructing metamodel for the Loss of strategy",
                        d))
          }
          if (length(sel.params) == 1) {
            if (n.params == 1) {
              lrm[[d]] <- gam(Loss[, d] ~ s(params))
            }
            else {
              lrm[[d]] <- gam(as.formula(paste("Loss[, d] ~ s(",
                                               colnames(params)[sel.params], ")")), data = params)
            }
          }
          else {
            lrm[[d]] <- gam(as.formula(paste("Loss[, d] ~ s(",
                                             paste(colnames(params[, sel.params]), collapse = ") + s("),
                                             ") + ti(", paste(colnames(params[, sel.params]),
                                                              collapse = ", "), ")")), data = params)
          }
          Lhatp[, d] <- lrm[[d]]$fitted
        }
      }
    } else{
      print(paste0("Spline with k = ", k, " basis functions selected"))
      for (d in 1:n.strategies) {
        if(d == d.star){
          if (verbose) {
            print(paste("Strategy", d, "is d*; loss of d* = 0"))
          }
          Lhatp[, d] <- 0
        } else {
          if (verbose) {
            print(paste("Constructing metamodel for the Loss of strategy",
                        d))
          }
          if (length(sel.params) == 1) {
            if (n.params == 1) {
              lrm[[d]] <- gam(Loss[, d] ~ s(params, k = k))
            }
            else {
              lrm[[d]] <- gam(as.formula(paste("Loss[, d] ~ s(",
                                               colnames(params)[sel.params], ", k = ", k,")")), data = params)
            }
          }
          else {
            lrm[[d]] <- gam(as.formula(paste("Loss[, d] ~ s(",
                                             paste(colnames(params[, sel.params]), collapse = paste0(", k = ", k,") + s(")),
                                             ",k=", k, ") + ti(", paste(colnames(params[, sel.params]),
                                                                        collapse = ", "), ", k=",3,")")), data = params)
          }
          Lhatp[, d] <- lrm[[d]]$fitted
        }
      }
    }
  } else { ## sel.gam != T
    if (is.null(k)){
      k <- 1
    }
    print(paste0("Polynomial of degree k = ", k, " selected"))
    for (d in 1:n.strategies) { # d <- 2
      if(d == d.star){
        if (verbose) {
          print(paste("Strategy", d, "is d*; loss of d* = 0"))
        }
        Lhatp[, d] <- 0
      } else {
        if (verbose) {
          print(paste("Constructing metamodel for the Loss of strategy",
                      d))
        }
        if (length(sel.params) == 1) {
          if (n.params == 1) {
            lrm[[d]] <- lm(Loss[, d] ~ polym(params, degree = k))
          }
          else {
            lrm[[d]] <- lm(as.formula(paste("Loss[, d] ~ polym(",
                                            paste(colnames(params)[sel.params], collapse = ","),
                                            ", degree = ", k, ", raw = TRUE)")), data = params)
          }
        }
        else {
          lrm[[d]] <- lm(as.formula(paste("Loss[, d] ~ polym(",
                                          paste(colnames(params[, sel.params]), collapse = ","),
                                          ", degree = ", k, ", raw = TRUE)")), data = params)
        }
        Lhatp[, d] <- lrm[[d]]$fitted
      }
    }
  }
  evppi <- mean(rowMaxs(Lhatp))
  return(c(evppi = round(evppi, 1)))
}


#' Expected Value of Perfect Information (EVPI)
#'
#' \code{evpi} is used to compute the expected value of perfect information 
#' (EVPI) from a probabilistic sensitivity analysis (PSA) dataset.
#' @param v.wtp Numeric vector with willingness-to-pay (WTP) thresholds
#' @param m.e Matrix of effectiveness. Each column corresponds to a vector of
#' effectiveness.
#' @param m.c Matrix of costs. Each column corresponds to a vector of
#' costs.
#' @param pop A scalar that corresponds to the total population
#' @keywords expected value of perfect information; net monetary benefit
#' @section Details:
#' \code{evpi} calculates the value of eliminating all the uncertainty of a 
#' cost-effectiveness analysis at each WTP threshold.
#' @return evpi A data frame with the EVPI at each WTP threshold. 
#'
evpi <- function(v.wtp, m.e, m.c, pop = 1){
  # Load required packages
  require(matrixStats)
  if(!(ncol(m.e) == ncol(m.c))){
    stop("Matrices of effectiveness and costs do not have same number of strategies.")
  }
  if(ncol(m.e)<2){
    stop("You need at least two different strategies to compute EVPI.")
  }
  # Create scalar with number of simulations
  n.sim <- nrow(m.e)
  # Create scalar with number of strategies (i.e. number of columns of 
  # effectiveness matrix)
  n.str <- ncol(m.e)
  # Data frame to store EVPI for each WTP threshold
  df.evpi <- as.data.frame(array(0, dim = c(length(v.wtp), 2)))
  # Name data frame's columns
  colnames(df.evpi) <- c("WTP", "EVPI")
  # Assign vector with WTP thresholds to first column of `evpi`
  df.evpi$WTP <- v.wtp
  # Estimate the Loss matrix and EVPI at each WTP threshold
  for(l in 1:length(v.wtp)){
    # Compute NMB with vector indexing
    nmb <-  v.wtp[l]*m.e - m.c
    ## Find the optimal strategy with current info
    d.star <- which.max(colMeans(nmb))
    ## Calculate the opportunity loss from choosing d.star for each strategy
    loss <- nmb - nmb[, d.star]
    ## Compute EVPI
    df.evpi$EVPI[l] <- mean(rowMaxs(as.matrix(loss))) * pop# needs to be a numeric matrix
  }
  #Return a data frame
  class(df.evpi) <- "evpi"
  return(df.evpi)
}


ceaf_changed <- function(v.wtp, strategies = NULL, m.e, m.c, currency = "$", effectunit = "LY", ceaf.out = FALSE){
  # Load required packages
  library(reshape2)
  library(ggplot2)
  library(scales)
  if(!ncol(m.e) == ncol(m.c)){
    stop("Matrices of effectiveness and costs do not have same number of strategies.")
  }
  if(ncol(m.e) < 2){
    stop("You need at least two different strategies to compute EVPI.")
  }
  # Create scalar with number of simulations
  n.sim <- nrow(m.e)
  # Create scalar with number of strategies (i.e. number of columns of 
  # effectiveness matrix)
  n.str <- ncol(m.e)
  # Matrix to store indicator of CEAC
  m.cea  <- array(0, dim = c(length(v.wtp), n.str))
  # Vector to store indicator of strategy at CEAF
  v.ceaf <- array(0, dim = c(length(v.wtp), 1))
  
  # If the name of the strategies is not provided, generate a generic vector
  # with strategy names
  if (is.null(strategies)){
    strategies <- paste(rep("Strategy_", n.strategies), seq(1, n.strategies), sep = "")
  }
  for(l in 1:length(v.wtp)){
    m.nmb <-  v.wtp[l] * m.e - m.c # Effectiveness minus Costs
    # Calculate point of CEAF, i.e., the strategy with the highest expected NMB
    v.ceaf[l, 1] <- which.max(colMeans(m.nmb))
    # Calculate points in CEAC, i.e, the probability that each strategy is cost-effective
    max.nmb <- max.col(m.nmb)
    opt <- table(max.nmb)
    m.cea[l, as.numeric(names(opt))] <- opt/n.sim
  }
  m.ceaf <- m.cea[cbind(1:length(v.wtp), v.ceaf)]
  
  df.cea <- data.frame(cbind(v.wtp, m.cea, m.ceaf))
  colnames(df.cea) <- c("WTP", strategies, "Frontier")
  
  df.ceac <- melt(df.cea, 
                  id.vars = "WTP", 
                  variable.name = "Strategy")
  
  ## Plot CEAC & CEAF
  # Format to plot frontier
  strats <- 1:(length(unique(df.ceac$Strategy))-1)
  point.shapes <- c(strats+14, 0) # Shapes: http://sape.inf.usi.ch/quick-reference/ggplot2/shape
  colors <- c(gg_color_hue(length(strats)), "#696969")
  point.size <- c(rep(2, length(strats)), 4) # Trick consists on firts define size as aes then manually change it
  # Plot CEAC & CEAF
  gg.ceaf <- ggplot(data = df.ceac, aes(x = WTP/1000, y = value)) +
    geom_point(aes(shape = Strategy, color = Strategy, size = Strategy)) +
    geom_line(aes(color = Strategy)) +
    ggtitle("Cost-Effectiveness Acceptability Curves and Frontier") + 
    # scale_colour_hue(l=50, values=colors) +
    scale_x_continuous(breaks = number_ticks(20))+
    scale_shape_manual(values = point.shapes) +
    #scale_shape(solid=TRUE) +
    scale_color_manual(values = colors) + 
    scale_size_manual(values = point.size) +
    #scale_alpha_manual(values=c(rep(0, length(strats)), 0.5)) + 
    xlab(paste("Willingness to pay (Thousand ", currency, "/", effectunit, ")", sep = "")) +
    ylab("Pr Cost-Effective") +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
  # Return a data frame of class ceaf
  class(df.ceac) <- "ceaf"
  if(ceaf.out){
    print(gg.ceaf)
  }
  out <- list(df.ceac = df.ceac,
              gg.ceaf = gg.ceaf)
  #return(gg.ceaf)
  return(out)
}


#gg_color_hue.R
#' \code{ggplot2}-like colour scale in HCL space
#'
#' Function for number of ticks in axis of \code{ggplot2} plots.
#' @param n Number of colours to return
#' @param hue_min Minimum hue value in the range [0, 360].
#' @param hue_max Maximum hue value in the range [0, 360].
#' @param l Luminance in the range [0, 100].
#' @param c Chroma of the colour.
#' @keywords ggplot2
#' @section Details:
#' @examples 
#' gg_color_hue(10)


# Cost of Research
CostRes <- function(fixed.cost = 0, 
                    samp.size, 
                    cost.per.patient, 
                    INMB, 
                    clin.trial = TRUE, n.arms = 2){
  # Computes the cost of collecting information (i.e., through a research study)
  #
  # Args:
  #   fixed.cost:       fixed cost of collecting information
  #                     (e.g., fixed cost of a clinical trial); default = 0
  #   samp.size:               vector with sample sizes
  #   cost.per.patient: cost per patient in research study
  #   INMB:             Incremental Net Monetary Benefit
  #   clin.trial:       indicator whether calculation is for a clinical trial;
  #                     default = TRUE
  #   n.arms:           Number of arms in research study design; default = 2
  #
  # Returns:
  #   cost.res: vector with the total cost of collecting information for each simple size
  #
  if (clin.trial){
    Cost.Res <- fixed.cost + n.arms*samp.size*cost.per.patient + samp.size*INMB
  } else { # E.g., cohort study
    Cost.Res <- fixed.cost + samp.size*cost.per.patient
  }
  return(Cost.Res)
}

## Obtain the NMB for all strategies 
NMB_function <- function(wtp, effect, costs){
  WTP <- wtp
  ly.s1 <- effect[,1]
  ly.s2 <- effect[,2]
  cost.s1 <- costs[,1]
  cost.s2 <- costs[,2]
  
  ## NMB for strategy 1-2
  #Strategy 1
  nmb.s1 <- (ly.s1 * WTP) - cost.s1
  #Strategy 2
  nmb.s2 <- (ly.s2 * WTP) - cost.s2
  
  # Save the results
  
  nmb <- matrix(c(nmb.s1, nmb.s2), ncol = 2, byrow = FALSE)
  nmb.l.names <- c("nmb.s1", "nmb.s2")
  colnames(nmb) <- nmb.l.names
  
  return(nmb)
}

NMB_function2 <- function(wtp, effect, costs){
  WTP <- wtp
  ly.s1 <- effect[,1]
  ly.s2 <- effect[,2]
  cost.s1 <- costs[,1]
  cost.s2 <- costs[,2]
  
  ## NMB for strategy 1-2
  #Strategy 1
  nmb.s1 <- ly.s1-(cost.s1/WTP)
  #Strategy 2
  nmb.s2 <- ly.s2-(cost.s2/WTP)
  
  # Save the results
  
  nmb<- matrix(c(nmb.s1, nmb.s2), ncol=2, byrow=FALSE)
  nmb.l.names <- c("nmb.s1", "nmb.s2")
  colnames(nmb)<- nmb.l.names
  
  return(nmb)
}


NMB_function3 <- function(wtp, effect){
  WTP <- wtp
  ly.s1 <- effect[,1]
  ly.s2 <- effect[,2]  
  
  ## NMB for strategy 1-2
  #Strategy 1
  nmb.s1 <- ly.s1
  #Strategy 2
  nmb.s2 <- ly.s2
  
  # Save the results
  
  nmb<- matrix(c(nmb.s1, nmb.s2), ncol=2, byrow=FALSE)
  nmb.l.names <- c("nmb.s1", "nmb.s2")
  colnames(nmb)<- nmb.l.names
  
  return(nmb)
}


