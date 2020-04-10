rm(list = ls())


library(tidyverse)
data <- read_csv("dati.csv")


ridge <- function(x, y, lambda) {
  n <- nrow(x)
  p <- ncol(x)
  beta <- solve(crossprod(x) + lambda * diag(p)) %*% crossprod(x, y)
  #beta <- solve(t(X) %*% X + lambda*diag(p)) %*% t(X) %*% y
  y_hat <- x %*% beta
  list(coef = beta,
       prevision = y_hat)
}


lars <- function(X, y, tol = 1e-10) {
  
  "DOCUMENTATION:
  This function project the y vector into the space L between
  the active variables (active); then it calculate the direction
  of the update (u) and its length (gamma).
  After max_iter iterations, the algorithm reaches convergence
  and than beta can be computed (using OLS method) between the
  active matrix Xa and the projection mu.
  Direction of the update is equiangular between each active
  variable at each step (and here come the name LAR).
  
  See _Least Angle Regression_ (Efron, Hastie, Johnstone,
  Tibshirami), January 2003, for details."
  
  least_squares <- function(x, y)
    solve(crossprod(x)) %*% crossprod(x, y)
  
  n <- nrow(X)
  p <- ncol(X)
  max_iter <- min(n - 1, p)
  
  beta_tot <- matrix(0, ncol = max_iter, nrow = p)
  rownames(beta_tot) <- colnames(X)
  
  # page 6
  ## mu is the projection of y on L
  mu <- matrix(rep(0, n))             # n × 1
  
  for (i in 1:max_iter) {
    # equation 2.1
    ## each c_j represent the correlation of the j-th variable
    ## between X and the projection on the sub-space L
    c_hat <- crossprod(X, y - mu)   # vector, p × 1
    
    # equation 2.9
    ## the "active" subset includes each variable that is as much
    ## correlate with y as the maximum-correlated x_j
    C <- as.double(max(abs(c_hat)))  # scalar
    active <- abs((abs(c_hat) - C)) <= tol  # due to R approximation
    alpha <- sum(active)            # scalar, value of a
    
    # equation 2.10
    ## a vector of signs of correlation, for multiply the X matrix
    ## to use only positive correlations
    s <- as.vector(sign(c_hat))     # vector, p × 1
    
    # equation 2.4
    ## the X matrix that includes only active variables (with
    ## positive correlation)
    Xa <- (X %*% diag(s))[, active]  # matrix, n × a
    
    # equation 2.5
    # (inverse is computed for performance)
    ## this part is quite complicated, see Paper for details
    Ga <- solve(crossprod(Xa))      # matrix, a × a
    ones <- matrix(rep(1, alpha))   # vector, a × 1
    A <- as.double(crossprod(ones, Ga) %*% ones) ^-0.5
    # scalar
    
    # equation 2.6
    ## u is the direction of the update
    w <- A * Ga %*% ones            # vector a × 1
    u <- Xa %*% w                   # vector n × 1
    
    # equation 2.11
    ## this part is not well described in the Paper
    a <- crossprod(X, u)            # vector p × 1
    
    # equation 2.13
    ## gamma is the intensity of the update: "We take the largest
    ## step possible in the direction of this predictor" (page 5)
    gamma <- Inf                    # scalar
    ## 2p passages: the for loop is not critical
    ## TODO? implement in R-Cpp
    for (j in 1:p) {  # functional "min+"
      cj <- c_hat[j, 1]
      aj <- a[j, 1]
      Ac <- c((C - cj) / (A - aj),
              (C + cj) / (A + aj))
      for (new_gamma in Ac) {
        if (!is.nan(new_gamma) & gamma > new_gamma & new_gamma > 0)
          gamma <- new_gamma
      }
    }
    
    # equation 2.12
    ## mu is now updated; the updated value is than used to
    ## compute the OLS solution and compute new beta vector
    mu <- mu + gamma * u
    beta_tot[active, i] <- least_squares(Xa, mu)
  }
  
  ## at the end of the iterative process, beta is returned as
  ## parameters of the model; the projection is used as prevision of
  ## the training data
  list(coef = beta_tot[, max_iter, drop = FALSE],
       prevision = mu,
       log = beta_tot)
}


train <- function(x_train, y_train, method, ...) {
  ## x_train
  ## y_train
  ## x_test
  ## y_test
  ## mod <- train(x_train, y_train, lars)
  ## mod(x_test) >>> y_hat
  n <- nrow(x_train)
  p <- ncol(x_train)
  
  standardizer <- function(x) {
    mu <- colMeans(x)
    sigma <- as.double(sqrt(diag(var(x))))
    function(y, reverse = FALSE) {
      if (reverse) {
        I <- diag(length(sigma)) * sigma
        t(t(y %*% I) + mu)
      } else {
        I <- diag(length(sigma)) / sigma
        t(t(y) - mu) %*% I
      }
    }
  }
  
  x_standardizer <- standardizer(x_train)
  y_standardizer <- standardizer(y_train)
  
  beta <- as.matrix(method(x_standardizer(x_train), y_standardizer(y_train), ...)$coef)
  
  predict <- function(x_new) {
    y_standardizer(x_standardizer(x_new) %*% beta, reverse = TRUE)
  }
}