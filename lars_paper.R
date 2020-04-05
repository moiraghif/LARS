rm(list = ls())


lars <- function(formula, tol = 1e-10) {

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

    y <- eval(formula[[2]], parent.frame())
    X <- eval(formula[[3]], parent.frame())
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

    beta <- as.matrix(method(y_standardizer(y_train) ~ x_standardizer(x_train), ...)$coef)

    predict <- function(x_new) {   
        y_standardizer(x_standardizer(x_new) %*% beta, reverse = TRUE)
    }
}



get_data <- function(n, p, random_fn) {
    "Generate some random data: p features mutually independent taken randomly"
    out <- matrix(nrow = n, ncol = p)
    for (i in 1:p) {
        out[, i] <- random_fn(n)
    }
    out
}

get_function <- function(beta, sigma) {
    "Generate values according to a ground truth function"
    function(x, noise = TRUE) {
        n <- nrow(x)
        epsilon <- rnorm(n, 0, ifelse(noise, sigma, 0))
        x %*% beta + epsilon
    }
}


simulation_fixed <- function(data_generator, beta, sigma, method, iterations, verbose = TRUE) {
    error <- function(y, y_hat) mean((y - y_hat)^2)
    phi <- get_function(beta, sigma)
    x <- data_generator()
    mu <- phi(x, noise = FALSE)
    n <- nrow(x)
    p <- ncol(x)

    err_total <- c()
    y_hat_total <- matrix(0, nrow = n, ncol = iterations)
    x <- data_generator()
    mu <- phi(x, noise = FALSE)
    start_time <- Sys.time()
    for (i in 1:iterations) {
        y <- phi(x)
        model <- train(x, y, method = method)
        y_hat <- model(x)
        y_hat_total[, i] <- y_hat
        err_total <- c(err_total, error(y, y_hat))
    }
    stop_time <- Sys.time()
    av_error = mean(err_total)
    variance = mean(apply(y_hat, 2, var))
    bias = mean((colMeans(y_hat) - mu)^2)
    ## TODO: il valore teorico e calcolato dell'MSE sono MOLTO diversi
    if (verbose) {
        print(paste0("Execution time: ", stop_time - start_time))
        print(paste0("Sigma ^2: ", sigma^2))
        print(paste0("Variance: ", variance))
        print(paste0("Bias ^2:  ", bias))
        print(paste0("MSE (theoretical): ", bias + variance))
        print(paste0("MSE (simulation):  ", av_error))
        print(paste0("       difference: ", round(abs(bias + variance - av_error), 4)))
        print(paste0("Prediction Error: ", bias + variance + sigma^2))
    }
    list("bias^2" = bias,
         variance = variance,
         sigma = sigma,
         "MSE_theoretical" = bias + variance,
         "MSE" = av_error,
         "Prediction error" = bias + variance + sigma^2,
         "execution time" = stop_time - start_time)
}



n <- 30
p <- 50
p_useless <- 15
sigma <- 0.001  # noise variance

beta_generator <- function(n) matrix(runif(n, +1, +2.5))
beta_true <- beta_generator(p)
beta_true[sample(p_useless), 1] <- 0

x_generator <- function() get_data(n, p, function(n) runif(n, -10, +10))

## strting simulation (fixed settings)
monte_carlo <- 1e3
lars_simulation <- simulation_fixed(x_generator, beta_true, sigma, lars, monte_carlo)
