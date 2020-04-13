rm(list = ls())


ridge <- function(x, y, lambda) {
  p <- ncol(x)
  beta <- solve(crossprod(x) + lambda * diag(p), tol = 0) %*% crossprod(x, y)
  ## same as: beta <- solve(t(X) %*% X + lambda*diag(p)) %*% t(X) %*% y
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
    ## execution example: mod <- train(x_train, y_train, lars)
    ## mod(x_test)$prediction >>> y_hat

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

    beta <- as.matrix(method(x_standardizer(x_train),
                          y_standardizer(y_train), ...)$coef)

    predict <- function(x_new) {
        x_standardized <- x_standardizer(x_new)
        y_hat <- x_standardized %*% beta
        predictions <- y_standardizer(y_hat, reverse = TRUE)
        list(prediction = predictions,
             y_standardized = y_hat,
             x_standardized = x_standardized,
             coef = beta)
    }
}


get_data <- function(n, p, random_fn) {
    "Generate some random data: p features mutually independent taken randomly"
    out <- matrix(nrow = n, ncol = p)
    for (i in 1:p)out[, i] <- random_fn(n)
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


simulation_fixed <- function(data_generator, beta, sigma, method,
                      iterations, verbose = TRUE, ...) {
    error <- function(y, y_hat) mean((y - y_hat)^2)
    phi <- get_function(beta, sigma)
    x <- data_generator()
    mu <- phi(x, noise = FALSE)
    n <- nrow(x)

    err_total <- c()
    bias_total <- c()
    variance_total <- c()
    sum_total <- c()
    y_hat_total <- matrix(NA, nrow = n, ncol = iterations)

    start_time <- Sys.time()
    for (i in 1:iterations) {
        y_train <- phi(x)
        y_test <- phi(x)
        model <- train(x, y_train, method = method, ...)
        y_hat <- model(x)$prediction
        y_hat_total[, i] <- y_hat
        err_total <- c(err_total, error(y_test, y_hat))
        bias_total <- c(bias_total,
                        mean((rowMeans(y_hat_total, na.rm = T) - mu)^2))
        variance_total <- c(variance_total,
                            mean(apply(y_hat_total, 1, var, na.rm = TRUE)))
        sum_total <- c(sum_total, bias_total[i] + variance_total[i])
        ## Print completion time
        cat("\r", "Iteration n.", i)
    }
    stop_time <- Sys.time()
    av_error <- mean(err_total)
    variance <- variance_total[iterations]
    bias <- bias_total[iterations]
    if (verbose) {
        print(paste0("Execution time: ", stop_time - start_time))
        print(paste0("Sigma²:   ", sigma^2))
        print(paste0("Variance: ", variance))
        print(paste0("Bias²:    ", bias))
        print(paste0("MSE (theoretical): ", bias + variance))
        print(paste0("MSE (simulation):  ", av_error))
        print(paste0("       difference: ", round(abs(bias + variance -
                                                      av_error), 4)))
        print(paste0("Prediction Error: ", bias + variance + sigma^2))
    }
    list("bias^2" = bias,
         variance = variance,
         sigma = sigma,
         "MSE_theoretical" = bias + variance,
         "MSE" = av_error,
         "Prediction error" = bias + variance + sigma^2,
         "execution time" = stop_time - start_time,
         "err_log" = err_total,
         "variance_log" = variance_total,
         "bias_log" = bias_total,
         "log" = sum_total)
}


simulation_random <- function(data_generator, beta, sigma,
                       method, iterations, verbose = TRUE, ...) {
    error <- function(y, y_hat) mean((y - y_hat)^2)
    phi <- get_function(beta, sigma)

    err_total <- c()
    bias_total <- c()
    variance_total <- c()
    sum_total <- c()
    x_test <- data_generator()
    mu <- phi(x_test, noise = FALSE)
    y_hat_total <- matrix(NA, nrow = n, ncol = iterations)

    start_time <- Sys.time()
    for (i in 1:iterations) {
        x_train <- data_generator()
        y_train <- phi(x_train)
        n <- nrow(x_train)
        y_test <- phi(x_test)
        model <- train(x_train, y_train, method = method, ...)
        y_hat <- model(x_test)$prediction
        y_hat_total[, i] <- y_hat
        err_total <- c(err_total, error(y_test, y_hat))
        bias_total <- c(bias_total,
                        mean((rowMeans(y_hat_total, na.rm = T) - mu)^2))
        variance_total <- c(variance_total,
                            mean(apply(y_hat_total, 1, var, na.rm = TRUE)))
        sum_total <- c(sum_total, bias_total[i] + variance_total[i])
        ## Print completion time
        cat("\r", "Iteration n.", i)
    }
    stop_time <- Sys.time()
    av_error <- mean(err_total)
    variance <- variance_total[iterations]
    bias <- bias_total[iterations]
    if (verbose) {
        print(paste0("Execution time: ", stop_time - start_time))
        print(paste0("Sigma²:   ", sigma^2))
        print(paste0("Variance: ", variance))
        print(paste0("Bias²:    ", bias))
        print(paste0("MSE (theoretical): ", bias + variance))
        print(paste0("MSE (simulation):  ", mean(err_total)))
        print(paste0("       difference: ", round(abs(bias + variance -
                                                      av_error), 4)))
        print(paste0("Prediction Error: ", bias + variance + sigma^2))
    }
    list("bias^2" = bias,
         variance = variance,
         sigma = sigma,
         "MSE_theoretical" = bias + variance,
         "MSE" = av_error,
         "Prediction error" = bias + variance + sigma^2,
         "execution time" = stop_time - start_time,
         "err_log" = err_total,
         "variance_log" = variance_total,
         "bias_log" = bias_total,
         "log" = sum_total)
}


n <- 30
p <- 50
p_useless <- 15
sigma <- 0.01  # noise variance
monte_carlo <- 1e3

beta_generator <- function(n) matrix(runif(n, +1, +2.5))
beta_true <- beta_generator(p)
beta_true[sample(1:p, p_useless), 1] <- 0

x_generator <- function() get_data(n, p, function(n) runif(n, -10, +10))


## LARS SIMULATIONS

## LARS fixed simulation
lars_fixed <- simulation_fixed(x_generator, beta_true, sigma,
                               lars, monte_carlo)

## LARS random simulation
lars_random <- simulation_random(x_generator, beta_true, sigma,
                                 lars, monte_carlo)


## RIDGE SIMULATIONS (1)

## RIDGE fixed simulation
lambdas <- seq(.5, 10, .5)
ridge_fixed <- list(lambda = c(),
                    MSE = c(),
                    MSE_t = c(),
                    bias = c(),
                    variance = c(),
                    pred_error = c())
for (lambda in lambdas) {
    sim <- simulation_fixed(x_generator, beta_true, sigma,
                            ridge, monte_carlo, lambda = lambda)
    ridge_fixed[["lambda"]] <- c(ridge_fixed[["lambda"]], lambda)
    ridge_fixed[["MSE_t"]] <- c(ridge_fixed[["MSE_t"]], sim$MSE_theoretical)
    ridge_fixed[["MSE"]]   <- c(ridge_fixed[["MSE"]], sim$MSE)
    ridge_fixed[["bias"]]     <- c(ridge_fixed[["bias"]], sim$"bias^2")
    ridge_fixed[["variance"]] <- c(ridge_fixed[["variance"]], sim$variance)
    ridge_fixed[["pred_error"]] <- c(ridge_fixed[["pred_error"]],
                                   sim$"Prediction error")
}

## RIDGE random simulation
lambdas <- seq(.5, 10, .5)
ridge_random <- list(lambda = c(),
                     MSE = c(),
                     MSE_t = c(),
                     bias = c(),
                     variance = c(),
                     pred_error = c())
for (lambda in lambdas) {
    sim <- simulation_random(x_generator, beta_true, sigma,
                             ridge, monte_carlo, lambda = lambda)
    ridge_random[["lambda"]]   <- c(ridge_random[["lambda"]], lambda)
    ridge_random[["MSE_t"]] <- c(ridge_random[["MSE_t"]], sim$MSE_theoretical)
    ridge_random[["MSE"]]   <- c(ridge_random[["MSE"]], sim$MSE)
    ridge_random[["bias"]]     <- c(ridge_random[["bias"]], sim$"bias^2")
    ridge_random[["variance"]] <- c(ridge_random[["variance"]], sim$variance)
    ridge_random[["pred_error"]] <- c(ridge_random[["pred_error"]],
                                    sim$"Prediction error")
}


## save workspace
save.image("simulation.RData")



## RIDGE SIMULATIONS (2)
tr <- function(x) sum(diag(x))
lambdas <- seq(0.5, 5, .5)
phi <- get_function(beta_true, sigma)

monte_carlo <- 1e3

## RIDGE fixed simulation
ridge_sim_fixed <- list(lambda = c(),
                        var = c(),
                        bias = c(),
                        mse = c())
for (lambda in lambdas) {
    print(lambda)
    x <- x_generator()
    mu <- phi(x, noise = FALSE)
    var_total <- c()
    bias_total <- c()
    mse_total <- c()
    for (i in 1:monte_carlo) {
        y_train <- phi(x)
        y_test <- phi(x)
        model <- train(x, y_train, ridge, lambda = lambda)
        results <- model(x)
        y_hat <- results$prediction
        x_std <- results$x_standardized
        W <- ridge(x_std, x_std, lambda)$coef
        s <- sigma^2
        var_total <- c(var_total,
                       tr(s * W %*% solve(crossprod(x_std), tol = 0) %*% t(W)))
        bias_total <- c(bias_total,
                        mean((y_hat - mu)^2))
        mse_total <- c(mse_total,
                       mean((y_test - y_hat)^2))
    }
    ridge_sim_fixed$lambda <- c(ridge_sim_fixed$lambda, lambda)
    ridge_sim_fixed$var    <- c(ridge_sim_fixed$var,    mean(var_total))
    ridge_sim_fixed$bias   <- c(ridge_sim_fixed$bias,   mean(bias_total))
    ridge_sim_fixed$mse    <- c(ridge_sim_fixed$mse,    mean(mse_total))
}
plot(ridge_sim_fixed$lambda,
     abs(ridge_sim_fixed$mse - (ridge_sim_fixed$bias + ridge_sim_fixed$var)),
     type = "h")



## RIDGE random simulation
ridge_sim_random <- list(lambda = c(),
                        var = c(),
                        bias = c(),
                        mse = c())
for (lambda in lambdas) {
    print(lambda)
    var_total <- c()
    bias_total <- c()
    mse_total <- c()
    for (i in 1:monte_carlo) {
        x_train <- x_generator()
        y_train <- phi(x_train)
        x_test <- x_generator()
        y_test <- phi(x_test)
        mu <- phi(x_test, noise = FALSE)
        model <- train(x_train, y_train, ridge, lambda = lambda)
        results <- model(x_test)
        y_hat <- results$prediction
        x_std <- results$x_standardized
        W <- ridge(x_std, x_std, lambda)$coef
        s <- sigma^2
        var_total <- c(var_total,
                       tr(s * W %*% solve(crossprod(x_std), tol = 0) %*% t(W)))
        bias_total <- c(bias_total,
                        mean((y_hat - mu)^2))
        mse_total <- c(mse_total,
                       mean((y_test - y_hat)^2))
    }
    ridge_sim_random$lambda <- c(ridge_sim_random$lambda, lambda)
    ridge_sim_random$var    <- c(ridge_sim_random$var,    mean(var_total))
    ridge_sim_random$bias   <- c(ridge_sim_random$bias,   mean(bias_total))
    ridge_sim_random$mse    <- c(ridge_sim_random$mse,    mean(mse_total))
}
plot(ridge_sim_random$lambda,
     abs(ridge_sim_random$mse - (ridge_sim_random$bias + ridge_sim_random$var)),
      type = "h")
