rm(list = ls())


library(tidyverse)

real2 <- read_csv("dati.csv")
rownames(real2) <- real2$X
real2 <- real2[,-1]


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
    mu <- matrix(rep(0, n))             # n Ã 1

    for (i in 1:max_iter) {
                                        # equation 2.1
        ## each c_j represent the correlation of the j-th variable
        ## between X and the projection on the sub-space L
        c_hat <- crossprod(X, y - mu)   # vector, p Ã 1

                                        # equation 2.9
        ## the "active" subset includes each variable that is as much
        ## correlate with y as the maximum-correlated x_j
        C <- as.double(max(abs(c_hat)))  # scalar
        active <- abs((abs(c_hat) - C)) <= tol  # due to R approximation
        alpha <- sum(active)            # scalar, value of a

                                        # equation 2.10
        ## a vector of signs of correlation, for multiply the X matrix
        ## to use only positive correlations
        s <- as.vector(sign(c_hat))     # vector, p Ã 1

                                        # equation 2.4
        ## the X matrix that includes only active variables (with
        ## positive correlation)
        Xa <- (X %*% diag(s))[, active]  # matrix, n Ã a

                                        # equation 2.5
                                        # (inverse is computed for performance)
        ## this part is quite complicated, see Paper for details
        Ga <- solve(crossprod(Xa))      # matrix, a Ã a
        ones <- matrix(rep(1, alpha))   # vector, a Ã 1
        A <- as.double(crossprod(ones, Ga) %*% ones) ^-0.5
                                        # scalar

                                        # equation 2.6
        ## u is the direction of the update
        w <- A * Ga %*% ones            # vector a Ã 1
        u <- Xa %*% w                   # vector n Ã 1

                                        # equation 2.11
        ## this part is not well described in the Paper
        a <- crossprod(X, u)            # vector p Ã 1

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



size_tr <- floor(0.75 * nrow(real2)) #  caso p > n
set.seed(123)

train_ind <- sample(seq_len(nrow(real2)), size = size_tr)
train_sample <- real2[train_ind, ]
test_sample <- real2[-train_ind, ]

x_train <- as.matrix(train_sample[,-1])
y_train <- as.matrix(train_sample[,1])
x_test <- as.matrix(test_sample[,-1])
y_test <- as.matrix(test_sample[,1])


mod_lars <- train(x_train, y_train, lars)
y_hat <- mod_lars(x_test)$prediction
mse.lars <- mean((y_test-y_hat)^2)
mse.lars



ridge_mse <- c()
lambdas <- seq(1,2000 , 50)

for (lambda in lambdas) {
  mod <- train(x_train, y_train, ridge, lambda=lambda)
  y_hat <- mod(x_test)$prediction
  mse <- mean((y_test-y_hat)^2)
  ridge_mse <- c(ridge_mse, mse)
}


plot(lambdas,ridge_mse,type="l",ylim=c(4,4.8),xlim=c(0,2000),
     xlab="Lambda",ylab="MSE",col="black")
abline(h=mse.lars,lty=2)
legend("topleft", lty=c(1),
       legend=c("Ridge Regression"))


#######

#Prove per trovare lambda ottimo nel RIDGE con glmnet


library(glmnet)
fit.ridge = glmnet(x_train,y_train,family="gaussian",alpha=0)

K <- 10
ridge.cv<-cv.glmnet(x_train,y_train,alpha=0, nfolds = K, grouped=FALSE)
plot(ridge.cv) # la seconda barra corrisponde al minimo errore di cross validation + una volta il suo standard error
hatlambda <-ridge.cv$lambda.min
hatlambda

predict(fit.ridge, s=hatlambda, type = "coefficients")
yhat.ridge = predict(fit.ridge, s=hatlambda, newx=x_test, exact=T)
mse_ridge <- mean((y_test-yhat.ridge)^2)
mse_ridge

