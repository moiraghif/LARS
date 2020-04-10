rm(list = ls())

library(tidyverse)
library(gganimate)
library(tibble)


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


simulation_fixed <- function(data_generator, beta, sigma, method, iterations, verbose = TRUE, ...) {
  error <- function(y, y_hat) mean((y - y_hat)^2)
  phi <- get_function(beta, sigma)
  x <- data_generator()
  mu <- phi(x, noise = FALSE)
  n <- nrow(x)
  p <- ncol(x)
  
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
    y_hat <- model(x)
    y_hat_total[, i] <- y_hat
    err_total <- c(err_total, error(y_test, y_hat))
    bias_total <- c(bias_total, mean((rowMeans(y_hat_total,na.rm = T) - mu)^2)) 
    variance_total <- c(variance_total, mean(apply(y_hat_total, 1, var,na.rm=TRUE)))
    sum_total <- c(sum_total,bias_total[i]+variance_total[i])
    #Print completion time
    Sys.sleep(0.1)
    cat("\r","Iteration n.",i)
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
    print(paste0("       difference: ", round(abs(bias + variance - av_error), 4)))
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
       "log" = sum_total
  )
}

simulation_random <- function(data_generator, beta, sigma, method, iterations, verbose = TRUE, ...) {
  error <- function(y, y_hat) mean((y - y_hat)^2)
  phi <- get_function(beta, sigma)
  
  err_total <- c()
  bias_total <- c()
  variance_total <- c()
  sum_total <- c()
  y_hat_total <- matrix(NA, nrow = n, ncol = iterations)

  start_time <- Sys.time()
  for (i in 1:iterations) {
    x_train <- data_generator()
    y_train <- phi(x_train)
    n <- nrow(x_train)
    p <- ncol(x_train)
    mu <- phi(x_train, noise = FALSE)
    x_test <- data_generator()
    y_test <- phi(x_test)
    model <- train(x_train, y_train, method = method, ...)
    y_hat <- model(x_test)
    y_hat_total[, i] <- y_hat
    err_total <- c(err_total, error(y_test, y_hat))
    bias_total <- c(bias_total, mean((y_hat - mu))^2)
    variance_total <- c(variance_total, mean(apply(y_hat_total, 1, var, na.rm=TRUE)))
    sum_total <- c(sum_total,bias_total[i]+variance_total[i])
    #Print completion time
    cat("\r","Iteration n.",i)
  }
  stop_time <- Sys.time()
  av_error <- mean(err_total)
  variance <- variance_total[iterations]
  bias <- mean(bias_total)
  if (verbose) {
    print(paste0("Execution time: ", stop_time - start_time))
    print(paste0("Sigma²:   ", sigma^2))
    print(paste0("Variance: ", variance))
    print(paste0("Bias²:    ", bias))
    print(paste0("MSE (theoretical): ", bias + variance))
    print(paste0("MSE (simulation):  ", mean(err_total)))
    print(paste0("       difference: ", round(abs(bias + variance - av_error), 4)))
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
       "log" = sum_total
  )
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

## lars fixed simulation
lars_simulation_fixed <- simulation_fixed(x_generator, beta_true, sigma, lars, monte_carlo, verbose = FALSE)

## lars random simulation
lars_simulation_random <- simulation_random(x_generator, beta_true, sigma, lars, monte_carlo, verbose = FALSE)

###########################################################################
##Measures

#Fixed
hatsigma2 = (n*MSEs.tr)/(n-p)
Cps = MSEs.tr + (2*hatsigma2*p)/n 
opt.fixed = 2 * sigma^2 * p / n

#Random
fixed_to_random = (sigma^2 * p / n) * (p + 1) / (n - p - 1)
opt.random = opt.fixed*(2+(p+1)/(n - p -1))
err.random = err.train + opt.random
Cps.random = Cps.fixed + fixed_to_random
###########################################################################
##Plot

#Function to create db
create_db <- function(x) {
  db <- data.frame(cbind(1:length(x$err_log),x$err_log,x$log,x$variance_log,x$bias_log))
  names(db) <- c("Iter","Errors","Bias² + Variance","Variance","Bias²")
  return(db)
}

#db <- create_db(lars_simulation_fixed)
db <- create_db(lars_simulation_random)

#Arrow plot for Variance
if(is.na(db$Variance[1])) {db$Variance[1]<-0}
arrow <- c()
for(k in 2:monte_carlo) {
  if(db$Variance[k]<db$Variance[k-1]) {
    arrow <- c(arrow,"\u2193")
  }
  else if (db$Variance[k]==db$Variance[k-1]) {
    arrow <- c(arrow,"\u2192")
  }
  else {
    arrow <- c(arrow,"\u2191")
  }
}

#Arrow plot for Bias
arrow_bias <- c()
for(k in 2:monte_carlo) {
  if(db$"Bias²"[k]<db$"Bias²"[k-1]) {
    arrow_bias <- c(arrow_bias,"\u2193")
  }
  else if (db$"Bias²"[k]==db$"Bias²"[k-1]) {
    arrow_bias <- c(arrow_bias,"\u2192")
  }
  else {
    arrow_bias <- c(arrow_bias,"\u2191")
  }
}

##Convergence animated Plot
z <- seq(1,nrow(db))
cols <- c("Bias² + Variance"="#00AFBB","Estimates Fluctuation"="#FC4E07","Expected Value"="#C4961A")
conv <- ggplot() +
  geom_hline(yintercept = mean(db$Errors),linetype="dashed")+
  #Theoretical Error
  geom_point(data = db, aes(x = db$Iter, y = db$`Bias² + Variance`, colour="Bias² + Variance"), size=6) +
  geom_line(data = db, aes(x = db$Iter, y = db$`Bias² + Variance`, colour="Bias² + Variance"), size=0.1) +
  #Calculated Error
  geom_point(data = db, aes(x = db$Iter, y = db$Errors, colour="Estimates Fluctuation"), size=6) +
  geom_line(data = db, aes(x = db$Iter, y = db$Errors, colour="Estimates Fluctuation"), size=0.1) +
  labs(title = "Convergence to the Expected MSE (Fixed-X Setting)")+
  ylim(0,5000)+
  xlim(0,length(z))+
  xlab("Iteration")+
  ylab("Value") +
  #Final points (Averages)
  geom_point(aes(x=length(z),y=mean(db$Errors),colour="Expected Value"), size=6)+
  geom_text(data=data.frame(z=z),
            mapping = aes(x = 0.8*length(z), y = 4700,
                          label = paste0("ITERATION N.  ",z,"\n","E(Bias²):       ",as.integer(db$"Bias²")," ",arrow_bias,"\n",
                                         "E(Variance):    ",as.integer(db$Variance)," ",arrow)))+
  theme(plot.title = element_text(size = 20, face = "bold"),
        legend.title=element_text(size=15, face = "bold"), 
        legend.text=element_text(size=15),
        plot.caption = element_text(face = "italic", hjust = 0.5, size = 15))+
  theme_bw()+
  labs(colour="MSE Values")+
  transition_reveal(z)+
  ease_aes("linear")+
  enter_appear()
plot1 <- animate(conv, fps=10, end_pause = 10)
plot1
#Save
anim_save("/Users/riccardocervero/Desktop/Plot1.gif",plot1)

##E(Bias) and E(Variance)
z <- seq(1,nrow(db))
cols <- c("E(Bias²)"="#00AFBB","E(Variance)"="#C4961A")
BiasVariance <- ggplot() +
  #Bias
  geom_line(data = db, aes(x = db$Iter, y = db$`Bias²`, colour="E(Bias²)"), size=0.5) +
  #Variance
  geom_line(data = db, aes(x = db$Iter, y = db$Variance, colour="E(Variance)"), size=0.5) +
  labs(title = "Bias² and Variance Trend")+
  ylim(0,1.01*max(db))+
  xlim(0,length(z))+
  xlab("Iteration")+
  ylab("Value") +
  geom_text(data=data.frame(z=z),
            mapping = aes(x = 0.8*length(z), y = 1.01*max(db),
                          label = paste0("ITERATION N.  ",z)))+
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.title=element_text(size=20, face = "bold"), 
        legend.text=element_text(size=20),
        plot.caption = element_text(face = "italic", size = 15))+
  theme_bw()+
  labs(colour="Value")+
  transition_reveal(z)+
  ease_aes("linear")+
  enter_appear()
plot2 <- animate(BiasVariance, fps=10, end_pause = 10)
plot2
#Save
anim_save("/Users/riccardocervero/Desktop/Plot2.gif",plot2)

###########################################################################
monte_carlo <- 2
ridge_simulations.mse <- c()
ridge_simulations.lambda <- c()
ridge_simulations.beta <- c()
ridge_simulations.variance <- c()
ridge_simulations.bias <- c()
lambdas <- seq(0.5, 5, .5)
for (lambda in lambdas) {
  print(paste0("Lambda: ",lambda))
  sim <- simulation_fixed(x_generator, beta_true, sigma, ridge, monte_carlo, verbose=FALSE, lambda=lambda)
  ridge_simulations.mse <- c(ridge_simulations.mse, sim$"MSE")
  ridge_simulations.lambda <- c(ridge_simulations.lambda, lambda)
  ridge_simulations.bias <- c(ridge_simulations.bias, sim$"bias^2")
  ridge_simulations.variance <- c(ridge_simulations.variance, sim$variance)
  #ridge_simulations.beta <- c()
  break
}


##Ridge Random-simulation
monte_carlo <- 1e3
ridge_simulations.mse <- c()
ridge_simulations.lambda <- c()
ridge_simulations.beta <- c()
ridge_simulations.variance <- c()
ridge_simulations.bias <- c()
lambdas <- seq(0.5, 5, .5)
lars_simulation <- simulation_random(x_generator, beta_true, sigma, lars, monte_carlo, verbose=FALSE)
for (lambda in lambdas) {
  print(paste0("Lambda: ",lambda,"\n"))
  sim <- simulation_random(x_generator, beta_true, sigma, ridge, monte_carlo, verbose=FALSE, lambda=lambda)
  ridge_simulations.mse <- c(ridge_simulations.mse, sim$"MSE")
  ridge_simulations.lambda <- c(ridge_simulations.lambda, lambda)
  ridge_simulations.bias <- c(ridge_simulations.bias, sim$"bias^2")
  ridge_simulations.variance <- c(ridge_simulations.variance, sim$variance)
  #ridge_simulations.beta <- c()
}

#Plot
bias_variance <- tibble(
  value = c(ridge_simulations.bias,     lars_simulation$"bias^2",
            ridge_simulations.variance, lars_simulation$variance),
  index_name = rep(c("bias", "variance"), each = length(names)),
  model_name = factor(rep(names, times = 2), levels = names))


ggplot(bias_variance, aes(x = model_name, fill = index_name)) +
  geom_bar(aes(y = value), stat = "identity", position = "stack") +
  xlab("") + ylab("") + ggtitle("Reducible Error scomposition") +
  guides(fill = guide_legend(title = "")) +
  scale_y_log10()
ggsave("random.png", dpi = 500,
       width = 16, height = 9)


phi <- get_function(beta_true, sigma)
X <- x_generator()
y <- phi(X)

ridge <- function(x, y, lambda) {
  n <- nrow(x)
  p <- ncol(x)
  beta <- solve(crossprod(x) + lambda * diag(p)) %*% crossprod(x, y)
  #beta <- solve(t(X) %*% X + lambda*diag(p)) %*% t(X) %*% y
  y_hat <- x %*% beta
  list(coef = beta,
       prevision = y_hat)
}
lambdas = seq(0.5, 5, .5)
nlambda = length(lambdas)
hatbetas <- matrix(NA,ncol=p, nrow=nlambda)
Vars <- matrix(NA,ncol=p, nrow=nlambda)

for (l in 1:nlambda){
  lambda = lambdas[l]
  hatbetas[l,] = ridge(X, y, lambda)$coef  # solve(t(X) %*% X + lambda*diag(p)) %*% t(X) %*% y 
  W = ridge(X, X, lambda)$coef  # solve(diag(p) + lambda * solve(t(X) %*% X))
  Vars[l,] = diag( W %*% solve(crossprod(X)) %*% t(W) )
  }
matplot(lambdas, hatbetas, type="l")
