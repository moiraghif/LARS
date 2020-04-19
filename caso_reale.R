rm(list = ls())

real2 <- read_csv("dati.csv")

real2 <- real2[real2$delta=='0',] # tengo solo i non censurati
real2 <- real2[,-1]
real2 <- real2[,-2]

sum(is.na(real2)) #  non ci sono missing




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
    mu <- matrix(rep(0, n))             # n ÃÂ 1

    for (i in 1:max_iter) {
                                        # equation 2.1
        ## each c_j represent the correlation of the j-th variable
        ## between X and the projection on the sub-space L
        c_hat <- crossprod(X, y - mu)   # vector, p ÃÂ 1

                                        # equation 2.9
        ## the "active" subset includes each variable that is as much
        ## correlate with y as the maximum-correlated x_j
        C <- as.double(max(abs(c_hat)))  # scalar
        active <- abs((abs(c_hat) - C)) <= tol  # due to R approximation
        alpha <- sum(active)            # scalar, value of a

                                        # equation 2.10
        ## a vector of signs of correlation, for multiply the X matrix
        ## to use only positive correlations
        s <- as.vector(sign(c_hat))     # vector, p ÃÂ 1

                                        # equation 2.4
        ## the X matrix that includes only active variables (with
        ## positive correlation)
        Xa <- (X %*% diag(s))[, active]  # matrix, n ÃÂ a

                                        # equation 2.5
                                        # (inverse is computed for performance)
        ## this part is quite complicated, see Paper for details
        Ga <- solve(crossprod(Xa))      # matrix, a ÃÂ a
        ones <- matrix(rep(1, alpha))   # vector, a ÃÂ 1
        A <- as.double(crossprod(ones, Ga) %*% ones) ^-0.5
                                        # scalar

                                        # equation 2.6
        ## u is the direction of the update
        w <- A * Ga %*% ones            # vector a ÃÂ 1
        u <- Xa %*% w                   # vector n ÃÂ 1

                                        # equation 2.11
        ## this part is not well described in the Paper
        a <- crossprod(X, u)            # vector p ÃÂ 1

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

size_tr <- floor(0.67 * nrow(real2)) #  caso p > n
set.seed(123)

train_ind <- sample(seq_len(nrow(real2)), size = size_tr)
train_sample <- real2[train_ind, ]
test_sample <- real2[-train_ind, ]

x_train <- as.matrix(train_sample[,-1])
y_train <- as.matrix(train_sample[,1])
x_test <- as.matrix(test_sample[,-1])
y_test <- as.matrix(test_sample[,1])

# MSE
mod_lars <- train(x_train, y_train, lars)
y_hat_l <- mod_lars(x_test)$prediction
mse.lars <- mean((y_test-y_hat_l)^2)
mse.lars

#R2ADJ
n_test <- nrow(y_test)
p <- ncol(x_test)
num <- n_test - 1
den <- n_test - p - 1
RSS_l <- sum((y_test - y_hat_l)^2)
TSS <- sum((y_test-mean(y_test))^2)

R2Adj_l <- 1-((num/den)*(RSS_l/TSS)) 
# 1-((RSS/den)/(TSS/num)) # forumula prof
R2Adj_l



# Cross Validation Lars
n <- nrow(real2)
K = 10
folds <- sample( rep(1:K,length=n) )
KCV <- vector()

for (k in 1:K){
  testIndexes <- which(folds==k,arr.ind=TRUE)
  testData <- real2[testIndexes, ]
  trainData <- real2[-testIndexes, ]
  mod_lcv <- train(as.matrix(trainData[-1]), as.matrix(trainData[1]), lars)
  y_hat_lcv <- mod_lcv(as.matrix(testData[-1]))$prediction
  KCV[k] <- mean((as.matrix(testData[1])-y_hat_lcv)^2)
}

mean(KCV) # mse




### Cross validation + feature selection LARS

n <- nrow(real2)
K = 10
folds <- sample( rep(1:K,length=n) )
KCV <- matrix(NA,K,100)


for (k in 1:K){
  testIndexes <- which(folds==k,arr.ind=TRUE)
  testData <- real2[testIndexes, ]
  trainData <- real2[-testIndexes, ]
  x_train <- as.matrix(trainData[,-1])
  y_train <- as.matrix(trainData[,1])
  x_test <- as.matrix(testData[,-1])
  y_test <- as.matrix(testData[,1])
  x_standardizer <- standardizer(x_train)
  y_standardizer <- standardizer(y_train)
  mod_lcv <- lars(x_standardizer(x_train), y_standardizer(y_train))
  beta_compl <- mod_lcv$log
  beta_sing <- mod_lcv$coef
  x_test_stand <- x_standardizer(x_test)
  I = ncol(beta_compl)
  
  for (i in 1:I) {
    yhat_k <- x_test_stand %*% beta_compl[,i]
    yhat_k <- y_standardizer(yhat_k, reverse = TRUE)
    KCV[k,i] <- mean((y_test-yhat_k)^2)
  }
  
  
}

KCVmean<-apply(KCV,2,mean)
KCVsd<-apply(KCV,2,sd)

KCVmean <-na.omit(KCVmean)
KCVsd <- na.omit(KCVsd)

plot(0:55,KCVmean, type="b",xlab="p", ylab="CV Error", pch=19,cex=.7)
lines(0:55,KCVmean+KCVsd/sqrt(K),lty=2,col="darkred")
lines(0:55,KCVmean-KCVsd/sqrt(K),lty=2,col="darkred")


#################




ridge_mse <- c()
ridge_r2adj <- c()
lambdas <- seq(1,2000,50)

for (lambda in lambdas) {
  mod_r <- train(x_train, y_train, ridge, lambda=lambda)
  y_hat_r <- mod_r(x_test)$prediction
  mse_r <- mean((y_test-y_hat_r)^2)
  ridge_mse <- c(ridge_mse, mse_r)
  RSS_r <- sum((y_test - y_hat_r)^2)
  R2Adj_r <- 1-((num/den)*(RSS_r/TSS)) 
  ridge_r2adj <- c(ridge_r2adj, R2Adj_r)
}


plot(lambdas,ridge_mse,type="l",ylim=c(4,4.8),xlim=c(0,2000),
     xlab="Lambda",ylab="MSE",col="black")
abline(h=mse.lars,lty=2)
legend("topleft", lty=c(1),
       legend=c("Ridge Regression"))


ggplot() + geom_line( aes(x = lambdas, y = ridge_mse), color='red', lwd=1.2) + theme_minimal() +
xlab("Lambda") + ylab("MSE") + xlim(0,10000) + ylim(1.9,3)


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


mod_r <- train(x_train, y_train, ridge, lambda=hatlambda)
y_hat_r <- mod_r(x_test)$prediction
mse_r <- mean((y_test-y_hat_r)^2)
mse_r

################# Grafici #################

# per vederli meglio su windows
trace(grDevices:::png, quote({
  if (missing(type) && missing(antialias)) {
    type <- "cairo-png"
    antialias <- "subpixel"
  }
}), print = FALSE)


##
## Corrplot
##

# da valutare, le variabili sono moltissime, qua ci sono quelle con la correlazione più alta con la y e quelle più correlate tra loro


corr <- cor(real2)

corr[lower.tri(corr,diag=TRUE)] <- NA 
corr[corr == 1] <- NA 
corr <- as.data.frame(as.table(corr))
corr <- na.omit(corr) 
corr2 <- corr[corr$Var1=='y' ,] 
corr2 <- subset(corr2, abs(Freq) > 0.3) 
corr3 <- subset(corr, abs(Freq) > 0.98) # la soglia Ã¨ cosÃ¬ alta per limitare il numero di variabili nel grafico
corrp <- rbind(corr2, corr3)
mtx_corrp <- t(reshape2::acast(corrp, Var1~Var2, value.var="Freq"))
ggcorrplot(mtx_corrp, lab = TRUE, ggtheme = ggplot2::theme_gray)


##
## Osservati vs predetti
##

db <- data.frame(seq(1:nrow(y_test)),cbind(y_test, y_hat_l))
names(db) <- c("index","y_test", "y_hat")

ggplot() + 
  geom_line(data = db, aes(x = index, y = y_hat_l, color = "#00AFBB"), lwd=1.1) +
  geom_line(data = db, aes(x = index, y = y_test, color = "#FC4E07"), lwd=1.1) +
  xlab("Index") +
  ylab('Values') +
  theme_minimal() +
  labs(title = "Observed values vs Predicted values") +
  scale_color_discrete(name = "Values", labels = c("Predicted", "Observed"))



db <- data.frame(seq(1:nrow(y_test)),cbind(y_test, y_hat_r))
names(db) <- c("index","y_test", "y_hat")

ggplot() + 
  geom_line(data = db, aes(x = index, y = y_hat_r, color = "#00AFBB"), lwd=1.1) +
  geom_line(data = db, aes(x = index, y = y_test, color = "#FC4E07"), lwd=1.1) +
  xlab("Index") +
  ylab('Values') +
  theme_minimal() +
  labs(title = "Observed values vs Predicted values") +
  scale_color_discrete(name = "Values", labels = c("Predicted", "Observed"))






##
## grafici LARS
##


## presa da train()
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


## usa tutti i dati
x_train <- as.matrix(select(real2, -y))
y_train <- as.matrix(select(real2,  y))

x_standardizer <- standardizer(x_train)
y_standardizer <- standardizer(y_train)

mod <- lars(x_standardizer(x_train),
            y_standardizer(y_train))

## rimuovi i coefficienti nulli
betas <- mod$log
betas[betas == 0] <- NA

selected_variables <- rowSums(is.na(betas)) < ncol(betas)
betas <- betas[selected_variables, ]
rownames(betas) <- colnames(x_train)[selected_variables]


## rigira i dati in un formato facilmente interpretabile da ggplot
data_plot <- tibble(
    variable = c(),
    value = c(),
    iteration = c())

p <- ncol(betas)
for (i in seq_len(nrow(betas))) {
    data_plot <- data_plot %>%
        add_row(variable = rownames(betas)[i],
                value = betas[i, ],
                iteration = 1:p)
}


## plot dei risultati
(ggplot(data_plot, aes(x = iteration, y = value, colour = variable)) +
    geom_line() + geom_hline(yintercept = 0) +
    theme_minimal() + theme(legend.position = "none") +
    xlab("Iteration") + ylab("coefficient") +
    ggtitle("LARS coefficient")) %>%

    ggsave(filename = "caso_reale_lars.png", device = png(),
           dpi = 500, width = 16, height = 9)
