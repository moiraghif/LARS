library(tidyverse)
library(dplyr)
library(ggcorrplot)
library(ggplot2)
set.seed(123)

# Per una maggior qualità dei grafici su Windows

trace(grDevices:::png, quote({
  if (missing(type) && missing(antialias)) {
    type <- "cairo-png"
    antialias <- "subpixel"
  }
}), print = FALSE)


# Funzioni prese da simulation.R

ridge <- function(x, y, lambda) {
  p <- ncol(x)
  beta <- solve(crossprod(x) + lambda * diag(p), tol = 0) %*% crossprod(x, y)
  y_hat <- x %*% beta
  list(coef = beta,
       prevision = y_hat)
}


lars <- function(X, y, tol = 1e-10) {

    least_squares <- function(x, y)
        solve(crossprod(x)) %*% crossprod(x, y)

    n <- nrow(X)
    p <- ncol(X)
    max_iter <- min(n - 1, p)

    beta_tot <- matrix(0, ncol = max_iter, nrow = p)
    rownames(beta_tot) <- colnames(X)

    mu <- matrix(rep(0, n))            

    for (i in 1:max_iter) {
                                        
        c_hat <- crossprod(X, y - mu)   

                                        
        C <- as.double(max(abs(c_hat)))  
        active <- abs((abs(c_hat) - C)) <= tol  
        alpha <- sum(active)            

                                        
        s <- as.vector(sign(c_hat))     
      
        Xa <- (X %*% diag(s))[, active]  
      
        Ga <- solve(crossprod(Xa))      
        ones <- matrix(rep(1, alpha))  
        A <- as.double(crossprod(ones, Ga) %*% ones) ^-0.5
                                       
        w <- A * Ga %*% ones            
        u <- Xa %*% w                  
      
        a <- crossprod(X, u)            
      
        gamma <- Inf                    
      
        for (j in 1:p) { 
            cj <- c_hat[j, 1]
            aj <- a[j, 1]
            Ac <- c((C - cj) / (A - aj),
            (C + cj) / (A + aj))
            for (new_gamma in Ac) {
                if (!is.nan(new_gamma) & gamma > new_gamma & new_gamma > 0)
                    gamma <- new_gamma
            }
        }

                                        
        mu <- mu + gamma * u
        beta_tot[active, i] <- least_squares(Xa, mu)
    }

    
    list(coef = beta_tot[, max_iter, drop = FALSE],
         prevision = mu,
         log = beta_tot)
}


train <- function(x_train, y_train, method, force_independent = FALSE, ...) {
  ## execution example: mod <- train(x_train, y_train, lars)
  ## mod(x_test)$prediction >>> y_hat

  make_independent <- function(x) {
    eg <- eigen(var(x))
    M <- t(tcrossprod(eg$vectors, diag(eg$values, nrow = length(eg$values))))
    M_inv <- solve(M, tol = 0)

    function(x_new, reverse)
      if (reverse) x_new %*% M
      else x_new %*% M_inv
  }

  standardizer <- function(x, ind = FALSE) {
    mu <- colMeans(x)
    sigma <- as.double(sqrt(diag(var(x))))
    mk_ind <- ifelse(ind, make_independent(x),
                          function(x, reverse) x)

    function(y, reverse = FALSE) {
      if (reverse) {
        I <- diag(length(sigma)) * sigma
        mk_ind(t(t(y %*% I) + mu), reverse)
      } else {
        I <- diag(length(sigma)) / sigma
        mk_ind(t(t(y) - mu) %*% I, reverse)
      }
    }
  }

  x_standardizer <- standardizer(x_train, force_independent)
  y_standardizer <- standardizer(y_train, FALSE)

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
                     
make_independent <- function(x) {
    eg <- eigen(var(x))
    M <- t(tcrossprod(eg$vectors, diag(eg$values, nrow = length(eg$values))))
    M_inv <- solve(M, tol = 0)

    function(x_new, reverse)
      if (reverse) x_new %*% M
      else x_new %*% M_inv
  }
                     
standardizer <- function(x, ind = FALSE) {
    mu <- colMeans(x)
    sigma <- as.double(sqrt(diag(var(x))))
    mk_ind <- ifelse(ind, make_independent(x),
                          function(x, reverse) x)

    function(y, reverse = FALSE) {
      if (reverse) {
        I <- diag(length(sigma)) * sigma
        mk_ind(t(t(y %*% I) + mu), reverse)
      } else {
        I <- diag(length(sigma)) / sigma
        mk_ind(t(t(y) - mu) %*% I, reverse)
      }
    }
  }
                     
###################################################################
                     

# Importazione dataset
real2 <- read_csv("dati.csv") 
real2 <- real2[,-1]                     

# Corrplot
corr <- cor(real2)
corr[lower.tri(corr,diag=TRUE)] <- NA 
corr[corr == 1] <- NA 
corr <- as.data.frame(as.table(corr))
corr <- na.omit(corr) 
corr2 <- corr[corr$Var1=='y' ,] 
corr2 <- subset(corr2, abs(Freq) > 0.3) 
corr3 <- subset(corr, abs(Freq) > 0.98) # la soglia è alta per limitare i risultati visibil nel grafico
corrp <- rbind(corr2, corr3)
mtx_corrp <- t(reshape2::acast(corrp, Var1~Var2, value.var="Freq"))
ggcorrplot(mtx_corrp, lab = TRUE, ggtheme = ggplot2::theme_gray)
                     
                     
real2 <- real2[real2$delta=='0',] # Si tengono solo i non censurati
real2 <- real2[,-2]

sum(is.na(real2)) #  Non ci sono missing
                    
              

# Suddivisione Train-Test
size_tr <- floor(0.67 * nrow(real2)) 
train_ind <- sample(seq_len(nrow(real2)), size = size_tr)
train_sample <- real2[train_ind, ]
test_sample <- real2[-train_ind, ]

x_train <- as.matrix(train_sample[,-1])
y_train <- as.matrix(train_sample[,1])
x_test <- as.matrix(test_sample[,-1])
y_test <- as.matrix(test_sample[,1])
                     
                     
# RIDGE  
n_test <- nrow(y_test)
p <- ncol(x_test)
num <- n_test - 1
den <- n_test - p - 1
TSS <- sum((y_test-mean(y_test))^2)
ridge_mse <- c()
ridge_r2adj <- c()
lambdas <- seq(0,15000,100)

for (lambda in lambdas) {
  mod_r <- train(x_train, y_train, ridge, lambda=lambda)
  y_hat_r <- mod_r(x_test)$prediction
  mse_r <- mean((y_test-y_hat_r)^2)
  ridge_mse <- c(ridge_mse, mse_r)
  RSS_r <- sum((y_test - y_hat_r)^2)
  R2Adj_r <- 1-((num/den)*(RSS_r/TSS)) 
  ridge_r2adj <- c(ridge_r2adj, R2Adj_r)
}

min(ridge_mse)
                     
ggplot() + geom_line( aes(x = lambdas, y = ridge_mse), color='red', lwd=1.2)  +  
theme_minimal() +
theme(axis.line = element_line(colour = "black", size = 1)) +
xlab("Lambda") + ylab("MSE") + xlim(0,15000) + ylim(1.9,3) + geom_vline(xintercept = 7800)


# Lambda ottimo tramite cross validation
library(glmnet)
fit.ridge = glmnet(as.matrix(select(real2, -y)),as.matrix(select(real2, y)),family="gaussian",alpha=0)
plot(fit.ridge, xvar = "lambda")

K <- 10
ridge.cv<-cv.glmnet(as.matrix(select(real2, -y)),as.matrix(select(real2, y)),alpha=0, nfolds = K, grouped=FALSE)
plot(ridge.cv) # la seconda barra corrisponde al minimo errore di cross validation + una volta il suo standard error
hatlambda <-ridge.cv$lambda.1se
hatlambda                     
                     
n <- nrow(real2)
folds <- sample( rep(1:K,length=n) )
KCV_r <- vector()

for (k in 1:K) {
  testIndexes <- which(folds==k,arr.ind=TRUE)
  testData <- real2[testIndexes, ]
  trainData <- real2[-testIndexes, ]
  mod_rcv <- train(as.matrix(trainData[-1]), as.matrix(trainData[1]), ridge, lambda = hatlambda)
  y_hat_rcv <- mod_rcv(as.matrix(testData[-1]))$prediction
  KCV_r[k] <- mean((as.matrix(testData[1])-y_hat_rcv)^2)
}
mean(KCV_r) 
                     
                     
                     
                     
# LARS
 
# Train-Test                     
mod_lars <- train(x_train, y_train, lars)
y_hat_l <- mod_lars(x_test)$prediction
mse_lars <- mean((y_test-y_hat_l)^2)
mse_lars
                     

x_standardizer <- standardizer(x_train)
y_standardizer <- standardizer(y_train)
mod_lars2 <- lars(x_standardizer(x_train), y_standardizer(y_train))
mse_lars2 <- c()
beta_compl <- mod_lars2$log
beta_sing <- mod_lars2$coef
x_test_stand <- x_standardizer(x_test)
I <- ncol(beta_compl)

for (i in 1:I) {
  yhat_l2 <- x_test_stand %*% beta_compl[,i]
  yhat_l2 <- y_standardizer(yhat_l2, reverse = TRUE)
  mse_lars2[i] <- mean((y_test-yhat_l2)^2)
}

mse_lars2 # come varia l'errore ad ogni iterazione
                     

K <- 5
folds <- sample( rep(1:K,length=n) )
KCV <- matrix(NA,K,100)

for (k in 1:K) {
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

i1 <- which.min(KCVmean) # modello minCV
i2 <- min(which(KCVmean<=KCVmean[i1]+KCVsd[i1]/sqrt(K))) # modello min one standard error rule



ggplot() + geom_point( aes(x = 1:50, y = KCVmean), color='red', lwd=2.3) +
  geom_line( aes(x = 1:50, y = KCVmean), color='red', lwd=0.2) +
  geom_line( aes(x = 1:50, y = KCVmean+KCVsd/sqrt(K)), color='black', lwd=1, linetype = "dashed") +
  geom_line( aes(x = 1:50, y = KCVmean-KCVsd/sqrt(K)), color='black', lwd=1, linetype = "dashed") +
  theme_minimal()   +
  xlab("Iteration") + ylab("MSE") +
  theme(axis.line = element_line(colour = "black", size = 1)) +
  geom_vline(xintercept = i2)         
                     
                     
                     
# Coefficient path
beta_sin <- order(abs(beta_sing), decreasing = TRUE)
beta_sin <- beta_sin[1:10]
beta_compl <- beta_compl[beta_sin,]

beta_compl[beta_compl == 0] <- NA


selected_variables <- rowSums(is.na(beta_compl)) < ncol(beta_compl) # se il numero di na è maggiore del numero di colonne
beta_compl <- beta_compl[selected_variables, ]
#rownames(betas) <- colnames(x_train)[selected_variables]


## rigira i dati in un formato facilmente interpretabile da ggplot
data_plot <- tibble(
    variable = numeric(),
    value = numeric(),
    iteration = numeric())

p <- ncol(beta_compl)
for (i in seq_len(nrow(beta_compl))) {
    data_plot <- data_plot %>%
        add_row(variable = i,
                value = beta_compl[i, ],
                iteration = 1:p)
}

data_plot$variable <- as.character(data_plot$variable)

## plot dei risultati
(ggplot(data_plot, aes(x = iteration, y = value, colour = variable)) +
    geom_line() + geom_hline(yintercept = 0) +
    theme_minimal() + theme(legend.position = "none") +
    xlab("Iteration") + ylab("Coefficient") +
    ggtitle("LARS coefficient"))                     
                     
                     
