library(matrixcalc)
library(MASS)
library(MESS)
library(glmnet)
library(glasso)
library(caret)
library(ggplot2)
library(stringr)
library(gridExtra)

generate_b <- function(p, prob=0.1) {
  
  B <- matrix(rep(0,p*p), nrow=p)
  sample <- 0.5 * rbinom(size=1, n=sum(upper.tri(B, diag=F)), p=prob)
  
  B[lower.tri(B)] <- sample
  B <- t(B)
  B[lower.tri(B)] <- sample
  
  return(B)
}

generate_theta <- function(B, delta=0, standardized=T, return_full=F) {
  #
  # Generates theta as described in the project instructions.
  #  - chooses delta to be lowest satisfactory increment of 0.1 above 0
  #  - takes optional parameter to denote whether output should be standardized
  #  - prints chosen lambda if desired
  #
  while(T) {
    theta <- B + delta * diag(nrow(B))
    
    if(is.positive.definite(theta)) break
    
    delta <- delta + 0.1
  }
  if(standardized) theta <- theta / theta[1][1]
  if(return_full) return(list(theta=theta, delta=delta))
  return(theta)
}

nwlasso <- function(X_train, X_test=NULL, lambdalist) {
  #
  # Implements a Node-wise LASSO approach to graphical estimation.
  # 
  # Accepts:
  #  - X_train: the input training data, (n x p)-array 
  #  - X_test: optional test data for parameter tuning, if this is not supplied then no error
  #            metric is output
  #  - lambdalist: the parameter search space
  #
  # Returns on object 'results' with attributes:
  #  - $E1: a list of the E1 estimates as described in the project instructions, each estimate
  #         is a matrix with non-zero entries where nodes i & j are connected by an edge
  #  - $E2: a list of the E2 estimates 
  #  - $error (conditional on supply of X_test): a list containing an error metric for each value
  #         in the parameter search space (the metric is descibred in the report)
  #  - $lambda: lambdalist as supplied
  #
  p <- ncol(X_train)
  
  results <- list()
  # iterate through the parameter search space
  for(i in 1:length(lambdalist)) {
    
    # instantiate sum_mse to store the error metric and an empty C hat matrix
    sum_mse <- 0
    Chat <- matrix(rep(0, p^2), nrow=p)
    
    for(j in 1:p) {
      # fit a LASSO on the covariates with variable j as the target, ensure no intercept!
      fit_j <- glmnet(X_train[,-j], X_train[,j], lambda=lambdalist[i], intercept=F, standardize=F)
      
      # if we want to output the error metric add the test MSE for 'regression j' to sum_mse
      if(!is.null(X_test)) sum_mse <- sum_mse + assess.glmnet(fit_j, X_test[,-j], X_test[,j])$mse
      
      # update C hat values to 1 for all non-zero coefficient estimates 
      Chat[j,] <- ifelse(append(as.vector(coef(fit_j))[-1], 0, after=j-1)==0, 0, 1)
    }
    
    # store our edge estimates in matrices E1 & E2 using definitions in instructions
    E1 <- E2 <- matrix(rep(1, p^2), nrow=p) - diag(p)
    E1[(Chat + t(Chat)==0)] <- 0
    E2[!(Chat + t(Chat)==2)] <- 0
    
    results$E1[[i]] <- E1
    results$E2[[i]] <- E2
    if(!is.null(X_test)) results$error[[i]] <- sum_mse
  }
  results$lambda <- lambdalist
  return(results)
}

glasso_ <- function(X_train, X_test=NULL, lambdalist){
  #
  # Implements a graphical LASSO approach to graphical estimation using the library glasso.
  # 
  # Accepts:
  #  - X_train: the input training data, (n x p)-array 
  #  - X_test: optional test data for parameter tuning, if this is not supplied then no error
  #            metric is output
  #  - lambdalist: the parameter search space
  #
  # Returns on object 'results' with attributes:
  #  - $E: a list of the estimates, as described in the project instructions, in matrix form
  #        where non-zero entries indicate nodes i & j are connected by an edge
  #  - $loglikelihood (conditional on supply of X_test): a list of the loglikelihood of the 
  #        test data given the fitted model 
  #
  p <- ncol(X_train)
  
  # glassopath doesn't terminate for singular S so we have to use glasso instead
  theta_hat <- glassopath(cov(X_train), lambdalist, trace=0)$wi
  results <- list()
  
  # if test data is supplied compute the loglikelihood and store in results
  if(!is.null(X_test)) {
    logllh <- apply(theta_hat, 3, function(P_train) {
      logdet <- determinant(P_train, logarithm=T)$modulus
      logllh <- logdet - sum(cov(X_test) * P_train) 
      # note this is equivalent to logdet - sum(diag(cov(X_test) %*% P_train))
    })
    results$loglikelihood <- logllh
  }
  
  # iterate through theta_hat (estimated precision matrices) and store edges in matrix form 
  for(i in 1:length(lambdalist)){
    Ehat <- matrix(rep(1, p^2), nrow=p) - diag(p)
    Ehat[!(theta_hat[,,i]!=0)] <- 0
    
    results$E[[i]] <- Ehat
  }
  return(results)
}

#################################################################################################
#                                                                                               #
#                                       ROC Curve Plotting                                      #
#                                                                                               #
#################################################################################################

plot_ROC <- function(n=100, p=10, delta=0, lambdalist=seq(0, 1, length.out=101), fig_no='') {
  #n=50; p=50; delta=0; lambdalist=seq(0, 1, length.out=1001); fig_no='1.2'
  #
  # Plots ROC curves for the 3 methods described in the instructions.
  #
  # Accepts:
  #  - n: the number of simulated points to generate.
  #  - p: the number of nodes in the simulated graph.
  #  - delta: the size of delta used in constructing theta
  #  - lambdalist: the lambda space over which to calculate the TPR and FPR. 
  #
  B <- generate_b(p, prob=0.25)
  theta <- generate_theta(B, delta=delta, standardized=T)
  X <- mvrnorm(n, mu=rep(0, p), Sigma=solve(theta))
  
  # store E1, E2 & E3 edge estimates in the list estimates before transforming into factor 
  # vectors of the upper triangular
  estimates <- list()
  nwl <- nwlasso(X, lambdalist=lambdalist)
  estimates$E1 <- nwl$E1
  estimates$E2 <- nwl$E2
  estimates$E3 <- glasso_(X, lambdalist=lambdalist)$E
  estimates <- lapply(estimates, function(x) 
    lapply(x, function(y) factor(y[upper.tri(y)], levels=c(0, 1), labels=c(0, 1))))
  
  
  # compute and store the ground truth factor vector
  E <- B*2
  E <- as.factor(E[upper.tri(E)])
  
  # use a confusion matrix to calculate TPR (Sensitivity) and FPR (1-Specificity) and store in ROC 
  results <- lapply(estimates, function(x) 
    lapply(x, function(y) confusionMatrix(y, E, positive='1')$byClass[1:2]))
  
  E1.tpr <- E2.tpr <- E3.tpr <- E1.fpr <- E2.fpr <- E3.fpr <- vector()
  for(i in 1:length(lambdalist)){
    E1.tpr[i] <- results$E1[[i]]['Sensitivity']
    E2.tpr[i] <- results$E2[[i]]['Sensitivity']
    E3.tpr[i] <- results$E3[[i]]['Sensitivity']
    E1.fpr[i] <- 1 - results$E1[[i]]['Specificity']
    E2.fpr[i] <- 1 - results$E2[[i]]['Specificity']
    E3.fpr[i] <- 1 - results$E3[[i]]['Specificity']
  }
  ROC <- cbind.data.frame(E1.tpr, E2.tpr, E3.tpr, E1.fpr, E2.fpr, E3.fpr)
  
  # plot the curves
  plot2 <- ggplot(data=ROC) + 
    geom_step(aes(E1.fpr, E1.tpr, color='blue')) +
    geom_step(aes(E2.fpr, E2.tpr, color='green')) +
    geom_step(aes(E3.fpr, E3.tpr, color='red')) +
    theme_classic() + xlab('FPR') + ylab('TPR') +
    theme(legend.position='none', 
          axis.title=element_text(size=15,face='bold'), 
          plot.title=element_text(size=15, face='bold', hjust=0.5)) + 
    annotate('text', 
             x=c(.75, rep(c(.67, .78), 3)), 
             y=c(.48, .43, .43, .38, .38, .33, .33), 
             label=c('   AUROC', 
                     'E1 ', round(auc(E1.fpr, E1.tpr), 4), 
                     'E2 ', round(auc(E2.fpr, E2.tpr), 4),
                     'E3 ', round(auc(E3.fpr, E3.tpr), 4)), 
             colour=c('black', 'black', 'blue', 'black', 'green', 'black', 'red'), 
             fontface=c(2, rep(c(2, 1), 3)), size=rep(5, 7)) +
    annotate('rect', xmin=.62, xmax=.88, ymin=.3, ymax=.52, alpha=.05, color='black') +
    ggtitle(paste('figure ', fig_no, ': n=', n, ', p=', p, ', prob=', 0.25, sep=''))
}

# this is very slow: don't run it unless you want the plots!
plot1.1 <- plot_ROC(n=50, p=25, delta=0, lambdalist=seq(0.001, 1, length.out=1000), fig_no=1.1)
plot1.3 <- plot_ROC(n=50, p=75, delta=0, lambdalist=seq(0.001, 1, length.out=1000), fig_no=1.3)
plot1.2 <- plot_ROC(n=50, p=50, delta=0, lambdalist=seq(0.001, 1, length.out=1000), fig_no=1.2)
plot1.4 <- plot_ROC(n=100, p=25, delta=0, lambdalist=seq(0.001, 1, length.out=1000), fig_no=1.4)
plot1.5 <- plot_ROC(n=100, p=50, delta=0, lambdalist=seq(0.001, 1, length.out=1000), fig_no=1.5)
plot1.6 <- plot_ROC(n=100, p=75, delta=0, lambdalist=seq(0.001, 1, length.out=1000), fig_no=1.6)

grid.arrange(plot1.1, plot1.2, plot1.3, plot1.4, plot1.5, plot1.6, ncol=3)

#################################################################################################
#                                                                                               #
#                           Parameter Tuning: plotting (1)                                      #
#                                                                                               #
#################################################################################################

# Let's take it easy with the number of nodes and simulations because it's quite 
# computationally heavy. Generate the simulation parameters and simulate some data.
# We will also assign fig_no as the figure number for plotting. 
p<-50; n<-100; fig_no<-'2.2'

# create an empty array in which to store the test metrics
results <- array(rep(0, 50*3*5), dim=c(50, 3, 5), 
                 dimnames=list(1:50, c('E1', 'E2', 'E3'), c('MMCE', 'TPR', 'FPR', 'Precision', 'F1')))

cat('Optimum lambda\n    GLASSO NWLASSO')
for(i in 1:50){
  B <- generate_b(p)
  theta <- generate_theta(B, standardized=T)
  X <- mvrnorm(n, mu=rep(0, p), Sigma=solve(theta))
  
  # we will do k-fold cross-validation but for ease of computation let's set k=3
  k <- 3 
  folds <- sample(rep(1:k, length=n))
  
  # define the tuning parameter search space; we don't include 0 as it can result in convergence
  # issues if S is singular (i.e. p<n)
  lambdalist <- seq(0.01, 1, length.out=100) 
  
  # compute glasso loglikelihoods for each lambda
  glasso.llh <- sapply(1:k, function(ki) {
    fit_ki <- glasso_(X[which(folds!=ki),], X[which(folds==ki),], lambdalist=lambdalist)
    fit_ki$loglikelihood
  })
  
  # compute nwlasso loglikelihoods for each lambda
  nwlasso.mse <- sapply(1:k, function(ki) {
    fit_ki <- nwlasso(X[which(folds!=ki),], X[which(folds==ki),], lambdalist=lambdalist)
    fit_ki$error
  })
  
  # print the results to the console (helps keep track of progress too)
  cat('\n', str_pad(i, 4, 'right'), 
      str_pad(glasso.lambda <- lambdalist[which.max(rowMeans(glasso.llh))], 7),
      nwlasso.lambda <- lambdalist[which.min(rowMeans(nwlasso.mse))])
  
  # fit the models on the training data again
  glasso.cv <- glasso_(X, lambdalist=glasso.lambda)
  nwlasso.cv <- nwlasso(X, lambdalist=nwlasso.lambda)
  
  # define E as the ground truth and store as a factor vector
  E <- B*2
  E <- as.factor(E[upper.tri(E)])
  
  # define a list of the estimates as factor vectors
  estimates <- lapply(list(E1=nwlasso.cv$E1, E2=nwlasso.cv$E2, E3=glasso.cv$E), function(x) 
    lapply(x, function(y) factor(y[upper.tri(y)], levels=c(0, 1), labels=c(0, 1))))
  estimates <- lapply(estimates, unlist)
  
  # update the test error metrics in results
  confusion <- unlist(lapply(estimates, function(x) confusionMatrix(x, E, positive='1')$byClass))
  results[i,,'TPR'] <- confusion[c('E1.Sensitivity', 'E2.Sensitivity', 'E3.Sensitivity')]
  results[i,,'FPR'] <- 1 - confusion[c('E1.Specificity', 'E2.Specificity', 'E3.Specificity')]
  results[i,,'Precision'] <- confusion[c('E1.Precision', 'E2.Precision', 'E3.Precision')]
  results[i,,'F1'] <- confusion[c('E1.F1', 'E2.F1', 'E3.F1')]
  results[i,,'MMCE'] <- unlist(lapply(estimates, function(x) mean(x!=E)))
}

# plot the test errors by method
ggplot(data=stack(as.data.frame(results[,,'MMCE'])), aes(x=ind, y=values)) +
  geom_boxplot(fill=c('blue', 'green', 'red'), alpha=0.3) + 
  theme_classic() + xlab('') + ylab('MMCE') +
  theme(axis.title=element_text(size=15, face='bold'), 
        plot.title=element_text(size=15, face='bold', hjust=0.5),
        axis.text.x=element_text(size=15, face='bold')) +
  ggtitle(paste('figure ', 2.1, ': n=', n, ', p=', p, sep=''))

# print the mean and sd to the console
cat('     MEAN  E1   E2   E3')
for(i in dimnames(results)[[3]]) {
  cat(str_pad(i, 10), round(colMeans(results[,,i]), 2), '\n')
}

cat('\n\n       SD  E1   E2   E3')
for(i in dimnames(results)[[3]]) {
  cat(str_pad(i, 10), round(apply(results[,,i], 2, sd), 2), '\n')
}

#################################################################################################
#                                                                                               #
#                           Parameter Tuning: plotting (2)                                      #
#                                                                                               #
#################################################################################################

# let's do something similar again but this time we will vary delta and only look at performance
# for model E3 to ease computation
results2 <- array(rep(0, 50*5*3), dim=c(50, 5, 3), 
                  dimnames=list(1:50, c('2', '4', '6', '8', '10'), c('MMCE', 'TPR', 'FPR')))

# lowest possible value of delta will vary simulation to simulation; store them for the average
low_delta <- as.vector(rep(0, 50))
cat('Optimum lambda\n    GLASSO NWLASSO')
for(i in 1:50){
  cat('\n', i)
  # choose our list of deltas, 2 will likely be below the threshold and will get increased
  for(delta in c(2, 4, 6, 8, 10)) {
    B <- generate_b(p, prob=0.5)
    
    if(delta==2) {
      temp <- generate_theta(B, delta=delta, standardized=T, return_full=T)
      theta <- temp$theta
      low_delta[i] <- temp$delta
    } else theta <- generate_theta(B, delta=delta, standardized=T)
    
    X <- mvrnorm(n, mu=rep(0, p), Sigma=solve(theta))
    
    # we will do k-fold cross-validation but for ease of computation let's set k=3
    k <- 3 
    folds <- sample(rep(1:k, length=n))
    
    # define the tuning parameter search space; we don't include 0 as it can result in convergence
    # issues if S is singular (i.e. p<n)
    lambdalist <- seq(0.01, 1, length.out=100) 
    
    # compute glasso loglikelihoods for each lambda
    glasso.llh <- sapply(1:k, function(ki) {
      fit_ki <- glasso_(X[which(folds!=ki),], X[which(folds==ki),], lambdalist=lambdalist)
      fit_ki$loglikelihood
    })
  
    # fit the models on the training data again
    glasso.cv <- glasso_(X, lambdalist=lambdalist[which.max(rowMeans(glasso.llh))])
    
    # define E as the ground truth and store as a factor vector
    E <- B*2
    E <- as.factor(E[upper.tri(E)])
    
    # define a list of the estimates as factor vectors
    E3 <- glasso.cv$E[[1]]
    E3 <- as.factor(E3[upper.tri(E3)])# levels=c(0, 1), labels=c(0, 1))
    
    results2[i,toString(delta),'MMCE'] <- mean(E3!=E)
    results2[i,toString(delta), 'TPR'] <- confusionMatrix(E3, E, positive='1')$byClass[1]
    results2[i,toString(delta), 'FPR'] <- 1- confusionMatrix(E3, E, positive='1')$byClass[2]
  }
}

colnames(results2)[1] <- round(mean(low_delta), 2)

# plot the test errors by delta
ggplot(data=stack(as.data.frame(results2[,,'MMCE'])), aes(x=ind, y=values)) +
  geom_boxplot() + 
  theme_classic() + xlab(expression(delta)) + ylab('MMCE') +
  theme(axis.title=element_text(size=15, face='bold'), 
        plot.title=element_text(size=15, face='bold', hjust=0.5),
        axis.text.x=element_text(size=15, face='bold')) +
  ggtitle(paste('figure ', 3.9, ': n=', n, ', p=', p, ', prob=', 0.5, sep='')) 
