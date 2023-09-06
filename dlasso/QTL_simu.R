## cvb - variable selection simulation

#########################################################################################################
library(MASS)
library(glmnet)
library(ggplot2)
library(ncvreg)
library(SIS)
library(EMVS)
library(BMS)
library(tictoc)

library(iterators)
# library(foreach)
library(doParallel)
setwd("/home/youchong/projects/collapsed_vb/QTL/")
source('vs_func.R')
source('epBVS.R')

#########################################################################################################
##### read data

response <- read.table("phe_simulat.csv", header = FALSE, sep = ",")
covariates <- read.table("gen_simulat.csv", header = FALSE, sep = ",")

n <- nrow(response)
p <- ncol(covariates)
p
vbeta <- c()
X <- matrix(0, n, 7381)
count <- 1
for (i in 1:p) {
    for (j in i:p) {
        if (i == j) {
            X[, count] <- covariates[, i]
        } else {
            X[, count] <- covariates[, i] * covariates[, j]
        }
        
        vbeta[count] <- 0
        
        if ((i == 1) & (j == 1))     { vbeta[count] <- 4.47; }
        if ((i == 21) & (j == 21))   { vbeta[count] <- 3.16; }
        if ((i == 31) & (j == 31))   { vbeta[count] <- 2.24; }
        if ((i == 51) & (j == 51))   { vbeta[count] <- 1.58; }
        if ((i == 71) & (j == 71))   { vbeta[count] <- 1.58; }
        if ((i == 91) & (j == 91))   { vbeta[count] <- 1.10; }
        if ((i == 101) & (j == 101)) { vbeta[count] <- 1.10; }
        if ((i == 111) & (j == 111)) { vbeta[count] <- 0.77; }
        if ((i == 121) & (j == 121)) { vbeta[count] <- 0.77; }
        if ((i == 1) & (j == 11))    { vbeta[count] <- 1.00; }
        if ((i == 2) & (j == 119))   { vbeta[count] <- 3.87; }
        if ((i == 10) & (j == 91))   { vbeta[count] <- 1.30; }
        if ((i == 15) & (j == 75))   { vbeta[count] <- 1.73; }
        if ((i == 20) & (j == 46))   { vbeta[count] <- 1.00; }
        if ((i == 21) & (j == 22))   { vbeta[count] <- 1.00; }
        if ((i == 26) & (j == 91))   { vbeta[count] <- 1.00; }
        if ((i == 41) & (j == 61))   { vbeta[count] <- 0.71; }
        if ((i == 56) & (j == 91))   { vbeta[count] <- 3.16; }
        if ((i == 65) & (j == 85))   { vbeta[count] <- 2.24; }
        if ((i == 86) & (j == 96))   { vbeta[count] <- 0.89; }
        if ((i == 101) & (j == 105)) { vbeta[count] <- 1.00; }
        if ((i == 111) & (j == 121)) { vbeta[count] <- 2.24; }
        
        count <- count + 1
    }
}

p <- ncol(X)
p

sigma2.true <- 20
beta_true <- vbeta
tag_true <- 1 * (beta_true != 0)
nzc <- sum(tag_true)
cat('model size is : ', nzc, '\n')


y_mu <- X %*% matrix(vbeta)
beta_c <- c(0, beta_true)
X_std <- scale(X, center = T, scale = T)
# X_center <- scale(X, center = T, scale = F)



#########################################################################################################

N <- 50  # repeat N times
res_all <- matrix(0, N, 32)


clnum <- detectCores()    # Get the number of CPU core
cat('cores', clnum, '\n')
cl <- makeCluster(60)
registerDoParallel(cl)


time_start <- Sys.time()

# for (k in 1:N) {

res_all <- foreach (
    k = 1:N,
    .combine = rbind,
    .inorder = TRUE,
    .export = ls(.GlobalEnv),
    .packages = c('MASS', 'ncvreg', 'EMVS', 'BMS', 'SIS', 'tictoc')
) %dopar% {
    
    
    
    #############################################################################
    
    set.seed(k)
    
    cat('doing the loop', k, '\n')
    
    y <- y_mu + rnorm(nrow(X), 0, sqrt(sigma2.true))
    
    y_center <- y - mean(y)
    
    
    
    #############################################################################
    
    Xy <- data.frame(x = X, y = y)
    
    Xc <- cbind(rep(1, n), X)
    
    #############################################################################
    #### lasso - ncvreg
    tic()
    # cv_lasso <- cv.ncvreg(X, y, penalty = 'lasso')
    model_lasso <- ebic.ncvreg(X, y, penalty = 'lasso')
    lasso_time <- toc()
    lasso_time <- lasso_time$toc - lasso_time$tic
    
    beta_hat <- model_lasso$beta
    # beta_hat <- coef(cv_lasso)
    tag_hat <- 1 * (as.vector(beta_hat[-1]) != 0)
    y_hat <- Xc %*% beta_hat
    
    f1_lasso <- cal_score(tag_true, tag_hat)$F1
    mse_lasso <- sum((y_mu - y_hat)^2)
    bias_lasso <- sum((beta_c - beta_hat)^2)
    
    cat('lasso finished', '\n')
    cat(f1_lasso, '\n')
    
    
    
    #############################################################################
    #### SCAD - ncvreg
    tic()
    model_scad <- ebic.ncvreg(X, y, penalty = 'SCAD')
    scad_time <- toc()
    scad_time <- scad_time$toc - scad_time$tic
    
    beta_hat <- model_scad$beta
    tag_hat <- 1 * (as.vector(beta_hat[-1]) != 0)
    y_hat <- Xc %*% beta_hat
    
    f1_scad <- cal_score(tag_true, tag_hat)$F1
    mse_scad <- sum((y_mu - y_hat)^2)
    bias_scad <- sum((beta_c - beta_hat)^2)
    
    cat('scad finished', '\n')
    cat(f1_scad, '\n')
    
    
    
    #############################################################################
    #### MCP - ncvreg
    # http://pbreheny.github.io/ncvreg/
    tic()
    # cv_mcp <- cv.ncvreg(X, y, penalty = 'MCP')
    model_mcp <- ebic.ncvreg(X, y, penalty = 'MCP')
    mcp_time <- toc()
    mcp_time <- mcp_time$toc - mcp_time$tic
    
    beta_hat <- model_mcp$beta
    tag_hat <- 1 * (as.vector(beta_hat[-1]) != 0)
    y_hat <- Xc %*% beta_hat
    
    f1_mcp <- cal_score(tag_true, tag_hat)$F1
    mse_mcp <- sum((y_mu - y_hat)^2)
    bias_mcp <- sum((beta_c - beta_hat)^2)
    
    cat('mcp finished', '\n')
    cat(f1_mcp, '\n')
    
    
    
    #############################################################################
    #### initial value chosen by ncvreg
    
    res.init.lasso   <- ncvreg(X, y, penalty = "lasso")
    mbeta.lasso      <- res.init.lasso$beta[-1, ]
    screening.lasso  <- as.vector(which(mbeta.lasso[, ncol(mbeta.lasso)] != 0))
    
    res.init.scad    <- ncvreg(X, y, penalty = "SCAD")
    mbeta.scad       <- res.init.scad$beta[-1, ]
    screening.scad   <- as.vector(which(mbeta.scad[, ncol(mbeta.scad)] != 0))
    
    res.init.mcp    <- ncvreg(X, y, penalty = "MCP")
    mbeta.mcp       <- res.init.mcp$beta[-1, ]
    screening.mcp   <- as.vector(which(mbeta.mcp[, ncol(mbeta.mcp)] != 0))
    
    screening.mix <- sort(unique(c(screening.mcp, screening.lasso, screening.scad)))
    # Xs            <- cbind(1, X[, screening.mix])
    
    
    
    #############################################################################
    ####### EMVS
    tic()
#    v0 = seq(0.1, 2, length.out = 20)
    v0 = exp(seq(log(0.001),log(1),,50))
    v1 = 1000
    beta_init = rep(0, p)
    beta_init[screening.mix] = 1
    sigma_init = 1
    a = b = 1
    epsilon = 10E-4
    model_emvs <- EMVS(y_center, X_std, v0 = v0, v1 = v1, type = 'betabinomial',
                       beta_init = beta_init, independent = F,
                       sigma_init = sigma_init, epsilon=epsilon,a = a, b = b)
    emvs_time <- toc()
    emvs_time <- emvs_time$toc - emvs_time$tic
    
    # emvs_best <- EMVSbest(model_emvs)
    emvs_best <- EMVSbest.my(model_emvs)
    
    emvs_index <- emvs_best$indices
    tag_hat <- rep(0, p)
    tag_hat[emvs_index] <- 1
    
    beta_hat <- emvs_best$beta.hat
    beta_hat_c <- transfer(X, emvs_best$beta.hat)
    beta_hat_c[1] <- beta_hat_c[1]+mean(y)
    y_hat <- Xc %*% beta_hat_c

    f1_emvs <- cal_score(tag_true, tag_hat)$F1
    mse_emvs <- sum((y_mu - y_hat)^2)
    bias_emvs <- sum((beta_c - beta_hat_c)^2)
    
    cat('EMVS finished', '\n')
    cat(f1_emvs, '\n')
    
    
    
    #############################################################################
    #### BMS, with screening.mix / sis.index as initial values 
    
    tic()
    model_bms <- bms(y ~ X, burn = 1000, iter = 1E6, nmodel = 1000, mcmc = "bd",
                     g = "hyper=3", mprior = "random", mprior.size = NA, user.int = TRUE,
                     start.value = screening.mcp, g.stats = TRUE, logfile = FALSE, logstep = 10000, 
                     force.full.ols = FALSE, fixed.reg = numeric(0))
    
    #model_bms <- bms(y ~ X, burn = 1000, iter = 1E6, nmodel = 500, mcmc = "bd",
    #                 g = "hyper=3", mprior = "random", mprior.size = NA, user.int = TRUE,
    #                 start.value = screening.mix, g.stats = TRUE, logfile = FALSE, logstep = 10000, 
    #                 force.full.ols = FALSE, fixed.reg = numeric(0))
    bms_time <- toc()
    bms_time <- bms_time$toc - bms_time$tic
    
    coef_bms   <- coef(model_bms, order.by.pip = FALSE)
    ind        <- coef_bms[, 5]
    tag_hat <- rep(0, p)
    tag_hat[ind] <- round(coef_bms[, 1])
    beta_hat <- coef_bms[, 2]
    y_hat <- predict(model_bms, exact = TRUE)
    
    f1_bms <- cal_score(tag_true, tag_hat)$F1
    mse_bms <- sum((y_mu - y_hat)^2)
    bias_bms <- sum((beta_true - beta_hat)^2)
    
    cat('BMS finished', '\n')
    cat(f1_bms, '\n')
    
    
    
    #############################################################################
    #### vb, path search
    tic()
    model_vb <- pathsearchVB(y, Xc)
    vb_time <- toc()
    vb_time <- vb_time$toc - vb_time$tic
    
    tag_hat <- 1 * (model_vb$best.gamma > 0.5)
    beta_hat <- model_vb$best.beta
    y_hat <- Xc %*% beta_hat
    
    f1_vb <- cal_score(tag_true, tag_hat)$F1
    mse_vb <- sum((y_mu - y_hat)^2)
    bias_vb <- sum((beta_c - beta_hat)^2)
    
    cat('vb finished', '\n')
    cat(f1_vb, '\n')
    
    
    
    #############################################################################
    #### cvb, with screening.mcp as initial values
    A <- 1
    B <- 1
    sigb <- 100000
    pi_v <- 0.5  # pi
    
    max_iter <- 100
    eps <- 1e-3
    
    tic()
    w_ini <- rep(0, p+1)  # null model
    w_ini[screening.mcp+1] <- 1
    model_cvb_ncvreg <- cvb_lr(w_ini, Xc, y, max_iter, eps, A, B, sigb, pi_v, 0.1)
    
    cvb_ncvreg_time <- toc()
    cvb_ncvreg_time <- cvb_ncvreg_time$toc - cvb_ncvreg_time$tic
    
    tag_hat_c <- round(model_cvb_ncvreg$w)
    tag_hat   <- tag_hat_c[-1]
    ind <- which(tag_hat_c == 1)
    lm_cvb <- lm(y ~ Xc[, ind]-1)
    beta_hat_c <- rep(0, p+1)
    beta_hat_c[ind] <- lm_cvb$coefficients
    y_hat <- Xc %*% beta_hat_c
    
    
    f1_cvb_ncvreg <- cal_score(tag_true, tag_hat)$F1
    mse_cvb_ncvreg <- sum((y_mu - y_hat)^2)
    bias_cvb_ncvreg <- sum((beta_c - beta_hat_c)^2)
    
    cat('cvb ncvreg finished', '\n')
    cat(f1_cvb_ncvreg, '\n')
    
    
    
    #############################################################################
    #### cvb, null as initial values
    A <- 1
    B <- 1
    sigb <- 100000
    pi_v <- 0.5  # pi
    
    max_iter <- 100
    eps <- 1e-3
    
    tic()
    w_ini <- rep(0, p+1)  # null model
    model_cvb0 <- cvb_lr(w_ini, Xc, y, max_iter, eps, A, B, sigb, pi_v,0.1)
    
    cvb0_time <- toc()
    cvb0_time <- cvb0_time$toc - cvb0_time$tic
    
    tag_hat_c <- round(model_cvb0$w)
    tag_hat   <- tag_hat_c[-1]
    ind <- which(tag_hat_c == 1)
    lm_cvb <- lm(y ~ Xc[, ind]-1)
    beta_hat_c <- rep(0, p+1)
    beta_hat_c[ind] <- lm_cvb$coefficients
    y_hat <- Xc %*% beta_hat_c
    
    f1_cvb0 <- cal_score(tag_true, tag_hat)$F1
    mse_cvb0 <- sum((y_mu - y_hat)^2)
    bias_cvb0 <- sum((beta_c - beta_hat_c)^2)
    
    cat('cvb0 finished', '\n')
    cat(f1_cvb0, '\n')
    
    
    #############################################################################
    
    cat(f1_lasso, f1_scad, f1_mcp, f1_emvs, f1_bms, f1_vb, f1_cvb_ncvreg, f1_cvb0, '\n')
    cat(mse_lasso, mse_scad, mse_mcp, mse_emvs, mse_bms, mse_vb, mse_cvb_ncvreg, mse_cvb0, '\n')
    cat(bias_lasso, bias_scad, bias_mcp, bias_emvs, bias_bms, bias_vb, bias_cvb_ncvreg, bias_cvb0, '\n')
    
    
    linkData <- c(f1_lasso, f1_scad, f1_mcp, f1_emvs, f1_bms, f1_vb, f1_cvb_ncvreg, f1_cvb0,  
                  mse_lasso, mse_scad, mse_mcp, mse_emvs, mse_bms, mse_vb, mse_cvb_ncvreg, mse_cvb0, 
                  bias_lasso, bias_scad, bias_mcp, bias_emvs, bias_bms, bias_vb, bias_cvb_ncvreg, bias_cvb0, 
                  lasso_time, scad_time, mcp_time, emvs_time, bms_time, vb_time, cvb_ncvreg_time, cvb0_time)
    
    # res_all[k, ] <- linkData
    
    return(linkData)
    
    
}


stopCluster(cl)


time_end <- Sys.time()
run_time <- time_end - time_start
print(run_time)





#########################################################################################################

method_names <- c('LASSO', 'SCAD', 'MCP', 'EMVS', 'BMS', 'VB', 'ACVB-MCP',"ACVB-null")


f1_all <- as.data.frame(res_all[, 1:8])
mse_all <- as.data.frame(log(res_all[, 9:16]/n))
bias_all <- as.data.frame(log(res_all[, 17:24]))
time_all <- as.data.frame(res_all[, 25:32])

colnames(f1_all) <- method_names
colnames(mse_all) <- method_names
colnames(bias_all) <- method_names
colnames(time_all) <- method_names

f1_mean <- colMeans(f1_all, na.rm = T)
mse_mean <- colMeans(mse_all, na.rm = T)
bias_mean <- colMeans(bias_all, na.rm = T)
time_mean <- colMeans(time_all, na.rm = T)

f1_sd <- apply(f1_all, 2, sd, na.rm = T)
time_sd <- apply(time_all, 2, sd, na.rm = T)

cat('mean of f1: ', '\n')
print(f1_mean)

cat('mean of mse: ', '\n')
print(mse_mean)

cat('mean of bias: ', '\n')
print(bias_mean)

cat('mean of time: ', '\n')
print(time_mean)
round(cbind(f1_mean,f1_sd,time_mean,mse_mean,apply(mse_all, 2,sd)),4)



save.image(file = 'QTL-MCP-NULL20220315.RData')



#########################################################################################################
##### box plot
my.method <-c("LASSO","SCAD","MCP","EMVS","BMS","VB","ACVB\nMCP","ACVB\nnull")

##### f1
fig_f1 <- ggplot(stack(f1_all), aes(x = ind, y = values, fill = ind)) + 
    stat_boxplot(geom ='errorbar', width = .45, linetype = 'solid') + 
    geom_boxplot(outlier.shape = NA) + xlab('') + ylab('F1 score') + 
    scale_x_discrete(labels = my.method)  + 
    theme(text = element_text(size = 18), legend.position = 'none', 
          axis.text.x = element_text(size = 15))
fig_f1

ggsave('QTL_f1.pdf', fig_f1, width = 6, height = 6, units = 'in', scale = 1)


##### mse
fig_mse <- ggplot(stack(mse_all), aes(x = ind, y = values, fill = ind)) + 
    stat_boxplot(geom ='errorbar', width = .45, linetype = 'solid') + 
    geom_boxplot(outlier.shape = NA) + xlab('') + ylab('log MSE') + 
    scale_x_discrete(labels = my.method)  +
    theme(text = element_text(size = 18), legend.position = 'none', 
          axis.text.x = element_text(size = 15)) 

fig_mse

ggsave('QTL_mse.pdf', fig_mse, width = 6, height = 6, units = 'in', scale = 1)

