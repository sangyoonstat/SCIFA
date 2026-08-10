rm(list = ls())

require(doMC) # load this package for parallel computing
ncores <- detectCores()
registerDoMC(cores = ncores - 1) # set the number of cores we are going to use for parallel computing

my_path <- "/Users/sayi/Library/CloudStorage/OneDrive-OklahomaAandMSystem/SCIFA-manuscript/" 

source(file.path(my_path, "source_scifa.R")) # load the source code we have
load(file.path(my_path, "all_rdat.RData"))
# all_rdat.RData file has the following objects:
# 1. aux_df: data frame of auxiliary covariates
# 2. names: taxon name in each view
# 3. X1,X2,X3: clr-transformed microbiome data
# Original source: https://github.com/bioFAM/MOFA2_tutorials/blob/master/R_tutorials/microbiome_vignette.html

X_list <- list(X1, X2, X3)
ss <- stdX(X_list, center = T)
X_centered <- list(X_list[[1]] - ss$mean_mat[[1]], X_list[[2]] - ss$mean_mat[[2]], 
                   X_list[[3]] - ss$mean_mat[[3]])
pvec <- unlist(lapply(X_list, ncol))
n <- nrow(X1)
clr_indc <- rep(T, 3)

scifa_out1 <- run_algo(X_list, clr_indc, pvec, rhat_vec = NULL, nlamb = 64, c_min = 10^-5, c_max = 0.4, 
                       niter_free = 20, max_iter = 500, eps = 10^-5, center = T, seednum = 2025)

######################################
## For the SCIFA results in Table 5 ##
######################################

sd_hat <- scifa_out1$shat

aic_M1 <- tcrossprod(scifa_out1$chosen_by_AIC$U0, sd_hat[1]*scifa_out1$chosen_by_AIC$Vd[[1]]) +
  tcrossprod(scifa_out1$chosen_by_AIC$Ud[[1]], sd_hat[1]*scifa_out1$chosen_by_AIC$Ad[[1]])
aic_M2 <- tcrossprod(scifa_out1$chosen_by_AIC$U0, sd_hat[2]*scifa_out1$chosen_by_AIC$Vd[[2]])
aic_M3 <- tcrossprod(scifa_out1$chosen_by_AIC$U0, sd_hat[3]*scifa_out1$chosen_by_AIC$Vd[[3]]) +
  tcrossprod(scifa_out1$chosen_by_AIC$Ud[[3]], sd_hat[3]*scifa_out1$chosen_by_AIC$Ad[[3]])

round(sum(aic_M1**2)/sum(X_centered[[1]]**2)*100, 3)
round(sum(aic_M2**2)/sum(X_centered[[2]]**2)*100, 3)
round(sum(aic_M3**2)/sum(X_centered[[3]]**2)*100, 3)

bic_M1 <- tcrossprod(scifa_out1$chosen_by_BIC$U0, sd_hat[1]*scifa_out1$chosen_by_BIC$Vd[[1]]) +
  tcrossprod(scifa_out1$chosen_by_BIC$Ud[[1]], sd_hat[1]*scifa_out1$chosen_by_BIC$Ad[[1]])

bic_M2 <- tcrossprod(scifa_out1$chosen_by_BIC$U0, sd_hat[2]*scifa_out1$chosen_by_BIC$Vd[[2]])

bic_M3 <- tcrossprod(scifa_out1$chosen_by_BIC$U0, sd_hat[3]*scifa_out1$chosen_by_BIC$Vd[[3]]) +
  tcrossprod(scifa_out1$chosen_by_BIC$Ud[[3]], sd_hat[3]*scifa_out1$chosen_by_BIC$Ad[[3]])

round(sum(bic_M1**2)/sum(X_centered[[1]]**2)*100, 3)
round(sum(bic_M2**2)/sum(X_centered[[2]]**2)*100, 3)
round(sum(bic_M3**2)/sum(X_centered[[3]]**2)*100, 3)

######################################
## For the SCIFA results in Table 6 ##
######################################

aic_resid1 <- X_centered[[1]] - aic_M1 ; aic_sv1 <- svd(aic_resid1)
aic_resid2 <- X_centered[[2]] - aic_M2 ; aic_sv2 <- svd(aic_resid2)
aic_resid3 <- X_centered[[3]] - aic_M3 ; aic_sv3 <- svd(aic_resid3)

proj_U0_aic <- (diag(n) - tcrossprod(scifa_out1$chosen_by_AIC$U0))  
aic_M1_extended <- aic_M1 + proj_U0_aic%*%aic_sv1$u[,1:3]%*%diag(aic_sv1$d[1:3])%*%t(aic_sv1$v[,1:3])
aic_M2_extended <- aic_M2 + proj_U0_aic%*%aic_sv2$u[,1:5]%*%diag(aic_sv2$d[1:5])%*%t(aic_sv2$v[,1:5])
aic_M3_extended <- aic_M3 + proj_U0_aic%*%aic_sv3$u[,1:3]%*%diag(aic_sv3$d[1:3])%*%t(aic_sv3$v[,1:3])

round(sum(aic_M1_extended**2)/sum(X_centered[[1]]**2)*100, 3)
round(sum(aic_M2_extended**2)/sum(X_centered[[2]]**2)*100, 3)
round(sum(aic_M3_extended**2)/sum(X_centered[[3]]**2)*100, 3)

bic_resid1 <- X_centered[[1]] - bic_M1 ; bic_sv1 <- svd(bic_resid1)
bic_resid2 <- X_centered[[2]] - bic_M2 ; bic_sv2 <- svd(bic_resid2)
bic_resid3 <- X_centered[[3]] - bic_M3 ; bic_sv3 <- svd(bic_resid3)

proj_U0_bic <- (diag(n) - tcrossprod(scifa_out1$chosen_by_BIC$U0))  
bic_M1_extended <- bic_M1 + proj_U0_bic%*%bic_sv1$u[,1:3]%*%diag(bic_sv1$d[1:3])%*%t(bic_sv1$v[,1:3])
bic_M2_extended <- bic_M2 + proj_U0_bic%*%bic_sv2$u[,1:5]%*%diag(bic_sv2$d[1:5])%*%t(bic_sv2$v[,1:5])
bic_M3_extended <- bic_M3 + proj_U0_bic%*%bic_sv3$u[,1:3]%*%diag(bic_sv3$d[1:3])%*%t(bic_sv3$v[,1:3])

round(sum(bic_M1_extended**2)/sum(X_centered[[1]]**2)*100, 3)
round(sum(bic_M2_extended**2)/sum(X_centered[[2]]**2)*100, 3)
round(sum(bic_M3_extended**2)/sum(X_centered[[3]]**2)*100, 3)



