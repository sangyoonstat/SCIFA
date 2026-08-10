
library(MASS)
library(glmnet)
library(Matrix)
library(mvtnorm)
library(CVXR)

########################################################
##### 0.1 function to standardize data as in SLIDE #####
########################################################

## Input
# X_list : data from each view in each element of X_list
# If center = F, do scaling only

## Output
# Z_list : list where each element is standardized data for each view
# mean_mat: list where each element is resulting matrix for each view after column centering
# fb_vec : vector of scaled Frobenius norms for scaling each view
# If center = F, each element of fb_vec is the squared Frobenius norm of the original data
# If center = T, each element of fb_vec is the squared Frobenius norm of the column-centered data

stdX <- function(X_list, center = F){
  
  D <- length(X_list) ; Z_list <- mean_mat <- list() ; fb_vec <- c() 
  
  if(center){
    for(d in 1:D){
      col_means <- colMeans(X_list[[d]])
      mean_mat[[d]] <- matrix(col_means, nrow = nrow(X_list[[d]]), ncol = ncol(X_list[[d]]), byrow = T)
      centered_X <- X_list[[d]] - mean_mat[[d]] 
      fb_vec[d] <- sqrt(sum(centered_X**2)) ; Z_list[[d]] <- centered_X/fb_vec[d] 
    }
  } else{
    for(d in 1:D){
      mean_mat[[d]] <- matrix(0, nrow = nrow(X_list[[d]]), ncol = ncol(X_list[[d]])) 
      fb_vec[d] <- sqrt(sum(X_list[[d]]**2)) ; Z_list[[d]] <- X_list[[d]]/fb_vec[d] 
    }
  }
  
  return(list(Z_list = Z_list, mean_mat = mean_mat, fb_vec = fb_vec))
  
}

##############################################################
### 0.2 Function to do refitting to induce sparse loadings ###
##############################################################

estimate_sparse_loading_cv <- function(X, score){
  p <- ncol(X)
  k <- ncol(score)
  L <- matrix(0, nrow = p, ncol = k)
  
  if(k == 1){
    for (j in 1:p) {
      cv_fit <- cv.glmnet(cbind(score, 0), X[, j], alpha = 1, intercept = FALSE)
      coef_j <- coef(cv_fit, s = "lambda.min")[2]  # intercept 제거
      L[j, ] <- as.numeric(coef_j)
    }
  } else{
    for (j in 1:p) {
      cv_fit <- cv.glmnet(score, X[, j], alpha = 1, intercept = FALSE)
      coef_j <- coef(cv_fit, s = "lambda.min")[-1]  # intercept 제거
      L[j, ] <- as.numeric(coef_j)
    } 
  }
  return(L)
}


###############################################################
#### 1. function to compute the MAD estimates for noise sd ####
###############################################################
# This is a modified version of denoiseR::estim_sigma
my_estim_sigma_eff <- function(X, n_eff = NULL, p_eff = NULL){
  n <- nrow(X); p <- ncol(X)
  #if (is.null(n_eff)) n_eff <- n - 1  
  #if (is.null(p_eff)) p_eff <- p      
  s <- svd(X, nu = 0, nv = 0)$d
  beta <- min(n_eff, p_eff) / max(n_eff, p_eff)
  lambdastar <- sqrt( 2*(beta + 1) + 8*beta/(beta + 1 + sqrt(beta^2 + 14*beta + 1)) )
  wbstar <- 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43
  median(s) / ( sqrt(max(n_eff, p_eff)) * (lambdastar / wbstar) )
}

#######################################################################################
#### 2. function to scale the data by the MAD estimates and get the rank estimates ####
#######################################################################################

process_dat <- function(X_list, pvec, clr_indc, center = T){
  
  D <- length(X_list)
  n <- nrow(X_list[[1]])
  rankMd <- dim_indiv <- shat <- rep(NA, D)
  Z_list <- centered_X <- mean_mat <- list()
  
  # set the effective sample size and number of variable
  #n = 100 ;  center = F
  n_eff <- n - ifelse(center, 1, 0)
  p_eff <- pvec - as.integer(clr_indc)
  
  for(d in 1:D){
    if(center){
      col_means <- colMeans(X_list[[d]])
      mean_mat[[d]] <- matrix(col_means, nrow = nrow(X_list[[d]]), ncol = ncol(X_list[[d]]), byrow = T)
      centered_X[[d]] <- X_list[[d]] - mean_mat[[d]] 
    } else{
      mean_mat[[d]] <- matrix(0, nrow = nrow(X_list[[d]]), ncol = ncol(X_list[[d]]))
      centered_X[[d]] <- X_list[[d]] - mean_mat[[d]] 
    }
    shat[d] <- my_estim_sigma_eff(centered_X[[d]], n_eff = n_eff, p_eff = p_eff[d])
    Z_list[[d]] <- centered_X[[d]]/shat[d]
    rankMd[d] <- sum(svd(Z_list[[d]])$d > (sqrt(n_eff) + sqrt(p_eff[d])))
  }
  
  Zmat <- do.call(cbind, Z_list)
  rankM <- sum(svd(Zmat)$d > (sqrt(n_eff) + sqrt(sum(p_eff))))
  numer <- sum(rankMd) - rankM
  
  dim_joint <- ifelse(numer < 0, 0, round(numer/(D - 1)))
  dim_indiv <- rankMd - dim_joint
  dim_indiv[dim_indiv < 0] <- 0
  
  return(list(Z_list = Z_list, rankM = rankM, rankMd = rankMd, 
              dim_joint = dim_joint, dim_indiv = dim_indiv, shat = shat, mean_mat = mean_mat))
}

#####################################
### 3. functions to simulate data ###
#####################################

###################################################################
#### 3.1 function to orthogonalize each column of given matrix ####
###################################################################

gs_wo_normalize <- function(X){
  n <- ncol(X)
  Q <- matrix(0, nrow = nrow(X), ncol = n)
  
  for (i in 1:n) {
    v <- X[, i]
    if (i > 1) {
      for (j in 1:(i - 1)) {
        proj <- sum(Q[, j] * X[, i]) / sum(Q[, j]^2) * Q[, j]
        v <- v - proj
      }
    }
    Q[, i] <- v  # Do NOT normalize
  }
  return(Q)
}

##############################################################
#### 3.2 function to generate a sparse, orthogonal matrix ####
##############################################################

# Inputs
# p: the total number of variables to be considered
# r: the rank of resulting matrix
# zero_prop: the desired proportion of 0's out of p times r many total elements 
# (note that the actual proportion of zero elements in the resulting matrix 
# might not be exactly the same with the specified zero_prop due to the orthonormalization step) 
# sum_zero: a logical to indicate whether the zero-sum constraint would be satisfied or not

# Output
# Vmat: the resulting p by k orthonormal matrix with zero elements

gen_sp_orth_load <- function(p, r, zero_prop, sum_zero = F, normalize = F){
  
  num_elem <- p * r
  Vtmp <- matrix(rnorm(num_elem), nrow = p, ncol = r)
  set_zeros <- sample(1:num_elem, size = ceiling(zero_prop*num_elem))
  Vtmp[set_zeros] <- 0
  if(sum_zero){
    for(k in 1:r){
      nonz_ind <- which(abs(Vtmp[,k]) > 10^-6)
      Vtmp[nonz_ind,k] <- Vtmp[nonz_ind,k] - mean(Vtmp[nonz_ind, k]) 
    }
  }
  
  Vtmp2 <- gs_wo_normalize(Vtmp)
  
  if(normalize){
    cnorm <- sqrt(colSums(Vtmp2^2))
    Vmat <- sweep(Vtmp2, 2, cnorm, "/")
  } else{
    Vmat <- Vtmp2
  }
  
  return(Vmat)
}


gen_sp_orth_load2 <- function(p, r, zero_prop, sum_zero = F, normalize = F){
  
  num_elem <- p * r
  Vtmp <- matrix(rnorm(num_elem), nrow = p, ncol = r)
  set_zeros <- sample(1:num_elem, size = ceiling(zero_prop*num_elem))
  Vtmp[set_zeros] <- 0
  if(sum_zero){
    for(k in 1:r){
      nonz_ind <- which(abs(Vtmp[,k]) > 10^-6)
      Vtmp[nonz_ind,k] <- Vtmp[nonz_ind,k] - mean(Vtmp[nonz_ind, k]) 
    }
  }
  
  if(normalize){
    cnorm <- sqrt(colSums(Vtmp^2))
    Vmat <- sweep(Vtmp2, 2, cnorm, "/")
  } else{
    Vmat <- Vtmp
  }
  
  return(Vmat)
}


#####################################################################################
### 3.3. function to simulate multi-view data matrices with clr-transformed views ###
#####################################################################################
# Inputs
# n: sample size
# pvec: the D by 1 vector where each element is the total number of variables for each view 
# clr_indc: the D by 1 logical vector to indicate whether the d-th view is assumed to be clr-transformed (if TRUE)
# sgv_joint: r0 by 1 vector that consist of singular values for joint structure
# sgv_indiv: the list of length D where each element is r_d by 1 vector of singular values for each individual structure
# zp_joint: the D by 1 vector where each element is the desired proportion of 0's out of pvec[d]*r0 many total elements
# zp_indiv: the D by 1 vector where each element is the desired proportion of 0's out of pvec[d]*r_d many total elements
# SNR_vec: D by 1 vector of signal-to-noise ratio for each view
# orth: a logical to indicate whether the each loading vector would be orthogonalized or not
# norm: a logical to indicate whether the each loading vector would be normalized or not

# Outputs
# X_list: the list of length D where each element is n by pvec[d] observed data matrix
# U0: n by r0 score matrix for joint structure
# Ud_list: the list of length D where each score matrix is for each individual structure
# V0_list: the list of length D where each element corresponds to the loading matrix of each view for joint structure
# Ad_list: the list of length D where each element corresponds to the loading matrix of each view for individual structure

dgen_fn <- function(n, pvec, clr_indc, sgv_joint, sgv_indiv, zp_joint, zp_indiv, snr_vec, 
                    gen_mu = F, min_mu = -0.5, max_mu = 0.5, norm_load = F){
  
  p_total <- sum(pvec) ; D <- length(pvec) 
  r0 <- length(sgv_joint)
  rvec <- sapply(c(1:D), FUN = function(d){length(sgv_indiv[[d]])})
  
  ##########################################
  ### Step 1: Compute the score matrices ###
  ##########################################
  
  tempU0 = matrix(rnorm(n*r0), n, r0)
  if(gen_mu){
    U0 <- qr.Q(qr((diag(1, n) - matrix(1/n, nrow = n, ncol = n))%*%tempU0))
  } else{
    U0 <- qr.Q(qr(tempU0)) 
  }
  
  # Individual scores
  Ud_list <- list()
  
  for(d in 1:D){
    tempUd <- matrix(rnorm(n*rvec[d]), nrow = n, ncol = rvec[d])
    if(gen_mu){
      U_0d <- qr.Q(qr(cbind(1, U0)))
      Ud_list[[d]] <-  qr.Q(qr((diag(1, n) - tcrossprod(U_0d))%*%tempUd))
    } else{
      Ud_list[[d]] <-  qr.Q(qr((diag(1, n) - tcrossprod(U0))%*%tempUd)) 
    }
  }
  
  ############################################
  ### Step 2: Compute the loading matrices ###
  ############################################
  
  # V: Joint loading, A: Individual loading
  V_list <- A_list <- list()
  for(d in 1:D){
    if(clr_indc[d]){
      # for the view with compositional data
      tempVd <- gen_sp_orth_load2(pvec[d], r0, zp_joint[d], sum_zero = T, normalize = norm_load)
      tempAd <- gen_sp_orth_load2(pvec[d], rvec[d], zp_indiv[d], sum_zero = T, normalize = norm_load) 
    } else{
      tempVd <- gen_sp_orth_load2(pvec[d], r0, zp_joint[d], sum_zero = F, normalize = norm_load)
      tempAd <- gen_sp_orth_load2(pvec[d], rvec[d], zp_indiv[d], sum_zero = F, normalize = norm_load)
    }
    V_list[[d]] <- tempVd%*%diag(sgv_joint)
    A_list[[d]] <- tempAd%*%diag(sgv_indiv[[d]])
  }
  
  ##################################################
  ### Step 3: Compute the observed data matrices ###
  ##################################################
  
  # Compute the signal matrices and generate the noise matrices to compute the observed data matrices
  X_list <- M_list <- M_ctr_list <- err_list <- mu_list <- vector("list", D)
  sig_vec <- rep(NA, D)
  # d <- 2
  for(d in 1:D){
    conca_load <- cbind(V_list[[d]], A_list[[d]]) 
    #M_list[[d]] <- cbind(U0, Ud_list[[d]])%*%diag(c(sgv_joint, sgv_indiv[[d]]))%*%t(conca_load)
    if(gen_mu){
      mu_list[[d]] <- runif(pvec[d], min_mu, max_mu)
      #mu_list[[d]] <- runif(pvec[d], -0.5, 0.5)
      M_ctr_list[[d]] <- tcrossprod(cbind(U0, Ud_list[[d]]), conca_load)
      M_list[[d]] <- outer(rep(1, n), mu_list[[d]]) + M_ctr_list[[d]]
    } else{
      mu_list[[d]] <- NULL
      M_list[[d]] <- M_ctr_list[[d]] <- tcrossprod(cbind(U0, Ud_list[[d]]), conca_load) 
    }
    
    if(clr_indc[d]){
      sig_vec[d] <- sqrt(sum(M_ctr_list[[d]]^2) / (n * (pvec[d] - 1) * snr_vec[d]))
      C_p <- diag(1, pvec[d]) - matrix(1/pvec[d], nrow = pvec[d], ncol = pvec[d])
      cov_mat <- C_p%*%diag(sig_vec[d]^2, pvec[d])%*%C_p
      err_list[[d]] <- rmvnorm(n, mean = rep(0, pvec[d]), sigma = cov_mat)
      X_list[[d]] <- M_list[[d]] + err_list[[d]]
    } else{
      sig_vec[d] <- sqrt(sum(M_ctr_list[[d]]^2) / (n * pvec[d] * snr_vec[d]))
      err_list[[d]] <- matrix(rnorm(n*pvec[d], sd = sig_vec[d]), nrow = n, ncol = pvec[d])
      X_list[[d]] <- M_list[[d]] + err_list[[d]]
    }
  }
  
  return(list(X_list = X_list, mu_list = mu_list, U0 = U0, Ud_list = Ud_list, V_list = V_list, 
              A_list = A_list, M_list = M_list, err_list = err_list, sig_vec = sig_vec))
}


########################################################
#### 4. function to evaluate the objective function ####
########################################################

obj_fn <- function(X_list, U0_mat, Ud_list, V_list, A_list, lvec_V, lvec_A){
  
  D <- length(X_list) 
  scaled_loss <- penV_vec <- penA_vec <- rep(NA, D)
  
  for(d in 1:D){
    fb_err <- sum((X_list[[d]] - tcrossprod(cbind(U0_mat, Ud_list[[d]]), cbind(V_list[[d]], A_list[[d]])))**2)
    scaled_loss[d] <- fb_err/(2*prod(dim(X_list[[d]])))
    penV_vec[d] <- lvec_V[d]*sum(abs(V_list[[d]])) 
    penA_vec[d] <-lvec_A[d]*sum(abs(A_list[[d]]))
  }
  
  sum <- sum(scaled_loss + penV_vec + penA_vec)
  
  return(list(sum = sum, scaled_loss = scaled_loss, penV_vec = penV_vec, penA_vec = penA_vec))
}


#####################################################################
#### 4. function to compute the initial values based on the data ####  
#####################################################################

get_init <- function(X_list, rhat_vec){
  
  D <- length(X_list) ; n <- nrow(X_list[[1]])
  conc_X <- do.call(cbind, X_list)
  svd_conc <- svd(conc_X)
  
  init_U0 <- svd_conc$u[,c(1:rhat_vec[1]), drop = F]
  
  init_Ud <- list()
  for(d in 1:D){
    svd_d <- svd((diag(1, n) - tcrossprod(init_U0))%*%X_list[[d]])
    if(rhat_vec[d+1] > 0){
      init_Ud[[d]] <- svd_d$u[,c(1:rhat_vec[d+1]), drop = F] 
    } else{
      init_Ud[[d]] <- matrix(numeric(0), nrow = n, ncol = 0)
    }
  }
  
  return(list(init_U0 = init_U0, init_Ud = init_Ud))
}

################################################################################
#### 5. function to update the score matrix as in the Orthogonal Procrustes ####
################################################################################

upd_Umat <- function(X, V){
  r <- ncol(V)
  sobj <- svd(X%*%V, nu = r, nv = r)
  return(tcrossprod(sobj$u, sobj$v))
}

###################################################
#### 6. functions to update the loading matrix ####
###################################################

# Separate functions depending on the loading for compositional or non-compositional below

################################
#### for compositional view ####
################################

## function below to update the loading matrix corresponding to compositional view
# subject to the constraint based on cvxr

upd_load_with_const_cvxr <- function(Xd, U0, Ud, lambda,
                                     solver = "ECOS",
                                     feastol = 1e-8, reltol = 1e-8, abstol = 1e-8){
  
  n  <- nrow(Xd); p <- ncol(Xd)
  r0 <- ncol(U0); rd <- ncol(Ud)
  
  V <- CVXR::Variable(shape = c(p, r0), name = "V")
  A <- CVXR::Variable(shape = c(p, rd), name = "A")
  
  X_hat <- U0 %*% t(V) + Ud %*% t(A)
  
  loss <- sum((Xd - X_hat)^2)
  penalty <- lambda * (sum(abs(V)) + sum(abs(A)))
  
  ones_p <- matrix(1, 1, p)
  constraints <- list(
    ones_p %*% V == matrix(0, 1, r0),
    ones_p %*% A == matrix(0, 1, rd)
  )
  
  problem <- CVXR::Problem(
    CVXR::Minimize((1/(2*n*p)) * loss + penalty),
    constraints
  )
  
  CVXR::psolve(
    problem,
    solver  = solver,
    feastol = feastol,
    reltol  = reltol,
    abstol  = abstol
  )
  
  problem_status <- CVXR::status(problem)
  
  if (!(problem_status %in%
        c("optimal", "optimal_inaccurate"))) {
    return(list(
      Vhat   = matrix(NA_real_, p, r0),
      Ahat   = matrix(NA_real_, p, rd),
      status = problem_status
    ))
  }
  
  Vhat <- matrix(
    as.numeric(CVXR::value(V)),
    nrow = p,
    ncol = r0
  )
  
  Ahat <- matrix(
    as.numeric(CVXR::value(A)),
    nrow = p,
    ncol = rd
  )
  
  if (!all(is.finite(Vhat)) ||
      !all(is.finite(Ahat))) {
    return(list(
      Vhat   = matrix(NA_real_, p, r0),
      Ahat   = matrix(NA_real_, p, rd),
      status = "nonfinite_solution"
    ))
  }
  
  list(Vhat = Vhat, Ahat = Ahat, status = problem_status)
}



upd_load_with_const_cvxr_ver2 <- function(Xd, U0, lambda, solver = "ECOS",
                                          feastol = 1e-8,
                                          reltol = 1e-8,
                                          abstol = 1e-8){
  
  n <- nrow(Xd) ; p <- ncol(Xd)
  r0 <- ncol(U0) 
  
  # Define variables
  V <- CVXR::Variable(shape = c(p, r0), name = "V")
  
  # Reconstruct Xd using U0 and Ud
  X_hat <- U0 %*% t(V) 
  
  # Objective: squared loss + lasso penalty
  loss <- CVXR::sum_squares(Xd - X_hat)
  penalty <- lambda * CVXR::norm1(V)
  
  # Constraints: column sums of V and A are zero
  ones_p <- matrix(1, nrow = 1, ncol = p)  # 1 x p
  constraints <- list(
    ones_p %*% V == matrix(0, nrow = 1, ncol = r0)
  )
  
  problem <- CVXR::Problem(
    CVXR::Minimize(
      loss / (2 * n * p) + penalty
    ),
    constraints
  )
  
  CVXR::psolve(
    problem,
    solver  = solver,
    feastol = feastol,
    reltol  = reltol,
    abstol  = abstol
  )
  
  problem_status <- CVXR::status(problem)
  
  if (!(problem_status %in%
        c("optimal", "optimal_inaccurate"))) {
    return(list(
      Vhat   = matrix(NA_real_, p, r0),
      status = problem_status
    ))
  }
  
  Vhat <- matrix(
    as.numeric(CVXR::value(V)),
    nrow = p,
    ncol = r0
  )
  
  if (!all(is.finite(Vhat))) {
    return(list(
      Vhat   = matrix(NA_real_, p, r0),
      status = "nonfinite_solution"
    ))
  }
  
  list(
    Vhat   = Vhat,
    status = problem_status
  )
}

#########################################################
#### for non-compositional view (w/o any constraint) ####
#########################################################

## (a) function to update the loading matrix based on glmnet

upd_load_wo_const_old <- function(Xd, U0, Ud, lambda){
  
  p <- ncol(Xd) ; r0 <- ncol(U0) ; rd <- ncol(Ud)
  pmat <- cbind(kronecker(diag(1, p), U0), kronecker(diag(1, p), Ud))
  resp_vec <- as.vector(Xd)
  fit <- glmnet(pmat, resp_vec, lambda = lambda, intercept = F, standardize = F)
  coef_est <- as.vector(coef(fit))[-1] # remove the intercept  
  V_trans <- matrix(coef_est[c(1:(r0*p))], nrow = r0, ncol = p)
  A_trans <- matrix(coef_est[-c(1:(r0*p))], nrow = rd, ncol = p)
  
  return(list(Vhat = t(V_trans), Ahat = t(A_trans)))
}

####################################################
#### 7. Algorithm based on the functions above  ####  
####################################################

algo_fn <- function(X_list, rhat_vec, clr_indc, lamb, curr_U0, curr_Ud, 
                    niter_free = 10, max_iter = 500, eps = 10^-6, seednum = NULL){
  
  D <- length(X_list) ; n <- nrow(X_list[[D]]) 
  
  obj_val <- c() ; obj_val[1] <- -Inf 
  obj_mat <- matrix(NA, nrow = max_iter, ncol = 3)
  diff <- 1 ; num_iter <- 0
  Ud_hat <- Vd_hat <- Ad_hat <- list()
  
  if(!is.null(seednum)){
    set.seed(seednum)
  }
  
  while(num_iter <= 2*niter_free || abs(diff) > eps){
    
    num_iter <- num_iter + 1
    #upd_load_with_const_cvxr_ver2 <- function(Xd, U0, lambda)    
    # Update loading 
    for(d in 1:D){
      if(rhat_vec[d + 1] > 0){
        if(clr_indc[d]){
          uobj <- upd_load_with_const_cvxr(X_list[[d]], curr_U0, curr_Ud[[d]], lamb)
        } else{
          uobj <- upd_load_wo_const_old(X_list[[d]], curr_U0, curr_Ud[[d]], lamb)
        }
        Vd_hat[[d]] <- uobj$Vhat ; Ad_hat[[d]] <- uobj$Ahat        
      } else{
        if(clr_indc[d]){
          uobj <- upd_load_with_const_cvxr_ver2(X_list[[d]], curr_U0, lamb)
        } else{
          uobj <- upd_load_wo_const_old(X_list[[d]], curr_U0, curr_Ud[[d]], lamb)
        }
        Vd_hat[[d]] <- uobj$Vhat 
        Ad_hat[[d]] <- matrix(0, nrow = ncol(X_list[[d]]), ncol = 0)
      }
    }
    
    obj_mat[num_iter, 1] <- obj_fn(X_list, curr_U0, curr_Ud, Vd_hat, Ad_hat, rep(lamb, D), rep(lamb, D))$sum
    
    if(num_iter <= niter_free){
      # Update individual score
      for(d in 1:D){
        if(rhat_vec[d+1] > 0){
          Ud_hat[[d]] <- upd_Umat(X_list[[d]] - tcrossprod(curr_U0, Vd_hat[[d]]), Ad_hat[[d]])   
        } else{
          Ud_hat[[d]] <- matrix(0, nrow = n, ncol = 0)
        }
      }
      
      obj_mat[num_iter, 2] <- obj_fn(X_list, curr_U0, Ud_hat, Vd_hat, Ad_hat, rep(lamb, D), rep(lamb, D))$sum
      
      # Update joint score
      curr_Id <- mapply(function(Ud, Ad){tcrossprod(Ud, Ad)}, Ud_hat, Ad_hat, SIMPLIFY = F)
      conc_R <- do.call(cbind, mapply(function(Xd, Id){Xd - Id}, X_list, curr_Id, SIMPLIFY = F))
      U0_hat <- upd_Umat(conc_R, do.call(rbind, Vd_hat))
      
      obj_mat[num_iter, 3] <- obj_fn(X_list, U0_hat, Ud_hat, Vd_hat, Ad_hat, rep(lamb, D), rep(lamb, D))$sum
      
    } else{
      # Update individual score
      proj_mat <- diag(n) - tcrossprod(curr_U0)
      for(d in 1:D){
        proj_resid <- proj_mat%*%(X_list[[d]] - tcrossprod(curr_U0, Vd_hat[[d]]))
        if(rhat_vec[d + 1] > 0){
          Ud_hat[[d]] <- upd_Umat(proj_resid, Ad_hat[[d]])   
          Ud_hat[[d]] <- Ud_hat[[d]] - curr_U0 %*% crossprod(curr_U0, Ud_hat[[d]])
          if (ncol(Ud_hat[[d]]) > 0) {
            Ud_hat[[d]] <- qr.Q(qr(Ud_hat[[d]]))[, 1:ncol(Ud_hat[[d]]), drop = FALSE]
          } else {
            Ud_hat[[d]] <- matrix(numeric(0), nrow = n, ncol = 0)
          }
        } else{
          Ud_hat[[d]] <- matrix(numeric(0), nrow = n, ncol = 0)
        }
      }
      
      obj_mat[num_iter, 2] <- obj_fn(X_list, curr_U0, Ud_hat, Vd_hat, Ad_hat, rep(lamb, D), rep(lamb, D))$sum
      
      # Update joint score
      curr_Id <- mapply(function(Ud, Ad){tcrossprod(Ud, Ad)}, Ud_hat, Ad_hat, SIMPLIFY = F)
      conc_R <- do.call(cbind, mapply(function(Xd, Id){Xd - Id}, X_list, curr_Id, SIMPLIFY = F))
      
      # span(Ud_hat[[1]],...,Ud_hat[[D]])의 정규직교 기저 Wmat 생성 (SVD로 수치 랭크 컷)
      Uall <- do.call(cbind, Ud_hat)                           # n x (sum r_d) 또는 n x 0
      if (is.null(Uall) || ncol(Uall) == 0) {
        Wmat <- matrix(numeric(0), nrow = n, ncol = 0)         # 공차원
      } else {
        sv <- svd(Uall)
        rk <- sum(sv$d > 1e-10)                                # 임계값은 데이터 스케일에 맞게 조정 가능
        Wmat <- if (rk == 0) matrix(numeric(0), nrow = n, ncol = 0) else sv$u[, 1:rk, drop = FALSE]
      }
      
      # U0는 Ud들의 여공간으로 투영된 잔차에서 업데이트
      proj_R <- (diag(n) - tcrossprod(Wmat)) %*% conc_R
      U0_hat <- upd_Umat(proj_R, do.call(rbind, Vd_hat))       # n x r0 (임시)
      
      # (3) U0를 Ud들의 공간에 대해 정확히 직교화 + 정규직교화
      if (!is.null(Uall) && ncol(Uall) > 0) {
        U0_hat <- U0_hat - Uall %*% crossprod(Uall, U0_hat)
      }
      U0_hat <- qr.Q(qr(U0_hat))[, 1:ncol(U0_hat), drop = FALSE]
      
      obj_mat[num_iter, 3] <- obj_fn(X_list, U0_hat, Ud_hat, Vd_hat, Ad_hat, rep(lamb, D), rep(lamb, D))$sum      
      #Wmat <- qr.Q(qr(do.call(cbind, Ud_hat))) # the orthonormal basis of span(Ud_hat[[1]],...,Ud_hat[[D]])
      #proj_R <- (diag(n) - tcrossprod(Wmat))%*%conc_R
      #U0_hat <- upd_Umat(proj_R, do.call(rbind, Vd_hat)) 
    }
    
    obj_val[num_iter + 1] <- obj_fn(X_list, U0_hat, Ud_hat, Vd_hat, Ad_hat, rep(lamb, D), rep(lamb, D))$sum
    diff <- obj_val[num_iter + 1] - obj_val[num_iter]  
    
    curr_U0 <- U0_hat ; curr_Ud <- Ud_hat
    
    if(num_iter==max_iter){
      break
    }
  }
  
  return(list(U0 = U0_hat, Vd = Vd_hat, Ud = Ud_hat, Ad = Ad_hat, 
              obj_val = obj_val, obj_mat = obj_mat))
}

#############################################################################
### wrapper to fit SCIFA over tuning parameter and print out final result ### 
#############################################################################

run_algo <- function(X_list, clr_indc, pvec, rhat_vec = NULL, nlamb, c_min = 10^-5, c_max = 0.4, 
                     niter_free = 15, max_iter = 200, eps = 10^-5, center = T, seednum){
  
  n <- nrow(X_list[[1]])
  D <- length(pvec)
  # preprocessing input data
  pp <- process_dat(X_list, pvec, clr_indc, center)
  if(is.null(rhat_vec)){
    rhat_vec <- c(pp$dim_joint, pp$dim_indiv)  
  }
  mean_mat <- do.call(cbind, pp$mean_mat)
  
  # initialization 
  init <- get_init(pp$Z_list, rhat_vec) 
  lmax_vec <- rep(NA, D)
  
  # setting up the grid for tuning parameters
  for(d in 1:D){
    cp_mat <- crossprod(cbind(init$init_U0, init$init_Ud[[d]]), pp$Z_list[[d]])
    lmax_vec[d] <- max(abs(cp_mat)/prod(dim(pp$Z_list[[d]])))
  }
  
  lmin <- min(lmax_vec)
  lamb_seq <- seq(lmin*c_min, lmin*c_max, len = nlamb)
  
  # run Alrogithm 1 for each element in lamb_seq
  res <- foreach(l = 1:nlamb, .errorhandling="pass") %dopar% {
    tryCatch(algo_fn(pp$Z_list, rhat_vec, clr_indc, lamb_seq[l],
                     init$init_U0, init$init_Ud,
                     niter_free, max_iter, eps, seednum + l), error=function(e) NA)
  }
  
  # compute the two terms appearing in the GIC-type formula over each lambda
  loss_vec <- df_vec <- rep(NA, nlamb)
  
  for(l in 1:nlamb){
    if(all(!is.na(res[[l]]))){
      # need to scale back both loadings by the noise level estimates
      scaled_V <- mapply(function(Vmat, shat){Vmat*shat}, 
                         res[[l]]$Vd, pp$shat, SIMPLIFY = FALSE)
      
      scaled_A <- mapply(function(Amat, shat){Amat*shat}, 
                         res[[l]]$Ad, pp$shat, SIMPLIFY = FALSE)
      
      joint_est <- do.call(cbind, lapply(scaled_V, 
                                         FUN = function(mat){tcrossprod(res[[l]]$U0, mat)}))
      
      indiv_est <- do.call(cbind, mapply(function(score, loading){tcrossprod(score, loading)}, 
                                         res[[l]]$Ud, scaled_A, SIMPLIFY = F))
      
      Mhat <- joint_est + indiv_est + mean_mat # estimated signal 
      loss_vec[l] <- sum((do.call(cbind, X_list) - Mhat)^2) # compute Frobenius norm loss 
      df_vec[l] <- sum(unlist(lapply(scaled_A, function(mat) sum(mat > 10^-10)))) + sum(unlist(lapply(scaled_V, function(mat) sum(mat > 10^-10))))
    }
  }
  
  # select lambda for each criterion 
  total <- n*sum(pvec)
  tau_vec <- c(2, log(total))
  aic_ind <- which.min(log(loss_vec/total) + tau_vec[1]*(df_vec - rhat_vec[1] - sum(rhat_vec[clr_indc]))/total)
  bic_ind <- which.min(log(loss_vec/total) + tau_vec[2]*(df_vec - rhat_vec[1] - sum(rhat_vec[clr_indc]))/total)
  
  return(list("chosen_by_AIC" = res[[aic_ind]], "chosen_by_BIC" = res[[bic_ind]],
              "shat" = pp$shat, "mean_mat" = mean_mat, "rhat_vec" = rhat_vec, "aic_ind" = aic_ind, "bic_ind" = bic_ind))
}
