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

algo_fn <- function(X_list, pvec, rhat_vec, clr_indc, lamb, curr_U0, curr_Ud, 
                    niter_free = 10, max_iter = 500, eps = 10^-6, seednum = NULL){
  
  D <- length(X_list) ; n <- nrow(X_list[[D]]) 
  
  orth_err_vec <- rep(NA_real_, max_iter)
  
  orth_err_by_view <- matrix(NA_real_, nrow = max_iter, ncol = D, 
                             dimnames = list(NULL, paste0("view_", seq_len(D))))

  calc_orth_err <- function(U0, Uk) {
    if (is.null(Uk) ||
        ncol(Uk) == 0L ||
        is.null(U0) ||
        ncol(U0) == 0L) {
      return(0)
    }
    
    cross_mat <- crossprod(Uk, U0)
    
    if (!all(is.finite(cross_mat))) {
      return(Inf)
    }
    
    max(abs(cross_mat))
  }
  
  obj_val <- c() ; obj_val[1] <- -Inf 
  obj_mat <- matrix(NA, nrow = max_iter, ncol = 3)
  diff <- 1 ; num_iter <- 0
  Ud_hat <- Vd_hat <- Ad_hat <- list()
  
  if(!is.null(seednum)){
    set.seed(seednum)
  }
  
  while(num_iter <= 2*niter_free || abs(diff) > eps){
    
    num_iter <- num_iter + 1
    
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
      R_list <- mapply(function(Xd, Id) Xd - Id, X_list, curr_Id, SIMPLIFY = F)
      weighted_R <- mapply(function(Rd, pk) Rd / sqrt(pk), R_list, pvec, SIMPLIFY = F)
      weighted_V <- mapply(function(Vd, pk) Vd / sqrt(pk), Vd_hat, pvec, SIMPLIFY = F)
      conc_R <- do.call(cbind, weighted_R) ; stack_V <- do.call(rbind, weighted_V)
      U0_hat <- upd_Umat(conc_R, stack_V)
      
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
      
      # Compute the bases of span(Ud_hat[[1]],...,Ud_hat[[D]]) 
      Uall <- do.call(cbind, Ud_hat)                           
      if (is.null(Uall) || ncol(Uall) == 0) {
        Wmat <- matrix(numeric(0), nrow = n, ncol = 0)         
      } else {
        sv <- svd(Uall)
        rk <- sum(sv$d > 1e-10)                               
        Wmat <- if (rk == 0) matrix(numeric(0), nrow = n, ncol = 0) else sv$u[, 1:rk, drop = FALSE]
      }

      # Update joint score
      curr_Id <- mapply(function(Ud, Ad){tcrossprod(Ud, Ad)}, Ud_hat, Ad_hat, SIMPLIFY = F)
      R_list <- mapply(function(Xd, Id) Xd - Id, X_list, curr_Id, SIMPLIFY = F)
      weighted_R <- mapply(function(Rd, pk) Rd / sqrt(pk), R_list, pvec, SIMPLIFY = F)
      weighted_V <- mapply(function(Vd, pk) Vd / sqrt(pk), Vd_hat, pvec, SIMPLIFY = F)
      conc_R <- do.call(cbind, weighted_R) 
      proj_R <- (diag(n) - tcrossprod(Wmat)) %*% conc_R      
      stack_V <- do.call(rbind, weighted_V)
      U0_hat <- upd_Umat(proj_R, stack_V)       

      obj_mat[num_iter, 3] <- obj_fn(X_list, U0_hat, Ud_hat, Vd_hat, Ad_hat, rep(lamb, D), rep(lamb, D))$sum      
    }
    
    ## Record joint-individual score orthogonality for each iteration
    orth_err_by_view[num_iter, ] <- vapply(
      seq_len(D),
      function(k) {
        calc_orth_err(
          U0 = U0_hat,
          Uk = Ud_hat[[k]]
        )
      },
      numeric(1)
    )
    
    orth_err_vec[num_iter] <- max(orth_err_by_view[num_iter, ], na.rm = TRUE)

    obj_val[num_iter + 1] <- obj_fn(X_list, U0_hat, Ud_hat, Vd_hat, Ad_hat, rep(lamb, D), rep(lamb, D))$sum
    diff <- obj_val[num_iter + 1] - obj_val[num_iter]  
    
    curr_U0 <- U0_hat ; curr_Ud <- Ud_hat
    
    if(num_iter==max_iter){
      break
    }
  }
  
  used_iter <- seq_len(num_iter)

  return(list(U0 = U0_hat, Vd = Vd_hat, Ud = Ud_hat, Ad = Ad_hat, obj_val = obj_val[seq_len(num_iter + 1)], 
              obj_mat = obj_mat[used_iter, , drop = FALSE], orth_err = orth_err_vec[used_iter]))
}

.scifa_iter_soft <- function(x, t) sign(x) * pmax(abs(x) - t, 0)

scifa_iter_zero_sum_soft <- function(c, kappa, clr = TRUE, root_rel_tol = 1e-12) {
  stopifnot(is.numeric(c), length(c) > 0L, all(is.finite(c)),
            length(kappa) == 1L, is.finite(kappa), kappa >= 0,
            is.logical(clr), length(clr) == 1L, !is.na(clr),
            length(root_rel_tol) == 1L, is.finite(root_rel_tol),
            root_rel_tol > 0, root_rel_tol < 1)
  c <- as.numeric(c)
  scale_c <- max(1, abs(c), kappa)
  y <- c / scale_c
  t <- kappa / scale_c
  alpha <- 0
  width <- 0
  iterations <- 0L
  zero_solution <- FALSE
  
  if (!clr) {
    b_scaled <- .scifa_iter_soft(y, t)
  } else if (kappa == 0) {
    alpha <- mean(y)
    b_scaled <- y - alpha
  } else if ((max(y) - min(y)) <= 2 * t) {
    # The all-zero solution has a nonempty dual interval.
    # Its midpoint is a valid multiplier, including ties and p=1.
    alpha <- min(y) / 2 + max(y) / 2
    b_scaled <- rep(0, length(y))
    zero_solution <- TRUE
  } else {
    lo <- min(y)
    hi <- max(y)
    for (it in seq_len(200L)) {
      mid <- lo + (hi - lo) / 2
      g <- sum(.scifa_iter_soft(y - mid, t))
      iterations <- it
      if (g == 0) {
        lo <- hi <- mid
        break
      }
      if (g > 0) lo <- mid else hi <- mid
      if (hi - lo <= root_rel_tol) break
    }
    width <- hi - lo
    alpha <- lo + width / 2
    # Given a support/sign pattern the root is linear. Accept a polishing step
    # only if it preserves that pattern (up to exact floating-point comparisons).
    for (it in seq_len(8L)) {
      u <- y - alpha
      active <- abs(u) > t
      if (!any(active)) break
      sg <- sign(u[active])
      a2 <- (sum(y[active]) - t * sum(sg)) / sum(active)
      active2 <- abs(y - a2) > t
      if (!identical(active, active2) ||
          !identical(sign((y - a2)[active2]), sg)) break
      old_alpha <- alpha
      alpha <- a2
      if (alpha == old_alpha) break
    }
    b_scaled <- .scifa_iter_soft(y - alpha, t)
  }
  b <- b_scaled * scale_c
  alpha_out <- alpha * scale_c
  # This is a numerical-boundary flag, NOT a rule for deleting coefficients.
  guard <- max(root_rel_tol, width,
               64 * .Machine$double.eps * (1 + max(abs(y)) + t + abs(alpha)))
  margin <- abs(abs(y - alpha) - t)
  list(b = b, alpha = alpha_out, near_boundary = sum(margin <= guard),
       min_boundary_margin_scaled = min(margin),
       root_width_scaled = width, root_iterations = iterations,
       zero_solution = zero_solution, scale = scale_c)
}

# Loading update when the CURRENT score matrix has Q'Q = I.
# Xd must already be on the original algo_fn fitting scale (pp$Z_list[[d]]).
# Solves RSS/(2*n*p) + lambda*sum(abs(B)); keeps lambda and CLR constraints.
# No coefficient thresholding, score rotation, or CVXR fallback is performed.
upd_load_orth_soft <- function(Xd, U0, Ud, lambda, clr = TRUE,
                               root_rel_tol = 1e-12,
                               orth_tol = 1e-10, check_tol = 1e-9) {
  is_mat <- function(x) is.matrix(x) && is.numeric(x) && all(is.finite(x))
  maxabs <- function(x) if (length(x)) max(abs(x)) else 0
  stopifnot(is_mat(Xd), is_mat(U0), is_mat(Ud),
            nrow(Xd) > 0L, ncol(Xd) > 0L,
            nrow(U0) == nrow(Xd), nrow(Ud) == nrow(Xd),
            is.logical(clr), length(clr) == 1L, !is.na(clr),
            is.numeric(lambda), length(lambda) == 1L,
            is.finite(lambda), lambda >= 0,
            length(root_rel_tol) == 1L, is.finite(root_rel_tol),
            root_rel_tol > 0, root_rel_tol < 1,
            length(orth_tol) == 1L, is.finite(orth_tol), orth_tol > 0,
            length(check_tol) == 1L, is.finite(check_tol), check_tol > 0)
  
  n <- nrow(Xd); p <- ncol(Xd)
  r0 <- ncol(U0); rd <- ncol(Ud)
  Q <- cbind(U0, Ud)
  q <- ncol(Q)
  G <- crossprod(Q)
  gram_error <- maxabs(G - diag(q))
  if (gram_error > orth_tol) {
    stop(sprintf(paste0("Score Gram error %.3g exceeds orth_tol = %.3g. ",
                        "Direct loading update requires Q'Q = I. ",
                        "Scores were not changed; no CVXR fallback."),
                 gram_error, orth_tol), call. = FALSE)
  }
  
  C <- crossprod(Xd, Q)           # p x (r0 + rd)
  kappa <- n * p * lambda         # from the ORIGINAL loss denominator 2*n*p
  if (!is.finite(kappa)) stop("Nonfinite n*p*lambda.", call. = FALSE)
  B <- matrix(0, nrow = p, ncol = q)
  alpha <- numeric(q)
  for (ell in seq_len(q)) {
    sol <- scifa_iter_zero_sum_soft(C[, ell], kappa, clr, root_rel_tol)
    B[, ell] <- sol$b
    alpha[ell] <- sol$alpha
  }
  
  # Check the KKT equations with the ACTUAL Gram matrix, not just I.
  if (q > 0L) {
    BG <- B %*% G
    station <- sweep(BG - C, 2L, alpha, "+")
    residual <- pmax(abs(station) - kappa, 0)
    active <- B != 0
    residual[active] <- abs(station[active] + kappa * sign(B[active]))
    kkt_scaled <- maxabs(residual) / max(1, maxabs(C), kappa, maxabs(BG))
    clr_sum <- if (clr) maxabs(colSums(B)) else 0
    clr_scaled <- if (clr) {
      max(abs(colSums(B)) / pmax(1, colSums(abs(B))))
    } else 0
  } else {
    kkt_scaled <- clr_sum <- clr_scaled <- 0
  }
  if (!is.finite(kkt_scaled) || !is.finite(clr_scaled) ||
      kkt_scaled > check_tol || clr_scaled > check_tol) {
    stop(sprintf("Loading check failed: scaled KKT=%.3g, scaled CLR=%.3g.",
                 kkt_scaled, clr_scaled), call. = FALSE)
  }
  
  # seq_len(0) safely preserves p x 0 matrices for zero-rank blocks.
  list(Vhat = B[, seq_len(r0), drop = FALSE],
       Ahat = B[, r0 + seq_len(rd), drop = FALSE],
       gram_error = gram_error, kkt_scaled = kkt_scaled, clr_sum = clr_sum)
}


algo_fn_orth_loading <- function(X_list, pvec, rhat_vec, clr_indc, lamb, curr_U0, curr_Ud, 
                                 niter_free = 10, max_iter = 500, eps = 10^-6, seednum = NULL,
                                 root_rel_tol = 1e-12, orth_tol = 1e-10,
                                 loading_check_tol = 1e-9){
  
  D <- length(X_list) ; n <- nrow(X_list[[D]]) 
  
  # NEW: diagnostics for the loading INPUT scores and direct loading solve.
  load_gram <- load_kkt <- load_clr <- matrix(NA_real_, max_iter, D)
  
  orth_err_vec <- rep(NA_real_, max_iter)
  
  orth_err_by_view <- matrix(NA_real_, nrow = max_iter, ncol = D, 
                             dimnames = list(NULL, paste0("view_", seq_len(D))))
  
  calc_orth_err <- function(U0, Uk) {
    if (is.null(Uk) ||
        ncol(Uk) == 0L ||
        is.null(U0) ||
        ncol(U0) == 0L) {
      return(0)
    }
    
    cross_mat <- crossprod(Uk, U0)
    
    if (!all(is.finite(cross_mat))) {
      return(Inf)
    }
    
    max(abs(cross_mat))
  }
  
  obj_val <- c() ; obj_val[1] <- -Inf 
  obj_mat <- matrix(NA, nrow = max_iter, ncol = 3)
  diff <- 1 ; num_iter <- 0
  Ud_hat <- Vd_hat <- Ad_hat <- list()
  
  if(!is.null(seednum)){
    set.seed(seednum)
  }
  
  while(num_iter <= 2*niter_free || abs(diff) > eps){
    
    num_iter <- num_iter + 1
    
    # Update loading: only this block changes the fitting calculation.
    for(d in 1:D){
      if(num_iter >= niter_free + 2L){
        # NEW: input scores come from a completed constrained score update.
        uobj <- tryCatch(
          upd_load_orth_soft(
            Xd = X_list[[d]], U0 = curr_U0, Ud = curr_Ud[[d]],
            lambda = lamb, clr = clr_indc[d],
            root_rel_tol = root_rel_tol, orth_tol = orth_tol,
            check_tol = loading_check_tol
          ),
          error = function(e) stop(sprintf(
            "Iteration %d, view %d: %s", num_iter, d, conditionMessage(e)
          ), call. = FALSE)
        )
        Vd_hat[[d]] <- uobj$Vhat ; Ad_hat[[d]] <- uobj$Ahat
        load_gram[num_iter, d] <- uobj$gram_error
        load_kkt[num_iter, d] <- uobj$kkt_scaled
        load_clr[num_iter, d] <- uobj$clr_sum
      } else{
        # UNCHANGED: original loading update through niter_free + 1.
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
      R_list <- mapply(function(Xd, Id) Xd - Id, X_list, curr_Id, SIMPLIFY = F)
      weighted_R <- mapply(function(Rd, pk) Rd / sqrt(pk), R_list, pvec, SIMPLIFY = F)
      weighted_V <- mapply(function(Vd, pk) Vd / sqrt(pk), Vd_hat, pvec, SIMPLIFY = F)
      conc_R <- do.call(cbind, weighted_R) ; stack_V <- do.call(rbind, weighted_V)
      U0_hat <- upd_Umat(conc_R, stack_V)
      
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
      
      # Compute the bases of span(Ud_hat[[1]],...,Ud_hat[[D]]) 
      Uall <- do.call(cbind, Ud_hat)                           
      if (is.null(Uall) || ncol(Uall) == 0) {
        Wmat <- matrix(numeric(0), nrow = n, ncol = 0)         
      } else {
        sv <- svd(Uall)
        rk <- sum(sv$d > 1e-10)                               
        Wmat <- if (rk == 0) matrix(numeric(0), nrow = n, ncol = 0) else sv$u[, 1:rk, drop = FALSE]
      }
      
      # Update joint score
      curr_Id <- mapply(function(Ud, Ad){tcrossprod(Ud, Ad)}, Ud_hat, Ad_hat, SIMPLIFY = F)
      R_list <- mapply(function(Xd, Id) Xd - Id, X_list, curr_Id, SIMPLIFY = F)
      weighted_R <- mapply(function(Rd, pk) Rd / sqrt(pk), R_list, pvec, SIMPLIFY = F)
      weighted_V <- mapply(function(Vd, pk) Vd / sqrt(pk), Vd_hat, pvec, SIMPLIFY = F)
      conc_R <- do.call(cbind, weighted_R) 
      proj_R <- (diag(n) - tcrossprod(Wmat)) %*% conc_R      
      stack_V <- do.call(rbind, weighted_V)
      U0_hat <- upd_Umat(proj_R, stack_V)       
      
      obj_mat[num_iter, 3] <- obj_fn(X_list, U0_hat, Ud_hat, Vd_hat, Ad_hat, rep(lamb, D), rep(lamb, D))$sum      
    }
    
    ## Record joint-individual score orthogonality for each iteration
    orth_err_by_view[num_iter, ] <- vapply(
      seq_len(D),
      function(k) {
        calc_orth_err(
          U0 = U0_hat,
          Uk = Ud_hat[[k]]
        )
      },
      numeric(1)
    )
    
    orth_err_vec[num_iter] <- max(orth_err_by_view[num_iter, ], na.rm = TRUE)
    
    obj_val[num_iter + 1] <- obj_fn(X_list, U0_hat, Ud_hat, Vd_hat, Ad_hat, rep(lamb, D), rep(lamb, D))$sum
    diff <- obj_val[num_iter + 1] - obj_val[num_iter]  
    
    curr_U0 <- U0_hat ; curr_Ud <- Ud_hat
    
    if(num_iter==max_iter){
      break
    }
  }
  
  used_iter <- seq_len(num_iter)
  
  # Diagnostic summaries only; NO extra final loading update is performed.
  max_recorded <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  loading_diagnostics <- data.frame(
    iteration = used_iter,
    method = ifelse(used_iter >= niter_free + 2L, "orth_soft", "original"),
    input_gram_error = apply(load_gram[used_iter, , drop = FALSE], 1L, max_recorded),
    KKT_scaled = apply(load_kkt[used_iter, , drop = FALSE], 1L, max_recorded),
    CLR_colsum_max = apply(load_clr[used_iter, , drop = FALSE], 1L, max_recorded)
  )
  
  return(list(U0 = U0_hat, Vd = Vd_hat, Ud = Ud_hat, Ad = Ad_hat, obj_val = obj_val[seq_len(num_iter + 1)],
              obj_mat = obj_mat[used_iter, , drop = FALSE], orth_err = orth_err_vec[used_iter],
              loading_diagnostics = loading_diagnostics,
              lambda = lamb,
              settings = list(niter_free = niter_free, max_iter = max_iter, eps = eps,
                              seednum = seednum, first_direct_iteration = niter_free + 2L,
                              root_rel_tol = root_rel_tol, orth_tol = orth_tol,
                              loading_check_tol = loading_check_tol,
                              final_loading_update = FALSE)))
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
    tryCatch(algo_fn(pp$Z_list, pvec, rhat_vec, clr_indc, lamb_seq[l],
                     init$init_U0, init$init_Ud,
                     niter_free, max_iter, eps, seednum + l), error=function(e) NA)
  }
  
  # Compute two versions of the GIC loss:
  # 1. standardized fitting scale
  # 2. original data scale after back-scaling
  loss_std_vec <- loss_orig_vec <- nonz_vec <- rep(NA_real_, nlamb)

  for(l in 1:nlamb){
    if(all(!is.na(res[[l]]))){
      loss_std_vec[l]  <- 0 ; loss_orig_vec[l] <- 0
      for(k in seq_len(D)){
        # Estimated signal on the standardized scale
        Mhat_std_k <- tcrossprod(res[[l]]$U0, res[[l]]$Vd[[k]]) + tcrossprod(res[[l]]$Ud[[k]],res[[l]]$Ad[[k]])
        # GIC loss evaluated on the standardized fitting scale
        loss_std_vec[l] <- loss_std_vec[l] + sum((pp$Z_list[[k]] - Mhat_std_k)^2) / (2 * n * pvec[k])
        
        # Estimated signal after returning to the original data scale
        Mhat_orig_k <- pp$mean_mat[[k]] + pp$shat[k] * Mhat_std_k
        # GIC loss evaluated on the original data scale
        loss_orig_vec[l] <- loss_orig_vec[l] + sum((X_list[[k]] - Mhat_orig_k)^2) / (2 * n * pvec[k])
      }
      # Count nonzero loading coefficients on the fitting scale
      nonz_vec[l] <- sum(sapply(res[[l]]$Vd, function(mat) sum(abs(mat) > 10^-10))) + sum(sapply(res[[l]]$Ad, function(mat) sum(abs(mat) > 10^-10)))
    }
  }
  
  # Each CLR view contributes r_0 + r_k zero-sum constraints
  num_clr_constraint <- sum((rhat_vec[1] + rhat_vec[-1])[clr_indc])
  effective_df_vec <- nonz_vec - num_clr_constraint
  
  # select lambda for each criterion 
  total <- n*sum(pvec)
  tau_vec <- c(2, log(total))
  
  # Standardized-scale GIC
  aic_score_std <- log(loss_std_vec) + tau_vec[1] * effective_df_vec / total
  bic_score_std <- log(loss_std_vec) + tau_vec[2] * effective_df_vec / total
  aic_ind_std  <- which.min(aic_score_std) ; bic_ind_std  <- which.min(bic_score_std)
  
  # Original-scale GIC
  aic_score_orig <- log(loss_orig_vec) + tau_vec[1] * effective_df_vec / total
  bic_score_orig <- log(loss_orig_vec) + tau_vec[2] * effective_df_vec / total
  aic_ind_orig <- which.min(aic_score_orig) ; bic_ind_orig <- which.min(bic_score_orig)

  return(list("chosen_by_AIC_std" = res[[aic_ind_std]], "chosen_by_BIC_std" = res[[bic_ind_std]],
              "chosen_by_AIC_orig" = res[[aic_ind_orig]], "chosen_by_BIC_orig" = res[[bic_ind_orig]],
              "shat" = pp$shat, "mean_mat" = mean_mat, "rhat_vec" = rhat_vec, 
              "aic_ind_std" = aic_ind_std, "bic_ind_std" = bic_ind_std,
              "aic_ind_orig" = aic_ind_orig, "bic_ind_orig" = bic_ind_orig))
}

run_algo_final <- function(X_list, clr_indc, pvec, rhat_vec = NULL, nlamb, c_min = 10^-5, c_max = 0.4, 
                           niter_free = 15, max_iter = 200, eps = 10^-5, center = T, seednum = NULL, nnz_tol = 10^-10){
  
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
  res <- foreach(l = seq_len(nlamb), .errorhandling = "pass") %dopar% {
    tryCatch(
      algo_fn_orth_loading(
        pp$Z_list, pvec, rhat_vec, clr_indc, lamb_seq[l],
        init$init_U0, init$init_Ud,
        niter_free, max_iter, eps,
        seednum = if (is.null(seednum)) NULL else seednum + l
      ),
      error = function(e) e
    )
  }
  
  # Validate only fitted matrices, NOT diagnostic fields that can contain NA
  # or the original objective trace whose first element is -Inf.
  is_mat <- function(x, nr, nc) {
    is.matrix(x) && is.numeric(x) &&
      all(dim(x) == c(nr, nc)) && all(is.finite(x))
  }
  valid_fit <- function(fit) {
    if (!is.list(fit) || inherits(fit, "error") ||
        !is_mat(fit$U0, n, rhat_vec[1]) ||
        !is.list(fit$Ud) || !is.list(fit$Vd) || !is.list(fit$Ad) ||
        length(fit$Ud) != D || length(fit$Vd) != D || length(fit$Ad) != D) {
      return(FALSE)
    }
    all(vapply(seq_len(D), function(k) {
      is_mat(fit$Ud[[k]], n, rhat_vec[k + 1L]) &&
        is_mat(fit$Vd[[k]], pvec[k], rhat_vec[1]) &&
        is_mat(fit$Ad[[k]], pvec[k], rhat_vec[k + 1L])
    }, logical(1)))
  }
  
  # Compute the GIC
  loss_std_vec <- nonz_vec <- rep(NA_real_, nlamb)
  fit_errors <- rep(NA_character_, nlamb)
  n_iter <- rep(NA_integer_, nlamb)
  
  for (l in seq_len(nlamb)) {
    fit <- res[[l]]
    if (!valid_fit(fit)) {
      fit_errors[l] <- if (inherits(fit, "error")) {
        conditionMessage(fit)
      } else {
        "Missing, nonfinite or incorrectly dimensioned fitted matrices."
      }
      next
    }
    
    loss_std_vec[l] <- 0
    for (k in seq_len(D)) {
      Mhat_std_k <- tcrossprod(fit$U0, fit$Vd[[k]]) +
        tcrossprod(fit$Ud[[k]], fit$Ad[[k]])
      loss_std_vec[l] <- loss_std_vec[l] +
        sum((pp$Z_list[[k]] - Mhat_std_k)^2)
    }
    nonz_vec[l] <-
      sum(vapply(fit$Vd, function(mat) sum(abs(mat) > nnz_tol), numeric(1))) +
      sum(vapply(fit$Ad, function(mat) sum(abs(mat) > nnz_tol), numeric(1)))
    
    if (is.matrix(fit$obj_mat)) n_iter[l] <- nrow(fit$obj_mat)
    if (!is.finite(loss_std_vec[l])) {
      fit_errors[l] <- "Pooled standardized RSS is nonfinite."
    }
  }
  
  # Each CLR view contributes r_0 + r_k zero-sum constraints
  num_clr_constraint <- sum((rhat_vec[1] + rhat_vec[-1])[clr_indc])
  effective_df_vec <- nonz_vec - num_clr_constraint
  if (any(effective_df_vec < 0, na.rm = TRUE)) {
    warning("Negative fixed-correction df values found; retained as specified. Inspect effective_df_vec.")
  }
  
  # select lambda for each criterion 
  total <- n*sum(pvec)
  tau_vec <- c(2, log(total))
  
  good <- which(is.finite(loss_std_vec) & loss_std_vec >= 0 & is.finite(effective_df_vec))
  
  if (!length(good)) {
    why <- unique(fit_errors[!is.na(fit_errors)])
    stop("No valid candidates for GIC. ", paste(head(why, 3L), collapse = " | "))
  }
  
  aic_score <- bic_score <- rep(NA_real_, nlamb)
  fit_term <- loss_std_vec[good] / total
  aic_score[good] <- fit_term + tau_vec[1] * effective_df_vec[good] / total
  bic_score[good] <- fit_term + tau_vec[2] * effective_df_vec[good] / total
  if (any(!is.finite(aic_score[good])) || any(!is.finite(bic_score[good]))) {
    stop("Nonfinite GIC values; inspect RSS and df before selection.")
  }
  aic_ind <- good[which.min(aic_score[good])]
  bic_ind <- good[which.min(bic_score[good])]
  
  if (length(good) < nlamb) {
    warning(sprintf("GIC evaluated on %d/%d valid candidates. See fit_errors.",
                    length(good), nlamb))
  }
  
  # Retain fits that reach max_iter, as in the original wrapper, but flag them.
  hit_max_iter <- n_iter >= max_iter
  if (any(hit_max_iter[good], na.rm = TRUE)) {
    warning("Some candidates reached max_iter; see hit_max_iter. They were not automatically excluded.")
  }
  
  return(list(chosen_by_AIC = res[[aic_ind]], chosen_by_BIC = res[[bic_ind]], shat = pp$shat, 
              mean_mat = mean_mat, rhat_vec = rhat_vec, aic_ind_std = aic_ind, 
              bic_ind_std = bic_ind, path = res, lamb_seq = lamb_seq, 
              loss_unscaled_std = loss_std_vec, nonz_vec = nonz_vec, 
              effective_df_vec = effective_df_vec, aic_score = aic_score, 
              bic_score = bic_score, valid_indices = good, fit_errors = fit_errors,
              n_iter = n_iter, hit_max_iter = hit_max_iter))
}

