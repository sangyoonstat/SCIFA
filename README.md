# SCIFA 

## 1. File description

**final_source_scifa.R** - R code for our main method and several auxiliary functions 

**fit_rdat_github.R** - R code to apply SCIFA for the real data

**all_rdat.RData** - file that includes real data we used. The original source can be found https://github.com/bioFAM/MOFA2_tutorials/blob/master/R_tutorials/microbiome_vignette.html and we saved and minimaly edited the original files to be consistent with our usage. 

## Supported view types

SCIFA supports both CLR-transformed compositional views and non-compositional views. The logical vector `clr_indc` specifies the type of each view, in the same order as the matrices in `X_list`:

- `TRUE`: the view is CLR-transformed, and its joint and individual loading columns are constrained to sum to zero.
- `FALSE`: the view is non-compositional, and its sparse loadings are estimated without zero-sum constraints.

Supply CLR-transformed data for compositional views and ordinary non-compositional measurements for the other views. The fitting function does not apply the CLR transformation automatically. Any subset of views can be marked `TRUE`; the two-view example below uses `clr_indc = c(TRUE, TRUE)`.

## 2. Example with two clr-transformed views

### 2.1 Generate data
```{r}

my_path <- "/Users/sayi/Library/CloudStorage/OneDrive-OklahomaAandMSystem/SCIFA-manuscript/" 
source(paste(my_path, sep = "", "final_source_scifa.R")) # load the source code we have

num_clr <- 2 # the number of clr-transformed microbial views
n = 100 # sample size
p1 = 100 ; p2 = 100 # the number of variables in each view
zpj1 = 0.8 ; zpj2 = 0.85 # target zero proportion for joint loading in each view
zpi1 = 0.7 ; zpi2 = 0.7 # target zero proportion for individual loading in each view
snr1 = 0.25 ; snr2 = 0.25 # Signal-to-noise ratio values for each view 

pvec <- c(p1, p2)
zp_joint <- c(zpj1, zpj2)
zp_indiv <- c(zpi1, zpi2)
snr_vec <- c(snr1, snr2)

# logical indicator for clr-transformed views 
if(num_clr < length(pvec)){
  clr_indc <- c(rep(F, length(pvec) - num_clr), rep(T, num_clr))  
} else{
  clr_indc <- rep(T, length(pvec))
}

# specify singular values 
cons = -0.3
sgv_joint <- c(2, 1.6) + cons
sgv_indiv <- list(c(1.7, 1.5) + cons, c(1.9, 1.3) + cons) 

# simulate data
seednum <- 2026
set.seed(seednum)
dat <- dgen_fn(n, pvec, clr_indc, sgv_joint, sgv_indiv, zp_joint, zp_indiv, snr_vec, gen_mu = F)

```

### 2.2 Run SCIFA
```{r}
require(doMC) # load this package for parallel computing
ncores <- detectCores()
registerDoMC(cores = ncores - 2) # set the number of cores we are going to use for parallel computing
# fit SCIFA 
scifa_out <- run_algo_final(dat$X_list, clr_indc, pvec, rhat_vec = NULL, nlamb = 80, c_min = 10^-5, c_max = 0.4,
                            niter_free = 20, max_iter = 200, eps = 10^-5, center = T, seednum = seednum)

# save the result
scifa_res <- list("chosen_by_AIC" = scifa_out$chosen_by_AIC, 
                  "chosen_by_BIC" = scifa_out$chosen_by_BIC,
                  "shat" = scifa_out$shat, "mean_mat" = scifa_out$mean_mat, 
                  "rhat" = scifa_out$rhat_vec)

# compute the estimate signal based on the BIC-type criterion
# for 1st view
s1_hat <- scifa_res$shat[1]
bic_M1 <- scifa_res$chosen_by_BIC$U0 %*% t(s1_hat*scifa_res$chosen_by_BIC$Vd[[1]]) + scifa_res$chosen_by_BIC$Ud[[1]] %*% t(s1_hat*scifa_res$chosen_by_BIC$Ad[[1]])
sum((dat$M_list[[1]] - bic_M1)^2) / sum(dat$M_list[[1]]^2) # reconstruction error
# for 2nd view
s2_hat <- scifa_res$shat[2]
bic_M2 <- scifa_res$chosen_by_BIC$U0 %*% t(s2_hat*scifa_res$chosen_by_BIC$Vd[[2]]) + scifa_res$chosen_by_BIC$Ud[[2]] %*% t(s2_hat*scifa_res$chosen_by_BIC$Ad[[2]])
sum((dat$M_list[[2]] - bic_M2)^2) / sum(dat$M_list[[2]]^2) # reconstruction error

```

## 3. Example with mixed view types

This example reuses the parameter values from Section 2.1 and the parallel backend registered in Section 2.2. The first view is non-compositional and the second is CLR-transformed. Pass the same indicator vector to both data generation and model fitting.

```r
clr_mixed <- c(FALSE, TRUE)

set.seed(seednum)
dat_mixed <- dgen_fn(
  n, pvec, clr_mixed, sgv_joint, sgv_indiv,
  zp_joint, zp_indiv, snr_vec, gen_mu = FALSE
)

fit_mixed <- run_algo_final(
  dat_mixed$X_list, clr_mixed, pvec, rhat_vec = NULL,
  nlamb = 80, c_min = 1e-5, c_max = 0.4,
  niter_free = 20, max_iter = 200, eps = 1e-5,
  center = TRUE, seednum = seednum
)
```

For your own data, set `clr_indc` according to the order and type of the views in `X_list`, then pass it to `run_algo_final()`.
