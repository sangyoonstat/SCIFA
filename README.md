# SCIFA 

## 1. File description

**source_scifa.R** - R code for our main method and several auxiliary functions 

**fit_rdat_github.R** - R code to apply SCIFA for the real data

**all_rdat.RData** - file that includes real data we used. The original source can be found https://github.com/bioFAM/MOFA2_tutorials/blob/master/R_tutorials/microbiome_vignette.html and we saved and minimaly edited the original files to be consistent with our usage. 

## 2. Example with two clr-transformed views

### 2.1 Generate data
```{r}

my_path <- "/Users/sayi/Library/CloudStorage/OneDrive-OklahomaAandMSystem/SCIFA-manuscript/" 
source(paste(my_path, sep = "", "source_scifa.R")) # load the source code we have

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
scifa_out <- run_algo(dat$X_list, clr_indc, pvec, rhat_vec = NULL, nlamb = 80, c_min = 10^-5, c_max = 0.4,
                      niter_free = 20, max_iter = 200, eps = 10^-5, center = T, seednum = seednum)

# save the result
scifa_res <- list("chosen_by_AIC" = scifa_out[[1]], "chosen_by_BIC" = scifa_out[[2]],
                  "shat" = scifa_out[[3]], "mean_mat" = scifa_out[[4]], "rhat" = scifa_out[[5]])

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
