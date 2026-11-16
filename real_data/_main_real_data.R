# Franzolini & Pozza (2026) - Laplace and skew-Laplace sampling for Dirichlet process mixtures posterior density

rm(list = ls())

library(rstudioapi)     # version 0.17.1
library(ggplot2)        # version 4.0.1 
library(optimx)         #
library(LaplacesDemon)  # version 16.1.6      
library(progress)       # version 1.2.3      
library(patchwork)      # version 1.3.2       
library(fossil)         # version 0.4.0       
library(salso)          # version 0.3.53  

#set working directory to Source file directory:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../utils_dplap.R")  #our code
source("../dplap.R")        
source("../dpvb.R")
source("../plot_top2_atoms.R")
source("../utils_slice.R")  #code from Franzolini & Gaffi, ICML, 2026
source("../dpslice.R")    
################################################################################
#real data
################################################################################
par(mfrow=c(1,5))
data(faithful, package = "datasets")
hist(faithful$eruptions, breaks = 30)
data(galaxies, package = "MASS")
hist(galaxies, breaks = 30)
data(iris, package = "datasets")
hist(iris$Petal.Length, breaks = 30)
data(rock, package = "datasets")
hist(rock$peri, breaks = 30)
data(quakes, package = "datasets")
hist(quakes$depth, breaks = 30)

data = vector(mode = "list", length = 5)

VB       = vector(mode = "list", length = 5)
LAP      = vector(mode = "list", length = 5)
Skew_LAP = vector(mode = "list", length = 5)
SS       = vector(mode = "list", length = 5)

MCMC_draws = 100000
Burn_in    = 20000
MC_samples = 2000

set.seed(0)

i_count = 1
y = scale(faithful$eruptions)[,1]
xgrid = seq(min(y)-0.5, max(y)+0.5, length.out = 400)
dx    = diff(xgrid)[1]
data[[i_count]] = y

VB[[i_count]] = dp_vb_mixmodel_normal_normal(
  y,
  K = floor(5*log(length(y))),
  n_samp = MC_samples,
  sigma = 0.5, s0 = 0.25, b = 3,
  skew_mu = FALSE,
  seed = 0
)
message(paste("VB for faithful completed in", 
              as.numeric(VB[[i_count]]$time$total, units = "secs"), "secs"))

LAP[[i_count]]      = dp_lap_mixmodel_normal_normal(y, K = floor(5*log(length(y))), n_samp = MC_samples, 
                                                    sigma = 0.5, s0 = 0.25, b = 3, seed = 0)  
message(paste("lap for faithful completed in", 
              as.numeric(LAP[[i_count]]$time$total, units = "secs"), "secs"))

Skew_LAP[[i_count]] = dp_lap_mixmodel_normal_normal(y, K = floor(5*log(length(y))), skew = TRUE, n_samp = MC_samples, 
                                                    sigma = 0.5, s0 = 0.25, b = 3, seed = 0)  
message(paste("skew_lap for faithful completed in", 
              as.numeric(Skew_LAP[[i_count]]$time$total, units = "secs"), "secs"))

SS[[i_count]] = 
  dp_slice_mixmodel_normal_normal(y, Tot = MCMC_draws, b = 3, seed = 0, hyper = c(0.25,0,0.0625))
message(paste("ss for faithful completed in", 
              as.numeric(SS[[i_count]]$time, units = "secs"), "secs"))

#compute density estimations 
VB[[i_count]]$density = posterior_density_dpm_vb(xgrid = xgrid,
                                                 fit   = VB[[i_count]],
                                                 sigma = .5)

LAP[[i_count]]$density = posterior_density_dpm(xgrid = xgrid, 
                                               draws = LAP[[i_count]]$samples, 
                                               sigma = .5)

Skew_LAP[[i_count]]$density = posterior_density_dpm(xgrid = xgrid, 
                                                    draws = Skew_LAP[[i_count]]$samples, 
                                                    sigma = .5)

SS[[i_count]]$density = estimate_dpmm_density_ss(fit    = SS[[i_count]],
                                                 grid   = xgrid,
                                                 sigma2 = .25,     # must match hyper[1]
                                                 burn   = Burn_in,
                                                 mu0    = 0,
                                                 tau20  = 0.0625)

#point estimate TV 


#compute TV from est to SS (mean)
fSS = SS[[i_count]]$density$mean_density
fSS = fSS / (sum(fSS) * dx)

f1 = VB[[i_count]]$density$mean_density
f1 = f1 / (sum(f1) * dx)
TV = 0.5 * sum(abs(f1 - fSS )) * dx
VB[[i_count]]$TVSS = TV
message(paste("TV(VB, SS) =", TV))

f1 = LAP[[i_count]]$density$mean_density
f1 = f1 / (sum(f1) * dx)
TV = 0.5 * sum(abs(f1 - fSS )) * dx
LAP[[i_count]]$TVSS = TV
message(paste("TV(Lap, SS) =", TV))

f1 = Skew_LAP[[i_count]]$density$mean_density
f1 = f1 / (sum(f1) * dx)
TV = 0.5 * sum(abs(f1 - fSS )) * dx
Skew_LAP[[i_count]]$TVSS = TV
message(paste("TV(Skew_LAP, SS) =", TV))

LAP[[i_count]]$TV_grid_vs_SS      = compute_pointwise_tv_vs_ss(
  LAP[[i_count]], SS[[i_count]])
Skew_LAP[[i_count]]$TV_grid_vs_SS = compute_pointwise_tv_vs_ss(
  Skew_LAP[[i_count]], SS[[i_count]])
VB[[i_count]]$TV_grid_vs_SS  = compute_pointwise_tv_vs_ss(
  VB[[i_count]], SS[[i_count]])

boxplot(LAP[[i_count]]$TV_grid_vs_SS, ylim=c(0,1))
boxplot(Skew_LAP[[i_count]]$TV_grid_vs_SS, ylim=c(0,1))
boxplot(VB[[i_count]]$TV_grid_vs_SS, ylim=c(0,1))

#########################################################################################################################
set.seed(0)
i_count = 2
y = scale(galaxies)[,1]
xgrid = seq(min(y)-0.5, max(y)+0.5, length.out = 400)
dx    = diff(xgrid)[1]
data[[i_count]] = y

VB[[i_count]] = dp_vb_mixmodel_normal_normal(
  y,
  K = floor(5*log(length(y))),
  n_samp = MC_samples,
  sigma = 0.5, s0 = 0.25, b = 3,
  skew_mu = FALSE,
  seed = 0
)
message(paste("VB for galaxies completed in", 
              as.numeric(VB[[i_count]]$time$total, units = "secs"), "secs"))

LAP[[i_count]]      = dp_lap_mixmodel_normal_normal(y, K = floor(5*log(length(y))), n_samp = MC_samples, 
                                                    sigma = 0.5, s0 = 0.25, b = 3, seed = 0)  
message(paste("lap for galaxies completed in", 
              as.numeric(LAP[[i_count]]$time$total, units = "secs"), "secs"))

Skew_LAP[[i_count]] = dp_lap_mixmodel_normal_normal(y, K = floor(5*log(length(y))), skew = TRUE, n_samp = MC_samples, 
                                                    sigma = 0.5, s0 = 0.25, b = 3, seed = 0)  
message(paste("skew_lap for galaxies completed in", 
              as.numeric(Skew_LAP[[i_count]]$time$total, units = "secs"), "secs"))

SS[[i_count]] = 
  dp_slice_mixmodel_normal_normal(y, Tot = MCMC_draws, b = 3, seed = 0, hyper = c(0.25,0,0.0625))
message(paste("ss for galaxies completed in", 
              as.numeric(SS[[i_count]]$time, units = "secs"), "secs"))

#compute density estimations 
VB[[i_count]]$density = posterior_density_dpm_vb(xgrid = xgrid,
                                                 fit   = VB[[i_count]],
                                                 sigma = .5)

LAP[[i_count]]$density = posterior_density_dpm(xgrid = xgrid, 
                                               draws = LAP[[i_count]]$samples, 
                                               sigma = .5)

Skew_LAP[[i_count]]$density = posterior_density_dpm(xgrid = xgrid, 
                                                    draws = Skew_LAP[[i_count]]$samples, 
                                                    sigma = .5)

SS[[i_count]]$density = estimate_dpmm_density_ss(fit    = SS[[i_count]],
                                                 grid   = xgrid,
                                                 sigma2 = .25,     # must match hyper[1]
                                                 burn   = Burn_in,
                                                 mu0    = 0,
                                                 tau20  = 0.0625)

#point estimate TV 



#compute TV from est to SS (mean)
fSS = SS[[i_count]]$density$mean_density
fSS = fSS / (sum(fSS) * dx)

f1 = VB[[i_count]]$density$mean_density
f1 = f1 / (sum(f1) * dx)
TV = 0.5 * sum(abs(f1 - fSS )) * dx
VB[[i_count]]$TVSS = TV
message(paste("TV(VB, SS) =", TV))

f1 = LAP[[i_count]]$density$mean_density
f1 = f1 / (sum(f1) * dx)
TV = 0.5 * sum(abs(f1 - fSS )) * dx
LAP[[i_count]]$TVSS = TV
message(paste("TV(Lap, SS) =", TV))

f1 = Skew_LAP[[i_count]]$density$mean_density
f1 = f1 / (sum(f1) * dx)
TV = 0.5 * sum(abs(f1 - fSS )) * dx
Skew_LAP[[i_count]]$TVSS = TV
message(paste("TV(Skew_LAP, SS) =", TV))

LAP[[i_count]]$TV_grid_vs_SS      = compute_pointwise_tv_vs_ss(
  LAP[[i_count]], SS[[i_count]])
Skew_LAP[[i_count]]$TV_grid_vs_SS = compute_pointwise_tv_vs_ss(
  Skew_LAP[[i_count]], SS[[i_count]])
VB[[i_count]]$TV_grid_vs_SS  = compute_pointwise_tv_vs_ss(
  VB[[i_count]], SS[[i_count]])

boxplot(LAP[[i_count]]$TV_grid_vs_SS, ylim=c(0,1))
boxplot(Skew_LAP[[i_count]]$TV_grid_vs_SS, ylim=c(0,1))
boxplot(VB[[i_count]]$TV_grid_vs_SS, ylim=c(0,1))

############################################################################################################################
set.seed(0)
i_count = 3
y = scale(iris$Petal.Length)[,1]
xgrid = seq(min(y)-0.5, max(y)+0.5, length.out = 400)
dx    = diff(xgrid)[1]
data[[i_count]] = y

VB[[i_count]] = dp_vb_mixmodel_normal_normal(
  y,
  K = floor(5*log(length(y))),
  n_samp = MC_samples,
  sigma = 0.5, s0 = 0.25, b = 3,
  skew_mu = FALSE,
  seed = 0
)
message(paste("VB for iris completed in", 
              as.numeric(VB[[i_count]]$time$total, units = "secs"), "secs"))

LAP[[i_count]]      = dp_lap_mixmodel_normal_normal(y, K = floor(5*log(length(y))), n_samp = MC_samples, 
                                                    sigma = 0.5, s0 = 0.25, b = 3, seed = 0)  
message(paste("lap for iris completed in", 
              as.numeric(LAP[[i_count]]$time$total, units = "secs"), "secs"))

Skew_LAP[[i_count]] = dp_lap_mixmodel_normal_normal(y, K = floor(5*log(length(y))), skew = TRUE, n_samp = MC_samples, 
                                                    sigma = 0.5, s0 = 0.25, b = 3, seed = 0)  
message(paste("skew_lap for iris completed in", 
              as.numeric(Skew_LAP[[i_count]]$time$total, units = "secs"), "secs"))

SS[[i_count]] = 
  dp_slice_mixmodel_normal_normal(y, Tot = MCMC_draws, b = 3, seed = 0, hyper = c(0.25,0,0.0625))
message(paste("ss for iris completed in", 
              as.numeric(SS[[i_count]]$time, units = "secs"), "secs"))

#compute density estimations 
VB[[i_count]]$density = posterior_density_dpm_vb(xgrid = xgrid,
                                                 fit   = VB[[i_count]],
                                                 sigma = .5)

LAP[[i_count]]$density = posterior_density_dpm(xgrid = xgrid, 
                                               draws = LAP[[i_count]]$samples, 
                                               sigma = .5)

Skew_LAP[[i_count]]$density = posterior_density_dpm(xgrid = xgrid, 
                                                    draws = Skew_LAP[[i_count]]$samples, 
                                                    sigma = .5)

SS[[i_count]]$density = estimate_dpmm_density_ss(fit    = SS[[i_count]],
                                                 grid   = xgrid,
                                                 sigma2 = .25,     # must match hyper[1]
                                                 burn   = Burn_in,
                                                 mu0    = 0,
                                                 tau20  = 0.0625)
#point estimate TV 



#compute TV from est to SS (mean)
fSS = SS[[i_count]]$density$mean_density
fSS = fSS / (sum(fSS) * dx)

f1 = VB[[i_count]]$density$mean_density
f1 = f1 / (sum(f1) * dx)
TV = 0.5 * sum(abs(f1 - fSS )) * dx
VB[[i_count]]$TVSS = TV
message(paste("TV(VB, SS) =", TV))

f1 = LAP[[i_count]]$density$mean_density
f1 = f1 / (sum(f1) * dx)
TV = 0.5 * sum(abs(f1 - fSS )) * dx
LAP[[i_count]]$TVSS = TV
message(paste("TV(Lap, SS) =", TV))

f1 = Skew_LAP[[i_count]]$density$mean_density
f1 = f1 / (sum(f1) * dx)
TV = 0.5 * sum(abs(f1 - fSS )) * dx
Skew_LAP[[i_count]]$TVSS = TV
message(paste("TV(Skew_LAP, SS) =", TV))


LAP[[i_count]]$TV_grid_vs_SS      = compute_pointwise_tv_vs_ss(
  LAP[[i_count]], SS[[i_count]])
Skew_LAP[[i_count]]$TV_grid_vs_SS = compute_pointwise_tv_vs_ss(
  Skew_LAP[[i_count]], SS[[i_count]])
VB[[i_count]]$TV_grid_vs_SS  = compute_pointwise_tv_vs_ss(
  VB[[i_count]], SS[[i_count]])

boxplot(LAP[[i_count]]$TV_grid_vs_SS, ylim=c(0,1))
boxplot(Skew_LAP[[i_count]]$TV_grid_vs_SS, ylim=c(0,1))
boxplot(VB[[i_count]]$TV_grid_vs_SS, ylim=c(0,1))

#################################################################################################################
set.seed(0)

i_count = 4
y = scale(rock$peri)[,1]
xgrid = seq(min(y)-0.5, max(y)+0.5, length.out = 400)
dx    = diff(xgrid)[1]
data[[i_count]] = y

VB[[i_count]] = dp_vb_mixmodel_normal_normal(
  y,
  K = floor(5*log(length(y))),
  n_samp = MC_samples,
  sigma = 0.5, s0 = 0.25, b = 3,
  skew_mu = FALSE,
  seed = 0
)
message(paste("VB for rock completed in", 
              as.numeric(VB[[i_count]]$time$total, units = "secs"), "secs"))

LAP[[i_count]]      = dp_lap_mixmodel_normal_normal(y, K = floor(5*log(length(y))), n_samp = MC_samples, 
                                                    sigma = 0.5, s0 = 0.25, b = 3, seed = 0)  
message(paste("lap for rock completed in", 
              as.numeric(LAP[[i_count]]$time$total, units = "secs"), "secs"))

Skew_LAP[[i_count]] = dp_lap_mixmodel_normal_normal(y, K = floor(5*log(length(y))), skew = TRUE, n_samp = MC_samples, 
                                                    sigma = 0.5, s0 = 0.25, b = 3, seed = 0)  
message(paste("skew_lap for rock completed in", 
              as.numeric(Skew_LAP[[i_count]]$time$total, units = "secs"), "secs"))

SS[[i_count]] = 
  dp_slice_mixmodel_normal_normal(y, Tot = MCMC_draws, b = 3, seed = 0, hyper = c(0.25,0,0.0625))
message(paste("ss for rock completed in", 
              as.numeric(SS[[i_count]]$time, units = "secs"), "secs"))

#compute density estimations 
VB[[i_count]]$density = posterior_density_dpm_vb(xgrid = xgrid,
                                                 fit   = VB[[i_count]],
                                                 sigma = .5)

LAP[[i_count]]$density = posterior_density_dpm(xgrid = xgrid, 
                                               draws = LAP[[i_count]]$samples, 
                                               sigma = .5)

Skew_LAP[[i_count]]$density = posterior_density_dpm(xgrid = xgrid, 
                                                    draws = Skew_LAP[[i_count]]$samples, 
                                                    sigma = .5)

SS[[i_count]]$density = estimate_dpmm_density_ss(fit    = SS[[i_count]],
                                                 grid   = xgrid,
                                                 sigma2 = .25,     # must match hyper[1]
                                                 burn   = Burn_in,
                                                 mu0    = 0,
                                                 tau20  = 0.0625)


#point estimate TV 
#compute TV from est to SS (mean)
fSS = SS[[i_count]]$density$mean_density
fSS = fSS / (sum(fSS) * dx)

f1 = VB[[i_count]]$density$mean_density
f1 = f1 / (sum(f1) * dx)
TV = 0.5 * sum(abs(f1 - fSS )) * dx
VB[[i_count]]$TVSS = TV
message(paste("TV(VB, SS) =", TV))

f1 = LAP[[i_count]]$density$mean_density
f1 = f1 / (sum(f1) * dx)
TV = 0.5 * sum(abs(f1 - fSS )) * dx
LAP[[i_count]]$TVSS = TV
message(paste("TV(Lap, SS) =", TV))

f1 = Skew_LAP[[i_count]]$density$mean_density
f1 = f1 / (sum(f1) * dx)
TV = 0.5 * sum(abs(f1 - fSS )) * dx
Skew_LAP[[i_count]]$TVSS = TV
message(paste("TV(Skew_LAP, SS) =", TV))


LAP[[i_count]]$TV_grid_vs_SS      = compute_pointwise_tv_vs_ss(
  LAP[[i_count]], SS[[i_count]])
Skew_LAP[[i_count]]$TV_grid_vs_SS = compute_pointwise_tv_vs_ss(
  Skew_LAP[[i_count]], SS[[i_count]])
VB[[i_count]]$TV_grid_vs_SS  = compute_pointwise_tv_vs_ss(
  VB[[i_count]], SS[[i_count]])

boxplot(LAP[[i_count]]$TV_grid_vs_SS, ylim=c(0,1))
boxplot(Skew_LAP[[i_count]]$TV_grid_vs_SS, ylim=c(0,1))
boxplot(VB[[i_count]]$TV_grid_vs_SS, ylim=c(0,1))

#########################################################################################################################
set.seed(0)
i_count = 5

y = scale(quakes$depth)[,1]
xgrid = seq(min(y)-0.5, max(y)+0.5, length.out = 400)
dx    = diff(xgrid)[1]
data[[i_count]] = y

VB[[i_count]] = dp_vb_mixmodel_normal_normal(
  y,
  K = floor(5*log(length(y))),
  n_samp = MC_samples,
  sigma = 0.5, s0 = 0.25, b = 3,
  skew_mu = FALSE,
  seed = 0
)
message(paste("VB for quakes completed in", 
              as.numeric(VB[[i_count]]$time$total, units = "secs"), "secs"))

LAP[[i_count]]      = dp_lap_mixmodel_normal_normal(y, K = floor(5*log(length(y))), n_samp = MC_samples, 
                                                    sigma = 0.5, s0 = 0.25, b = 3, seed = 0)  
message(paste("lap for quakes completed in", 
              as.numeric(LAP[[i_count]]$time$total, units = "secs"), "secs"))

Skew_LAP[[i_count]] = dp_lap_mixmodel_normal_normal(y, K = floor(5*log(length(y))), skew = TRUE, n_samp = MC_samples, 
                                                    sigma = 0.5, s0 = 0.25, b = 3, seed = 0)  
message(paste("skew_lap for quakes completed in", 
              as.numeric(Skew_LAP[[i_count]]$time$total, units = "secs"), "secs"))

SS[[i_count]] = 
  dp_slice_mixmodel_normal_normal(y, Tot = MCMC_draws, b = 3, seed = 0, hyper = c(0.25,0,0.0625))
message(paste("ss for quakes completed in", 
              as.numeric(SS[[i_count]]$time, units = "secs"), "secs"))

#compute density estimations 
VB[[i_count]]$density = posterior_density_dpm_vb(xgrid = xgrid,
                                                 fit   = VB[[i_count]],
                                                 sigma = .5)

LAP[[i_count]]$density = posterior_density_dpm(xgrid = xgrid, 
                                               draws = LAP[[i_count]]$samples, 
                                               sigma = .5)

Skew_LAP[[i_count]]$density = posterior_density_dpm(xgrid = xgrid, 
                                                    draws = Skew_LAP[[i_count]]$samples, 
                                                    sigma = .5)

SS[[i_count]]$density = estimate_dpmm_density_ss(fit    = SS[[i_count]],
                                                 grid   = xgrid,
                                                 sigma2 = .25,     # must match hyper[1]
                                                 burn   = Burn_in,
                                                 mu0    = 0,
                                                 tau20  = 0.0625)

#point estimate TV 

#compute TV from est to SS (mean)
fSS = SS[[i_count]]$density$mean_density
fSS = fSS / (sum(fSS) * dx)

f1 = VB[[i_count]]$density$mean_density
f1 = f1 / (sum(f1) * dx)
TV = 0.5 * sum(abs(f1 - fSS )) * dx
VB[[i_count]]$TVSS = TV
message(paste("TV(VB, SS) =", TV))

f1 = LAP[[i_count]]$density$mean_density
f1 = f1 / (sum(f1) * dx)
TV = 0.5 * sum(abs(f1 - fSS )) * dx
LAP[[i_count]]$TVSS = TV
message(paste("TV(Lap, SS) =", TV))

f1 = Skew_LAP[[i_count]]$density$mean_density
f1 = f1 / (sum(f1) * dx)
TV = 0.5 * sum(abs(f1 - fSS )) * dx
Skew_LAP[[i_count]]$TVSS = TV
message(paste("TV(Skew_LAP, SS) =", TV))


LAP[[i_count]]$TV_grid_vs_SS      = compute_pointwise_tv_vs_ss(
  LAP[[i_count]], SS[[i_count]])
Skew_LAP[[i_count]]$TV_grid_vs_SS = compute_pointwise_tv_vs_ss(
  Skew_LAP[[i_count]], SS[[i_count]])
VB[[i_count]]$TV_grid_vs_SS  = compute_pointwise_tv_vs_ss(
  VB[[i_count]], SS[[i_count]])

boxplot(LAP[[i_count]]$TV_grid_vs_SS, ylim=c(0,1))
boxplot(Skew_LAP[[i_count]]$TV_grid_vs_SS, ylim=c(0,1))
boxplot(VB[[i_count]]$TV_grid_vs_SS, ylim=c(0,1))

############################################################################################################################


make_time_summary_table <- function(VB, LAP, Skew_LAP, SS, digits = 4) {
  n = length(LAP)
  
  tab = data.frame(
    dataset    = c("faithful", "galaxy", "iris", "rock", "quakes"),
    
    VB_sec      = vapply(VB, function(x) 
      as.numeric(x$time$total, units = "secs"), 
      numeric(1)),
    
    LAP_sec      = vapply(LAP, function(x) 
      as.numeric(x$time$total, units = "secs"), 
      numeric(1)),
    
    Skew_LAP_sec = vapply(Skew_LAP, function(x) 
      as.numeric(x$time$total, units = "secs"), 
      numeric(1)),
    
    SS_sec       = vapply(SS, function(x) 
      as.numeric(x$time, units = "secs"), 
      numeric(1))
  )
  
  tab[-1] <- lapply(tab[-1], round, digits)
  tab
}

time_summary = make_time_summary_table(VB, LAP, Skew_LAP, SS)
print(time_summary)


make_tv_summary_table <- function(VB, LAP, Skew_LAP, SS, digits = 4) {
  
  tab = data.frame(
    dataset    = c("faithful", "galaxy", "iris", "rock", "quakes"),
    VB_SS          = vapply(VB,       function(x) x$TVSS, numeric(1)),
    LAP_SS         = vapply(LAP,      function(x) x$TVSS, numeric(1)),
    Skew_LAP_SS    = vapply(Skew_LAP, function(x) x$TVSS, numeric(1))
  )
  
  tab[-1] = lapply(tab[-1], round, digits)
  tab
}

tv_summary = make_tv_summary_table(VB, LAP, Skew_LAP, SS)
print(tv_summary)

make_tv_boxplot_df <- function(VB, LAP, Skew_LAP, dataset_names) {
  do.call(rbind, lapply(seq_along(dataset_names), function(i) {
    xgrid = VB[[i]]$density$xgrid
    rbind(
      data.frame(
        i_count = i,
        n = dataset_names[i],
        method = "VB",
        grid_point = seq_along(xgrid),
        x = xgrid,
        TV = VB[[i]]$TV_grid_vs_SS
      ),
      data.frame(
        i_count = i,
        n = dataset_names[i],
        method = "Lap",
        grid_point = seq_along(xgrid),
        x = xgrid,
        TV = LAP[[i]]$TV_grid_vs_SS
      ),
      data.frame(
        i_count = i,
        n = dataset_names[i],
        method = "Skew-Lap",
        grid_point = seq_along(xgrid),
        x = xgrid,
        TV = Skew_LAP[[i]]$TV_grid_vs_SS
      )
    )
  }))
}

dataset_names = c("faithful", "galaxy", "iris", "rock", "quakes")

tv_box_df = make_tv_boxplot_df(VB, LAP, Skew_LAP,
                               dataset_names)
tv_box_df$n = factor(tv_box_df$n, levels = dataset_names)

ggplot(tv_box_df, aes(x = method, y = TV, fill = method)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  facet_wrap(~ n, scales = "free_y", nrow = 1) +
  labs(
    x = NULL,
    y = "TV distance to Slice posterior",
    title = "Density ordinates TV to Slice",
  ) +
  theme_bw() +
  theme(legend.position = "none") + coord_cartesian(ylim = c(0, 0.4))

set.seed(1)
op <- par(
  mfrow = c(1, 5),
  mar = c(4, 4, 3, 1),
  mgp = c(2.2, 0.7, 0),
  tcl = -0.3
)

i_count = 0

for (n in dataset_names) {
  i_count = i_count + 1
  y       = data[[i_count]]
  xgrid   = seq(min(y) - 0.5, max(y) + 0.5, length.out = 400)
  
  hist(
    y,
    breaks = 20,
    prob = TRUE,
    border = "white",
    col = "grey70",
    xlim = range(xgrid),
    main = n,
    xlab = "",
    ylab = "",
    cex.main = 2.2,
    cex.axis = 1.9,
    cex.lab = 1,
    las = 1
  )
  
  lines(xgrid, SS[[i_count]]$density$mean_density,
        lwd = 2, col = "black", lty = 2)
  lines(xgrid, LAP[[i_count]]$density$mean_density,
        lwd = 2, col = "tomato2")
  lines(xgrid, Skew_LAP[[i_count]]$density$mean_density,
        lwd = 2, col = "turquoise4")
  lines(xgrid, VB[[i_count]]$density$mean_density,
        lwd = 2, col = "blue")
  
  if (i_count == 1) {
    legend(
      "topleft",
      legend = c("Slice sampler", "VB", "Laplace", "Skew-Laplace"),
      col = c("black", "blue", "tomato2", "turquoise4"),
      lwd = 2,
      lty = c(2, 1, 1, 1),
      bty = "n",
      cex = 1.85
    )
  }
}
print(time_summary)
print(tv_summary)
