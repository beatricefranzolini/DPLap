# Franzolini & Pozza (2026) - Laplace and skew-Laplace sampling for Dirichlet process mixtures posterior density
# This code reproduces the results of Simulation Study 4.
# IMPORTANT: run this script directly; do not source it.

rm(list = ls())

library(rstudioapi) # version 0.17.1
library(ggplot2)    # version 4.0.1 

#set working directory to Source file directory:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("dplap.R")        #our code
source("DPalg/DPalg.R")  #code from Franzolini & Gaffi (2026)


################################################################################
#Simulation setting 4 (100 components + student t)
################################################################################
set.seed(1)
K_0 = 100
mus = rnorm(100, 0, 1.5)                      # components' centers
sig = 1                                     # components  scale
W   = (1:100)^(-2)
prs = W / sum(W)

#true density of a grid 
xgrid = seq(min(mus)-3*sig, max(mus)+3*sig, length.out = 400)
dx    = diff(xgrid)[1]
truth = colSums(t(dnorm(xgrid, matrix(rep(mus, length(xgrid)), 
                                      ncol = length(mus), byrow = TRUE)))*prs)
truth_dens_norm = truth / (sum(truth) * dx)

sample_sizes = c(20, 50, 100, 200, 500, 1000, 1500, 2000)
num_simul    = length(sample_sizes)  

data = vector(mode = "list", length = num_simul)

LAP      = vector(mode = "list", length = num_simul)
Skew_LAP = vector(mode = "list", length = num_simul)
SS       = vector(mode = "list", length = num_simul)

MCMC_draws = 10000

par(mfrow=c(2,4))

i_count = 0
for(n in sample_sizes){
  set.seed(n)
  i_count = i_count + 1
  y = sample(mus, size = n, prob = prs, replace = TRUE) + rt(n, 5) 
  hist(y, # histogram
       border="gray",
       prob = TRUE, # show densities instead of frequencies
       xlim = c(min(y),max(y)),
       main = paste("n =", n), 
       xlab = "",
       breaks = 100)
  lines(xgrid, truth, lwd = 2, col = 1)
  
  data[[i_count]] = y
  
  LAP[[i_count]]      = dp_lap_mixmodel_normal_normal(y, K = 20, n_samp = 2000, 
                                                      b = 3*log(length(y)), 
                                                      seed = 0)  
  message(paste("lap for n =", n, "completed in", 
                as.numeric(LAP[[i_count]]$time$total, units = "secs"), "secs"))
  
  Skew_LAP[[i_count]] = dp_lap_mixmodel_normal_normal(y, K = 20, skew = TRUE, 
                                                      n_samp = 2000, 
                                                      b = 3*log(length(y)), 
                                                      seed = 0) 
  message(paste("skew_lap for n =", n, "completed in", 
                as.numeric(Skew_LAP[[i_count]]$time$total, units = "secs"), "secs"))
  
  SS[[i_count]] = 
    dp_slice_mixmodel_normal_normal(y, Tot = MCMC_draws, b = 3*log(length(y)), 
                                    seed = 0)
  message(paste("ss for n =", n, "completed in", 
                as.numeric(SS[[i_count]]$time, units = "secs"), "secs"))
}

#compute density estimations 
i_count = 0
for(n in sample_sizes){
  i_count = i_count + 1
  LAP[[i_count]]$density = posterior_density_dpm(xgrid = xgrid, 
                                                 draws = LAP[[i_count]]$samples, 
                                                 K     = 20,
                                                 sigma = 1)
  
  Skew_LAP[[i_count]]$density = posterior_density_dpm(xgrid =xgrid, 
                                                      draws = Skew_LAP[[i_count]]$samples, 
                                                      K     = 20,
                                                      sigma = 1)
  
  SS[[i_count]]$density = estimate_dpmm_density_ss(fit    = SS[[i_count]],
                                                   grid   = xgrid,
                                                   sigma2 = 1,     # must match hyper[1]
                                                   burn   = 2000)
  
}

make_time_summary_table <- function(LAP, Skew_LAP, SS, digits = 4) {
  n = length(LAP)
  
  tab = data.frame(
    sample_size    = sample_sizes,
    
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

time_summary = make_time_summary_table(LAP, Skew_LAP, SS)
print(time_summary)

#compute point estimate metrics 
i_count = 0
for(n in sample_sizes){
  i_count = i_count + 1
  
  #compute TV from est to truth
  f1 = LAP[[i_count]]$density$mean_density
  f1 = f1 / (sum(f1) * dx)
  TV = 0.5 * sum(abs(f1 - truth_dens_norm)) * dx
  LAP[[i_count]]$TV = TV
  #message(paste("TV(Lap, Truth) =", TV))
  
  f1 = Skew_LAP[[i_count]]$density$mean_density
  f1 = f1 / (sum(f1) * dx)
  TV = 0.5 * sum(abs(f1 - truth_dens_norm)) * dx
  Skew_LAP[[i_count]]$TV = TV
  #message(paste("TV(Skew_LAP, Truth) =", TV))
  
  f1 = SS[[i_count]]$density$mean_density 
  f1 = f1 / (sum(f1) * dx)
  TV = 0.5 * sum(abs(f1 - truth_dens_norm)) * dx
  SS[[i_count]]$TV = TV
  #message(paste("TV(SS, Truth) =", TV))
  
  
  #compute TV from est to SS (mean)
  f1 = LAP[[i_count]]$density$mean_density
  f1 = f1 / (sum(f1) * dx)
  TV = 0.5 * sum(abs(f1 - SS[[i_count]]$density$mean_density )) * dx
  LAP[[i_count]]$TVSS = TV
  message(paste("TV(Lap, SS) =", TV))
  
  f1 = Skew_LAP[[i_count]]$density$mean_density
  f1 = f1 / (sum(f1) * dx)
  TV = 0.5 * sum(abs(f1 - SS[[i_count]]$density$mean_density )) * dx
  Skew_LAP[[i_count]]$TVSS = TV
  message(paste("TV(Skew_LAP, SS) =", TV))
  
}

make_tv_summary_table <- function(LAP, Skew_LAP, SS, digits = 4) {
  
  tab = data.frame(
    sample_size    = sample_sizes,
    LAP_Truth      = vapply(LAP,      function(x) x$TV,   numeric(1)),
    Skew_LAP_Truth = vapply(Skew_LAP, function(x) x$TV,   numeric(1)),
    SS_Truth       = vapply(SS,       function(x) x$TV,   numeric(1)),
    LAP_SS         = vapply(LAP,      function(x) x$TVSS, numeric(1)),
    Skew_LAP_SS    = vapply(Skew_LAP, function(x) x$TVSS, numeric(1))
  )
  
  tab[-1] = lapply(tab[-1], round, digits)
  tab
}

tv_summary = make_tv_summary_table(LAP, Skew_LAP, SS)


#compute all posterior metrics 
i_count = 0
for(n in sample_sizes){
  i_count = i_count + 1
  LAP[[i_count]]$TV_grid_vs_SS      = compute_pointwise_tv_vs_ss(
    LAP[[i_count]], SS[[i_count]])
  Skew_LAP[[i_count]]$TV_grid_vs_SS = compute_pointwise_tv_vs_ss(
    Skew_LAP[[i_count]], SS[[i_count]])
}

tv_box_df = make_tv_boxplot_df(LAP, Skew_LAP, sample_sizes, xgrid)

ggplot(tv_box_df, aes(x = method, y = TV, fill = method)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  facet_wrap(~ n, nrow = 2, ncol = 4, scales = "free_y") +
  labs(
    x = NULL,
    y = "TV distance to Slice posterior",
    title = "Density ordinates TV to Slice",
  ) +
  theme_bw() +
  theme(legend.position = "none")+ coord_cartesian(ylim = c(0, 1))

print(time_summary)
print(tv_summary)
1 - tv_summary[6]/tv_summary[5]


# i_count = 0
# for (n in sample_sizes) {
#   i_count = i_count + 1
#   y       = data[[i_count]]
#   
#   hist(
#     y,
#     breaks = 20,
#     prob = TRUE,
#     border = "white",
#     col = "grey70",
#     xlim = range(xgrid),
#     main = n,
#     xlab = "",
#     ylab = "Density",
#     cex.main = 1.2,
#     cex.axis = 0.9,
#     cex.lab = 1,
#     las = 1
#   )
#   
#   lines(xgrid, SS[[i_count]]$density$mean_density,
#         lwd = 2, col = "black", lty = 2)
#   lines(xgrid, LAP[[i_count]]$density$mean_density,
#         lwd = 2, col = "tomato2")
#   lines(xgrid, Skew_LAP[[i_count]]$density$mean_density,
#         lwd = 2, col = "turquoise4")
#   
#   if (i_count == 1) {
#     legend(
#       "topleft",
#       legend = c("Slice sampler", "Laplace", "Skew-Laplace"),
#       col = c("black", "tomato2", "turquoise4"),
#       lwd = 2,
#       lty = c(2, 1, 1),
#       bty = "n",
#       cex = 0.85
#     )
#   }
# }

