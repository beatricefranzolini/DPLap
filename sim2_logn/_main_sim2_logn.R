# Franzolini & Pozza (2026) - Laplace and skew-Laplace sampling for Dirichlet process mixtures posterior density
# This code reproduces the results of Simulation Study 2.

rm(list = ls())

library(rstudioapi)     # version 0.17.1
library(ggplot2)        # version 4.0.1 
library(optimx)         # version 2025-4.9
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
#Simulation setting 2 (100 components + Normal)
################################################################################
set.seed(0)
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

VB       = vector(mode = "list", length = num_simul)
LAP      = vector(mode = "list", length = num_simul)
Skew_LAP = vector(mode = "list", length = num_simul)
SS       = vector(mode = "list", length = num_simul)

MCMC_draws = 100000
Burn_in    = 20000
MC_samples = 2000

par(mfrow=c(2,4))

seed = 2   #seed for data
i_count = 0
y = NULL
for(n in sample_sizes){
  i_count = i_count + 1
  y = c(y, rnorm(n - length(y),
           mean = sample(mus, size = n - length(y), prob = prs, replace = TRUE),
           sd   = sig))
  
  hist(y, # histogram
       border="gray",
       prob = TRUE, # show densities instead of frequencies
       xlim = range(xgrid),
       main = paste("n =", n), 
       xlab = "",
       breaks = 20)
  
  data[[i_count]] = y
  
  VB[[i_count]] = dp_vb_mixmodel_normal_normal(
    y,
    K = floor(5*log(length(y))),,
    n_samp = MC_samples,
    b = 3 * log(length(y)),
    skew_mu = FALSE,
    seed = 0
  )
  
  message(paste("VB for n =", n, "completed in",
                as.numeric(VB[[i_count]]$time$total, units = "secs"), "secs"))
  
  LAP[[i_count]] = dp_lap_mixmodel_normal_normal(
    y, 
    K = floor(5*log(length(y))), 
    n_samp = MC_samples, 
    b = 3*log(length(y)), 
    seed = 0,
    init_par = as.vector(colMeans(VB[[i_count]]$samples))
  )
  
  message(paste("lap for n =", n, "completed in", 
                as.numeric(LAP[[i_count]]$time$total, units = "secs"), "secs"))
  
  Skew_LAP[[i_count]] = dp_lap_mixmodel_normal_normal(
    y, 
    K = floor(5*log(length(y))),
    skew = TRUE, 
    n_samp = MC_samples, 
    b = 3*log(length(y)), 
    seed = 0,
    init_par = as.vector(colMeans(VB[[i_count]]$samples))
  ) 
  
  message(paste("skew_lap for n =", n, "completed in", 
                as.numeric(Skew_LAP[[i_count]]$time$total, units = "secs"), "secs"))
  
  SS[[i_count]] = dp_slice_mixmodel_normal_normal(
    y, 
    Tot = MCMC_draws, 
    b = 3*log(length(y)), 
    seed = 0
  )
  
  message(paste("ss for n =", n, "completed in", 
                as.numeric(SS[[i_count]]$time, units = "secs"), "secs"))
}

#compute density estimations 
i_count = 0
for(n in sample_sizes){
  i_count = i_count + 1
  VB[[i_count]]$density = posterior_density_dpm_vb(xgrid = xgrid,
                                                   fit   = VB[[i_count]],
                                                   sigma = 1)
  
  LAP[[i_count]]$density = posterior_density_dpm(xgrid = xgrid, 
                                                 draws = LAP[[i_count]]$samples, 
                                                 sigma = 1)
  
  Skew_LAP[[i_count]]$density = posterior_density_dpm(xgrid = xgrid, 
                                                      draws = Skew_LAP[[i_count]]$samples, 
                                                      sigma = 1)
  
  SS[[i_count]]$density = estimate_dpmm_density_ss(fit    = SS[[i_count]],
                                                   grid   = xgrid,
                                                   sigma2 = 1,     # must match hyper[1]
                                                   burn   = Burn_in)
  print(n)
}
#check convergence of LAP and Skew_Lap
# sapply(LAP, function(x) x$optimization$convergence) #should be all zeros 
# sapply(Skew_LAP, function(x) x$optimization$convergence)
# 
# LAP_diagnostics = do.call(
#   rbind,
#   lapply(seq_along(LAP), function(i) {
#     check_lap_fit(
#       fit = LAP[[i]],
#       y = data[[i]],
#       K = 20,
#       sigma = 1,
#       m0 = 0,
#       s0 = 1,
#       a = 3,
#       b = 3 * log(length(data[[i]]))
#     )
#   })
# )
# 
# rownames(LAP_diagnostics) = sample_sizes
# LAP_diagnostics
make_time_summary_table <- function(VB, LAP, Skew_LAP, SS, digits = 4) {
  n = length(LAP)
  
  tab = data.frame(
    sample_size    = sample_sizes,
    
    VB_sec      = vapply(VB, function(x) 
      as.numeric(x$time$total, units = "secs"), 
      numeric(1)),
    
    # CVB_sec      = vapply(CVB, function(x) 
    #   as.numeric(x$time$total, units = "secs"), 
    #   numeric(1)),
    
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

time_summary[-1] = lapply(time_summary[-1] , round, digits = 4)
print(time_summary)

write.csv(time_summary, file ="time.csv")

#compute point estimate metrics 
i_count = 0 
for(n in sample_sizes){
  i_count = i_count + 1
  
  #compute VB from est to truth
  f1 = VB[[i_count]]$density$mean_density
  f1 = f1 / (sum(f1) * dx)
  TV = 0.5 * sum(abs(f1 - truth_dens_norm)) * dx
  VB[[i_count]]$TV = TV
  #message(paste("TV(VB, Truth) =", TV))
  
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
  
  ##############################################
  
  #compute TV from est to SS (mean) Skew_VB
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
  
}

make_tv_summary_table <- function(VB, LAP, Skew_LAP, SS, digits = 4) {
  
  tab = data.frame( 
    sample_size    = sample_sizes,
    VB_Truth       = vapply(VB,       function(x) x$TV,   numeric(1)),
    LAP_Truth      = vapply(LAP,      function(x) x$TV,   numeric(1)),
    Skew_LAP_Truth = vapply(Skew_LAP, function(x) x$TV,   numeric(1)),
    SS_Truth       = vapply(SS,       function(x) x$TV,   numeric(1)),
    VB_SS          = vapply(VB,       function(x) x$TVSS, numeric(1)),
    LAP_SS         = vapply(LAP,      function(x) x$TVSS, numeric(1)),
    Skew_LAP_SS    = vapply(Skew_LAP, function(x) x$TVSS, numeric(1))
  )
  
  tab[-1] = lapply(tab[-1], round, digits)
  tab
}

tv_summary = make_tv_summary_table(VB, LAP, Skew_LAP, SS)
print(tv_summary[,c(1:5)])
print(tv_summary[,c(1,6:8)])

write.csv(tv_summary, file ="tv_point_estimate.csv")

i_count = 0
for(n in sample_sizes){
  i_count = i_count + 1
  VB[[i_count]]$TV_grid_vs_SS      = compute_pointwise_tv_vs_ss(
    VB[[i_count]], SS[[i_count]])
  LAP[[i_count]]$TV_grid_vs_SS      = compute_pointwise_tv_vs_ss(
    LAP[[i_count]], SS[[i_count]])
  Skew_LAP[[i_count]]$TV_grid_vs_SS = compute_pointwise_tv_vs_ss(
    Skew_LAP[[i_count]], SS[[i_count]])
}

make_tv_boxplot_df <- function(VB, LAP, Skew_LAP, sample_sizes, xgrid) {
  do.call(rbind, lapply(seq_along(sample_sizes), function(i) {
    rbind(
      data.frame(
        i_count = i,
        n = sample_sizes[i],
        method = "VB",
        grid_point = seq_along(xgrid),
        x = xgrid,
        TV = VB[[i]]$TV_grid_vs_SS
      ),
      data.frame(
        i_count = i,
        n = sample_sizes[i],
        method = "Lap",
        grid_point = seq_along(xgrid),
        x = xgrid,
        TV = LAP[[i]]$TV_grid_vs_SS
      ),
      data.frame(
        i_count = i,
        n = sample_sizes[i],
        method = "Skew-Lap",
        grid_point = seq_along(xgrid),
        x = xgrid,
        TV = Skew_LAP[[i]]$TV_grid_vs_SS
      )
    )
  }))
}

tv_box_df = make_tv_boxplot_df(VB, LAP, Skew_LAP, sample_sizes, xgrid)

pdf(file = "Boxplots_TVerror_densityordinatesposteriors.pdf", 
    width = 8.25, height = 6) 
p1 = ggplot(tv_box_df, aes(x = method, y = TV, fill = method)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  facet_wrap(~ n, nrow = 2, ncol = 4, scales = "free_y") +
  labs(
    x = NULL,
    y = "TV distance to Slice posterior",
    title = "Density ordinates TV to Slice",
  ) +
  theme_bw() +
  theme(legend.position = "none")+ coord_cartesian(ylim = c(0, 1))
print(p1) 
dev.off()
p1 


thin = floor((MCMC_draws - Burn_in) / MC_samples)
pdf(file = "bivariate_atom_n100.pdf", width = 5, height = 4)  
i <- 3
set.seed(123)

p2 <- plot_top2_atoms(
  fits = list(
    VB = VB[[i]],
    Laplace = LAP[[i]],
    `Skew-Laplace` = Skew_LAP[[i]],
    Slice = SS[[i]]
  ),
  burn_slice = Burn_in,
  thin_slice = thin,
  type = "both"
)
print(p2)
dev.off()
p2

pdf(file = "bivariate_atom_n1500.pdf", width = 5, height = 4) 
i <- 7
set.seed(123)

p3 <- plot_top2_atoms(
  fits = list(
    VB = VB[[i]],
    Laplace = LAP[[i]],
    `Skew-Laplace` = Skew_LAP[[i]],
    Slice = SS[[i]]
  ),
  burn_slice = Burn_in,
  thin_slice = thin,
  type = "both"
)
print(p3)
dev.off()
p3

i_count = 0
for(n in sample_sizes){
  i_count = i_count + 1
  SS[[i_count]]$ESS = ESS(dp_loglik_trace_mixmodel_normal_normal(data[[i_count]], SS[[i_count]]))
}
make_ESS_slice_table <- function(SS, digits = 4) {
  n = length(LAP)
  
  tab = data.frame(
    sample_size    = sample_sizes,
    
    ESS      = vapply(SS, function(x) 
      as.numeric(x$ESS), 
      numeric(1))
  )
  
  tab[-1] <- lapply(tab[-1], round, digits)
  tab
}

ESS_summary = make_ESS_slice_table(SS)

print(ESS_summary)

write.csv(ESS_summary, file ="ESS.csv")

i_count = 0
for (n in sample_sizes) {
  i_count = i_count + 1
  y       = data[[i_count]]
  
  hist(
    y,
    breaks = 20,
    prob = TRUE,
    border = "white",
    col = "grey70",
    xlim = range(xgrid),
    main = n,
    xlab = "",
    ylab = "Density",
    cex.main = 1.2,
    cex.axis = 0.9,
    cex.lab = 1,
    las = 1
  )
  
  lines(xgrid, SS[[i_count]]$density$mean_density,
        lwd = 2, col = "black", lty = 2)
  lines(xgrid, VB[[i_count]]$density$mean_density,
        lwd = 2, col = "blue")
  lines(xgrid, LAP[[i_count]]$density$mean_density,
        lwd = 2, col = "tomato2")
  lines(xgrid, Skew_LAP[[i_count]]$density$mean_density,
        lwd = 2, col = "turquoise4")
  
  
  if (i_count == 1) {
    legend(
      "topleft",
      legend = c("Slice sampler", "Laplace", "Skew-Laplace", "VB"),
      col = c("black", "blue", "tomato2", "turquoise4"),
      lwd = 2,
      lty = c(2, 1, 1, 1),
      bty = "n",
      cex = 0.85
    )
  }
}

save.image("sim2.RData")
