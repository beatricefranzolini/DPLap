################################################################################
## this code sources the required functions and packages 
################################################################################

library(LaplacesDemon)  # version 16.1.6      to sample from finite Dirichlet
library(progress)       # version 1.2.3       to draw the progress bar
#library(rstudioapi)    # version 0.17.1      (called by _main_to_run)
library(ggplot2)        # version 4.0.1       to plot
library(patchwork)      # version 1.3.2       to plot
library(fossil)         # version 0.4.0       to compute rand index
library(salso)          # version 0.3.53      to compute psm and point estimate

this_dir = dirname(normalizePath(sys.frame(1)$ofile))

source(file.path(this_dir, "Slice.R"))        #posterior Slice Sampler MCMC 
source(file.path(this_dir, "SlicenoAtoms.R")) #posterior Slice Sampler MCMC with marginalized atoms
source(file.path(this_dir, "CRPwithAtoms.R")) #posterior marginal MCMC with atom sampling 
source(file.path(this_dir, "CRPnoAtoms.R"))   #posterior marginal MCMC with marginalized atoms 
source(file.path(this_dir, "BGS.R"))          #Block Gibbs sampler
source(file.path(this_dir, "utils.R"))        #Likelihoods eval & atoms sampling 
