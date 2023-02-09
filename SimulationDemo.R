# load data
library(spatstat)
library(lattice)
library(ppmlasso)
library(data.table)

source("ppmlasso_functions.R")
source("New functions.R")
source("DiagnosticFunctions.R")

load("FinalEnv_Fr_4km.RData") # loads env.4km.FR, our grid
env = env.4km.Fr # rename grid
env$X = env$X/1000 # put on km scale
env$Y = env$Y/1000 # put on km scale

alreadyrun = TRUE

if (alreadyrun == FALSE)
{
  
  load("Simulated Lucanus points.RData")
  sp.dat = data.frame(X = simX, Y = simY) # coordinates as data frame
 
  N = dim(sp.dat)[1]
  
  form = ~ poly(temp, prec, degree = 2, raw = TRUE) + perch10 + perc311 + perc312 + 
    percagri + log(hp)
  
  n_var = 11 # for creating matrix dimensions for simulations (10 variables plus 1 intercept)
  
  n.sim = 1000
  
  subset_sizes = c(1000, 500, 200, 100, 50)
  
  beta.sims = array(NA, dim = c(length(subset_sizes), n.sim, n_var))
  s.means.sims = s.sds.sims = array(NA, dim = c(length(subset_sizes), n.sim, (n_var - 1))) # n_var minus the "intercept" (-1)
  dimnames(beta.sims)[[1]] = subset_sizes
  dimnames(beta.sims)[[3]] = names(fitN$beta)
  
  data_file = "LucanusData.RData"
  progress_file = "LucanusProgress.RData"
  
  for (subsetsize in 1:length(subset_sizes))
  {
    n.pts = subset_sizes[subsetsize]
    for (sim in 1:n.sim)
    {
      set.seed(sim)
      pts_keep = sample(1:dim(sp.dat)[1], n.pts)
      sp.sim = sp.dat[pts_keep,]
      data.po  = ppmdat(sp.xy = sp.sim, sp.scale = 4, back.xy = env, coord = c("X","Y"), sp.file = NA, quad.file = NA, file.name = "TestPPM") #sp.scale is the resolution of my pixel (here 4 km)
      
      fit.sim = ppmlasso(form, sp.xy = sp.sim, env.grid = env, sp.scale = 4, data = data.po, n.fits = 200, criterion = "bic", max.it = 10, lamb = NA) 
      
      beta.sims[subsetsize,sim,] = fit.sim$beta
      s.means.sims[subsetsize,sim,] = fit.sim$s.means
      s.sds.sims[subsetsize,sim,] = fit.sim$s.sds
      save(beta.sims, s.means.sims, s.sds.sims, file = data_file) 
      save(subsetsize, sim, file = progress_file)
      cat(paste(subsetsize, sim, "\n"))
      flush.console()
    }
  }
  data.po  = ppmdat(sp.xy = sp.dat, sp.scale = 4, back.xy = env, coord = c("X","Y"), sp.file = NA, quad.file = NA, file.name = "TestPPM") #sp.scale is the resolution of my pixel (here 4 km)
  
  fitN = ppmlasso(form, sp.xy = sp.dat, env.grid = env, sp.scale = 4, data = data.po, n.fits = 200, criterion = "bic", max.it = 10, lamb = NA)
  
}
if (alreadyrun == TRUE)
{
  load("LucanusDatalocal.RData")
  data.po  = ppmdat(sp.xy = sp.dat, sp.scale = 4, back.xy = env, coord = c("X","Y"), sp.file = NA, quad.file = NA, file.name = "TestPPM") #sp.scale is the resolution of my pixel (here 4 km)
  
  fitN = ppmlasso(form, sp.xy = sp.dat, env.grid = env, sp.scale = 4, data = data.po, n.fits = 200, criterion = "bic", max.it = 10, lamb = NA)
  
  dimnames(beta.sims)[[1]] = subset_sizes # our sharing code will do this already, but we are doing this after the fact
  dimnames(beta.sims)[[3]] = names(fitN$beta) # here, too
}


########## Intensity measures

# view predicted intensity for model using all points

pred.data = env
bias.var = c("hp")
bias.col = match(bias.var, names(pred.data))

for (v in bias.col)
{
  pred.data[,v] = mean(pred.data[,v])
}

correctedmu = predict.ppmlasso(fitN, newdata = pred.data)

library(lattice)
load("ColorScheme.RData")
levelplot(sqrt(correctedmu) ~ env$X + env$Y, cuts = 100, asp = "iso", 
          col.regions = col2, xlab = "", ylab = "")


all_mus = compute_intensity(fitN, beta.sims, s.means.sims, s.sds.sims, pred.data)
# all_mus calculates the intensity surfaces and rescaled intensity surfaces for the subset models

enveffect = ZeroEnvEffect(all_mus)

f = identity

zmax = rep(NA, 5)
for (i in 1:5)
{
  zmax[i] = max(apply(all_mus$rescale_mus[i,,], 2, mean))
}
zmax = c(zmax, max(correctedmu))

zplot = max(zmax)

avg_mu_plot(all_mus, 50, env[,1:2], mu.min = 1.e-5, intensity = "rescaled",
            col.regions = col2, main = "n = 50", xlab = "", ylab = "", cuts = 100, f = f, at = seq(0, zplot, length.out = 99))

avg_mu_plot(all_mus, 100, env[,1:2], mu.min = 1.e-5, intensity = "rescaled",
            col.regions = col2, main = "n = 100", xlab = "", ylab = "", cuts = 100, f = f, at = seq(0, zplot, length.out = 99))

avg_mu_plot(all_mus, 200, env[,1:2], mu.min = 1.e-5, intensity = "rescaled",
            col.regions = col2, main = "n = 200", xlab = "", ylab = "", cuts = 100, f = f, at = seq(0, zplot, length.out = 99))

avg_mu_plot(all_mus, 500, env[,1:2], mu.min = 1.e-5, intensity = "rescaled",
            col.regions = col2, main = "n = 500", xlab = "", ylab = "", cuts = 100, f = f, at = seq(0, zplot, length.out = 99))

avg_mu_plot(all_mus, 1000, env[,1:2], mu.min = 1.e-5, intensity = "rescaled",
            col.regions = col2, main = "n = 1000", xlab = "", ylab = "", cuts = 100, f = f, at = seq(0, zplot, length.out = 99))

levelplot(correctedmu ~ env$X + env$Y, cuts = 100, asp = "iso", 
          col.regions = col2, main = "n = 2576", xlab = "", ylab = "", f = f, at = seq(0, zplot, length.out = 99)) # all 2577 points



sd_plot(all_mus, 500, env[,1:2], col.regions = col2, main = "n = 500", xlab = "", ylab = "", cuts = 100)



# compare maps from different random subsets of size n = 50 with n = 500

#n = 50

levelplot(all_mus$rescale_mus[5, 1 ,] ~ env$X + env$Y, cuts = 100, xlab = "", 
          ylab = "", asp = "iso", cex = 1.8, col.regions = col2)
levelplot(all_mus$rescale_mus[5, 9 ,] ~ env$X + env$Y, cuts = 100, xlab = "", 
          ylab = "", asp = "iso", cex = 1.8, col.regions = col2)
levelplot(all_mus$rescale_mus[5, 10 ,] ~ env$X + env$Y, cuts = 100, xlab = "", 
          ylab = "", asp = "iso", cex = 1.8, col.regions = col2)

#n = 500

levelplot(all_mus$rescale_mus[2, 1 ,] ~ env$X + env$Y, cuts = 100, xlab = "", 
          ylab = "", asp = "iso", cex = 1.8, col.regions = col2)
levelplot(all_mus$rescale_mus[2, 2 ,] ~ env$X + env$Y, cuts = 100, xlab = "", 
          ylab = "", asp = "iso", cex = 1.8, col.regions = col2)
levelplot(all_mus$rescale_mus[2, 3 ,] ~ env$X + env$Y, cuts = 100, xlab = "", 
          ylab = "", asp = "iso", cex = 1.8, col.regions = col2)

# compute and plot IMSE

mu.min = 1.e-5 # truncate minimum intensity
muN = correctedmu
mean(correctedmu < mu.min)
muN[muN < mu.min] = mu.min 

plotIMSE = IMSE_plot(muN, all_mus, mu.min = 1.e-5, f = "log")

# compute and plot correlation

pearson = Corr_plot(muN, all_mus, mu.min = 1.e-5)
pearson = Corr_plot(muN, all_mus, mu.min = 1.e-5, f_intensity = log)
sum(pearson$Corrs[1,] < 0.9) # subset 1000

spearman = Corr_plot(muN, all_mus, mu.min = 1.e-5, method = "spearman")
kendall = Corr_plot(muN, all_mus, mu.min = 1.e-5, method = "kendall")

########## Fitted coefficient measures

mean.betas = data.frame(apply(beta.sims[1,,], 2, mean),
                        apply(beta.sims[2,,], 2, mean),
                        apply(beta.sims[3,,], 2, mean),
                        apply(beta.sims[4,,], 2, mean),
                        apply(beta.sims[5,,], 2, mean))

names(mean.betas) = c("1000", "500", "200", "100", "50")

sd.betas = data.frame(apply(beta.sims[1,,], 2, sd),
                      apply(beta.sims[2,,], 2, sd),
                      apply(beta.sims[3,,], 2, sd),
                      apply(beta.sims[4,,], 2, sd),
                      apply(beta.sims[5,,], 2, sd))

names(sd.betas) = c("1000", "500", "200", "100", "50")

# Preliminary exploration

miny = min(sd.betas)
maxy = max(sd.betas)

# Plot of standard deviation of fitted coefficients

scales = c(1000, 500, 200, 100, 50)
plot(scales, sd.betas[1,], log = "x", type = "l", ylim = c(miny, maxy))
for (i in 2:dim(mean.betas)[1])
{
  points(scales, sd.betas[i,], type = "l")
}

apply(sd.betas, 2, mean) # Table 2

# Plot of mean of fitted coefficients

miny = min(mean.betas)
maxy = max(mean.betas)

scales = c(N, 1000, 500, 200, 100, 50)
plot(scales[-1], mean.betas[1,], log = "x", type = "l", ylim = c(miny, maxy))
for (i in 2:dim(mean.betas)[1])
{
  points(scales[-1], mean.betas[i,], type = "l")
}

# using temperature, the 2nd covariate, as an example:

coef_plot(2, fitN, beta.sims)
coef_plot(3, fitN, beta.sims)
coef_plot(4, fitN, beta.sims)
coef_plot(5, fitN, beta.sims)
coef_plot(6, fitN, beta.sims)
coef_plot(7, fitN, beta.sims)
coef_plot(8, fitN, beta.sims)
coef_plot(9, fitN, beta.sims)
coef_plot(10, fitN, beta.sims)
coef_plot(11, fitN, beta.sims)

# Signs of coefficients

signs = signcoefs(beta.sims)





