IMSE = function(mu1, mu2, f_intensity = log, mu.min = 1.e-5)
{
  mu1.use = mu1
  mu1.use[mu1.use < mu.min] = mu.min
  mu2.use = mu2
  mu2.use[mu2.use < mu.min] = mu.min
  mu1.use = f_intensity(mu1.use)
  mu2.use = f_intensity(mu2.use)
  imse = sum((mu1.use - mu2.use)^2)
  imse
}

dotscale = function(x, f, main = NULL, xlab = NULL, ylab = NULL, col = c(228, 155, 15)/255)
{
  plot(x, f, log = "x", main = main, xlab = xlab, ylab = ylab,
       col = rgb(col[1], col[2], col[3], alpha = 0.2), pch = 19)
  x.unique = sort(unique(x))
  scale.means = rep(NA, length(x.unique))
  for (s in 1:length(x.unique))
  {
    scale.means[s] = mean(f[x == x.unique[s]], na.rm = TRUE)
  }
  points(x.unique, scale.means, type = "l", col = "black")
}

compute_intensity = function(fitN, simbetas, simmeans, simsds, PredData)
{
  N = sum(fitN$pres)
  scales = c(N, as.numeric(dimnames(simbetas)[[1]]))
  n.sims = dim(simbetas)[2]
  mus.scale = rescale_mus.scale = array(NA, dim = c(length(scales) - 1, n.sims, dim(PredData)[1]))
  dimnames(mus.scale)[[1]] = dimnames(rescale_mus.scale)[[1]] = scales[-1]
  for (scale_i in 2:length(scales))
  {
    s.sds.scale = simsds[(scale_i - 1),,]
    s.means.scale = simmeans[(scale_i - 1),,]
    betas.scale = simbetas[(scale_i - 1),,]
    for (sim in 1:n.sims)
    {
      fit.sub = fitN
      fit.sub$s.sds = s.sds.scale[sim,]
      fit.sub$s.means = s.means.scale[sim,]
      fit.sub$beta = betas.scale[sim,]
      mus.scale[(scale_i - 1), sim,] = predict.ppmlasso(fit.sub, newdata = PredData)
      rescale_mus.scale[(scale_i - 1), sim,] = mus.scale[(scale_i - 1), sim,]*scales[1]/scales[scale_i]
      cat(paste(scales[scale_i], " ", sim, "\n", sep = ""))
      flush.console()
    }
  }
  return(list(mus = mus.scale, rescale_mus = rescale_mus.scale))
}


coef_plot = function(v, fitN, simbetas, col = c(228, 155, 15)/255)
{
  N = sum(fitN$pres)
  scales = c(N, as.numeric(dimnames(simbetas)[[1]]))
  n.sim = dim(simbetas)[2]
  v.beta = c(fitN$beta[v], as.vector(simbetas[,, v]))
  sd_beta = apply(simbetas[,,v], 1, sd)
  x.sims = c(scales[1], rep(scales[-1], n.sim))
  plot(x.sims, v.beta, log = "x", main = names(fitN$beta)[v], xlab = "Number of points",
       ylab = expression(hat(beta)), col = rgb(col[1], col[2], col[3], 
                                               alpha = 0.2), pch = 19)
  mean_beta = apply(simbetas[,,v], 1, mean)
  points(scales, as.numeric(c(fitN$beta[v], mean_beta)), type = "l", col = "black")
  return(list(means = mean_beta, sds = sd_beta))
}

coef_se_plot = function(v, fitN, simbetas, col = "blue")
{
  N = sum(fitN$pres)
  scales = c(N, as.numeric(dimnames(simbetas)[[1]]))
  n.sim = dim(simbetas)[2]
  v.beta = c(fitN$beta[v], as.vector(simbetas[,, v]))
  sd_beta = apply(simbetas[,,v], 1, sd)
  x.sims = scales[-1]
  plot(x.sims, sd_beta, log = "x", main = names(fitN$beta)[v], xlab = "Number of points",
       ylab = expression(SE(hat(beta))), col = col, type = "o", lwd = 2)
  mean_beta = apply(simbetas[,,v], 1, mean)
  return(list(means = mean_beta, sds = sd_beta))
}


IMSEsims = function(muN, all_mus, mu.min = 1.e-5, ...)
{
  n_scales = dim(all_mus$mus)[1]
  n.sims = dim(all_mus$mus)[2]
  IMSE_out = num.inf = matrix(NA, n_scales, n.sims)
  dimnames(IMSE_out)[[1]] = dimnames(all_mus$mus)[[1]]
  for (scale_i in 1:n_scales)
  {
    rescale_mus.scale = all_mus$rescale_mus[scale_i,,]
    for (sim in 1:n.sims)
    {
      inf.pts = which(is.infinite(log(rescale_mus.scale[sim,])))
      IMSE_out[scale_i, sim] = if (length(inf.pts) > 0) IMSE(muN[-inf.pts], rescale_mus.scale[sim, -inf.pts], mu.min = mu.min, ...) else IMSE(muN, rescale_mus.scale[sim,], mu.min = mu.min, ...)
      num.inf[scale_i, sim] = length(inf.pts)
    }
  }
  return(list(IMSE = IMSE_out, num.inf = num.inf))
}

IMSE_plot = function(muN, mus, mu.min = 1.e-5, f = NULL, col = c(228, 155, 15)/255, ...)
{
  scales = as.numeric(dimnames(mus$mus)[[1]])
  IMSEs = IMSEsims(muN, mus, mu.min = mu.min, ...)
  meanIMSEs = apply(IMSEs$IMSE, 1, mean)
  n.sims = dim(mus$mus)[2]
  y_plot = as.vector(IMSEs$IMSE)
  x_plot = rep(scales, n.sims)
  main_paste = ""
  if (is.null(f) == FALSE)
  {
    if (f == "log")
    {
      y_plot = log(y_plot)
      main_paste = "log"
      meanIMSEs = apply(log(IMSEs$IMSE), 1, mean)
    }
    if (f == "sqrt")
    {
      y_plot = sqrt(y_plot)
      main_paste = "sqrt"
      meanIMSEs = apply(sqrt(IMSEs$IMSE), 1, mean)
    }
  }
  dotscale(x_plot, y_plot, main = paste(main_paste, "Integrated mean square error"), 
           xlab = "Number of points", ylab = paste(main_paste, "Integrated mean square error"), col = col)
  print(meanIMSEs)
  return(list(means = meanIMSEs, IMSEs = IMSEs))
}


CorrSims = function(muN, all_mus, mu.min = 1.e-5, f_intensity = identity, method = "pearson")
{
  mu_truncated = muN
  mu_truncated[mu_truncated < mu.min] = mu.min
  mu_truncated = f_intensity(mu_truncated)
  n_scales = dim(all_mus$mus)[1]
  n.sims = dim(all_mus$mus)[2]
  Corr_out = matrix(NA, n_scales, n.sims)
  dimnames(Corr_out)[[1]] = dimnames(all_mus$mus)[[1]]
  for (scale_i in 1:n_scales)
  {
    rescale_mus.scale = all_mus$rescale_mus[scale_i,,]
    for (sim in 1:n.sims)
    {
      inf.pts = which(is.infinite(log(rescale_mus.scale[sim,])))
      muscale.sim = rescale_mus.scale[sim,]
      muscale.sim[muscale.sim < mu.min] = mu.min
      muscale.sim = f_intensity(muscale.sim)
      Corr_out[scale_i, sim] = if (length(inf.pts) > 0) cor(mu_truncated[-inf.pts], muscale.sim[-inf.pts], method = method)
      else cor(mu_truncated, muscale.sim, method = method)
    }
  }
  Corr_out
}

Corr_plot = function(muN, mus, mu.min = 1.e-5, method = "pearson", f = NULL, col = c(228, 155, 15)/255, threshold = NULL, ...)
{
  scales = as.numeric(dimnames(mus$mus)[[1]])
  Corrs = CorrSims(muN, mus, mu.min = mu.min, method = method, ...)
  meanCorrs = apply(Corrs, 1, mean, na.rm = TRUE)
  n.sims = dim(mus$mus)[2]
  y_plot = as.vector(Corrs)
  x_plot = rep(scales, n.sims)
  main_paste = method
  if (is.null(f) == FALSE)
  {
    if (f == "sqrt")
    {
      y_plot = sqrt(y_plot)
      main_paste = paste("sqrt", method)
      meanCorrs = apply(sqrt(Corrs), 1, mean, na.rm = TRUE)
    }
  }
  dotscale(x_plot, y_plot, main = paste(main_paste, "correlation"), 
           xlab = "Number of points", ylab = paste(main_paste, "correlation"), col = col)
  if (is.null(threshold) == FALSE)
  {
    abline(h = threshold, lty = 2, lwd = 2, col = "grey70")
  }
  print(meanCorrs)
  return(list(means = meanCorrs, Corrs = Corrs))
}


signcoefs = function(simbetas)
{
  n_subsets = dim(simbetas)[1]
  n_cov = dim(simbetas)[3]
  signs = array(NA, dim = c(n_subsets, n_cov, 3))
  dimnames(signs)[[1]] = dimnames(simbetas)[[1]]
  dimnames(signs)[[2]] = dimnames(simbetas)[[3]]
  dimnames(signs)[[3]] = c("-", "0", "+")
  for (subset in 1:n_subsets)
  {
    for (v in 1:n_cov)
    {
      signs[subset, v, 1] = sum(simbetas[subset,,v] < 0)
      signs[subset, v, 2] = sum(simbetas[subset,,v] == 0)
      signs[subset, v, 3] = sum(simbetas[subset,,v] > 0)
    }
  }
  signs
}

avg_mu_plot = function(mus, subsetsize, XY, f = identity, mu.min = 1.e-5, intensity = "rescaled", ...)
{
  subset = match(as.character(subsetsize), dimnames(mus$mus)[[1]])
  muplot = if (intensity == "rescaled") apply(mus$rescale_mus[subset,,], 2, mean) else apply(mus$mus[subset,,], 2, mean)
  muplot[muplot < mu.min] = mu.min
  z = f(muplot)
  levelplot(z ~ XY[,1] + XY[,2], asp = "iso", ...)
}

sd_plot = function(mus, subsetsize, XY, f = identity, intensity = "rescaled", ...)
{
  subset = match(as.character(subsetsize), dimnames(mus$mus)[[1]])
  sdplot = if (intensity == "rescaled") apply(mus$rescale_mus[subset,,], 2, sd) else apply(mus$mus[subset,,], 2, sd)
  z = f(sdplot)
  levelplot(z ~ XY[,1] + XY[,2], asp = "iso", ...)
}

ZeroEnvEffect = function(mus)
{
  n_subsets = dim(mus$mus)[1]
  n_sims = dim(mus$mus)[2]
  zero_out = matrix(NA, n_subsets, 2)
  dimnames(zero_out)[[1]] = dimnames(mus$mus)[[1]]
  dimnames(zero_out)[[2]] = c("Env Effect", "No Env Effect")
  for (i in 1:n_subsets)
  {
    i_flat = apply(mus$mus[i,,], 1, var)
    zero_out[i, 2] = sum(i_flat == 0)
    zero_out[i, 1] = n_sims - zero_out[i, 2]
  }
  
  zero_out
}

signplot = function(simbetas, v, ...)
{
  signs = signcoefs(simbetas = simbetas)
  subsetsizes = as.numeric(rownames(signs[,2,]))
  varnames = names(signs[1,,1])
  barplot(t(signs[rank(subsetsizes),v,]), col = c("blue", "white", "red"), main = varnames[v])
}

quantilematch = function(muN, mus, subsetsize, XY, quantiles = c(0.2, 0.4, 0.6, 0.8), ...)
{
  all_qs = c(0, quantiles, 1)
  
  subset = match(as.character(subsetsize), dimnames(mus$mus)[[1]])
  ref_quantiles = .bincode(muN, breaks = quantile(muN, all_qs), 
                           include.lowest = TRUE)
  sim_quantiles = matrix(NA, dim(mus$mus)[2], dim(mus$mus)[3])
  for (i in 1:dim(sim_quantiles)[1])
  {
    if (var(mus$mus[subset, i, ]) > 0)
    {
      sim_quantiles[i, ] = .bincode(mus$mus[subset, i, ], breaks = quantile(mus$mus[subset, i, ], all_qs),
                                    include.lowest = TRUE)
    }
    if (var(mus$mus[subset, i, ]) == 0)
    {
      sim_quantiles[i, ] = -1
    }
  }
  propmatch = rep(NA, length(ref_quantiles))
  for (i in 1:length(propmatch))
  {
    propmatch[i] = sum(sim_quantiles[, i] == ref_quantiles[i])/dim(sim_quantiles)[1]
  }
  quantileplot = levelplot(1 - propmatch ~ XY[,1] + XY[,2], asp = "iso", ...)
  print(quantileplot)
  return(list(muN_q = ref_quantiles, musub_q = sim_quantiles, prop_no_match = 1 - propmatch))
}

makeraster = function(v_frame, proj = NULL, toproj = NULL, saveTIFF = FALSE, TIFFname = "myTIFF.tif", ...)
{
  x_min = min(v_frame$X)
  x_max = max(v_frame$X)
  x_step = min(diff(sort(unique(v_frame$X))))
  x_seq = seq(x_min, x_max, x_step)
  
  y_min = min(v_frame$Y)
  y_max = max(v_frame$Y)
  y_step = min(diff(sort(unique(v_frame$Y))))
  y_seq = seq(y_max, y_min, -y_step)
  
  my_matrix = matrix(NA, nrow = length(y_seq), ncol = length(x_seq))
  rownames(my_matrix) = y_seq
  colnames(my_matrix) = x_seq
  
  for (i in 1:dim(v_frame)[1])
  {
    #my_matrix[match(v_frame$Y, y_seq), match(v_frame$X, x_seq)] = v_frame$v
    col_i = match(v_frame$X[i], x_seq)
    row_i = match(v_frame$Y[i], y_seq)
    my_matrix[row_i, col_i] = v_frame$v[i]
  }
  
  if (is.null(proj))
  {
    myraster = raster(my_matrix, xmn = min(as.numeric(colnames(my_matrix))),
                      xmx = max(as.numeric(colnames(my_matrix))),
                      ymn = min(as.numeric(rownames(my_matrix))),
                      ymx = max(as.numeric(rownames(my_matrix))), ...)
  }
  
  if (is.null(proj) == FALSE)
  {
    myraster = raster(my_matrix, xmn = min(as.numeric(colnames(my_matrix))),
                      xmx = max(as.numeric(colnames(my_matrix))),
                      ymn = min(as.numeric(rownames(my_matrix))),
                      ymx = max(as.numeric(rownames(my_matrix))), 
                      crs = proj, ...)
    if (is.null(toproj) == FALSE)
    {
      if (identical(proj, toproj) == FALSE)
      {
        myraster = projectRaster(myraster, crs = toproj)
      }
    }
  }
  
  if (saveTIFF == TRUE)
  {
    writeRaster(myraster, filename = TIFFname, format = "GTiff", overwrite = TRUE)
  }
  myraster
}
