library(data.table)

DecimalCount.old = function(x, max.dec = max(10, max(nchar(x))))
{
	x.mult  = x*10^(0:max.dec)
	num.dec = min(which(x.mult %% 1 == 0)) - 1
	num.dec
}

DecimalCount = function(x, max.dec = max(10, max(nchar(x))), tol = 1.e-1)
{
	digits  = 0:max.dec
	x.diff = (x - round(x, digits))/(10^(-1*(digits + 3)))
	num.dec = digits[min(which(abs(x.diff) < tol))]
	num.dec
}

SpatRes = function(env.grid, coord = c("X", "Y"))
{
	x.col   = which(names(env.grid) == coord[1])
	y.col   = which(names(env.grid) == coord[2])
	x.uq    = sort(unique(env.grid[, x.col]))
	y.uq    = sort(unique(env.grid[, y.col]))
	n.dec   = max(unlist(lapply(x.uq, DecimalCount)))
	x.diff  = diff(x.uq)
	y.diff  = diff(y.uq)
	x.dec   = unlist(lapply(x.diff, DecimalCount))
	y.dec   = unlist(lapply(y.diff, DecimalCount))
	x.step  = min(floor(x.diff*10^max(x.dec) + 0.1))/(10^max(x.dec))
	y.step  = min(floor(y.diff*10^max(y.dec) + 0.1))/(10^max(y.dec))
	return(list(x.step = x.step, y.step = y.step))
}

ZapCoord = function(x, numdec = DecimalCount(x))
{
	x.out  = floor(x*10^numdec + 0.1)/(10^numdec)
	x.out
}

sample.quad.old = function(env.grid, sp.scale, coord = c("X", "Y"), file = "Quad")
{
	convert = FALSE
	#if (any(lapply(env.grid, class) == "factor"))
	#{
	#	is.cat    = which(lapply(env.grid, class) == "factor")
	#	out.grid  = CatConvert(env.grid)
	#	cat.names = out.grid$cat.names
	#	env.grid  = out.grid$X
	#	convert   = TRUE
	#}
	f.name = list()
	for(i in 1:length(sp.scale))
	{
		i.scale = sp.scale[i]
		x.col   = which(names(env.grid) == coord[1])
		y.col   = which(names(env.grid) == coord[2])
		x.uq    = sort(unique(env.grid[, x.col]))
		y.uq    = sort(unique(env.grid[, y.col]))
		
		scale.mult = max(c(unlist(lapply(x.uq, DecimalCount)), unlist(lapply(y.uq, DecimalCount))))
		
		i.scale = i.scale*10^scale.mult
		env.grid[, x.col] = env.grid[, x.col]*10^scale.mult
		env.grid[, y.col] = env.grid[, y.col]*10^scale.mult

		#x.step  = diff((x.uq*10^scale.mult)[1:2])
		x.step	= min(diff(x.uq*10^scale.mult))
		#y.step  = diff((x.uq*10^scale.mult)[1:2])
		y.step	= min(diff(y.uq*10^scale.mult))

		x.o = min(env.grid[,x.col]) - floor(min(env.grid[,x.col])/x.step)*x.step
		y.o = min(env.grid[,y.col]) - floor(min(env.grid[,y.col])/y.step)*y.step
		if (x.o/x.step > 0.5)
		{
			x.o = x.o - x.step
		}	
		if (y.o/y.step > 0.5)
		{
			y.o = y.o - y.step
		}

		is.on.scale   = abs((env.grid[,x.col]/i.scale) - round(env.grid[,x.col]/i.scale) - x.o/i.scale) + abs((env.grid[,y.col]/i.scale) - round(env.grid[,y.col]/i.scale) - y.o/i.scale) < 1.e-8
	  	dat.quad      = env.grid[is.on.scale,]
		dat.quad[, x.col] = dat.quad[, x.col]/(10^scale.mult)
		dat.quad[, y.col] = dat.quad[, y.col]/(10^scale.mult)
		if (is.na(file) == FALSE)
		{
			f.name[[i]] = paste(file, sp.scale[i], ".RData", sep = "")
			save(dat.quad, file = f.name[[i]])
			print(paste("Output saved in the file", f.name[[i]]))
		}
	}
	if (convert == TRUE)
	{
		dat.quad = list(X = dat.quad, cat.names = cat.names)
	}
	if (length(sp.scale) == 1)
		dat.quad
	else
		f.name
}

sample.quad = function(env.grid, sp.scale, coord = c("X", "Y"), file = "Quad")
{
	convert = FALSE
	x.col   = which(names(env.grid) == coord[1])
	y.col   = which(names(env.grid) == coord[2])
	res     = SpatRes(env.grid, coord = coord)
	x.step  = res$x.step
	y.step  = res$y.step
	f.name = list()
	for(i in 1:length(sp.scale))
	{
		i.scale = sp.scale[i]
		x.o = min(env.grid[,x.col]) - floor(min(env.grid[,x.col])/x.step)*x.step
		y.o = min(env.grid[,y.col]) - floor(min(env.grid[,y.col])/y.step)*y.step
		if (x.o/x.step > 0.5)
		{
			x.o = x.o - x.step
		}	
		if (y.o/y.step > 0.5)
		{
			y.o = y.o - y.step
		}
		
		
    # x.on.scale = abs(((env.grid[,x.col] - x.o)/i.scale) %% 1) < 1.e-4
    # y.on.scale = abs(((env.grid[,y.col] - y.o)/i.scale) %% 1) < 1.e-4
    is.on.scale = abs(((env.grid[,x.col] - x.o)/i.scale) %% 1 + ((env.grid[,y.col] - y.o)/i.scale) %% 1) < 1.e-4
    
		# is.on.scale   = abs((env.grid[,x.col]/i.scale) - round(env.grid[,x.col]/i.scale) - x.o/i.scale) + abs((env.grid[,y.col]/i.scale) - round(env.grid[,y.col]/i.scale) - y.o/i.scale) < 1.e-8
		# is.on.scale   = abs(((env.grid[,x.col] - x.o)/i.scale) - round((env.grid[,x.col] - x.o)/i.scale)) + abs(((env.grid[,y.col] - y.o)/i.scale) - round((env.grid[,x.col] - x.o)/i.scale)) < 1.e-8
		# table(is.on.scale)
		# 
		# 
			dat.quad      = env.grid[is.on.scale,]
		### Get rid of machine error in coordinates
		dec.x = max(unlist(lapply(dat.quad[,x.col], DecimalCount)))
		dec.y = max(unlist(lapply(dat.quad[,y.col], DecimalCount)))
		dat.quad[,x.col] = unlist(lapply(dat.quad[, x.col], ZapCoord, dec.x))
		dat.quad[,y.col] = unlist(lapply(dat.quad[, y.col], ZapCoord, dec.y))
		#n.x           = floor((dat.quad[,x.col] - min(dat.quad[,x.col]))/x.step + 0.1)
		#n.y           = floor((dat.quad[,y.col] - min(dat.quad[,y.col]))/y.step + 0.1)
		#dat.quad[,x.col] = min(dat.quad[,x.col]) + n.x*x.step
		#dat.quad[,y.col] = min(dat.quad[,y.col]) + n.y*y.step

		if (is.na(file) == FALSE)
		{
			f.name[[i]] = paste(file, sp.scale[i], ".RData", sep = "")
			save(dat.quad, file = f.name[[i]])
			print(paste("Output saved in the file", f.name[[i]]))
		}
	}
	if (convert == TRUE)
	{
		dat.quad = list(X = dat.quad, cat.names = cat.names)
	}
	if (length(sp.scale) == 1)
		dat.quad
	else
		f.name
}

ppmdat.old = function(sp.xy, sp.scale, back.xy, coord = c("X","Y"), sp.dat = env.var(sp.xy = sp.xy, env.scale = sp.scale, env.grid = back.xy, coord = coord, file.name = "SpEnvData"), sp.file = NA, quad.file = NA, file.name = NA)
{
	convert = FALSE
	if (class(sp.dat) == "list")
	{
		convert   = TRUE
		cat.names = sp.dat$cat.names
		sp.dat    = sp.dat$X
	}
	if (is.character(sp.xy) == TRUE)
	{
		sp.file = paste(sp.xy, "Env.RData")
		load(sp.file)
	}
	if (is.na(sp.file) == FALSE)
	{
		load(sp.file)
	}
	if (is.na(quad.file) == TRUE)
	{
		save(sp.scale, file = "Test.RData")
		dat.quad = sample.quad(env.grid = back.xy, sp.scale = sp.scale, coord = coord, file = NA)
	}
	if (is.na(quad.file) != TRUE)
	{
		load(paste(quad.file, sp.scale, ".RData", sep = ""))
	}

	dat.quad$Pres = 0
	sp.dat$Pres   = 1
	
	quad.x.col = which(names(dat.quad) == coord[1])
	quad.y.col = which(names(dat.quad) == coord[2])
	sp.x.col   = which(names(sp.dat) == coord[1])
	sp.y.col   = which(names(sp.dat) == coord[2])
	sp.var     = setdiff(1:dim(sp.dat)[2], c(sp.x.col, sp.y.col))

	sp.data         = data.frame(sp.dat[,sp.x.col], sp.dat[,sp.y.col], sp.dat[,sp.var])
	save(sp.data, file = "TestSP.RData")
	names(sp.data)  = c("X", "Y", names(sp.dat)[sp.var])
	
	if (any(lapply(dat.quad, class) == "factor"))
	{
		dat.quad = CatConvert(dat.quad)$X
	}
	add.var = setdiff(names(sp.dat), names(dat.quad))
	if (length(add.var) > 0)
	{
		quad.names = names(dat.quad)
		for (add in 1:length(add.var))
		{
			dat.quad = data.frame(dat.quad, 0)
		}
		names(dat.quad) = c(quad.names, add.var)

		dat.quad = dat.quad[,match(names(dat.quad), names(sp.dat))]
	}

	quad.var        = setdiff(1:dim(dat.quad)[2], c(quad.x.col, quad.y.col))
	quad.dat        = data.frame(dat.quad[,quad.x.col], dat.quad[,quad.y.col], dat.quad[,quad.var])
	names(quad.dat) = c("X", "Y", names(dat.quad)[quad.var])
	quad.dat        = quad.dat[,match(names(sp.dat), names(quad.dat))]
	dat.ppm         = rbind(sp.data, quad.dat)
	
	dat.ppm$wt = ppmlasso:::weights(sp.data, quad.dat, coord)	

	dimnames(dat.ppm)[[1]] = 1:dim(dat.ppm)[1]
	if (is.na(file.name) == FALSE)
	{
		save.name = paste(file.name, ".RData", sep = "")
		save(dat.ppm, file = save.name)
		print(paste("Output saved in the file", save.name))
	}
	if (convert == TRUE)
	{
		dat.ppm = list(X = dat.ppm, cat.names = cat.names)
	}
	dat.ppm
}

ppmdat = function(sp.xy, sp.scale, back.xy, coord = c("X","Y"), sp.dat = newenv.var(sp.xy = sp.xy, env.scale = sp.scale, env.grid = back.xy, coord = coord, file.name = "SpEnvData"), sp.file = NA, quad.file = NA, file.name = NA)
{
	convert = FALSE
	if (class(sp.dat) == "list")
	{
		convert   = TRUE
		cat.names = sp.dat$cat.names
		sp.dat    = sp.dat$X
	}
	if (is.character(sp.xy) == TRUE)
	{
		sp.file = paste(sp.xy, "Env.RData")
		load(sp.file)
	}
	if (is.na(sp.file) == FALSE)
	{
		load(sp.file)
	}
	if (is.na(quad.file) == TRUE)
	{
		save(sp.scale, file = "Test.RData")
		dat.quad = sample.quad(env.grid = back.xy, sp.scale = sp.scale, coord = coord, file = NA)
	}
	if (is.na(quad.file) != TRUE)
	{
		load(paste(quad.file, sp.scale, ".RData", sep = ""))
	}

	dat.quad$Pres = 0
	sp.dat$Pres   = 1
	
	quad.x.col = which(names(dat.quad) == coord[1])
	quad.y.col = which(names(dat.quad) == coord[2])
	sp.x.col   = which(names(sp.dat) == coord[1])
	sp.y.col   = which(names(sp.dat) == coord[2])
	sp.var     = setdiff(1:dim(sp.dat)[2], c(sp.x.col, sp.y.col))

	sp.data         = data.frame(sp.dat[,sp.x.col], sp.dat[,sp.y.col], sp.dat[,sp.var])
	save(sp.data, file = "TestSP.RData")
	names(sp.data)  = c("X", "Y", names(sp.dat)[sp.var])
	
	if (any(lapply(dat.quad, class) == "factor"))
	{
		dat.quad = CatConvert(dat.quad)$X
	}
	add.var = setdiff(names(sp.dat), names(dat.quad))
	if (length(add.var) > 0)
	{
		quad.names = names(dat.quad)
		for (add in 1:length(add.var))
		{
			dat.quad = data.frame(dat.quad, 0)
		}
		names(dat.quad) = c(quad.names, add.var)

		dat.quad = dat.quad[,match(names(dat.quad), names(sp.dat))]
	}

	quad.var        = setdiff(1:dim(dat.quad)[2], c(quad.x.col, quad.y.col))
	quad.dat        = data.frame(dat.quad[,quad.x.col], dat.quad[,quad.y.col], dat.quad[,quad.var])
	names(quad.dat) = c("X", "Y", names(dat.quad)[quad.var])
	quad.dat        = quad.dat[,match(names(sp.dat), names(quad.dat))]
	dat.ppm         = rbind(sp.data, quad.dat)
	
	dat.ppm$wt = ppmlasso:::weights(sp.data, quad.dat, coord)	

	dimnames(dat.ppm)[[1]] = 1:dim(dat.ppm)[1]
	if (is.na(file.name) == FALSE)
	{
		save.name = paste(file.name, ".RData", sep = "")
		save(dat.ppm, file = save.name)
		print(paste("Output saved in the file", save.name))
	}
	if (convert == TRUE)
	{
		dat.ppm = list(X = dat.ppm, cat.names = cat.names)
	}
	dat.ppm
}

interp = function(sp.xy, sp.scale, f, back.xy, coord = c("X","Y"))
{
	options(scipen = 999)
	x.dat   = sp.xy[,which(names(sp.xy) == coord[1])]
	y.dat   = sp.xy[,which(names(sp.xy) == coord[2])]
	x.back  = back.xy[,which(names(back.xy) == coord[1])]
	y.back  = back.xy[,which(names(back.xy) == coord[2])]

	ux    = sort(unique(x.back))
	uy    = sort(unique(y.back))

	x.col   = which(names(back.xy) == coord[1])
	y.col   = which(names(back.xy) == coord[2])

	scale.mult = max(c(unlist(lapply(ux, DecimalCount)), unlist(lapply(uy, DecimalCount))))

	sp.scale = sp.scale*10^scale.mult
	back.xy[, x.col] = back.xy[, x.col]*10^scale.mult
	back.xy[, y.col] = back.xy[, y.col]*10^scale.mult


	#grain = ux[2] - ux[1]
	#step  = sp.scale/grain

	x.step  = diff((ux*10^scale.mult)[1:2])
	y.step  = diff((uy*10^scale.mult)[1:2])

	x.o = min(back.xy[,x.col]) - floor(min(back.xy[,x.col])/x.step)*x.step
	y.o = min(back.xy[,y.col]) - floor(min(back.xy[,y.col])/y.step)*y.step

	#if (x.o/x.step > 0.5)
	#{
	#	x.o = x.o - x.step
	#}	
	#if (y.o/y.step > 0.5)
	#{
	#	y.o = y.o - y.step
	#}

	x.back = x.back*10^scale.mult
	y.back = y.back*10^scale.mult
	x.dat  = x.dat*10^scale.mult
	y.dat  = y.dat*10^scale.mult
	ux     = ux*10^scale.mult
	uy     = uy*10^scale.mult

	col.ref   = match(x.back, ux)
	row.ref   = match(y.back, uy)

	env.grid         = matrix(NA, length(uy), length(ux), dimnames = list(uy*10^scale.mult, ux*10^scale.mult))
	all.vec          = rep(NA, max(row.ref)*max(col.ref))
	vec.ref          = (col.ref - 1)*max(row.ref) + row.ref
	all.vec[vec.ref] = f
	f.grid           = matrix(all.vec, max(row.ref), max(col.ref), dimnames = list(uy*10^scale.mult, ux*10^scale.mult))

	x.1   = round(floor((x.dat - x.o)/sp.scale)*sp.scale + x.o, 1)
	y.1   = round(floor((y.dat - y.o)/sp.scale)*sp.scale + y.o, 1)
	x.2   = pmin(x.1 + sp.scale, max(ux))
	y.2   = pmin(y.1 + sp.scale, max(uy))

	#x.1   = round(floor(x.dat/sp.scale)*sp.scale + x.o, 1)
	#y.1   = round(floor(y.dat/sp.scale)*sp.scale + y.o, 1)
	#x.2   = pmin(x.1 + sp.scale, max(ux))
	#y.2   = pmin(y.1 + sp.scale, max(uy))

	w11   = (x.2 - x.dat)*(y.2 - y.dat)/((x.2 - x.1)*(y.2 - y.1))
	w12   = (x.2 - x.dat)*(y.dat - y.1)/((x.2 - x.1)*(y.2 - y.1))
	w21   = (x.dat - x.1)*(y.2 - y.dat)/((x.2 - x.1)*(y.2 - y.1))
	w22   = (x.dat - x.1)*(y.dat - y.1)/((x.2 - x.1)*(y.2 - y.1))

	x.1.id = match(x.1, ux)
	y.1.id = match(y.1, uy)
	x.2.id = match(x.2, ux)
	y.2.id = match(y.2, uy)

	if (length(x.1.id) != 1)
	{
		f11 = diag(f.grid[y.1.id, x.1.id])
		f12 = diag(f.grid[y.2.id, x.1.id])
		f21 = diag(f.grid[y.1.id, x.2.id])
		f22 = diag(f.grid[y.2.id, x.2.id])
	}

	if (length(x.1.id) == 1)
	{
		f11 = f.grid[y.1.id, x.1.id]
		f12 = f.grid[y.2.id, x.1.id]
		f21 = f.grid[y.1.id, x.2.id]
		f22 = f.grid[y.2.id, x.2.id]
	}

	c11 = 1 - is.na(f11)
	c12 = 1 - is.na(f12)
	c21 = 1 - is.na(f21)
	c22 = 1 - is.na(f22)

	env.wt.mat = cbind(f11*w11*c11, f12*w12*c12, f21*w21*c21, f22*w22*c22) 
	
	f.interp = apply(env.wt.mat, 1, sum, na.rm = TRUE)/(w11*c11 + w12*c12 + w21*c21 + w22*c22)
	f.interp
}

findres = function(scales, lambda = 0, coord = c("X", "Y"), sp.xy, env.grid, formula, ...)
{
	form.1 = formula
	likelihoods = rep(NA, length(scales))
	#sp.data = env.var(sp.xy = sp.xy, env.scale = min(scales), env.grid = env.grid, coord = coord, file.name = "SpEnvData")
	sp.data = newenv.var(sp.xy = sp.xy, env.scale = min(scales), env.grid = env.grid, coord = coord, file.name = "SpEnvData")
	for (sc in 1:length(scales))
	{
		formula = form.1
		if (lambda == 0)
		{
			data            = ppmdat(sp.xy = sp.xy, sp.scale = scales[sc], back.xy = env.grid, sp.dat = sp.data, sp.file = NA, quad.file = NA)
			if (class(data) == "list")
			{
				use.form = as.character(formula)[2]
				cat.names = setdiff(unique(data$cat.names), NA)
				for (i in 1:length(cat.names))
				{
					use.form = gsub(cat.names[i], paste(names(data$X)[which(data$cat.names == cat.names[i])], collapse = " + "), use.form)
				}
				formula = as.formula(paste("~", use.form))
				data = data$X
			}
			glm.form        = as.formula(paste("Pres/wt ~ ", as.character(formula)[2], sep = ""))
			glm.fit         = glm(glm.form, data = data, weights = data$wt, family = poisson())
			eps   = 1.e-9
			while (glm.fit$deviance > glm.fit$null.dev)
			{
				glm.fit = glm(glm.form, data = data, weights = data$wt, family = poisson(), control = list(epsilon = eps))
				eps     = eps/10
			}
			#glm.fit         = glm(glm.form, data = data, weights = data$wt, family = poisson(), control = list(epsilon = 1.e-12))
			likelihoods[sc] = sum(data$wt*(data$Pres/data$wt*log(glm.fit$fitted) - glm.fit$fitted)) - sum(log(1:sum(data$Pres > 0)))
		}
		if (lambda != 0)
		{
			sc.fit = ppmlasso(ppm.form, sp.scale = scales[sc], lamb = lambda, data = ppmdat(sp.xy = sp.xy, sp.scale = scales[sc], back.xy = env.grid, sp.dat = sp.data, sp.file = NA, quad.file = NA), ...)
			likelihoods[sc] = sc.fit$pen.likelihood[1]
		}
	}
	plot(scales, likelihoods, log = "x", type = "o", pch = 16, xlab = "Spatial Resolution", ylab = "Likelihood")
}

polynames = function(X)
{
	coef.names = dimnames(X)[[2]]
	which.poly = grep("poly\\(", coef.names)
	if (length(which.poly) > 0)
	{
	  split1     = strsplit(coef.names[which.poly], "\\)")
	  polyframe  =  plyr::ldply(split1, rbind)
	  varframe   = plyr::ldply(strsplit(unlist(strsplit(as.character(polyframe[,1]), "poly\\("))[2*(1:(length(which.poly)))], ","), rbind)
	  varframe   = varframe[,-dim(varframe)[2]]
	  varframe   = data.frame(lapply(varframe, as.character), stringsAsFactors = FALSE)
	  if (dim(varframe)[1] == 1)
	  {
	    varframe = t(varframe)
	  }
	  expframe   = plyr::ldply(strsplit(as.character(polyframe[,2]), "\\."), rbind)
	  expframe   = data.frame(lapply(expframe, as.character), stringsAsFactors = FALSE)
	  expframe[is.na(expframe)] = 0
	  vframe = varframe[,1:dim(expframe)[2]]
	  nameframe  = matrix(paste(as.matrix(vframe), "^", as.matrix(expframe), sep = ""), nrow(vframe), ncol(vframe))
	  nameframe[expframe == "0"] = ""
	  nameframe[expframe == "1"] = as.character(vframe[expframe == 1])
	  nameframe = gsub(" ", "", nameframe)
	  nameframe  = as.data.frame(nameframe)
	  if (is.null(dim(nameframe)) == FALSE)
	  {
	    names.out = apply(nameframe, 1, function(row) paste(row[nzchar(row)], collapse = "*"))
	  }
	  if (is.null(dim(nameframe)) == TRUE)
	  {
	    names.out = nameframe
	  }
	  id.out    = which.poly
	}
	else
	{
	  names.out = coef.names
	  id.out = NULL
	}
	return(list(names = names.out, ids = id.out))
}

oldpolynames = function(X) # delete
{
	coef.names = dimnames(X)[[2]]
	which.poly = grep("poly\\(", coef.names)
	if (length(which.poly) > 0)
	{
		poly.exp = unlist(strsplit(coef.names[which.poly], "\\)"))[2*(1:(length(which.poly)))]
		poly.var = unlist(strsplit(coef.names[which.poly], "\\)"))[1]
		var.list = unlist(strsplit(unlist(strsplit(poly.var, "poly\\("))[2], ","))
		var.list = gsub(" $", "", gsub("^ ","", var.list, perl = T), perl = T)[-length(var.list)]
		exp.mat  = t(matrix(as.numeric(unlist(strsplit(poly.exp, "\\."))), length(var.list), length(poly.exp)))
		max.deg  = max(exp.mat)

		name.mat = exp.mat
		for (v in 1:length(var.list))
		{
			name.mat[,v] = gsub("0", "", name.mat[,v], perl = T)
			name.mat[,v] = gsub("1", var.list[v], name.mat[,v], perl = T)
			for (deg in 2:max.deg)
			{
				name.mat[,v] = gsub(as.character(deg), paste(var.list[v], "^", deg, sep = ""), name.mat[,v], perl = T)
			}
		}

		names.out = apply(name.mat, 1, function(row) paste(row[nzchar(row)], collapse = "*"))
		id.out    = which.poly
		return(list(names = names.out, ids = id.out))
	}
}


ppmlasso = function(formula, sp.xy, env.grid, sp.scale, coord = c("X", "Y"), data = ppmdat(sp.xy = sp.xy, sp.scale = sp.scale, back.xy = env.grid, coord = c("X","Y"), sp.file = NA, quad.file = NA, file.name = "TestPPM"), lamb = NA, n.fits = 200, ob.wt = NA, criterion = "bic", alpha = 1, family = "poisson", tol = 1.e-9, gamma = 0, init.coef = NA, mu.min = 1.e-16, mu.max = 1/mu.min, r = NA, interactions = NA, availability = NA, max.it = 25, min.lamb = -10, standardise = TRUE, n.blocks = NA, block.size = sp.scale*100, seed = 1)
{
	error.flag = FALSE
	formula.out = formula
	if (class(data) == "list")
	{
		use.form = as.character(formula)[2]
		cat.names = setdiff(unique(data$cat.names), NA)
		for (i in 1:length(cat.names))
		{
			use.form = gsub(cat.names[i], paste(names(data$X)[which(data$cat.names == cat.names[i])], collapse = " + "), use.form)
		}
		formula = as.formula(paste("~", use.form))
		data = data$X
	}
	wt.calc = FALSE
	if (is.na(ob.wt) == TRUE)
	{
		wt.calc = TRUE
		ob.wt   = data$wt
	}

	call = match.call()
	mf   = model.frame(formula, data = data)
	mt   = attr(mf, "terms")
   	y    = data$Pres/data$wt
	X    = if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
   	 else matrix(, NROW(y), 0L)
	
	if (any(X[,1] != 1))
	{
		X = as.matrix(cbind(1, X))
	}
	cut.var = which(apply(X, 2, max) == 0)
	if (length(cut.var) > 0)
	{
		X       = X[, -cut.var]
	}
	Xnames  = polynames(X)
	if (is.null(Xnames) == FALSE)
	{
		dimnames(X)[[2]][Xnames$ids] = Xnames$names
	}
	dimnames(X)[[2]][1] = "Intercept"

	area.int = FALSE
	raw.int  = NA
	if (family == "area.inter")
	{
		family   = "poisson"
		area.int = TRUE
		if (is.na(interactions) == TRUE)
		{
			interactions = point.interactions(data, r, availability)
		}
		raw.int = interactions
	}

	s.means = NULL
	s.sds   = NULL

	if (standardise == TRUE)
	{
		stand.X = standardise.X(X[,-1])
		X       = as.matrix(cbind(1, stand.X$X))
		dimnames(X)[[2]][1] = "Intercept"
		s.means = stand.X$dat.means
		s.sds   = stand.X$dat.sds
		if (area.int == TRUE)
		{
			interactions = standardise.X(interactions)$X
		}
	}

	if (family == "poisson")
	{
        	family = get(family, mode = "function", envir = parent.frame())
	}
   	if (is.function(family))
	{
        	family = family()
	}
	if (family$family == "binomial")
	{
		mu.max  = 1 - mu.min
	}
	score.0    = score.int(y, X, ob.wt = ob.wt, area.int = area.int, int = interactions, family)

	if (is.na(init.coef))
	{
		gamma = 0
	}
	if (gamma == 0)
	{
		adapt.weights = rep(1, dim(X)[2])
	}
	if (gamma != 0)
	{
		adapt.weights = 1/abs(init.coef)^gamma
	}

	cv = rep(0, dim(data)[1])

	if (criterion == "blockCV")
	{
		cv      = blocks(n.blocks, block.size, data, seed = seed)
		pred.mu = matrix(NA, dim(data)[1], n.fits)
	}
	data.all = data
	y.all    = y
	X.all    = X
	wt.all   = ob.wt
	if (area.int == TRUE)
	{
		interactions.all = interactions
	}

	for (cv.i in 1:length(unique(cv)))
	{
	dat.test = data[cv == cv.i,]
	data     = data[cv != cv.i,]
	data$wt  = weights(data[data$Pres == 1,], data[data$Pres == 0,], coord)
	y        = data$Pres/data$wt
	X        = X[cv != cv.i,]
	ob.wt    = ob.wt[cv != cv.i]
	if (wt.calc == TRUE)
	{
		ob.wt = data$wt
	}
	if (area.int == TRUE)
	{
		interactions = interactions[cv != cv.i]
	}

	if (is.na(lamb) == TRUE)
	{
		new.score  = abs(score.0/adapt.weights)
		sub.score  = new.score[-1]
		max.lambda = max(sub.score[is.infinite(sub.score) == FALSE])
		if (is.na(max.lambda) == FALSE)
		{
			lambs = sort(exp(seq(min.lamb, log(max.lambda + 1.e-5), length.out = n.fits)), decreasing = TRUE)
		}
	}
	if (is.na(lamb) == FALSE)
	{
		lambs = lamb
	}

	if (criterion == "blockCV")
	{
		lambda.max = max(abs(score.int(data$Pres/data$wt, X, ob.wt = data$wt, area.int = area.int, int = interactions, family = poisson())[-1]))	
		lambs      = exp(seq(0, -12, length.out = n.fits))*lambda.max
	}

	n.pres     = sum(y > 1.e-8)
	n.var      = dim(X)[2] - 1
	
	if (criterion == "msi")
	{
		lambs = max.lambda/sqrt(n.pres)
	}

	if (area.int == TRUE)
	{
		X.0 = as.matrix(cbind(X, interactions))
	}
	if (area.int != TRUE)
	{
		X.0 = X
	}

	mod.0     = glm(y ~ X.0[,-1], family = family, weights = ob.wt)
	#coefs     = matrix(NA, (dim(X)[2] + area.int), (length(lambs) + 1), dimnames = list(dimnames(X)[[2]]))
	coefs     = matrix(NA, (dim(X)[2]), (length(lambs) + 1), dimnames = list(dimnames(X)[[2]]))
	if (area.int == 1)
	{
		coefs     = matrix(NA, (dim(X)[2] + area.int), (length(lambs) + 1), dimnames = list(c(dimnames(X)[[2]], "Interaction")))
	}
	num.param = rep(NA, (length(lambs) + 1))
	gcvs      = rep(NA, (length(lambs) + 1))
	aics      = rep(NA, (length(lambs) + 1))
	hqcs      = rep(NA, (length(lambs) + 1))
	bics      = rep(NA, (length(lambs) + 1))
	devs      = rep(NA, (length(lambs) + 1))
	ll        = rep(NA, (length(lambs) + 1))
	pll       = rep(NA, (length(lambs) + 1))
	nlgcvs    = rep(NA, (length(lambs) + 1))
	offset    = c(log(mean(y)), rep(0, n.var), rep(1, area.int))

	if (any(is.na(adapt.weights) == FALSE))
	{
		if (sum(is.infinite(adapt.weights[-1])) != (length(adapt.weights) - 1 - area.int))
		{
			it.max = 100
			for (reg.path in 1:length(lambs))
			{
				mod  = try(single.lasso(y, X, lamb = lambs[reg.path], ob.wt = ob.wt, alpha = alpha, b.init = offset, family = family, tol = 1.e-9, gamma = gamma, init.coef = init.coef, area.int = area.int, interactions = interactions, max.it = it.max, standardise = FALSE), TRUE)
				if(class(mod) == "try-error")
				{
					break
				}
				if (any(is.na(mod$b)))
				{
					break
				}
				coefs[,reg.path]    = mod$b
				gcvs[reg.path]      = mod$GCV
				aics[reg.path]      = mod$AIC
				hqcs[reg.path]      = mod$HQC
				bics[reg.path]      = mod$BIC
				devs[reg.path]      = mod$dev
				ll[reg.path]        = mod$like
				pll[reg.path]       = mod$pen
				num.param[reg.path] = sum(abs(mod$b) > 1.e-7) - 1 - area.int
				offset = mod$b
				it.max = max.it
				cat(paste("Fitting Models:", reg.path, "of", length(lambs), "\n"))
				flush.console()
			}

			num.done = length(lambs) + 1 - sum(is.na(aics))
			s.denom  = sum(abs(mod.0$coef[-c(1, (n.var + 2))]))
	
			if (n.var > n.pres)
			{
				s.denom = sum(abs(coefs[-c(1, (n.var + 2)), num.done]))
			}

			if (any(is.na(mod.0$coef)) == TRUE)
			{
				s.denom = sum(abs(coefs[-c(1, (n.var + 2)), num.done]))
			}

			ss     = rep(NA, (length(lambs) + 1))
			nlgcvs = rep(NA, (length(lambs) + 1))

			for (i.nlgcv in 1:num.done)
			{
				s.num       = sum(abs(coefs[-c(1, (n.var + 2)), i.nlgcv]))
				s           = s.num/s.denom
				ss[i.nlgcv] = s
				nlgcv       = devs[i.nlgcv]/(n.pres*(1 - (area.int + n.var*s)/n.pres)^2)
				if (area.int + n.var*s > n.pres)
				{
					nlgcv = NA
				}
				nlgcvs[i.nlgcv] = nlgcv
			}
		}
		
		if (sum(is.infinite(adapt.weights[-1])) == (length(adapt.weights) - 1 - area.int))
		{
			coefs     = matrix(rep(init.coef, (length(lambs) + 1)), length(adapt.weights), (length(lambs) + 1))
			gcvs      = rep(1, (length(lambs) + 1))
			aics      = rep(1, (length(lambs) + 1))
			bics      = rep(1, (length(lambs) + 1))
			hqcs      = rep(1, (length(lambs) + 1))
			nlgcvs    = rep(1, (length(lambs) + 1))
			devs      = rep(1, (length(lambs) + 1))
			ll        = rep(1, (length(lambs) + 1))
			num.param = rep(0, (length(lambs) + 1))
			num.done  = length(lambs) + 1 - sum(is.na(aics))
		}

		criterion.matrix        = data.frame(aics, bics, hqcs, gcvs, nlgcvs)
		names(criterion.matrix) = c("AIC", "BIC", "HQC", "GCV", "NLGCV")

		lambs[(length(lambs) + 1)] = 0
		coefs[,length(lambs)] = mod.0$coef
		if (gamma != 0)
		{
			coefs[,length(lambs)] = init.coef
		}
		num.param[(length(lambs) + 1)] = length(adapt.weights) - sum(is.infinite(adapt.weights[-c(1, (n.var + 2))])) - 1 - area.int
		
		meth.id  = paste(criterion, "s", sep = "")

		if (criterion == "msi" | criterion == "blockCV")
		{
			choice.id = 1
		}
		if (criterion != "msi" & criterion != "blockCV")
		{
			choice.id = max(which.min(get(meth.id)))
		}

		lambda.hat = lambs[choice.id]
		beta.hat   = coefs[,choice.id]
		eta.hat    = X.0 %*% beta.hat
		mu.hat     = family$linkinv(eta.hat)
		like.hat   = ll[choice.id]

		assign(paste("coefs.", cv.i, sep = ""), coefs)
		

		if (criterion == "blockCV")
		{
			for (i in 1:length(lambs) - 1)
			{ #Fix this for predicting to test locations
				fam.fit = poisson()
				if (area.int == TRUE)
				{
					fam.fit = "area.inter"
				}
				cv.fit = list(beta = coefs[,i], s.means = s.means, s.sds = s.sds, family = fam.fit, pt.interactions = interactions, formula = formula, mu = rep(0, dim(data)[1]))
				if (area.int == TRUE)
				{
					pred.mu[cv == cv.i, i] = predict.ppmlasso(cv.fit, newdata = dat.test, interactions = interactions.all[cv == cv.i])
				}
				if (area.int != TRUE)
				{
					pred.mu[cv == cv.i, i] = predict.ppmlasso(cv.fit, newdata = dat.test)
				}
			}
		}
	}

	data  = data.all
	y     = y.all
	X     = X.all
	ob.wt = wt.all
	if (area.int == TRUE)
	{
		interactions = interactions.all
	}

	} # cv.i

	if (criterion == "blockCV")
	{
		l.vec       = apply(pred.mu, 2, unp.likelihood, ob.wt = data$wt, y = data$Pres/data$wt)
		coef.mat    = matrix(NA, n.blocks, dim(coefs)[1])
		for (i in 1:n.blocks)
		{
			coefs.i       = get(paste("coefs.", i, sep = ""))
			coef.mat[i, ] = coefs.i[, which.max(l.vec)]
		}
		coef.av = apply(coef.mat, 2, mean, na.rm = TRUE)
		lambda.mult = exp(seq(0, -12, length.out = n.fits))[which.max(l.vec)]
		if (area.int != TRUE)
		{
			lambda.max = max(abs(score.int(data$Pres/data$wt, X, ob.wt = data$wt, family = poisson())[-1]))
		}	

		if (area.int == TRUE)
		{
			lambda.max = max(abs(score.int(data$Pres/data$wt, X, ob.wt = data$wt, family = poisson(), area.int = TRUE, int = scale(interactions))[-1]))
		}
		final.fit = try(single.lasso(y.all, X.all, lamb = lambda.mult*lambda.max, ob.wt = wt.all, alpha = alpha, b.init = coef.av, family = family, tol = 1.e-9, gamma = gamma, init.coef = init.coef, area.int = area.int, interactions = interactions, max.it = it.max, standardise = FALSE), TRUE)
		coefs     = final.fit$b
		lambs     = lambda.mult*lambda.max
		gcvs      = NA
		aics      = NA
		hqcs      = NA
		bics      = NA
		nlgcvs    = NA
		devs      = final.fit$dev
		ll        = final.fit$like
		pll       = final.fit$pen
		num.param = sum(abs(final.fit$b) > 1.e-7) - 1 - area.int

		if (area.int == TRUE)
		{
			X.0 = as.matrix(cbind(X, interactions))
		}
		if (area.int != TRUE)
		{
			X.0 = X
		}
		
		lambda.hat = lambs
		beta.hat   = coefs
		eta.hat    = X.0 %*% beta.hat
		mu.hat     = family$linkinv(eta.hat)
		like.hat   = ll
		criterion.matrix = data.frame(aics, bics, hqcs, gcvs, nlgcvs)
	}

	family.out = family$family
	if (area.int == TRUE)
	{	
		family.out = "area.inter"
	}
	output = list(betas = coefs, lambdas = lambs, likelihoods = ll, pen.likelihoods = pll, lambda = lambda.hat, beta = beta.hat, mu = mu.hat, likelihood = like.hat, criterion = criterion, family = family.out, gamma = gamma, alpha = alpha, init.coef = init.coef, criterion.matrix = criterion.matrix, data = X.0, pt.interactions = raw.int, wt = ob.wt, pres = data$Pres, x = data$X, y = data$Y, r = r, call = call, formula = formula.out, s.means = s.means, s.sds = s.sds, cv.group = cv, n.blocks = n.blocks)
	class(output) = c("ppmlasso", "list")
	output
}

newinterp = function(sp.xy, sp.scale, f, back.xy, coord = c("X","Y"))
{
	options(scipen = 999)
	x.dat   = sp.xy[,which(names(sp.xy) == coord[1])]
	y.dat   = sp.xy[,which(names(sp.xy) == coord[2])]
	x.back  = back.xy[,which(names(back.xy) == coord[1])]
	y.back  = back.xy[,which(names(back.xy) == coord[2])]
	res     = SpatRes(back.xy, coord = coord)

	grid    = data.table(x.back, y.back, f, key = c("x.back", "y.back"))

	ux    = sort(unique(x.back))
	uy    = sort(unique(y.back))

	x.col   = which(names(back.xy) == coord[1])
	y.col   = which(names(back.xy) == coord[2])

	x.step  = res$x.step
	y.step  = res$y.step

	x.o = min(back.xy[,x.col]) - floor(min(back.xy[,x.col])/x.step)*x.step
	y.o = min(back.xy[,y.col]) - floor(min(back.xy[,y.col])/y.step)*y.step

	x.1   = floor((x.dat - x.o)/sp.scale)*sp.scale + x.o
	y.1   = floor((y.dat - y.o)/sp.scale)*sp.scale + y.o
	x.2   = pmin(x.1 + sp.scale, max(ux))
	y.2   = pmin(y.1 + sp.scale, max(uy))

	w11   = (x.2 - x.dat)*(y.2 - y.dat)/((x.2 - x.1)*(y.2 - y.1))
	w12   = (x.2 - x.dat)*(y.dat - y.1)/((x.2 - x.1)*(y.2 - y.1))
	w21   = (x.dat - x.1)*(y.2 - y.dat)/((x.2 - x.1)*(y.2 - y.1))
	w22   = (x.dat - x.1)*(y.dat - y.1)/((x.2 - x.1)*(y.2 - y.1))

	f11   = grid[list(x.1, y.1)]$f
	f12   = grid[list(x.1, y.2)]$f
	f21   = grid[list(x.2, y.1)]$f
	f22   = grid[list(x.2, y.2)]$f

	c11 = 1 - is.na(f11)
	c12 = 1 - is.na(f12)
	c21 = 1 - is.na(f21)
	c22 = 1 - is.na(f22)

	env.wt.mat = cbind(f11*w11*c11, f12*w12*c12, f21*w21*c21, f22*w22*c22) 
	
	f.interp = apply(env.wt.mat, 1, sum, na.rm = TRUE)/(w11*c11 + w12*c12 + w21*c21 + w22*c22)
	f.interp
}


newenv.var = function(sp.xy, env.grid, env.scale, coord = c("X","Y"), file.name = NA)
{
	convert = FALSE
	if (any(lapply(env.grid, class) == "factor"))
	{
		convert  = TRUE
		out.grid = CatConvert(env.grid)
		env.grid = out.grid$X
	}
	x.dat   = sp.xy[,which(names(sp.xy) == coord[1])]
	y.dat   = sp.xy[,which(names(sp.xy) == coord[2])]
	x.back  = env.grid[,which(names(env.grid) == coord[1])]
	y.back  = env.grid[,which(names(env.grid) == coord[2])]
	x.col   = which(names(env.grid) == coord[1])
	y.col   = which(names(env.grid) == coord[2])
	var.col = setdiff(1:dim(env.grid)[2], c(x.col, y.col))
	s.res   = SpatRes(env.grid)
	
	sp.dat        = as.data.frame(matrix(NA, length(x.dat), length(var.col)))
	names(sp.dat) = names(env.grid[var.col])

	for (var in 1:length(var.col))
	{
		loop.scale = min(c(s.res$x.step, s.res$y.step))
		loc        = which(is.na(sp.dat[,var]))
		while(sum(is.na(sp.dat[,var])) > 0)
		{
			loc = which(is.na(sp.dat[,var]))
			sp.dat[loc, var] = newinterp(sp.xy[loc,], loop.scale, env.grid[,var.col[var]], env.grid, coord = c("X","Y"))
			loop.scale = loop.scale*2
		}
		cat(paste("Calculating species environmental data for variable:", names(sp.dat)[var], "\n"))
		flush.console()
	}
	
	sp.dat = data.frame(x.dat, y.dat, sp.dat)
	names(sp.dat)[1:2] = c("X", "Y")
	if (is.na(file.name) == FALSE)
	{
		save.name = paste(file.name, ".RData", sep = "")
		save(sp.dat, file = save.name)
		print(paste("Output saved in the file", save.name))
	}
	if (convert == TRUE)
	{
		sp.dat = list(X = sp.dat, cat.names = out.grid$cat.names)
	}
	sp.dat
}

makeMask = function(dat.ppm)
{
	if (class(dat.ppm) == "list")
	{
		dat.ppm = dat.ppm$X
	}
	dat.q = dat.ppm[dat.ppm$Pres == 0,]
	ux = sort(unique(dat.q$X))
	uy = sort(unique(dat.q$Y))
	nx = length(ux)
	ny = length(uy)

	col.ref = match(dat.q$X, ux)
	row.ref = match(dat.q$Y, uy)

	all.vec          = rep(0, max(row.ref)*max(col.ref))
	vec.ref          = (col.ref - 1)*max(row.ref) + row.ref
	all.vec[vec.ref] = 1
	mask.out         = matrix(all.vec, max(row.ref), max(col.ref), dimnames = list(uy, ux))
	mask.out
}


point.interactions = function(dat.ppm, r, availability = NA)
{
	if (class(dat.ppm) == "list")
	{
		dat.ppm = dat.ppm$X
	}
	if (any(is.na(availability)))
	{
		#grain = min(0.5, r/2)
		#min.x = min(dat.ppm$X) - r - grain
		#max.x = max(dat.ppm$X) + r + grain
		#min.y = min(dat.ppm$Y) - r - grain
		#max.y = max(dat.ppm$Y) + r + grain

		#availability = matrix(1, 1 + (max.y - min.y)/grain, 1 + (max.x - min.x)/grain)
		#rownames(availability) = seq(min.y, max.y, grain)
		#colnames(availability) = seq(min.x, max.x, grain)
		availability = makeMask(dat.ppm)
	}	

	cat(paste("Calculating point interactions", "\n"))
	flush.console()
	occupied = matrix(0, dim(availability)[1], dim(availability)[2])
	rownames(occupied) = rownames(availability)
	colnames(occupied) = colnames(availability)

	x.mat = availability
	y.mat = availability

	for (i in 1:dim(x.mat)[1])
	{
		x.mat[i,] = as.numeric(colnames(availability))
	}
	for (i in 1:dim(y.mat)[2])
	{
		y.mat[,i] = as.numeric(rownames(availability))
	}

	grain = as.numeric(colnames(availability)[2]) - as.numeric(colnames(availability)[1])

	pres.x = dat.ppm$X[dat.ppm$Pres > 0]
	pres.y = dat.ppm$Y[dat.ppm$Pres > 0]

	quad.x = dat.ppm$X[dat.ppm$Pres == 0]
	quad.y = dat.ppm$Y[dat.ppm$Pres == 0]

	quad.int = rep(0, length(quad.x))

	for (i in 1:length(pres.x))
	{
		sub.col = which(as.numeric(colnames(occupied)) >= pres.x[i] - (r + grain) & as.numeric(colnames(occupied)) <= pres.x[i] + (r + grain))
		sub.row = which(as.numeric(rownames(occupied)) >= pres.y[i] - (r + grain) & as.numeric(rownames(occupied)) <= pres.y[i] + (r + grain))

		sub.occ = occupied[sub.row, sub.col]
		sub.x   = x.mat[sub.row, sub.col]
		sub.y   = y.mat[sub.row, sub.col]

		sub.occ[(sub.x - pres.x[i])^2 + (sub.y - pres.y[i])^2 < r^2] = sub.occ[(sub.x - pres.x[i])^2 + (sub.y - pres.y[i])^2 < r^2] + 1
		occupied[sub.row, sub.col] = sub.occ
		quad.cells           = which((quad.x - pres.x[i])^2 + (quad.y - pres.y[i])^2 <= (2*r)^2)
		quad.int[quad.cells] = quad.int[quad.cells] + 1
	}

	int.q = rep(0, length(quad.int))

	for (quad.i in which(quad.int > 0))
	{
		sub.col = which(as.numeric(colnames(occupied)) >= quad.x[quad.i] - (r + grain) & as.numeric(colnames(occupied)) <= quad.x[quad.i] + (r + grain))
		sub.row = which(as.numeric(rownames(occupied)) >= quad.y[quad.i] - (r + grain) & as.numeric(rownames(occupied)) <= quad.y[quad.i] + (r + grain))

		sub.occ  = occupied[sub.row, sub.col]
		sub.availability = availability[sub.row, sub.col]
		sub.x    = x.mat[sub.row, sub.col]
		sub.y    = y.mat[sub.row, sub.col]

		sub.cell = (sub.x - quad.x[quad.i])^2 + (sub.y - quad.y[quad.i])^2 <= r^2 & sub.availability > 0

		#int.q[quad.i] = sum(sub.occ[sub.cell] > 0)/sum(sub.cell)
		int.q[quad.i] = sum(sub.occ[sub.cell] > 0, na.rm = TRUE)/sum(sub.cell, na.rm = TRUE)
	}

	int.p = rep(0, length(pres.x))

	for (pres.i in 1:length(pres.x))
	{
		sub.col = which(as.numeric(colnames(occupied)) >= pres.x[pres.i] - (r + grain) & as.numeric(colnames(occupied)) <= pres.x[pres.i] + (r + grain))
		sub.row = which(as.numeric(rownames(occupied)) >= pres.y[pres.i] - (r + grain) & as.numeric(rownames(occupied)) <= pres.y[pres.i] + (r + grain))

		sub.occ  = occupied[sub.row, sub.col]
		sub.availability = availability[sub.row, sub.col]
		sub.x    = x.mat[sub.row, sub.col]
		sub.y    = y.mat[sub.row, sub.col]

		sub.cell = (sub.x - pres.x[pres.i])^2 + (sub.y - pres.y[pres.i])^2 <= r^2 & sub.availability > 0

		#int.p[pres.i] = sum(sub.occ[sub.cell] > 1)/sum(sub.cell)
		int.p[pres.i] = sum(sub.occ[sub.cell] > 1, na.rm = TRUE)/sum(sub.cell, na.rm = TRUE)
	}

	interactions = c(int.p, int.q)
	interactions
}

plotpath = function(fit, colors = c("gold", "green3", "blue", "brown", "pink"))
{
	best.models = apply(fit$criterion.matrix, 2, which.min) #See which fit optimises other criteria
	unique.beta = unique(best.models)
	names = rep(NA, length(unique.beta))
	for (i in 1:length(unique.beta))
	{
		names[i] = paste(names(best.models[best.models == unique.beta[i]]), collapse = "/")
	}

	min.y = min(apply(fit$betas[,-dim(fit$betas)[2]], 1, min))
	max.y = max(apply(fit$betas[,-dim(fit$betas)[2]], 1, max))

	plot(fit$lambdas, fit$betas[1,], log = "x", type = "l", ylim = c(min.y, max.y), xlab = "LASSO penalty", ylab = "Coefficients")
	for (i in 2:dim(fit$betas)[1])
	{
		points(fit$lambdas, fit$betas[i,], type = "l")
	}

	for (i in 1:length(unique.beta))
	{
		abline(v = fit$lambdas[unique.beta[i]], lwd = 3, col = colors[i])
	}	

	legend("topright", names, lwd = rep(3, length(unique.beta)), col = colors[1:length(unique.beta)])
}

plotfit = function(fit, pred.data = data.frame(X = fit$x[fit$pres == 0], Y = fit$y[fit$pres == 0], 1, scale(fit$data[fit$pres == 0, -1], center = -fit$s.means/fit$s.sds, scale = 1/fit$s.sds)), 
	coord = c("X", "Y"), asp = "iso", ylab = "", xlab = "", col.regions = heat.colors(1024)[900:1], 
	cuts = length(col.regions), cex = 1.4, main.text = paste(toupper(fit$criterion), "fit"), cex.color = 1.4,
	log.z = FALSE, sqrt.z = FALSE, cutquantile = NULL)
{
	x.col    = which(names(pred.data) == coord[1])
	y.col    = which(names(pred.data) == coord[2])
	pred.int = predict(fit, newdata = pred.data[,-c(x.col, y.col)])
	if (log.z == TRUE)
	{
	  pred.int = log(pred.int)
	}
	if (sqrt.z == TRUE)
	{
	  pred.int = sqrt(pred.int)
	}
	if (is.null(cutquantile) == FALSE)
	{
	  max.z = quantile(pred.int, cutquantile)
	 pred.int[pred.int > max.z] = max.z
	}
	levelplot(pred.int ~ pred.data[,x.col] + pred.data[,y.col], asp = asp, 
		ylab = ylab, xlab = xlab, col.regions = col.regions, cuts = cuts, 
		main = list(main.text, cex = cex), at = seq(min(pred.int), max(pred.int), length.out = (cuts + 1)),
		scales = list(y = list(draw = FALSE), x = list(draw = FALSE)), 
		cex = cex, colorkey = list(labels = list(cex = cex.color)))

}

