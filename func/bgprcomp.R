bgprcomp <- function(x, g) {
	x <- as.matrix(x)
	x <- x[!is.na(g),]
	g <- factor(na.omit(g))
	gm <- do.call(rbind, by(x, g, colMeans))
	pca <- prcomp(gm)
	centered <- x - rep(1, nrow(x)) %*% t(pca$center)
	x_ind <- centered %*% pca$rotation
	sdev_ind <- apply(x_ind, 2, sd)
	result <- c(pca, list(group=g, gmeans=gm, x_ind=x_ind, sdev_ind=sdev_ind, data=x))
	class(result) <- "bgprcomp"
	return(result)
}


plot.bgprcomp <- function(bgpca, i=1, j=2, col=1, bg=8, legend=NA, asp=NA, cex=1, cex.group=2, cex.legend=c(1.25, 0.9), font.legend=1, ord.legend, ...) {
	i <- i[1]
	j <- j[1]
	group <- bgpca$group
	K <- nlevels(group)
	L <- levels(group)
	if (all(L %in% names(col))) {
		col <- col[L]
	} else {
		col <- setNames(rep_len(col, K), L)
	}
	if (all(L %in% names(bg))) {
		bg <- bg[L]
	} else {
		bg <- setNames(rep_len(bg, K), L)
	}
	plot(bgpca$x_ind[,c(i,j)], xlab=paste("bgPC", i), ylab=paste("bgPC", j), type="n", asp=asp, ...)
	points(bgpca$x_ind[,c(i,j)], pch=21, col=bg[group], bg=NA, lwd=1.5, cex=cex)
	points(bgpca$x[,c(i,j)], pch=21, col=col, bg=bg, lwd=2, cex=cex.group)
	if (!is.na(legend)) {
		if (missing(ord.legend)) {
			ord.legend <- seq(K)
		}
		legend(legend, legend=L[ord.legend], pch=21, col=col[ord.legend], pt.bg=bg[ord.legend], pt.cex=cex.legend[1], cex=cex.legend[2], text.font=font.legend, bty="n")
	}
}


predict.bgprcomp <- function(bgpca, newdata) {
	new <- as.matrix(newdata) - rep(1, nrow(newdata)) %*% t(bgpca$center)
	return(new %*% bgpca$rotation)
}


# bgpca: bgprcomp object or a list with components data and group
# nfold: specify 'n' of n-fold cross validation, defaults to leave-one-out cross validation
# ncomp: how many compomnents to use, default is no. of groups - 1

cv.bgprcomp <- function(bgpca, nfold, ncomp) {
	data <- bgpca$data
	group <- factor(bgpca$group)
	N <- nrow(data)
	K <- nlevels(group)
	L <- levels(group)
	if (missing(ncomp)) {
		ncomp <- K - 1
	}
	if (missing(nfold)) {
		nfold <- N
	} else if (nfold > N) {
		nfold <- N
		warning(paste("'nfold' was changed to", N))
	}	
	nmin <- min(table(group))
	if (nfold < N & nfold > nmin) {
		nfold <- nmin
		warning(paste("'nfold' was changed to", nmin))
	}
	if (nfold < N) {
		segments <- vector("list", nfold)
		for (j in seq(K)) {
			k <- which(group == L[j])
			ni <- length(k)
			intervals <- findInterval(seq(ni), seq(1, ni, length=nfold + 1), rightmost.closed=TRUE)
			segments <- Map(c, segments, split(sample(k, ni), intervals))
		}
	} else {
		segments <- as.list(seq(N))
	}
	confmatrix <- matrix(0, K, K, dimnames=list(L,L))
	cvpred <- group
	cvdata <- data
	for (i in seq(nfold)) {
		test <- segments[[i]]
		bgpca <- bgprcomp(data[-test,], group[-test])
		new <- predict(bgpca, data[test,,drop=FALSE])
		dst <- as.matrix(dist(rbind(bgpca$x, new)[,seq(ncomp),drop=FALSE]))[-seq(K),seq(K),drop=FALSE]
		grp <- colnames(dst)[apply(dst, 1, which.min)]
		grp <- group[match(grp, group)]
		cvpred[test] <- grp
		cvdata[test,] <- new %*% t(bgpca$rotation) + rep(1, nrow(new)) %*% t(bgpca$center) 
	}
	confmatrix <- table(group, cvpred)
	cv <- diag(confmatrix)
	corr <- sum((1 / K) * table(group))
	tau <- (sum(cv) - corr) / (N - corr)		#Cardini & Elton (2011), Int J Primatol, 32: 377–389
	result <- list(confmatrix=confmatrix, success=tau, prediction=cvpred, group=group, cvdata=cvdata)
	return(result)
}



cvscatter <- function(cv) {
	gmeans <- as.matrix(do.call(rbind, by(cv$cvdata, cv$group, colMeans)))
	pca <- prcomp(gmeans)
	x_ind <- cv$cvdata - rep(1, nrow(cv$cvdata)) %*% t(pca$center)
	x_ind <- x_ind %*% pca$rotation
	result <- list(sdev=pca$sdev, rotation=pca$rotation, center=pca$center, scale=pca$scale, x=pca$x, group=cv$group, gmeans=gmeans, x_ind=x_ind, sdev_ind=apply(x_ind, 2, sd), data=cv$cvdata)
	class(result) <- "bgprcomp"
	return(result)
}




coherence <- function(bgpca, d) {
	K <- nrow(bgpca$gmeans)
	if (missing(d)) {
		d <- K - 1
	} else {
		d <- min(c(K - 1, d))
	}
	dst <- as.matrix(dist(rbind(bgpca$x[,seq(d),drop=FALSE], bgpca$x_ind[,seq(d),drop=FALSE])))[-seq(K), seq(K)]
	diag(table(bgpca$group, apply(dst, 1, which.min))) / table(bgpca$group)
}


# bgpca: 'bgprcomp' object
# type: 'residual', 'gaussian' or 'isotropic'
# reorder: reordering of data rows

random.bgprcomp <- function(bgpca, type="residual", reorder=NULL) {
	data <- bgpca$data
	if (type[1] == "residual") {
		data <- data - bgpca$gmeans[bgpca$group,]
	} else if (type[1] == "gaussian") {
		covmat <- cov(data - bgpca$gmeans[bgpca$group,])
		data <- mvtnorm::rmvnorm(nrow(data), mean=bgpca$center, sigma=covmat)
	} else if (type[1] == "isotropic") {
		data <- mvtnorm::rmvnorm(nrow(data), mean=bgpca$center, sigma=diag(ncol(data)))
	}
	if (is.null(reorder)) {
		reorder <- sample(nrow(data))
	}
	data <- data[reorder,]
	return(bgprcomp(x=data, g=bgpca$group))
}



isotropic <- function(n, p) {
	return(replicate(p, rnorm(n, 0, 1)))
}


rotate.bgprcomp <- function(bgpca, which) {
	bgpca$rotation[,which] <- -1 * bgpca$rotation[,which]
	bgpca$x[,which] <- -1 * bgpca$x[,which]
	bgpca$x_ind[,which] <- -1 * bgpca$x_ind[,which]
	class(bgpca) <- "bgprcomp"
	return(bgpca)
}

