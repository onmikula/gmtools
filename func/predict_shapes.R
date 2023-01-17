### PREDICTED SHAPES
predict_shapes <- function(score, vector, mshape, d=2, mag=1, part=NULL, eigrot=NULL) {
	
	mat2arr <- function(mat, k) {
		n <- nrow(mat)
		p <- ncol(mat) / k
		arr <- aperm(array(t(mat), dim=c(k, p, n)), c(2,1,3))
		lmlabels <- paste0("L", formatC(seq(p), format="d", flag=0, width=nchar(p)))
		dimnames(arr) <- list(lmlabels, c("X","Y","Y")[seq(k)], rownames(mat))
		return(arr)
	}

	eigenrotation <- function(M, cons=NULL, lm=NULL) {
		if (is.null(cons)) {
		cons <- M
		}
		if (is.null(lm)) {
			rm <- eigen(cov(cons))$vectors
		} else {
			rm <- eigen(cov(cons[lm,]))$vectors
		}
		center <- rep(1, nrow(M)) %*% t(colMeans(cons))
		Mr <-  (M - center) %*% rm + center
		return(Mr)
	}

	score <- as.matrix(score)
	vector <- as.matrix(vector)
	shapes <- mat2arr(score %*% t(mag * vector) + rep(1, nrow(score)) %*% t(mshape), k=d)

	if (!is.null(part)) {
		if (sum(part) != nrow(shapes)) {
			warning("sum(part) != no. of landmarks, partitioning not performed")
		} else {
			if (is.null(names(part))) {
				part <- setNames(part, seq_along(part))
			}
			part <- rep(names(part), part)
			shapes <- setNames(lapply(unique(part), function(p) shapes[part == p,,,drop=FALSE]), unique(part))
		}
	}

	if (!is.null(eigrot) & d == 2) {
		if (!is.list(eigrot)) {
			eigrot <- list(eigrot)
		}
		if (!is.null(part) & length(eigrot) != length(part)) {
			warning("length(eigrot) != length(part), eigenrotation not performed")
		} else {
			if (is.list(shapes)) {
				for (i in seq_along(shapes)) {
					cons <- apply(shapes[[i]], c(1,2), mean)
					for (j in 1:dim(shapes[[1]])[3]) {
						shapes[[i]][,,j] <- eigenrotation(M=shapes[[i]][,,j], cons=cons, lm=eigrot[[i]])
					}
				}
			} else {
				cons <- apply(shapes, c(1,2), mean)
				for (j in 1:dim(shapes)[3]) {
					shapes[,,j] <- eigenrotation(M=shapes[,,j], cons=cons, lm=eigrot)	
				}
			}
		}
	}

	return(shapes)

}


predict_sizes <- function(score, vector, msize, mag=1, which) {

	score <- as.matrix(score)
	vector <- as.matrix(vector)
	sizes <- score %*% t(mag * vector)[,which]
	size <- exp(sizes + msize)

	return(sizes)
}



### PREDICTED SHAPES CORRESPONDING TO SHIFTS ALONG PCs

predict_pca_shapes <- function(x, d=2, part=NULL) {
	
	mat2arr <- function(mat, k) {
		n <- nrow(mat)
		p <- ncol(mat) / k
		arr <- aperm(array(t(mat), dim=c(k, p, n)), c(2,1,3))
		lmlabels <- paste0("L", formatC(seq(p), format="d", flag=0, width=nchar(p)))
		dimnames(arr) <- list(lmlabels, c("X","Y","Y")[seq(k)], rownames(mat))
		return(arr)
	}

	if (inherits(x, "gpagen")) {
		n <- dim(x$coords)[3]
		p <- x$p
		d <- x$k
		pca <- prcomp(t(matrix(aperm(x$coords, c(2,1,3)), p * d, n)))	
		r <- sum(pca$sdev > (pca$sdev[1] * sqrt(.Machine$double.eps)))
	} else if (inherits(x, "prcomp")) {
		pca <- x
		n <- nrow(pca$x)
		p <- as.integer(ncol(pca$x) / d)
		r <- sum(pca$sdev > (pca$sdev[1] * pca$tol))
	} else if (inherits(x, "bgprcomp")) {
		pca <- x
		n <- nrow(pca$gmeans)
		p <- as.integer(ncol(pca$gmeans) / d)
		r <- sum(pca$sdev > (pca$sdev[1] * sqrt(.Machine$double.eps)))
	} else {
		stop("'x' must be of class 'gpagen', 'prcomp' or 'bgprcomp'")
	}

	center <- rep(1, r) %*% t(pca$center)
	shift <- diag(pca$sdev[seq(r)]) %*% t(pca$rotation[,seq(r)])
	neg <- center - shift
	pos <- center + shift
	shapes <- lapply(seq(r), function(i) mat2arr(rbind(neg[i,], pos[i,]), k=d))
	type <- c("PC", "bgPC")[inherits(pca, "bgprcomp") + 1]
	names(shapes) <- paste0(type, formatC(seq(r), format="d", flag=0, width=nchar(r)))

	if (!is.null(part)) {
		if (sum(part) != nrow(shapes[[1]])) {
			warning("sum(part) != no. of landmarks, partitioning not performed")
		} else {
			if (is.null(names(part))) {
				part <- setNames(part, seq_along(part))
			}
			part <- rep(names(part), part)
			for (i in seq(r)) {
				shapes[[i]] <- setNames(lapply(unique(part), function(p) shapes[[i]][part == p,,,drop=FALSE]), names(part))
			}
		}
	}

	return(shapes)

}


plot_pca_shapes <- function(shapes, pc=1, part=NULL, mag=1, type=c("wireframe","tps")) {
	if (!is.null(part)) {
		shapes <- lapply(shapes, "[[", part)
	}
	if (mag != 1) {
		center <- apply(shapes[[pc]], c(1,2), mean)
		shapes[[pc]][,,1] <- center - mag * (center - shapes[[pc]][,,1])
		shapes[[pc]][,,2] <- center + mag * (shapes[[pc]][,,2] - center)
	}
	if (type[1] == "wireframe") {
		wireframe2d(shapes[[pc]][,,1], shapes[[pc]][,,2], col_points=c(4,1))
		legend("topleft", legend=c("-","+"), pch=21, col=1, pt.bg=c(4,1), pt.cex=1.15, cex=1.15, title=names(shapes)[pc], adj=c(0,0.4), bty="n")
	} else if (type[1] == "tps") {
		geomorph::plotRefToTarget(shapes[[pc]][,,1], shapes[[pc]][,,2])
		title(main=names(shapes)[pc], cex.main=1.15)
	}
}


