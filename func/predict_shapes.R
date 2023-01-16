### PREDICTED SHAPES
predict_shapes <- function(score, vector, mshape, msize=NULL, d=2, mag=1, part=NULL, eigrot=NULL) {

	score <- as.matrix(score)
	vector <- as.matrix(vector)
	N <- nrow(score)
	pred <- score %*% t(mag * vector)

	prmean <- function(pred, mshape, msize, d) { 
		N <- nrow(pred)
		if (!is.null(dim(mshape))) {
			mshape <- as.vector(t(mshape))
		}
		if (!is.null(msize)) {
			size <- exp(pred[,1] + msize)
			shape <- pred[,-1] + rep(1,N) %*% t(mshape)
			form <- diag(size) %*% shape
		} else {
			size <- NULL
			shape <- pred + rep(1,N) %*% t(mshape)
			form <- NULL
		}
		shape <- mat2arr(shape, k=d)
		if (!is.null(form)) {
			form <- mat2arr(form, k=d)
		}
		return(list(size=size, shape=shape, form=form))
	}

	mat2arr <- function(mat, k=2) {
		n <- nrow(mat)
		p <- ncol(mat) / k
		arr <- aperm(array(t(mat), dim=c(k, p, n)), c(2,1,3))
		lmlabels <- paste0("L", formatC(seq(p), format="d", flag=0, width=nchar(p)))
		dimnames(arr) <- list(lmlabels, c("X","Y","Y")[seq(k)], rownames(mat))
		return(arr)
	}

	Eigenrotation <- function(M, cons=NULL, lm=NULL) {
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

	center <- function(X) X - rep(1, nrow(X)) %*% t(colMeans(X))
	csize <- function(X) sqrt(sum(center(X)^2))

	if (!is.null(part)) {
		i <- 1; j <- 0
		size <- list()
		shape <- list()
		form <- list()
		while (i <= length(part)) {
			idx <- (j+1):sum(part[1:i])
			mshapei <- mshape[[i]]
			msizei <- msize[[i]]
			predi <- prmean(pred[,idx], mshapei, msizei, k=d)
			size <- c(size, predi[1])
			shape <- c(shape, predi[2])
			form <- c(form, predi[3])
			j <- sum(part[1:i])
			i <- i + 1
		}
	} else {
		predi <- prmean(pred, mshape, msize, k=d)
		size <- predi[1]
		shape <- predi[2]
		form <- predi[3]
	}
	if (!is.null(part)) {
		units <- sapply(shape, function(x) max(apply(x, 3, dist)))
		for (i in seq_along(part)) {
			for (j in seq(nrow(score))) {
				shape[[i]][,,j] <- center(shape[[i]][,,j]) / units[i]
			}
		}
	}

	if (!is.null(eigrot) & d == 2) {
		if (!is.list(eigrot)) {
			eigrot <- list(eigrot)
		}
		for (i in 1:length(shape)) {
			cons <- apply(shape[[i]], c(1,2), mean)
			for (j in 1:dim(shape[[1]])[3]) {
				shape[[i]][,,j] <- Eigenrotation(M=shape[[i]][,,j], cons=cons, lm=eigrot[[i]])
			}
		}
		if (!is.null(msize)) {
			for (i in 1:length(form)) {
				cons <- apply(form[[i]], c(1,2), mean)
				for (j in 1:dim(shape[[1]])[3]) {
					form[[i]][,,j] <- Eigenrotation(M=form[[i]][,,j], cons=cons, lm=eigrot[[i]])
				}
			}
		}
	}

	return(list(size=size, shape=shape, form=form))
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

	neg <- diag(pca$sdev[seq(r)]) %*% t(-pca$rotation[,seq(r)]) + rep(1, r) %*% t(pca$center)
	pos <- diag(pca$sdev[seq(r)]) %*% t(pca$rotation[,seq(r)]) + rep(1, r) %*% t(pca$center)
	shapes <- lapply(seq(r), function(i) mat2arr(rbind(neg[i,], pos[i,]), k=d))
	type <- c("PC", "bgPC")[inherits(pca, "bgprcomp") + 1]
	names(shapes) <- paste0(type, formatC(seq(r), format="d", flag=0, width=nchar(r)))

	if (!is.null(part)) {
		if (sum(part) != nrow(shapes)) {
			warning("partitioning not performed because the sum of 'part' is not equal to the no. of rows in the matrix of predicted shapes")
		} else {
			if (is.null(names(part))) {
				part <- setNames(part, seq_along(part))
			}
			part <- rep(names(part), diff(cumsum(c(0, part))))
			for (i in seq(r)) {
				shapes[[i]] <- setNames(lapply(unique(part), function(p) shapes[[i]][part == p,,,drop=FALSE]), names(part))
			}
		}
	}

	return(shapes)

}


plot_pca_shapes <- function(shapes, i=1, part=NULL, mag=1) {
	if (!is.null(part)) {
		shapes <- lapply(shapes, "[[", part)
	}
	geomorph::plotRefToTarget(shapes[[i]][,,1], shapes[[i]][,,2], mag=mag)
	title(main=names(shapes)[i], cex.main=2)
}


