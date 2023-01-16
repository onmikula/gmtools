get_scale_factor <- function(tps, scale=1, rm=NULL) {
	x <- readTPS(tps, asList=TRUE, rm=rm)
	nz <- sapply(x, length) > 0
	nz[nz] <- sapply(x[nz], nrow) > 1
	x[nz] <- lapply(x[nz], function(m) m[1:2,]) 
	d <- rep(NA, length(x))
	d[nz] <- sapply(x[nz], dist)
	f <- setNames(scale / d, names(x))
	return(f)
}


set_scale_factor <- function(X, scale=NULL, Y=NULL, lm=NULL, rmx=NULL, rmy=NULL, asList=FALSE) {

	arr2list <- function(a) {
		setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, 3]), array, dim = dim(a)[-3], dimnames(a)[-3]), dimnames(a)[[3]])
	}
	list2arr <- function(l) {
		array(unlist(l), dim=c(nrow(l[[1]]), ncol(l[[1]]), length(l)), dimnames=list(NULL, NULL, names(l)))
	}

	if (is.character(X) & length(X) == 1) {
		x <- readTPS(X, asList=TRUE, rm=rmx)
		semilm <- attr(x, "semilm")
	} else if (is.array(X)) {
		semilm <- attr(X, "semilm")
		x <- arr2list(X)
	} else if (is.list(X)) {
		semilm <- attr(X, "semilm")
		x <- X
	} else {
		stop("'X' has to be a .TPS file name or an array or a list with landmark configurations")
	}

	nz <- sapply(x, length) > 0
	if (length(scale) == 1) {
		scale <- rep_len(scale, sum(nz))
	} else if (length(scale) >= length(x)) {
		if (!is.null(names(x)) & !is.null(names(scale))) {
			m <- match(names(x), names(scale))
			if (any(is.na(m))) {
				stop("if provided, names of 'scale' must match those of 'X'")
			} else {
				scale <- scale[m]
			}
		}
		scale <- scale[nz]
	} else if (length(scale) > 0) {
		stop("if specified, 'scale' must be of length one or of the same length as the no. of configurations in 'X'")
	}
	if (!is.null(scale)) {
		if(any(is.na(scale))) {
			stop("'scale' must not contain any NAs")
		}
		x[nz] <- Map("*", x[nz], as.list(scale))
	}
	
	if (is.null(Y)) {
		if (isFALSE(asList)) {
			x <- list2arr(x)
		}
		attr(x, "semilm") <- semilm
		return(x)
	} else if (!is.null(lm)) {
		if (is.character(Y) & length(Y) == 1) {
			y <- readTPS(Y, asList=TRUE, rm=rmy)
			semilm <- attr(y, "semilm")
		} else if (is.array(Y)) {
			semilm <- attr(Y, "semilm")
			y <- arr2list(Y)
		} else if (is.list(Y)) {
			semilm <- attr(Y, "semilm")
			y <- Y
		} else {
			stop("'Y' has to be a .TPS file name or an array or a list with landmark configurations")
		}
		subs <- intersect(names(x), names(y))
		if (length(subs) < length(y)) {
			y <- y[subs]
			warning("subsetting done, not all configurations of 'Y' are present in 'X'")
		}
		xnz <- sapply(x, length) > 0
		ynz <- sapply(y, length) > 0
		if (all(xnz[ynz])) {
			xlm <- x
			xlm[ynz] <- lapply(x[ynz], function(p, ii) p[ii,], ii=lm[[1]])
			ylm <- y
			ylm[ynz] <- lapply(y[ynz], function(p, ii) p[ii,], ii=lm[[2]])
			dx <- rep(NA, length(xlm))
			dx[ynz] <- sapply(xlm[ynz], dist)
			dy <- rep(NA, length(ylm))
			dy[ynz] <- sapply(ylm[ynz], dist)
			yscale <- setNames(dx / dy, names(y))
			y[ynz] <- Map("*", y[ynz], as.list(yscale))
			if (isFALSE(asList)) {
				y <- list2arr(y)
			}
			attr(y, "semilm") <- semilm
			return(y)
		} else {
			stop("non-empty configurations of 'Y' must have non-empty counterparts in 'X'")
		}
	} else {
		stop("if 'Y' is provided, 'lm' must be provided as well")
	}

}



