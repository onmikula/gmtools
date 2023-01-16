# utils for conversion between file formats usal in morphometric analyses
# n ~ no. of configurations
# p ~ no. of landmarks
# k ~ no. of dimensions
# config ~ new data matrix


# array to list along arbitrary dimension (original name "split.along.dim")
# contributed by "flodel" at
# http://stackoverflow.com/questions/20198751/three-dimensional-array-to-list
arr2list <- function(a, n=3) {
	setNames(
		lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]), array, dim = dim(a)[-n], dimnames(a)[-n]),
	dimnames(a)[[n]])
}


# list to array
list2arr <- function(l) {
	stopifnot(all(diff(t(sapply(l, dim))) == 0))
	n <- length(l)
	p <- nrow(l[[1]])
	k <- ncol(l[[1]])
	arr <- array(unlist(l), dim=c(p, k, n))
	lmlabels <- rownames(l[[1]])
	if (is.null(lmlabels)) {
		lmlabels <- paste0(L, formatC(seq(p), format="d", flag=0, width=nchar(p)))
	}
	dimnames(arr) <- list(lmlabels, c("X","Y","Y")[seq(k)], names(l))
	return(arr)
}


# matrix to array
mat2arr <- function(mat, k=2) {
	n <- nrow(mat)
	p <- ncol(mat) / k
	arr <- aperm(array(t(mat), dim=c(k, p, n)), c(2,1,3))
	lmlabels <- paste0("L", formatC(seq(p), format="d", flag=0, width=nchar(p)))
	dimnames(arr) <- list(lmlabels, c("X","Y","Y")[seq(k)], rownames(mat))
	return(arr)
}


# array to matrix
arr2mat <- function(arr) {
	p <- nrow(arr)
	k <- ncol(arr)
	n <- dim(arr)[3]
	mat <- t(matrix(aperm(arr, c(2,1,3)), p*k, n))
	rownames(mat) <- dimnames(arr)[[3]]
	lmlab <- dimnames(arr)[[1]]
	if (is.null(lmlab)) {
		lmlab <- paste0(L, formatC(seq(p), format="d", flag=0, width=nchar(p)))
	}
	axislab <- c("X","Y","Y")[seq(k)]
	colnames(mat) <- paste0(rep(lmlab, each=k), axislab)
	return(mat)
}
