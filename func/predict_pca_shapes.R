predict_pca_shapes <- function(x, k=2) {
	if (inherits(x, "gpagen")) {
		n <- dim(x$coords)[3]
		p <- x$p
		k <- x$k
		d <- min(c(p * k, n))
		pca <- prcomp(t(matrix(aperm(x$coords, c(2,1,3)), p * k, n)))	
	} else if (inherits(x, "prcomp")) {
		pca <- x
		n <- nrow(pca$x)
		p <- as.integer(ncol(pca$x) / k)
		d <- sum(pca$sdev > (pca$sdev[1] * pca$tol))
	} else if (inherits(x, "bgprcomp")) {
		pca <- x
		n <- nrow(pca$gmeans)
		p <- as.integer(ncol(pca$gmeans) / k)
		d <- qr(pca$x)$r
	} else {
		stop("'x' must be of class 'gpagen', 'prcomp' or 'bgprcomp'")
	}
	shapes <- diag(pca$sdev) %*% t(pca$rotation) + rep(1, d) %*% t(pca$center)
	return(aperm(array(t(shapes), dim=c(k, p, n)), c(2,1,3)))
}
