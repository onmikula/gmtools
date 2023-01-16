plot_lm_outliers <- function(gpa) {
	mahalanobis <- function(x, s) t(x) %*% s %*% x
	resid <- lapply(seq(dim(gpa$coords)[3]), function(i) gpa$coords[,,i] - gpa$consensus)
	resid <- array(unlist(resid), dim=c(nrow(resid[[1]]), ncol(resid[[1]]), length(resid)), dimnames=list(NULL, NULL, names(resid)))
	icovm <- lapply(apply(resid, 1, function(x) cov(t(x)), simplify=FALSE), solve)
	dst <- sapply(seq(nrow(resid)), function(i) apply(resid[i,,], 2, mahalanobis, s=icovm[[i]]))
	maxdist <- apply(dst, 2, which.max)
	outliers <- t(sapply(seq_along(maxdist), function(i) gpa$coords[i,, maxdist[i],drop=FALSE]))
	plot(gpa)
	for (i in seq(gpa$p)) {
		lines(rbind(gpa$consensus[i,],outliers[i,]), col=2)
	}
	points(outliers, col=2, pch=16)
	points(gpa$consensus, pch=21, col=1, bg="white", cex=1.9, lwd=0.25)
	text(gpa$consensus, labels=seq(nrow(gpa$consensus)), cex=0.5)
}


print_lm_outliers <- function(gpa, lm=NULL, order=NULL) {
	mahalanobis <- function(x, s) t(x) %*% s %*% x
	resid <- lapply(seq(dim(gpa$coords)[3]), function(i) gpa$coords[,,i] - gpa$consensus)
	resid <- array(unlist(resid), dim=c(nrow(resid[[1]]), ncol(resid[[1]]), length(resid)), dimnames=list(NULL, NULL, names(resid)))
	icovm <- lapply(apply(resid, 1, function(x) cov(t(x)), simplify=FALSE), solve)
	dst <- sapply(seq(nrow(resid)), function(i) apply(resid[i,,], 2, mahalanobis, s=icovm[[i]]))
	maxdist <- apply(dst, 2, which.max)
	if (is.null(lm)) {
		lm <- seq(nrow(gpa$coords))
	}
	df <- data.frame(lm=as.integer(lm), id=dimnames(gpa$coords)[[3]][maxdist[lm]])
	if (!is.null(order)) {
		df$order <- unname(order[maxdist[lm]])
	}
	print(df)
}

