unify_reflection <- function(A) {
	d <- apply(!is.na(A),1,all)
	r <- sample(dim(A)[3],1)
	s <- t(scale(A[d,,r],scale=FALSE))
	m <- t(colMeans(A[d,,r]))
	for (i in seq(dim(A)[3])) {
		SVD <- svd(s %*% (A[d,,i] - rep(1,sum(d)) %*% m))
		RM <- SVD$v %*% t(SVD$u)
		A[,,i] <- (A[,,i] - rep(1,nrow(A)) %*% m) %*% RM
	}
	return(A)
}
