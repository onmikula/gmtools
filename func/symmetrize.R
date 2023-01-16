symmetrize <- function(X, midplane, rightside, leftside, matching=FALSE) {
	sym_one_config <- function(X, midplane, rightside, leftside, matching) {
		k <- ncol(X)
		if (isTRUE(matching)) {
			if (is.list(X) & length(X) == 2) {
				if (!is.null(names(X))) {
					X <- X[order(names(X))]	
				}
				L <- X[[1]]
				R <- X[[2]]
			} else {
				L <- X[leftside,]
				R <- X[rightside,]
			}
		} else {
			L <- X
			R <- X
			L[rightside,] <- R[leftside,]	
			L[leftside,] <- R[rightside,]
			L <- L %*% diag(c(1,-1,-1)[seq(k)])
		}
		p <- nrow(R)
		rshift <- colMeans(R)
		rrot <- eigen(var(R))$vectors
		rrot[,k] <- sign(det(rrot)) * rrot[,k]
		CR <- R - rep(1, p) %*% t(rshift)
		RR <- CR %*% rrot		
		lshift <- colMeans(L)
		CL <- L - rep(1, p) %*% t(lshift)
		CL <- CL %*% diag(c(1,-1,-1)[seq(k)])
		sv <- svd(t(RR) %*% CL)
		lrot <- sv$v %*% t(sv$u)
		#lrot[,k] <- sign(det(lrot)) * lrot[,k]
		RL <- CL %*% lrot
		RX <- 0.5 * (RR + RL)
		CX <- RX %*% solve(rrot)
		X <- CX + rep(1, p) %*% t(rshift)
		return(X)		
	}

	semilm <- attr(X, "semilm")
	arr <- length(dim(X))
	if (is.null(arr)) {
		stop("X must be an array of dimension 2 or 3")
	} else if (arr == 2) {
		S <- sym_one_config(X=X, midplane=midplane, rightside=rightside, leftside=leftside, matching=matching)
	} else if (arr == 3) {
		S <- apply(X, 3, sym_one_config, midplane=midplane, rightside=rightside, leftside=leftside, matching=matching, simplify=FALSE)
		S <- array(unlist(S), dim=c(nrow(S[[1]]), ncol(S[[1]]), length(S)), dimnames=dimnames(X))
	}
	attr(S, "semilm") <- semilm

	return(S)

}

