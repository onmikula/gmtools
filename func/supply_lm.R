supply_lm <- function(X, midplane, rightside, leftside, matching=FALSE) {
	supply_lm_core <- function(X, midplane, rightside, leftside, matching) {
		k <- ncol(X)
		d <- diag(c(rep(1, k-1), -1))
		if (sum(!is.na(X[midplane,1])) < 2) {
			return(X)
		} else if (isTRUE(matching)) {
			p <- nrow(X) / 2
			rna <- is.na(X[rightside,1])
			lna <- is.na(X[leftside,1])
			common <- !(rna | lna)
			rrot <- eigen(var(X[rightside[common],]))$vectors
			rrot[,k] <- sign(det(rrot)) * rrot[,k]
			Xr <- X[rightside,]
			rshift <- colMeans(Xr[common,])
			CXr <- Xr - rep(1, p) %*% t(rshift)
			RXr <- CXr %*% rrot		
			Xl <- X[leftside,]	%*% d	
			lshift <- colMeans(Xl[common,])
			CXl <- Xl - rep(1, p) %*% t(lshift)
			sv <- svd(t(RXr[common,]) %*% CXl[common,])
			lrot <- sv$v %*% t(sv$u)
			lrot[,k] <- sign(det(lrot)) * lrot[,k]
			RXl <- CXl %*% lrot
			if (length(rna) > 0) {
				RXr[rna,] <- RXl[rna,]
			}
			if (length(lna) > 0) {
				RXl[lna,] <- RXr[lna,]
			}
			CXr <- RXr %*% solve(rrot)
			Xr <- CXr + rep(1, p) %*% t(rshift)
			CXl <- RXl %*% solve(lrot)
			Xl <- (CXl + rep(1, p) %*% t(lshift)) %*% d
			return(list(right=Xr, left=Xl))
		} else {
			p <- nrow(X)
			shift <- colMeans(X[midplane,], na.rm=TRUE)
			rrot <- eigen(var(X[midplane,], na.rm=TRUE))$vectors
			rrot[,k] <- sign(det(rrot)) * rrot[,k]
			CX <- X - rep(1, p) %*% t(shift)
			RX <- CX %*% rrot
			rna <- which(is.na(RX[rightside,1]))
			lna <- which(is.na(RX[leftside,1]))
			if (length(rna) > 0) {
				RX[rightside[rna],] <- RX[leftside[rna],] %*% d
			}
			if (length(lna) > 0) {
				RX[leftside[lna],] <- RX[rightside[lna],] %*% d		
			}
			CX <- RX %*% solve(rrot)
			X <- CX + rep(1, p) %*% t(shift)
			return(X)	
		}
	}
	
	semilm <- attr(X, "semilm")
	if (is.array(X) & length(dim(X)) == 3) {
		if (isFALSE(matching)) {
			X <- apply(X, 3, supply_lm_core, midplane=midplane, rightside=rightside, leftside=leftside, matching=matching, simplify=FALSE)
			X <- array(unlist(X), dim=c(nrow(X[[1]]), ncol(X[[1]]), length(X)), dimnames=list(NULL, NULL, names(X)))
		} else {
			warning("not implemented for an array input with the option matching=TRUE")
		}
	} else {
		if (isFALSE(matching)) {
			X <- supply_lm_core(X=X, midplane=midplane, rightside=rightside, leftside=leftside, matching=FALSE)
		} else {
			X <- supply_lm_core(X=X, rightside=rightside, leftside=leftside, matching=TRUE)
		}
	}
	attr(X, "semilm") <- semilm

	return(X)
}

