filter_ind <- function(A, filter, subset) {
	semilm <- attr(A, "semilm")
	if (missing(filter)) {
		filter <- function(x) apply(!is.na(x), 3, all)
	}
	if (missing(subset)) {
		subset <- filter(A)
	}
	if (is.logical(subset)) {
		ord <- which(subset)
	} else if (is.character(subset)) {
		ord <- match(subset, dimnames(A)[[3]])
	}
	A <- A[,,subset]
	attr(A, "order") <- ord
	attr(A, "semilm") <- semilm
	return(A)
}

filter_lm <- function(A, mask) {
	if (missing(mask)) {
		return(A)
	} else if (!any(is.na(mask))) {
		return(A)
	} else {
		out <- which(apply(is.na(mask), 1, any))
		semilm <- attr(A, "semilm")
		if (!is.null(semilm)) {
			semilm[semilm %in% out] <- NA
			semilm_out <- which(apply(is.na(semilm), 1, any))
			if (length(semilm_out) > 0) {
				out <- sort(unique(c(out, na.omit(semilm[semilm_out, 2]))))
				semilm <- semilm[-semilm_out,,drop=FALSE]
			}
			semilm <- matrix(match(semilm, setdiff(seq(nrow(A)), out)), nrow(semilm), 3)
		}
		if (length(dim(A)) == 2) {
			A <- array(as.numeric(A), dim=c(dim(A),1))
		}
		A <- A[-out,,,drop=FALSE]
		attr(A, "semilm") <- semilm
		return(A)
	}
}

filter_sliders <- function(semilm, mask) {
	out <- which(apply(is.na(mask), 1, any))
	semilm[semilm %in% out] <- NA
	new <- setdiff(unique(na.omit(semilm[apply(is.na(semilm), 1, any),2])), out)
	while (length(new) > 0) {
		out <- sort(c(out, new))
		semilm[semilm %in% out] <- NA
		new <- setdiff(unique(na.omit(semilm[apply(is.na(semilm), 1, any),2])), out)
	}
	semilm_out <- which(apply(is.na(semilm), 1, any))
	if (length(semilm_out) > 0) {
#		out <- sort(unique(c(out, na.omit(semilm[semilm_out, 2]))))
		semilm <- semilm[-semilm_out,,drop=FALSE]
	}
	semilm <- matrix(match(semilm, setdiff(seq(nrow(mask)), out)), nrow(semilm), 3)	
	return(semilm)
}

filter_links <- function(links, mask) {
	out <- which(apply(is.na(mask), 1, any))
	links[links %in% out] <- NA
	links <- na.omit(links)
	links <- matrix(match(links, setdiff(seq(nrow(mask)), out)), nrow(links), 2)
	return(links)
}
