### ARGUMENTS
# X - array or list of landmark configurations
# file - name of .TPS file
# sliders - three-column matrix specifying semilandmarks, if supplied overrides sliders specified in attribute 'semilm' of X
# missing - a numeric indicating missing data
# imagename - function modifying image names
# extension - file extension to be added to umage name


writeTPS <- function(X, file, sliders=NULL, missing=-1, imagename=identity, extension=NULL) {
	split_sliders <- function(s) {
		if (is.null(s)) {
			return(NULL)
		} else {
			breaks <- s[-nrow(s),2] != s[-1,1]
			if (any(breaks)) breaks <- which(breaks) else breaks <- NULL
			curves <- Map(seq, c(1,breaks+1), c(breaks, nrow(s)))
			return(lapply(curves, function(x) unique(as.numeric(t(s[x,])))))			
		}
	}
	writeInd <- function(m, semi, missing) {
		m[is.na(m)] <- missing
		m <- apply(formatC(m, format="f", digit=5), 1, paste, collapse=" ")
		nc <- length(semi)
		if (nc > 0) {
			curves <- lapply(seq(nc), function(i) c(paste0("POINTS=", length(semi[[i]])), m[semi[[i]]]))
			m <- c(m[-unlist(semi)], paste0("CURVES=", nc), unlist(curves))
		}
		return(m)
	}
  	if (inherits(X, "array")) {
  		if (is.null(sliders)) {
  			sliders <- attr(X, "semilm")
  		}
  		semilm <- list(split_sliders(sliders))
		X <- setNames(
			lapply(split(X, arrayInd(seq_along(X), dim(X))[,3]), array, dim=dim(X)[1:2], dimnames(X)[1:2]),
			dimnames(X)[[3]])
	} else if (is.list(X)) {
  		if (is.null(sliders)) {
  			semilm <- lapply(lapply(X, attr, "semilm"), split_sliders)
  		} else {
  			semilm <- list(split_sliders(sliders))
  		}
	}
	semilm <- rep_len(semilm, length(X))
	images <- paste0("IMAGE=", imagename(names(X)))
	if (!is.null(extension)) {
		images <- paste(images, extension, sep=".")
	}
	X <- Map(writeInd, X, semilm, missing)
	lm <- sapply(X, function(x) ifelse(any(grepl("CURVES",x)), grep("CURVES", x) - 1, length(x)))
	lm <- paste0("LM=", lm)
	X <- unlist(Map(c, Map(c, lm, X), images))
	writeLines(X, file)
}

