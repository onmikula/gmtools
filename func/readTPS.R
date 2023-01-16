### ARGUMENTS
# file - name of .TPS file, if it contains CURVES, their boundary points are interpreted as landmarks, while the others as semilandmarks 
# asList - specifies whether the output should be a list in any case; if FALSE (default) and the number of landmarks & semilandmarks is the same in all objects, the output is an array of dimension c(p, k, n), where p=no. of points, k=dimensionality (2 or 3) and n=no. of objects
# missing - numeric, coordinate values <= missing are converted to NA; use missing=-Inf to disable the option
# rm - what to remove from the image file names to get specimen names
# extension - whether to retain extension of the image file names

readTPS <- function(file, asList=FALSE, missing=0, rm=NULL, extension=FALSE) {
	makeSliders <- function(x) {
		makeMatrix <- function(p) cbind(p[1:(length(p)-2)], p[-c(1,length(p))], p[-c(1:2)])
		coords <- grep("^[[:digit:]]|-", x)
		points <- grep("^POINTS", x)
		npoints <- as.numeric(gsub("^POINTS=", "", x[points]))
		semilm <- Map("+", as.list(points), lapply(npoints, seq))
		semilm <- lapply(semilm, match, table=coords)
		semilm <- semilm[sapply(semilm, length) >= 3]
		return(do.call(rbind, lapply(semilm, makeMatrix)))
	}
	readIndiv <- function(x, missing) {
		coords <- grep("^[[:digit:]]|-", x)
		if (length(coords) > 0) {
			lm <- do.call(rbind, lapply(strsplit(x[coords], "\\s+"), as.numeric))
			lm[which(lm <= missing)] <- NA
			attr(lm, "semilm") <- makeSliders(x)
			return(lm)
		} else {
			return(NULL)
		}
	}
	list2arr <- function(x) {
		nz <- which(sapply(x, length) > 0)[1]
		array(unlist(x), dim=c(nrow(x[[nz]]), ncol(x[[nz]]), length(x)), dimnames=list(NULL, NULL, names(x)))
	}
	make_names <- function(x) {
		list(paste0("L", formatC(seq(nrow(x)), format="d", flag=0, width=nchar(nrow(x)))), c("X","Y","Y")[seq(ncol(x))])
	}

	lines <- readLines(file)
	lines <- lines[nchar(lines) > 0]
	start <- grep("^LM", lines)
	image <- grep("^IMAGE", lines)
	names <- gsub("^IMAGE=", "", lines[image])
	if (isFALSE(extension)) {
		rm <- paste0(rm, "|\\.[[:alpha:]]+$")
	}
	if (!is.null(rm)) {
		names <- gsub(rm, "", names)
	}
	ids <- grep("^ID", lines)
	scale <- grep("^SCALE", lines)
	units <- as.numeric(gsub("^SCALE=", "", lines[scale]))
	which <- sapply(scale, function(x) sum(x > start))
	units <- units[match(seq_along(start), which)]
	units[is.na(units)] <- 1
	out <- c(image, scale, ids)
	lines <- lines[-out]
	start <- start - sapply(start, function(x) sum(out < x))
	
	indiv <- setNames(split(lines, rep(start, diff(c(start, length(lines)+1)))), names)
	config <- lapply(indiv, readIndiv, missing=missing)
	lm <- unlist(ifelse(sapply(config, length) == 0, 0, lapply(config, nrow)))
	semilm <- lapply(config, attr, which="semilm")
	if (length(scale) > 0) {
		config <- Map("*", config, units)
	}
	if (length(setdiff(unique(lm), 0)) == 1 & length(unique(semilm)) == 1 & isFALSE(asList)) {
		if (any(lm == 0)) {
			nr <- ncol(config[[which(lm > 0)[1]]])
			nc <- ncol(config[[which(lm > 0)[1]]])
			for (i in which(lm == 0)) {
				config[[i]] <- matrix(NA, nr, nc) 
			}	
		}
		config <- list2arr(config)
		dimnames(config) <- c(make_names(config), dimnames(config)[3])
		attr(config, "semilm") <- semilm[[1]]
	} else {
		if (any(lm == 0)) {
			for (i in which(lm > 0)) {
				dimnames(config[[i]]) <- make_names(config[[i]]) 
			}
		}
	}

	return(config)
}

