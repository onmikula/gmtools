# tps - the name of 'tps' file
# missing - specifies which coordinates are treated as missing data, with missing=0, every coordintae <= 0 is considered missing; use missing=-Inf to disable the option

check_lm_counts <- function(tps, missing=0) {
	LM <- readTPS(file=tps, missing=missing, asList=FALSE)
	if (!is.array(LM)) {
		nlm <- unname(lapply(LM, nrow))
		nlm[sapply(nlm, is.null)] <- 0
		nlm <- unlist(nlm) 
		tab <- table(nlm)
		common <- as.numeric(names(which.max(tab)))
		different <- nlm != common
		df <- data.frame(id=names(LM)[different], nlm=nlm[different], order=which(different))
		if (nrow(df) > 0) {
			cat(paste0("The most common count of landmarks is ", common, ".\n"))
			plural <- (sum(different) > 1) + 1
			specimen <- c("The specimen with a different count of landmarks is:\n",
			"The specimens with different counts of landmarks are:\n")[plural]
			cat(specimen)
			print(df)	
		} else {
			cat(paste("All configurations have", nrow(LM[[1]]), "landmarks,\nbut differ in the definition of semilandmarks.\n"))
		}
	} else {	 
		cat(paste("All configurations have", nrow(LM), "landmarks.\n"))
	}	
} 



# !todo - cases with some specimens without curves

check_semilm_curves <- function(tps, missing=0) {
	sort_semilm <- function(semi) {
		return(semi[order(semi[,2]),])
	}
	split_curves <- function(semi) {
		start <- c(1, which(semi[-nrow(semi),2] != semi[-1,1] & semi[-nrow(semi),2] != semi[-1,3]) + 1)
		return(by(semi, rep(seq_along(start), diff(c(start, nrow(semi) + 1))), unique))
	}

	LM <- readTPS(file=tps, missing=missing, asList=TRUE)
	N <- length(LM)
	semilm <- lapply(LM, attr, which="semilm")
	nsemilm <- sapply(semilm, nrow)
	curves <- lapply(lapply(semilm, sort_semilm), split_curves)
	ncurves <- sapply(curves, length)

	tsemilm <- table(nsemilm)
	comsemilm <- as.numeric(names(which.max(tsemilm)))
	difsemilm <- nsemilm != comsemilm
	if (any(difsemilm)) {
		dfsemilm <- data.frame(id=names(LM)[difsemilm], nsemilm=nsemilm[difsemilm], order=which(difsemilm))
		cat(paste0("The most common count of semilandmarks is ", comsemilm, ".\n"))
		plural <- (sum(difsemilm) > 1) + 1
		specimen <- c("The specimen with a different count of semilandmarks is:\n",
			"The specimens with different counts of semilandmarks are:\n")[plural]
		cat(specimen)
		print(dfsemilm)	
	}
		
	tcurves <- table(ncurves)
	comcurves <- as.numeric(names(which.max(tcurves)))
	difcurves <- ncurves != comcurves
	if (any(difcurves)) {
		dfcurves <- data.frame(id=names(LM)[difcurves], ncurves=ncurves[difcurves], order=which(difcurves))
		cat(paste0("The most common count of curves is ", comcurves, ".\n"))
		plural <- (sum(difcurves) > 1) + 1
		specimen <- c("The specimen with a different count of curves is:\n",
			"The specimens with different counts of curves are:\n")[plural]
		cat(specimen)
		print(dfcurves)	
	}

	if (!any(difsemilm) & !any(difcurves)) {
		curve_lengths <- do.call(rbind, lapply(curves, function(x) unname(sapply(x, nrow))))
		curvelen_patterns <- apply(curve_lengths, 1, paste, collapse="-")
		curvelen_patterns_tab <- table(curvelen_patterns)
		if (length(curvelen_patterns_tab) > 1) {
			comcurvelens <- names(curvelen_patterns_tab)[which.max(curvelen_patterns_tab)]
			difcurvelens <- curvelen_patterns != comcurvelens
			dfcurvelens <- data.frame(id=names(LM)[difcurvelens], pattern=curvelen_patterns[difcurvelens], order=which(difcurvelens))
			cat(paste0("The most common curve lengths are: ", comcurvelens, ".\n"))
			plural <- (sum(difcurvelens) > 1) + 1
			specimen <- c("The specimen with a different pattern of curve lengths is:\n",
			"The specimens with different different patterns of curve lengths are:\n")[plural]
			cat(specimen)
			print(dfcurvelens)	
		}
		
# different definition of curves
		curvedefinition <- function(x) paste(apply(x, 1, paste, collapse="-"), collapse="===")
		curve_definitions <- do.call(rbind, lapply(curves, function(x) sapply(x, curvedefinition)))
		curvedef_patterns <- unique(curve_definitions)
		if (nrow(curvedef_patterns) > 1) {
			curvedef_patterns_tab <- vector("list", nrow(curvedef_patterns))
			for (i in seq_along(curvedef_patterns_tab)) {
				curvedef_patterns_tab[[i]] <- which(apply(curve_definitions == curvedef_patterns[rep(i, N),], 1, all))
			}
			com <- which.max(sapply(curvedef_patterns_tab, length))
			messages <- vector("list", nrow(curvedef_patterns))
			for (i in setdiff(seq(nrow(curvedef_patterns)), com)) {
				specimens <- curvedef_patterns_tab[[i]]
				plural <- (length(specimens) > 1) + 1
				if (plural == 2) {
					specimens <- split(specimens, rep(c(1,2), c(length(specimens) - 1, 1)))
					specimens[[1]] <- paste(specimens[[1]], collapse=", ")
					specimens <- paste(specimens[[1]], "and", specimens[[2]])
				} else {
					specimens <- specimens
				}
				messages[[i]] <- paste(c("Specimen", "Specimens")[plural], specimens, c("differs", "differ")[plural])
				messages[[i]] <- paste(messages[[i]], "from the most common semilandmark definition in", c("its", "their")[plural])
				difcurves <- which(curvedef_patterns[i,] != curvedef_patterns[com,])
				plural <- (length(difcurves) > 1) + 1
				if (plural == 2) {
					difcurves <- split(difcurves, rep(c(1,2), c(length(difcurves) - 1, 1)))
					difcurves[[1]] <- paste(difcurves[[1]], collapse=", ")
					difcurves <- paste(difcurves[[1]], "and", difcurves[[2]])
				}
				difcurves <- paste0(difcurves, ".\n") 
				messages[[i]] <- paste(messages[[i]], c("curve", "curves")[plural], "no.", difcurves)
			}
			messages <- paste(messages[-com], collapse="")
			cat(messages)
		}
		
		if (length(curvelen_patterns_tab) == 1 & nrow(curvedef_patterns) == 1) {
			plural <- (ncurves[1] > 1) + 1
			cat(paste("All configurations share the same semilandmark", c("curves.\n", "curves.\n")[plural]))
		}

	}

}

