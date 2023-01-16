export_procrustes <- function(P, file=NULL, return=TRUE, size=TRUE, prefix=NULL) {
	arr2mat <- function(arr) {
		t(matrix(aperm(arr, c(2,1,3)), nrow(arr)*ncol(arr), dim(arr)[3]))
	}
	form <- cbind(log(P$Csize), arr2mat(P$coords))
	rownames(form) <- dimnames(P$coords)[[3]]
	lmlab <- formatC(seq(nrow(P$coords)), format="d", flag=0, width=nchar(nrow(P$coords)))
	lmlab <- paste0(ifelse(is.null(prefix), "X", prefix), rep(lmlab, each=ncol(P$coords)), c("X","Y","Y")[seq(ncol(P$coords))])
	cslab <- paste0(ifelse(is.null(prefix), "", prefix), "CS")
	colnames(form) <- c(cslab, lmlab)
	if (!isTRUE(size)) {
		form <- form[,-1]
	} 
	if (!is.null(file)) {
		write.table(form, file=file, row.names=TRUE, quote=FALSE, sep="\t")
	}
	if (isTRUE(return)) {
		return(form)
	}
}



import_procrustes <- function(file) {
	as.matrix(read.delim(file))
}
