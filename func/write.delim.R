write.delim <- function(x, file) {
	write.table(x, file, row.names=FALSE, quote=FALSE, sep="\t")
}
