make_links <- function(links) {
	do.call(rbind, lapply(links, function(x) cbind(x[-length(x)],x[-1])))
}

write_links <- function(links, file) {
	write.table(links, file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}

read_links <- function(file) {
	as.matrix(unname(read.table(file)))
}

