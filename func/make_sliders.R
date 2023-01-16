make_sliders <- function(x) {
	do.call(rbind, lapply(x, function(p) cbind(p[1:(length(p)-2)], p[-c(1,length(p))], p[-c(1:2)])))
}

write_sliders <- function(x, file) {
	write.table(x, file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}

read_sliders <- function(file) {
	as.matrix(unname(read.table(file)))
}
