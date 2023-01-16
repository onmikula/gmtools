plotlm <- function(X, col=1, cex=0.7, font=1, ...) {
	par(...)
	plot(X, asp=1, type="n", xlab="X", ylab="Y")
	text(X, labels=seq(nrow(X)), col=col, cex=cex, font=font)
}
