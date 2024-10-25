# ARGUMENTS
# M1, M2: matrices with landmark configurations to be compared, M1 can be an array of dimension p x k x 2
# rotat: rotation matrix
# points: whether to show landmarks
# links: matrix indicating landmark pairs to be joined by links
# texts: either logical (whether to show landmark names, i.e. rownames of M1, M2) or character (annotations of landmarks)
# subset: which landmarks should be displayed
# bg: which landmarks should serve as a background
# col_points: colors of landmarks
# bg_points: background colors of landmarks
# lwd_points: line width of landmarks
# cex_points: size of landmarks
# col_links: color(s) of links
# lwd_links: width(s) of links
# col_bg: color(s) of background links
# lwd_bg: width(s) of background links
# col_texts: color(s) of text
# cex_texts: size of text
# font_texts: size of text
# box: whether to display an enclosing box
# device: "quartz", "x11" (or "X11") or "pdf", if NULL, the objects are plotted into the current device
# file: name of pdf file if device == "pdf"
# width: width of the device
# height: height of the device, if NULL (default), it is set automatically
# mai: margins in inches
# expr: an expression whose evoluation adds a new element (e.g. an arrow) to the wireframe plot


wireframe2d <- function(M1, M2=NULL, rotat=c(0,0,0), points=TRUE, links=NULL, texts=FALSE, subset=NULL, bg=NULL, col_points=c(8,1), bg_points="white", lwd_points=1, cex_points=1, col_links=1, lwd_links=1, col_bg="lightgrey", lwd_bg=1, col_texts=1, cex_texts=1, font_texts=1, box=FALSE, device, file, add=FALSE, width=7, height=NULL, mai=rep(0,4), expr=NULL) {
	if (length(dim(M1)) == 3) {
		if (dim(M1)[3] > 1) {
			M2 <- M1[,,2]			
		}
		M1 <- M1[,,1]
	}
	single <- is.null(M2)
	k <- ncol(M1)
	if (is.vector(rotat)) {
		if (k == 3) {
			rx <- matrix(c(1, 0, 0, 0, cos(rotat[1]), sin(rotat[1]), 0, -sin(rotat[1]), cos(rotat[1])), 3, 3)
			ry <- matrix(c(cos(rotat[2]), 0, -sin(rotat[2]), 0, 1, 0, sin(rotat[2]), 0, cos(rotat[2])), 3, 3)
			rz <- matrix(c(cos(rotat[3]), sin(rotat[3]), 0, -sin(rotat[3]), cos(rotat[3]), 0, 0, 0, 1), 3, 3)
			rotat <- rz %*% ry %*% rx
		} else if (k == 2) {
			rotat <- matrix(c(cos(rotat[1]), sin(rotat[1]), - sin(rotat[1]), cos(rotat[1])), 2, 2)
		} else {
			rotat <- diag(k)
			warning("rotation in more than three dimensions is not implemented")
		}
	}
	M1c <- scale(M1, scale=FALSE)
	M1 <- M1c %*% rotat + rep(1,nrow(M1)) %*% t(colMeans(M1))
	M1 <- M1[,1:2]
	if (isFALSE(single)) {
		M2c <- scale(M2, scale=FALSE)
		M2 <- M2c %*% rotat + rep(1,nrow(M2)) %*% t(colMeans(M2))
		M2 <- M2[,1:2]
	}

	if (missing(device)) {
		if (.Platform$OS.type == "unix") {
			device <- "quartz"
		} else {
			device <- "x11"
		}
	} else if (isTRUE(device == "pdf") & missing(file)) {
		file <- "wireframe.pdf"
	}

	if (!is.null(device)) {
		w <- range(M1[,1])
		h <- range(M1[,2])
		if (isFALSE(single)) {
			w <- range(c(w, M2[,1]))
			h <- range(c(h, M2[,2]))
		}
		if (is.null(height)) {
			height <- width * diff(h) / diff(w)
		}
		if (device == "pdf") {
			pdf(file, width=width, height=height)
		} else {
			match.fun(tolower(device))(width=width, height=height)
		}
		par(mai=mai)
		plot(cbind(w, h), asp=1, type="n", axes=FALSE, ann=FALSE, bty=c("n","o")[box + 1])
	} else if (!isTRUE(add)) {
		w <- range(M1[,1])
		h <- range(M1[,2])
		if (isFALSE(single)) {
			w <- range(c(w, M2[,1]))
			h <- range(c(h, M2[,2]))
		}
		par(mai=mai)
		plot(cbind(w, h), asp=1, type="n", axes=FALSE, ann=FALSE, bty=c("n","o")[box + 1])
	}
	
	if (!is.null(subset)) {
		if (!is.list(subset)) {
			subset <- list(subset)
		}
	} else {
		subset <- list(seq(nrow(M1)))
	}
	subset <- rep_len(subset, (!single) + 1)
	subset <- lapply(subset, sort)
	if (!is.null(bg)) {
		if (!is.list(bg)) {
			bg <- list(bg)
		}
		col_bg <- rep_len(col_bg, (!single) + 1)
		lwd_bg <- rep_len(lwd_bg, (!single) + 1)
		bg <- rep_len(bg, (!single) + 1)
		bg <- Map(intersect, bg, subset)
	}
	if (!is.null(links)) {
		links <- rep(list(links), (!single) + 1)
		for (i in seq_along(links)) {
			links[[i]] <- links[[i]][apply(apply(links[[i]], 2, "%in%", subset[[i]]), 1, all),,drop=FALSE]
		}	
		col_links <- rep_len(as.list(col_links), (!single) + 1)
		lwd_links <- rep_len(as.list(lwd_links), (!single) + 1)
		for (i in seq_along(col_links)) {
			if (!is.null(bg)) {
				col_links[[i]] <- ifelse(apply(apply(links[[i]], 2, "%in%", bg[[i]]), 1, any), col_bg[i], col_links[[i]])
				lwd_links[[i]] <- ifelse(apply(apply(links[[i]], 2, "%in%", bg[[i]]), 1, any), lwd_bg[i], lwd_links[[i]])
			} else {
				col_links[[i]] <- rep(col_links[[i]], nrow(links[[i]]))
				lwd_links[[i]] <- rep(lwd_links[[i]], nrow(links[[i]]))
			}
		}
	}
	if (isTRUE(points)) {
		points_subs <- list(setdiff(subset[[1]], bg[[1]]))
		landmarks <- M1[points_subs[[1]],,drop=FALSE]
		if (isFALSE(single)) {
			points_subs <- c(points_subs, list(setdiff(subset[[2]], bg[[2]])))
			landmarks <- rbind(landmarks, M2[points_subs[[2]],,drop=FALSE])
		}
		col_points <- rep(rep_len(col_points, (!single) + 1), sapply(points_subs, length))
		bg_points <- rep(rep_len(bg_points, (!single) + 1), sapply(points_subs, length))
		lwd_points <- rep(rep_len(lwd_points, (!single) + 1), sapply(points_subs, length))
		cex_points <- rep(rep_len(cex_points, (!single) + 1), sapply(points_subs, length))
		points(landmarks, pch=21, col=col_points, bg=bg_points, lwd=lwd_points, cex=cex_points)
	}
	if (!is.null(links)) {
		if (isFALSE(single)) {
			ii <- setdiff(seq(nrow(links[[2]])), which(col_links[[2]] %in% c("NA","transparent","white")))
			for (i in ii) {
				lines(rbind(M2[links[[2]][i, 1], ], M2[links[[2]][i,2], ]), col=col_links[[2]][i], lwd=lwd_links[[2]][i])
			}
		}
		ii <- setdiff(seq(nrow(links[[1]])), which(col_links[[1]] %in% c("NA","transparent","white")))
		for (i in ii) {
			lines(rbind(M1[links[[1]][i, 1], ], M1[links[[1]][i,2], ]), col=col_links[[1]][i], lwd=lwd_links[[1]][i])
		}
	}
	if (isTRUE(texts[1]) | mode(texts) != "logical") {
		texts <- list(texts, rownames(M1))[[(mode(texts) == "logical") + 1]]
		text(M1, texts=texts, col=col_texts[seq(nrow(M1))], cex=cex_texts, font=font_texts)
	}
	if (!is.null(expr)) {
		if (is.expression(expr)) {
			eval(expr)
		} else {
			warning("expr' argument cannot be evaluated as it is not of class 'expression'")
		}
	}

	if (isTRUE(device == "pdf")) {
		invisible(dev.off())
	}

}

