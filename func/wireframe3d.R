# ARGUMENTS
# M1, M2: matrices with landmark configurations to be compared, M1 may be an array of dimension p x k x 2
# mag: magnification factor for actual differences
# colors: indicates colors of landmarks 
# cex: size of landmarks and text
# links: matrix indicating landmark pairs to be joined by links
# points: whether to show landmarks
# texts: either logical (whether to show landmark names, i.e. rownames of M1, M2) or character (annotations of landmarks)
# subset: which landmarks should be displayed
# bg: which landmarks should serve as a background
# add: whether to add the specified into an opened interactive plot
# rotat: specifies rotation of configurations along x-, y-, z- axes 
# view: specifies point of view
# userMatrix: userMatrix specification of a point of view (overrides the first two elements of view)
# windowRect: size of window
# box: whether to display an enclosing box
# col.links: color(s) of links
# col.bg: color(s) of background links
# lwd.links: width(s) of links
# lwd.bg: width(s) of background links

wireframe3d <- function (M1, M2=NULL, mag=1, colors=c(8,1), cex=1, links=NULL, points=TRUE, texts=FALSE, 
subset=NULL, bg=NULL, add=FALSE, rotat=NULL, view=c(0,15,60,1), userMatrix=NULL, windowRect=NULL, box=FALSE,
col.links=1, col.bg="lightgrey", lwd.links=1, lwd.bg=1, ...) {
	if (length(dim(M1)) == 3) {
		if (dim(M1)[3] > 1) {
			M2 <- M1[,,2]			
		}
		M1 <- M1[,,1]
	}
	single <- is.null(M2)
	k <- ncol(M1)
	mag <- mag - 1
	if (is.null(rotat)) {
		rotat <- rep(1,3) 
	}
	if (is.vector(rotat)) {
		rotat <- diag(rotat)
	}
	M1 <- M1 %*% rotat
	if (isFALSE(single)) {
		M2 <- M2 + (M2 - M1) * mag
		M2 <- M2 %*% rotat
	}
	
	if (add == FALSE) {
		if (!is.null(windowRect)) {
			rgl::open3d(windowRect=windowRect)			
		} else {
			rgl::open3d()			
		}
		rgl::plot3d(rbind(M1, M2), type="n", size=1.25, aspect=FALSE, xlab="", ylab="", zlab="", box=box, axes=box, ...)
		if (!is.null(userMatrix)) {
			rgl::rgl.viewpoint(userMatrix=userMatrix, fov=view[3], zoom=view[4])
		} else {
			rgl::rgl.viewpoint(theta=view[1], phi=view[2], fov=view[3], zoom=view[4])
		}
	}
	
	r <- mean(sqrt(rowSums(scale(M1, scale=FALSE)^2)), na.rm=TRUE) / 20
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
		col.bg <- rep_len(col.bg, (!single) + 1)
		lwd.bg <- rep_len(lwd.bg, (!single) + 1)
		bg <- rep_len(bg, (!single) + 1)
		bg <- Map(intersect, bg, subset)
	}
	if (!is.null(links)) {
		links <- rep(list(links), (!single) + 1)
		for (i in seq_along(links)) {
			links[[i]] <- links[[i]][apply(apply(links[[i]], 2, "%in%", subset[[i]]), 1, all),,drop=FALSE]
		}	
		col.links <- rep_len(as.list(col.links), (!single) + 1)
		lwd.links <- rep_len(as.list(lwd.links), (!single) + 1)
		for (i in seq_along(col.links)) {
			if (!is.null(bg)) {
				col.links[[i]] <- ifelse(apply(apply(links[[i]], 2, "%in%", bg[[i]]), 1, any), col.bg[i], col.links[[i]])
				lwd.links[[i]] <- ifelse(apply(apply(links[[i]], 2, "%in%", bg[[i]]), 1, any), lwd.bg[i], lwd.links[[i]])
			} else {
				col.links[[i]] <- rep(col.links[[i]], nrow(links[[i]]))
				lwd.links[[i]] <- rep(lwd.links[[i]], nrow(links[[i]]))
			}
		}
	}
	spheres_subs <- list(setdiff(subset[[1]], bg[[1]]))
	spheres <- M1[spheres_subs[[1]],,drop=FALSE]
	if (isFALSE(single)) {
		spheres_subs <- c(spheres_subs, list(setdiff(subset[[2]], bg[[2]])))
		spheres <- rbind(spheres, M2[spheres_subs[[2]],,drop=FALSE])
	}
	colors <- rep(rep_len(colors, (!single) + 1), sapply(spheres_subs, length))
	if (isTRUE(points)) {
		rgl::spheres3d(spheres, radius=r*cex, col=colors, ...)
	}
	if (!is.null(links)) {
		if (isFALSE(single)) {
			ii <- setdiff(seq(nrow(links[[2]])), which(col.links[[2]] %in% c("NA","transparent","white")))
			for (i in ii) {
				rgl::lines3d(rbind(M2[links[[2]][i, 1], ], M2[links[[2]][i,2], ]), col=col.links[[2]][i], lwd=lwd.links[[2]][i], ...)
			}
		}
		ii <- setdiff(seq(nrow(links[[1]])), which(col.links[[1]] %in% c("NA","transparent","white")))
		for (i in ii) {
			rgl::lines3d(rbind(M1[links[[1]][i, 1], ], M1[links[[1]][i,2], ]), col=col.links[[1]][i], lwd=lwd.links[[1]][i], ...)
		}
	}
	if (texts[1] == TRUE | mode(texts) != "logical") {
		texts <- list(texts, rownames(M1))[[(mode(texts) == "logical") + 1]]
		rgl::text3d(M1, texts=texts, col=colors[seq(nrow(M1))], cex=cex, ...)
	}
}

