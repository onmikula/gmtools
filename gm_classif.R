### TOOLS
library(geomorph)
source("func/bgprcomp.R")
source("func/predict_shapes.R")
source("func/make_links.R")
source("func/wireframe2d.R")
source("func/plot_lm.R")


### DATA
dshape <- read.delim("data/D_Beamys_shape.txt")
lshape <- read.delim("data/L_Beamys_shape.txt")
vshape <- read.delim("data/V_Beamys_shape.txt")
shapes <- cbind(dshape, lshape, vshape)
d <- length(unique(unlist(regmatches(colnames(shapes), gregexpr("[XYZ]", colnames(shapes))))))
views <- c(D=ncol(dshape), L=ncol(lshape), V=ncol(vshape)) / d
info <- read.delim("data/Beamys_info.txt", stringsAsFactors=FALSE)
info <- info[match(rownames(shapes), info$ID),]

color <- c(B1="#DF536B", B2="#61D04F", B3="#2297E6", B4="#28E2E5", B5="#CD0BBC", B6="#9E9E9E", B7="#F5C710")


### PCA
pca <- prcomp(shapes)
par(mai=c(1.02, 1.02, 0.42, 0.42))
plot(pca$x, pch=21, col=color[info$Lineage], cex=1.5, lwd=1.5, cex.lab=1.5, cex.axis=1.25)
legend("topleft", legend=names(color), pch=15, col=color, pt.cex=2.15, cex=1.1, bty="n")


### BETWEEN-GROUP PCA
subset <- info$Lineage %in% c("B1","B5","B7")
bgpca <- bgprcomp(shapes[subset,], info$Lineage[subset])
predicted <- predict(bgpca, newdata=shapes[!subset,])

# plot of between-group princinap components
# with three group means defining the axies, and individual shapes fron all groups projected on them
par(mai=c(1.02, 1.02, 0.42, 0.42))
plot(bgpca, i=1, j=2, col=1, bg=color, cex=1.5, cex.group=3, cex.lab=1.5, cex.axis=1.25)
points(predicted[,1:2], pch=24, col=color[info$Lineage[!subset]], cex=1.5, lwd=1.5)
legend("topright", legend=names(color), pch=15, col=color, pt.cex=2.15, cex=1.1, bty="n")

# cross-validated classification success
cv <- cv.bgprcomp(bgpca)
cv$success
cv$confmatrix

# explained variation of group means and individual shapes
100 * round(apply(bgpca$x, 2, var) / sum(diag(cov(bgpca$gmeans))), 3)
100 * round(apply(bgpca$x_ind, 2, var) / sum(diag(cov(bgpca$data))), 3)


# wireframe plot showing shape differences associated with unit shift along bgPC1
pcashapes <- predict_pca_shapes(bgpca, part=views, d=2)
plot_pca_shapes(pcashapes, pc=1, part=3, mag=3)

# wireframe plot showing shape differences between the most different group means
source("data/Beamys_links.R")
links <- lapply(list(D=d_links, L=l_links, V=v_links), make_links)

dst <- as.matrix(dist(bgpca$x))
dif <- sort(which(dst == max(dst), arr.ind=TRUE)[1,])
otus <- rownames(bgpca$x)[dif]
predshapes <- predict_shapes(score=bgpca$x[dif,], vector=bgpca$rotation[,1:ncol(bgpca$x)], mshape=bgpca$center, part=views)

npart <- length(predshapes)
h <- 1.05 * sapply(predshapes, function(x) diff(range(x[,2,])))
w <- 1.05 * max(sapply(predshapes, function(x) diff(range(x[,1,]))))
r <- sum(h) / w
layout(matrix(1:npart, npart, 1), heights=h/sum(h))
for (j in 1:npart) {
	wireframe2d(predshapes[[j]][,,1], predshapes[[j]][,,2], points=FALSE, links=links[[j]], col_links=color[otus], lwd_links=3, device=NULL)	
}




