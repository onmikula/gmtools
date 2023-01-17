### TOOLS
library(geomorph)
source("func/readTPS.R")
source("func/writeTPS.R")
source("func/format_conversion.R")
source("func/supply_lm.R")
source("func/unify_reflection.R")
source("func/symmetrize.R")
source("func/filtering.R")
source("func/make_scaling.R")
source("func/make_sliders.R")
source("func/export_procrustes.R")
source("func/plot_lm.R")
source("func/write.delim.R")



### IMPUTATION OF MISSING LANDMARKS
D <- readTPS("data/D_Beamys.TPS", missing=20, asList=FALSE)
Dsuppl <- supply_lm(D, midplane=1:5, rightside=6:37, leftside=38:69)
writeTPS(Dsuppl, "data/D_Beamys_imputed.TPS", extension="JPG")
V <- readTPS("data/V_Beamys.TPS", missing=20, asList=FALSE)
Vsuppl <- supply_lm(V, midplane=1:5, rightside=6:51, leftside=52:97)
writeTPS(Vsuppl, "data/V_Beamys_imputed.TPS", extension="JPG")

# plot highlighting the imputed data
par(mfrow=c(2,1), mai=c(0.82, 0.82, 0.22, 0.22))
plotlm(D[,,1])
plotlm(Dsuppl[,,1], col=apply(is.na(D[,,1]),1,any)+1)


### FILTERING OUT MISSING DATA
# filtering out all individuals with any missing data
D <- readTPS("data/D_Beamys_imputed.TPS", rm="^D_")
D <- filter_ind(D)
V <- readTPS("data/V_Beamys_imputed.TPS", rm="^V_")
V <- filter_ind(V)
L <- readTPS("data/L_Beamys.TPS", rm="^L_")
L <- filter_ind(L)

# finding individuals with no missing data in all views of the skull (individual IDs in 'complete') & retaining just these individuals in all data sets
complete <- Reduce(intersect, list(dimnames(D)[[3]], dimnames(L)[[3]], dimnames(V)[[3]]))
D <- filter_ind(D, subset=complete)
V <- filter_ind(V, subset=complete)
L <- filter_ind(L, subset=complete)


### SCALING
# re-scales the landmark configurations from pixels to units implicit in 'scale' argument of 'get_scale_factor' function
# the scaling is based on either two-landmark configurations setting the scale ("D" and "V") or on the shared landmarks ("L" is scaled according to the already scaled "D")
dscale <- get_scale_factor(tps="data/D_Beamys_scaling.TPS", rm="D_", scale=10)
D <- set_scale_factor(X=D, scale=dscale)
vscale <- get_scale_factor(tps="data/V_Beamys_scaling.TPS", rm="V_", scale=10)
V <- set_scale_factor(X=V, scale=vscale)
L <- set_scale_factor(X=D, Y=L, lm=list(c(1,5),c(1,53)))


### SYMMETRIZATION & REFLECTION
D <- symmetrize(D, midplane=1:5, rightside=6:37, leftside=38:69)
V <- symmetrize(V, midplane=1:5, rightside=6:51, leftside=52:97)
L <- unify_reflection(L)


### SUBSETTING
# for instance, retain just individuals of adult age
info <- read.delim("data/Beamys_info.txt", stringsAsFactors=FALSE)
info <- info[info$ID %in% dimnames(D)[[3]],]
age <- info$Age == "adult"
D <- filter_ind(D, subset=info$ID[age])
V <- filter_ind(V, subset=info$ID[age])
L <- filter_ind(L, subset=info$ID[age])


### PROCRUSTES SUPERIMPOSITION
# prepare objects defining semilandmark curves
d_semilm <- list(c(6:11,1), c(8,12:14), 15:18, 18:20, 25:29, c(38:43,1), c(40,44:46), 47:50, 50:52, 57:61)
d_semilm <- make_sliders(d_semilm)
v_semilm <- list(6:8, 8:11, 14:16, 16:21, 22:27, 37:39, c(43:46,37), 47:49, 52:54, 54:57, 60:62, 62:67, 68:73, 83:85, c(89:92,83), 93:95)
v_semilm <- make_sliders(v_semilm)     
l_semilm <- list(16:23, 32:46, 47:57)
l_semilm <- make_sliders(l_semilm)

d <- geomorph::gpagen(D, curves=d_semilm, print.progress=FALSE)
l <- geomorph::gpagen(L, curves=l_semilm, print.progress=FALSE)
v <- geomorph::gpagen(V, curves=v_semilm, print.progress=FALSE)

# exports Procrustes shape coordinates in tab-delimited text file,
# possibly with with log(centroid size) in the first colmun if size=TRUE
dshape <- export_procrustes(d, "data/D_Beamys_shape.txt", size=FALSE, prefix="D")
lshape <- export_procrustes(l, "data/L_Beamys_shape.txt", size=FALSE, prefix="L")
vshape <- export_procrustes(v, "data/V_Beamys_shape.txt", size=FALSE, prefix="V")

dform <- export_procrustes(d, "data/D_Beamys_form.txt", size=TRUE, prefix="D")
lform <- export_procrustes(l, "data/L_Beamys_form.txt", size=TRUE, prefix="L")
vform <- export_procrustes(v, "data/V_Beamys_form.txt", size=TRUE, prefix="V")


