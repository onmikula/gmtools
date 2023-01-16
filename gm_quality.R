### TOOLS
library(geomorph)
source("func/readTPS.R")
source("func/format_conversion.R")
source("func/supply_lm.R")
source("func/filtering.R")
source("func/check_lm_counts.R")
source("func/check_lm_outliers.R")
source("func/predict_pca_shapes.R")
source("func/unify_reflection.R")
source("func/make_sliders.R")
source("func/plot_lm.R")



### TPS FILE NAME
tps <- "data/L_Beamys.TPS"
#tps <- "data/V_Beamys.TPS"


### NO. OF LANDMARKS
# First, it is tested, whether the number of landmarks is identical in all specimens.
# If it is not, the most common number of landmarks is printed
# along with the list of specimens deviating from it.
# 'check_lm_counts' checks points put as landmarks in tpsDIG
# 'check_semilm_curves' does the same with curves created in the same software

check_lm_counts(tps, missing=20)
check_semilm_curves(tps, missing=20)


### OUTLIER DETECTION
# Three methods of outlier detection are applied, .
# all analyzing distribution of Procrustes shape coordinates.

# Procrustes superimposition requires filtering out incomplete configurations (by 'filter_ind').
# Some data sets may require also unify reflection of the objects (by 'unify_reflection'), namely those showing so called matching symmetry, e.g., lateral views of the skulls or mandibles.

# Optionally, you can apply 'supply_lm' or 'filter_lm' prior to filtering.
# - 'supply_lm' supplies missing landmarks as the mirror image of their left-right counterparts, which allows inclusion of some incomplete specimens, but requires a case-specific specification of midplane and left-/right-side landmarks (the latter in the corresponding order).
# - 'filter_lm' allows inspection of a subset of landmarks, which also means some incomplete specimens can be checked.
# You can also introduce semilandmark sliding at this stage (argument 'curves' in 'gpagen'),
# which may (depending on circumstances) highlight, but also obscure landmarking errors. The matrix specifying curves can be taken from the 'semi' attribute of the landmark array or created using 'make_sliders' function

LM <- readTPS(tps, missing=20, asList=FALSE)
#LM <- supply_lm(LM, midplane=1:5, rightside=6:50, leftside=51:95)
#LM <- filter_lm(LM, mask=LM[,,1])
LM <- filter_ind(LM)
LM <- unify_reflection(LM)

curves <- attr(LM, "semi")
#semilm <- list(6:8, 8:11, 14:16, 16:21, 22:27, 37:39, c(43:46,37), 47:49, 52:54, 54:57, 60:62, 62:67, 68:73, 83:85, c(89:92,83), 93:95)
#curves <- make_sliders(semilm)     


### PROCRUSTES OUTLIERS
# the Procrustes-superimposed configurations can be shown as clouds of points, scattered around the average landmark positions
# landmarking-errors manifest themselves here as points lying out of the clouds
# 'plot_lm_outliers shows the clouds along with landmark numbers (at the avergae positions) and highlits the most remote outlier from each of the clouds
# based on this plot landmarks with possible errors can be chosen and supplied to 'lm' argument of 'print_lm_outliers', which lists the number of landmark, id of the specimen the outlier belongs to and, optionally, the order of the specimen in the original .TPS file 

P <- geomorph::gpagen(LM, print.progress=FALSE, curves=curves)
plot_lm_outliers(P)
print_lm_outliers(P, lm=c(27,85), order=attr(LM, "order"))



### PCA ANALYSIS
# Some small-scale mismatches may not be obvious in the gpagen plot,
# although they entirely distort the ordination of specimens in PCA plots
# The other strategy is thus to run PCA and look for PC score outliers.
# pcs indicates the components whose pairwise plots are to be created
# ij indicates the selected pair of components for more detailed look 
# labels are either specimen IDs or their order in the original .TPS file.

P <- geomorph::gpagen(LM, print.progress=FALSE, curves=curves)
pca <- prcomp(arr2mat(P$coords))
labels <- rownames(pca$x)
# labels <- attr(LM, "order")

pcs <- 1:7
pairs(pca$x[,pcs])

ij <- 5:6
plot(pca$x[,ij], asp=1, type="n");
text(pca$x[,ij], labels=labels, cex=0.7)



### TPS GRIDS
# Sometimes the ordination in PC space shows moderate outliers,
# so it is doubtful if it relects a landmarking error or just an extreme of natural variation.
# In that case, TPS deformation grid can reveal the nature of underlying shape variation.
# The landmarking error would likely result in some oddity,
# e.g., a suspicious asymmetry or a grid 'warped into itself'.
# 'i' indicates i-th PC; 'mag' is an optional magnification of shape differences

P <- geomorph::gpagen(LM, print.progress=FALSE, curves=curves)
pcshapes <- predict_pca_shapes(P)
plot_pca_shapes(pcshapes, i=2, mag=3)


