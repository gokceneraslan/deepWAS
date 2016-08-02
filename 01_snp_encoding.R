
library(snpStats)
library(rhdf5)
library(SNPRelate)
library(ggplot2)
library(dplyr)
library(readr)
library(bigmemory)

#significant snps
funsig <- read.table('gsk_17Msnps_significant_overlap.tsv', header = F, stringsAsFactors = F)$V1
nfunsig <- length(funsig)

# Read plink files --------------------------------------------------------
pl <- read.plink('../../GSK_janine_imp07/data/all.GSK.best_guessed_07.qc_geno0')
inds <- pl$fam$member
ninds <- length(inds)

inds.outliers <- read.table('../../GSK_janine_imp07/mds/all_outliers.txt', header = F, stringsAsFactors = F)$V1
inds <- setdiff(inds, inds.outliers)
ninds <- length(inds)

#affection is missing in new imputed fam file, so use the old one instead
#response <- read.table('../all.GSK.best_guessed_remove.qc2.fam', stringsAsFactors = F)
#response <- response[inds,,drop=F]
#response <- response[complete.cases(response),,drop=F]
#response.affected <- response$affected
#inds <- rownames(response)
#ninds <- length(inds)


# Plot densities ----------------------------------------------------------
features <- read_csv('allfeatures_uniq.txt', col_names = F)$X1
nfeatures <- length(features)

deepsea.hdf <-"/localscratch/gokcen/GSK/deepsea_features.hdf5"
deepsea.snps <- read_csv(paste0(deepsea.hdf, '.snps'))$snp
nsnps <- length(deepsea.snps)

stopifnot(all(funsig %in% deepsea.snps))

#sort plink genotypes according to deepsea snp order
snp.matrix <- pl$genotypes[inds,funsig]


# Amazing closure to cache features ---------------------------------------
# get.feature.func <- function() {
#   feat.list <- list()
#   function(feat) {
#     if (is.null(feat.list[[feat]])) {
#       H5close()
#       tmpfeat <- h5read(deepsea.hdf, feat)
#       colnames(tmpfeat) <- deepsea.snps
#       feat.list[[feat]] <<- tmpfeat[, funsig]
#     }
#     feat.list[[feat]]
#   }
# }
# get.feature <- get.feature.func()
# gc()
#
#
# get.individual.data <- function(ind, feat) {
#   f <- get.feature(feat)
#   ind.profile <- cbind(as(snp.matrix[ind,], 'numeric')[1,]+1, seq_len(ncol(snp.matrix)))
#   f[ind.profile]
# }
#
# get.individual.data.multi <- function(inds, feat) {
#   f <- get.feature(feat)
#   ind.genotypes <- t(as(snp.matrix[inds,], 'numeric')+1)
#   dim(ind.genotypes) <- prod(dim(ind.genotypes)) #vectorize
#   ind.profile <- cbind(ind.genotypes, rep(seq_len(ncol(snp.matrix)), length(inds)))
#   f <- f[ind.profile]
#   dim(f) <- c(ncol(snp.matrix), length(inds))
#   t(f)
# }

get.all <- function(inds, feat) {
  H5close()
  f <- h5read(deepsea.hdf, feat)
  colnames(f) <- deepsea.snps
  f <- f[, funsig]
  ind.genotypes <- t(as(snp.matrix[inds,], 'numeric')+1)
  dim(ind.genotypes) <- prod(dim(ind.genotypes)) #vectorize
  ind.profile <- cbind(ind.genotypes, rep(seq_len(ncol(snp.matrix)), length(inds)))
  f <- f[ind.profile]
  dim(f) <- c(ncol(snp.matrix), length(inds))
  t(f)
}

#big.guy <- get.individual.data.multi(1:ninds, features[760])
#library(gplots)
#library(Rclusterpp)
#heatmap.2(big.guy, distfun=function(x)x, hclustfun=Rclusterpp.hclust, trace='none')

# Generate the big matrix -------------------------------------------------

big.fella <- filebacked.big.matrix(ninds, length(funsig)*length(features),
                                   backingpath = '/localscratch/gokcen/GSK/',
                                   backingfile = 'big_guy_fullimputation_1500ind.bigmat',
                                   descriptorfile = 'big_guy_fullimputation_1500ind.desc')

progress <- dplyr::progress_estimated(nfeatures)

for (feat.i in seq_along(features)) {
  progress$tick()$print()
  print(feat.i/nfeatures)

  feat <- features[feat.i]
  print(feat)
  mat <- get.all(seq_len(ninds), feat)

  col.start <- (feat.i-1)*ncol(mat) + 1
  col.end <- col.start + ncol(mat) - 1

  print(c(col.start, col.end))

  big.fella[, col.start:col.end] <- mat

  if (feat.i %% 20 == 0) gc()
}

print('Done!')

write.table(inds, '/localscratch/gokcen/GSK/big_guy_fullimputation_1500ind.individuals',
            row.names = F, col.names = F)

write.table(features, '/localscratch/gokcen/GSK/big_guy_fullimputation_1500ind.features',
            row.names = F, col.names = F)

write.table(funsig, '/localscratch/gokcen/GSK/big_guy_fullimputation_1500ind.snps',
            row.names = F, col.names = F)

# Generate binned final data matrix ---------------------------------------

# library(KernSmooth)
# nbins <- 50
# ind.splits <- split(1:ninds, ceiling(seq_len(ninds)/20))
#
# final.matrix <- do.call(rbind, lapply(ind.splits, function(ind.index){
#   feat <- get.individual.data.multi(ind.index, features[160])
#   t(apply(feat, 1, function(x)bkde(na.omit(x), gridsize = nbins, range.x = c(0,1))$y))
# }))
#
# pheatmap(final.matrix, cluster_cols = F)

