library(bigmemory)
library(glmnet)
library(readr)
library(snpStats)
library(ggplot2)
library(parallel)

#significant snps
funsig <- read.table('gsk_17Msnps_significant_overlap.tsv', header = F, stringsAsFactors = F)$V1
nfunsig <- length(funsig)

top.eval <- readRDS('top_evalues_allpairs.Rds')
top.eval.overlap <- top.eval[top.eval$snp %in% funsig,]
top.eval.overlap.split <- split(top.eval.overlap$snp, top.eval.overlap$mineval.feature)

big.fella <- attach.big.matrix('/localscratch/gokcen/GSK/big_guy_fullimputation_1500ind.desc')

source('analysis_toolkit.R')

# read genotypes
pl <- read.plink('../../GSK_janine_imp07/data/all.GSK.best_guessed_07.qc_geno0')


features <- read_csv('allfeatures_uniq.txt', col_names = F)$X1
features.uniq <- make.names(features, unique = T)
nfeatures <- length(features)

# find individuals --------------------------------------------------------

#pca <- readRDS('../../GSK_janine_imp07/data/PCA_snprelate_17Msnps.Rds')
inds <- read.table('/localscratch/gokcen/GSK/big_guy_fullimputation_1500ind.individuals')$V1
ninds <- length(inds)

#qplot(pca$eigenvect[,1], pca$eigenvect[,2], color=as.factor(fam$V6)) + theme_minimal()
#pairs(pca$eigenvect[,1:5], col=pl$fam$affected+1)

#inds <- inds[pca$eigenvect[,1] > 0]
#ninds <- length(inds)

response <- read.table('../../GSK_janine_imp07/gwas/all.GSK.best_guessed_07.qc_geno0_deepSEA_prune.fam', header = F, stringsAsFactors = F)
response.affected <- response[match(inds, response$V1), 'V6']
response.affected[response.affected == -9] = NA
response.affected <- as.factor(response.affected)
names(response.affected) <- response[match(inds, response$V1), 'V1']

covar <- read.table('../../GSK_janine_imp07/gwas/GSK.covar', header = T, stringsAsFactors = F)
covar$AGE[covar$AGE == -9] <- NA
covar$K_SEX[covar$K_SEX == -9] <- NA
covar <- data.frame(age=covar$AGE,
                    sex=as.factor(covar$K_SEX),
                    row.names = covar$IID)
covar <- covar[inds,]
covar.age <- covar$age
covar.sex <- covar$sex

rm(pl)
gc()

fits <- mclapply(features.uniq, function(feature) {
  print(feature)
  big.fella.feat <- big.fella[,feat2idx(feature)]
  colnames(big.fella.feat) <- funsig
  big.fella.feat <- big.fella.feat[, top.eval.overlap.split[[feature]]]
  rownames(big.fella.feat) <- inds

  #big.fella.feat <- big.fella.feat[inds, ]

  nona.idx <- rowSums(is.na(big.fella.feat)) == 0 &
    !is.na(response.affected) &
    complete.cases(covar.age) &
    complete.cases(covar.sex)

  mat <- cbind(big.fella.feat[nona.idx,], covar.age[nona.idx], as.numeric(covar.sex[nona.idx])-1)

  colnames(mat) <- c(colnames(big.fella.feat), 'age', 'sex')
  fit <- cv.glmnet(mat,
                   response.affected[nona.idx],
                   nfolds = 100,
                   alpha=1,
                   family='binomial')
  fit
}, mc.cores = 12)

print('Done')

names(fits) <- features.uniq
saveRDS(fits, '../../GSK_janine_imp07/lasso_fits_17M_1500inds.Rds')


