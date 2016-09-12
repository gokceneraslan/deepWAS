#!/usr/bin/env Rscript

library(snpStats)
library(ggplot2)
library(glmnet)
library(parallel)

source('rundeepsea.R')

DEEPSEA.DIR <- '~/Deepsea'


# Run DeepSEA -------------------------------------------------------------

snps <- read.csv('data/snps.txt', stringsAsFactors=F, col.names='snps')$snps

# Run deepsea or use existing file
if (file.exists('data/evalues.tsv')) {
  snp.eval <- read.delim('data/evalues.tsv', stringsAsFactors = F)
} else {
  outdir <- rundeepsea(snps, DEEPSEA.DIR)
  snp.eval <- get.significant.snps(snps, outdir)
}

summary(snp.eval)
snp.eval.split <- split(snp.eval$snp, 
                        snp.eval$mineval.feature)

# Visualize number of SNPs we have for each deepse feature
qplot(sapply(snp.eval.split, length), bins=100, color=I('white')) + 
  labs(x='Model feature sizes') +
  theme_minimal()

# Save names of DeepSEA features
features.uniq <- unique(snp.eval$mineval.feature)


# Read plink or dosage ----------------------------------------------------

pl <- read.plink('data/gwas')
ind.ids <- rownames(pl$genotypes)
covar <- read.delim('data/gwas.covar', stringsAsFactors = F)
pl$fam <- merge(pl$fam, covar) #merge extra covariates
covar <- cbind(age=pl$fam$age, sex=pl$fam$sex-1)
response <- as.factor(pl$fam$affected - 1)

# Remove missing cases, shouldn't be a problem with dosage data
cc.index <- complete.cases(as(pl$genotypes, 'numeric'), covar, response)
pl$genotypes <- pl$genotypes[cc.index,]
covar <- covar[cc.index,]
response <- response[cc.index]
ind.ids <- ind.ids[cc.index]


# Regression --------------------------------------------------------------

if (file.exists('data/models.Rds')) {
  cv.fits <- readRDS('data/models.Rds')
} else {
  set.seed(42)
  cv.fits <- mclapply(features.uniq, function(feature) {
    print(feature)
    mat <- cbind(as(pl$genotypes[, snp.eval.split[[feature]]], 'numeric'),
                 covar)
    
    fit <- cv.glmnet(mat,
                     response,
                     nfolds = 100,
                     alpha=1,
                     family='binomial')
    fit
  }, mc.cores = 12)
  
  names(cv.fits) <- features.uniq
  saveRDS(cv.fits, 'data/models.Rds')
}


# Model selection ---------------------------------------------------------

# Let's fit models on shuffled response
nperm <- 1000
perm.index <- sapply(seq_len(nperm), function(x)sample(ind.ids))
stopifnot(ncol(perm.index) == nperm)

if (file.exists('data/permutation_models.Rds')) {
  
  U <- readRDS('data/permutation_models.Rds')
  
} else {
  U <- mclapply(features.uniq, function(feature) {
    
    print(feature)
    mat <- cbind(mat,
                 covar)
    
    fit <- glmnet(mat,
                  response,
                  alpha=1,
                  lambda = cv.fits[[feature]]$lambda.1se,
                  family='binomial')
    
    devi <- lapply(seq_len(nperm), function(ix){
      
      # Use shuffled response here
      f <- cv.glmnet(mat,
                     response[perm.index[,ix]],
                     alpha=1,
                     #nfolds = 100, # too slow
                     family='binomial')
      glmnet(mat,
             response[perm.index[,ix]],
             alpha=1,
             lambda = f$lambda.1se,
             family='binomial')
    })
    
    list(fit=fit, devi=devi)
  }, mc.cores = 12)

  names(U) <- features.uniq
  saveRDS(U, 'data/permutation_models.Rds')
}


main.model.dev.ratio <- sapply(U, function(model)model$fit$dev.ratio)
model.order <- order(main.model.dev.ratio, decreasing = T)
main.model.dev.ratio.sorted <- main.model.dev.ratio[model.order]


D <- sapply(U, function(model)sapply(model[-1], function(model)model$dev.ratio))
D <- t(D)
D <- D[model.order,] #sort models based on the dev.ratio of the main model
#pheatmap(D, cluster_rows = F)

