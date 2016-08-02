library(bigmemory)
library(glmnet)
library(readr)
library(snpStats)
library(ggplot2)
library(parallel)
library(pheatmap)
library(gplots)
library(dplyr)

#significant snps
funsig <- read.table('gsk_17Msnps_significant_overlap.tsv', header = F, stringsAsFactors = F)$V1
nfunsig <- length(funsig)

top.eval <- readRDS('top_evalues_allpairs.Rds')
top.eval.overlap <- top.eval[top.eval$snp %in% funsig,]
top.eval.overlap.split <- split(top.eval.overlap$snp, top.eval.overlap$mineval.feature)

big.fella <- attach.big.matrix('/localscratch/gokcen/GSK/big_guy_fullimputation_1500ind.desc')

source('analysis_toolkit.R')

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
                    # c1=covar$C1,
                    # c2=covar$C2,
                    # c3=covar$C3,
                    row.names = covar$IID)
covar <- covar[inds,]
covar.age <- covar$age
covar.sex <- covar$sex
#covar.c1 <- covar$c1
#covar.c2 <- covar$c2
#covar.c3 <- covar$c3


# load lasso models -------------------------------------------------------
cv.lasso.fits <- readRDS('../../GSK_janine_imp07/results_fixedTF_1540ind/lasso_fits_17M_1500inds.Rds')

rm(pl)
gc()


# exclude missing values --------------------------------------------------

nona.idx <- !is.na(response.affected) &
  complete.cases(covar.age) &
  complete.cases(covar.sex)

nperm <- 1000
perm.index <- sapply(seq_len(nperm), function(x)sample(unname(which(nona.idx))))
stopifnot(ncol(perm.index) == nperm)

U <- mclapply(features.uniq, function(feature) {
  print(feature)
  big.fella.feat <- big.fella[,feat2idx(feature)]
  colnames(big.fella.feat) <- funsig
  big.fella.feat <- big.fella.feat[, top.eval.overlap.split[[feature]]]
  rownames(big.fella.feat) <- inds

  mat <- cbind(big.fella.feat[nona.idx,], covar.age[nona.idx], covar.sex[nona.idx])
  colnames(mat) <- c(colnames(big.fella.feat), 'age', 'sex')
  main.fit <- glmnet(mat,
                     response.affected[nona.idx],
                     alpha=1,
                     lambda = cv.lasso.fits[[feature]]$lambda.1se,
                     family='binomial')

  devi <- lapply(seq_len(nperm), function(ix){

    f <- cv.glmnet(mat,
                   response.affected[perm.index[,ix]],
                   alpha=1,
                   #nfolds = 100,
                   family='binomial')
    glmnet(mat,
           response.affected[perm.index[,ix]],
           alpha=1,
           lambda = f$lambda.1se,
           family='binomial')
  })

  c(list(main=main.fit), devi)

}, mc.cores = 12)

print('Done')

names(U) <- features.uniq
saveRDS(U, '../../GSK_janine_imp07/results_fixedTF_1540ind/permutation_models.Rds')
U <- readRDS('../../GSK_janine_imp07/results_fixedTF_1540ind/permutation_models.Rds')


lasso.fits.coef <- readRDS('../../GSK_janine_imp07/results_fixedTF_1540ind/lasso_fits_17M_1500inds_coef.Rds')
good.lasso.coef <- lasso.fits.coef[ sapply(lasso.fits.coef, length) > 0]
U <- U[names(good.lasso.coef)]

main.model.dev.ratio <- sapply(U, function(model)model[[1]]$dev.ratio)
model.order <- order(main.model.dev.ratio, decreasing = T)
main.model.dev.ratio.sorted <- main.model.dev.ratio[model.order]


D <- sapply(U, function(model)sapply(model[-1], function(model)model$dev.ratio))
D <- t(D)
D <- D[model.order,] #sort models based on the dev.ratio of the main model
#pheatmap(D, cluster_rows = F)

colMax <- function(data) apply(data, 2, max, na.rm = TRUE)

nD <- nrow(D)

#MaxT procedure
# U.final <- sapply(seq_len(nD), function(start){colMax(D[seq(start,nD),,drop=F])})
# U.final <- t(U.final)
# dim(U.final)
# p.values <- sapply(seq_len(nD), function(row)mean(U.final[row,] >= main.model.dev.ratio.sorted[row]))

p.values <- sapply(seq_len(nD), function(row)mean(D[row,] >= main.model.dev.ratio.sorted[row]))
names(p.values) <- names(main.model.dev.ratio.sorted)
p.values.adj <- p.adjust(p.values, 'fdr')
p.values.adj[p.values.adj<0.1]
unique(unlist(good.lasso.coef[names(p.values.adj[p.values.adj<0.1])]))


# let's also add rest of 919 models ---------------------------------------
failed.lasso.names <- setdiff(features.uniq, names(p.values.adj))
failed.lasso.pvals <- rep(NA_real_, length(failed.lasso.names))
names(failed.lasso.pvals) <- failed.lasso.names

#combine failed lasso models and observed pvalues
sigmodel.pvalues <- c(p.values.adj, failed.lasso.pvals)
sigmodel.pvalues.orig <- c(p.values, failed.lasso.pvals)
sigmodel.names <- features[match(names(sigmodel.pvalues), features.uniq)]
sigmodel.names.uniq <- features.uniq[match(names(sigmodel.pvalues), features.uniq)]
sigmodel.df <- as.data.frame(do.call(rbind, strsplit(sigmodel.names, '|', fixed = T)), stringsAsFactors = F)
colnames(sigmodel.df) <- c('Cell.line', 'Chromatin.feature', 'Treatment')
sigmodel.df$pval.adj <- sigmodel.pvalues
sigmodel.df$pval <- sigmodel.pvalues.orig
sigmodel.df$codename <- sigmodel.names.uniq
sigmodel.df$deepsea.name <- sigmodel.names

tidy.treatments <- gsub('None', '', gsub('\\.[1-9]$', '', sigmodel.df$Treatment))
sigmodel.df$Treatment <- tidy.treatments

sigmodel.df$devratio <- main.model.dev.ratio[match(sigmodel.df$codename, names(main.model.dev.ratio))]
sigmodel.df <- sigmodel.df[order(sigmodel.df$pval.adj),]
rownames(sigmodel.df) <- NULL
write.table(sigmodel.df, '../../GSK_janine_imp07/results_fixedTF_1540ind/permutation_pvalues.tsv', sep='\t', quote = F, row.names = F)



# now, plot the pvalue heatmap --------------------------------------------
tidy.treatments[tidy.treatments!=''] <- paste0('(', tidy.treatments[tidy.treatments!=''], ')')
sigmodel.df$Treatment <- tidy.treatments
sigmodel.df <- dplyr::select(sigmodel.df, Cell.line, Chromatin.feature, Treatment, pval.adj)
sigmodel.df = group_by(sigmodel.df, Cell.line, Chromatin.feature, Treatment) %>%
  mutate(pval.adj=min(pval.adj, na.rm=T)) %>% distinct %>% as.data.frame(., stringsAsFactors=F)
sigmodel.df$pval.adj[is.na(sigmodel.df$pval.adj)] <- -2


#fill `not measured` with NA so that we can color it differently
sigmodel.mat <- tidyr::spread(sigmodel.df, "Chromatin.feature", pval.adj)
rownames(sigmodel.mat) <- paste0(sigmodel.mat$Cell.line, sigmodel.mat$Treatment)
sigmodel.mat$Cell.line <- NULL
sigmodel.mat$Treatment <- NULL

f <- function(x){(!is.na(x)) & (x != -2) & (x < 0.1)}
sigmodel.mat <- sigmodel.mat[rowSums((!f(sigmodel.mat)), na.rm = T) != ncol(sigmodel.mat),]
sigmodel.mat <- sigmodel.mat[,colSums((!f(sigmodel.mat)), na.rm = T) != nrow(sigmodel.mat)]

#sort the rows based on the total number of significant models
ord <- order(rowSums(f(sigmodel.mat), na.rm = T), decreasing = T)

#the models that have no non-zero lasso coef. and the ones that failed the significance test are all
#represented by the same color
#set the pvalue of the lasso-dropped models (previously set to 1) somewhere between 0.1 and max so that it'll be painted in white
sigmodel.mat[sigmodel.mat == -2] <- mean(c(0.1, max(sigmodel.mat, na.rm = T)))

pdf('../../GSK_janine_imp07/Manuscript/Figures/snp_funcunit_pval_heatmap_landscape_ALLBLACK_stripped.pdf', 15, 11)
heatmap.2(((as.matrix(sigmodel.mat[ord,]))),
          Rowv = F, Colv = F, trace = 'none', dendrogram = 'none', na.color = 'black',
          col=c(colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'YlOrRd')[-1]),
                                 interpolate='spline')(20), 'white'),
          breaks=c(seq(0, 0.1, length.out = 21), max(sigmodel.mat, na.rm = T)),
          cexRow = 0.8,
          cexCol = 0.8,
          margins = c(7,13)
)
dev.off()


png('../../GSK_janine_imp07/Manuscript/Figures/snp_funcunit_pval_heatmap_landscape_ALLBLACK_stripped.png', 1500, 1100)
heatmap.2(((as.matrix(sigmodel.mat[ord,]))),
          Rowv = F, Colv = F, trace = 'none', dendrogram = 'none', na.color = 'black',
          col=c(colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'YlOrRd')[-1]),
                                 interpolate='spline')(20), 'white'),
          breaks=c(seq(0, 0.1, length.out = 21), max(sigmodel.mat, na.rm = T)),
          cexRow = 0.8,
          cexCol = 0.8,
          margins = c(7,13)
)
dev.off()
