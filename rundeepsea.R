
source('rsid2vcf.R')

library(tidyr)
library(dplyr)
library(readr)

rundeepsea <- function(rsid, deepsea.dir, temp.dir, read.all=F) {

  #TODO: implement paging i.e. run deepsea with at most 1M snps
  if(missing(temp.dir))
    temp.dir <- tempdir()

  vcffile <- tempfile(fileext = '.vcf', tmpdir=temp.dir)

  # get position and allele info and save to vcf file
  rsid2vcf(rsid, vcffile)

  setwd(deepsea.dir)

  message('Running DeepSEA...')

  #run deepsea
  ret <- system2('python', args=c('rundeepsea.py', vcffile, temp.dir))
  stopifnot(ret==0)

  if (read.all) {
    #filenames
    basefile <- 'infile.vcf.out.'
    filetypes <- c('alt', 'diff', 'evalue', 'logfoldchange',
                   'ref', 'snpclass', 'funsig')

    filenames <- paste0(basefile, filetypes)
    names(filenames) <- filetypes

    results <- lapply(filetypes, function(ft){
                      f <- filenames[ft]
                      f <- read.csv(paste0(outdir, '/', f), stringsAsFactors = F, check.names = F)
                      f
                   })
    names(results) <- filetypes
    return(results)
  } else {
    return(temp.dir)
  }
}

get.significant.snps <- function(outdir, cutoff=5e-5) {
  message('Pre-selecting variants based on deepsea predictions...')

  file <- paste0(outdir, '/infile.vcf.out.evalue')
  eval <- read_csv(file)
  mins <- which(eval[,-(1:6)] < cutoff, arr.ind=T)
  mins.val <- as.matrix(eval[,-(1:6)])[mins]

  df <- data.frame(snp=eval$name[mins[,1]],
                   mineval.feature=make.names(colnames(eval), unique = T)[mins[,2]+6],
                   mineval=mins.val,
                   stringsAsFactors=F)
  df
}

preprocess.results <- function(result.df, snp, col.subset) {

  #process ref and alt files
  allele.df.names <- c('ref', 'alt')

  #check if column names are same
  stopifnot(all(head(colnames(result.df$alt)) == head(colnames(result.df$ref))))

  #deduplicate colnames
  colnames(result.df$ref)[-(1:6)] <- colnames(result.df$ref[,-(1:6)])
  colnames(result.df$alt)[-(1:6)] <- colnames(result.df$alt[,-(1:6)])

  if(!missing(snp)) {
    result.df$ref <- result.df$ref[result.df$ref$name %in% snp,]
    result.df$alt <- result.df$alt[result.df$alt$name %in% snp,]
  }

  if(!missing(col.subset)) {
    result.df$ref <- result.df$ref[,c(1:6, which(colnames(result.df$ref) %in% col.subset))]
    result.df$alt <- result.df$alt[,c(1:6, which(colnames(result.df$alt) %in% col.subset))]
  }

  colnames(result.df$alt)[1] <- 'X'
  colnames(result.df$ref)[1] <- 'X'

  nref <- nrow(result.df$ref)
  nalt <- nrow(result.df$alt)

  result.df$ref <- tbl_df(result.df$ref) %>% select(-alt) %>% dplyr::rename(allele=ref)
  result.df$alt <- tbl_df(result.df$alt) %>% select(-ref) %>% dplyr::rename(allele=alt)

  #rbind all
  ds <- bind_rows(result.df$alt, result.df$ref) %>%
    mutate(ref.alt=c(rep('alt', nalt), rep('ref', nref))) %>%
    select(1:5, ref.alt, everything()) #move ref.alt to the beginning

  #melt prob & chr features
  ds <- gather(ds, chr.feature, prob, -(1:6), convert=T)

  #split cell type, chr feat and treatment into different columns
  ds <-  mutate(ds, cell.type=sapply(strsplit(ds$chr.feature, '|', fixed=T), `[`, 1),
               chr.feat=sapply(strsplit(ds$chr.feature, '|', fixed=T), `[`, 2),
               treatment=sapply(strsplit(ds$chr.feature, '|', fixed=T), `[`, 3)) %>%
    select(-chr.feature) #drop old col.
  ds

  #TODO: preprocess following thigs as well
  #   snp.df.names <- c('diff', 'evalue', 'logfoldchange')
  #   snp.class.name <- c('snpclass', 'funsig')
}