
library(SNPlocs.Hsapiens.dbSNP144.GRCh37) #Bioconductor annotation package for dbSNP
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome)
library(Biostrings)

rsid2vcf <- function(snps, outfile) {
  message('Generating vcf file...')

  snp.db <- SNPlocs.Hsapiens.dbSNP144.GRCh37
  gen <- BSgenome.Hsapiens.UCSC.hg19

  not.rs <- substr(snps, 1, 2) != 'rs'

  if(any(not.rs)) {
    warning(sprintf('Some SNP ids do not start with rs. Excluding %d snps...', sum(not.rs)))
    snps <- snps[!not.rs]
  }

  message('Fetching GRCh37 positions of variants...')

  #most-memory-and-time-consuming-function-of-all-times
  snp.loc <- snpsById(snp.db, snps, ifnotfound='warning')

  #convert dbSNP style ('ch1') to UCSC style ('chr1')
  seqlevelsStyle(snp.loc) <- 'UCSC'

  #it seems there is a mismatch btw grch37 and hg19, hmm
  snp.loc <- dropSeqlevels(snp.loc, 'chrM')

  #workaround to make grch37 assembly work
  genome(snp.loc) <- 'hg19'

  #get alleles and split into char vector
  als <- IUPAC_CODE_MAP[snp.loc$alleles_as_ambig]
  als <- strsplit(als, '')
  names(als) <- snp.loc$RefSNP_id

  multi <- sapply(als, length)
  if(any(multi!=2))
  {
    warning(paste0('Multiallelic SNPs identified, taking the first alt. allele:',
                   snp.loc$RefSNP_id[multi!=2], '\n'))
  }

  #get ref allele from reference sequence
  ref <- as.character(as.matrix(getSeq(gen, snp.loc)))

  #alt allele is the one that's not ref. allele :)
  alt <- unname(mapply(function(alleles, ref)alleles[alleles!=ref][1], als, ref))

  if(is.matrix(alt))
    stop('Ref. allele could not be identified')

  #generate a df
  df <- data.frame(chr=as.character(seqnames(snp.loc)),
                   pos=start(snp.loc),
                   rsid=snp.loc$RefSNP_id,
                   ref=ref,
                   alt=alt, stringsAsFactors = F)

  if(!missing(outfile))
    write.table(df, outfile, quote = F, sep = '\t', row.names = F, col.names = F)

  df

}