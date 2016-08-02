library(readr)
library(dplyr)

f = commandArgs(trailingOnly=TRUE)

for (file in f) {
  print(file)

  eval <- read_csv(file)
  mins = which(eval[,-(1:6)] < 5e-5, arr.ind=T)
  mins.val = as.matrix(eval[,-(1:6)])[mins]

  df <- data.frame(snp=eval$name[mins[,1]],
                    mineval.feature=make.names(colnames(eval), unique=T)[mins[,2]+6],
                    mineval=mins.val,
                    stringsAsFactors=F)

  write_tsv(df, 'evalues.tsv', append=T)
}
