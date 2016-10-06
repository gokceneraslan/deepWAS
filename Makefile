all:
	Rscript -e 'library(rmarkdown);render("deepwas.Rmd", html_document(toc=T,fig_caption=F,code_folding="show",code_download=T,toc_float=list(smooth_scroll=F,collapsed=F)))'
	mv deepwas.html index.html

