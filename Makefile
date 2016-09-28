all:
	Rscript -e 'library(knitr);library(markdown);knit("deepwas.Rmd","deepwas.md");markdownToHTML("deepwas.md","deepwas.html")'

