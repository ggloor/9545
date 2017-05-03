#! /opt/local/bin/Rscript --vanilla --default-packages=base,stats,utils
library(knitr)
library(rmarkdown)
file <- list.files(pattern='Correlation_compositions.Rmd')
rmarkdown::render(file)
