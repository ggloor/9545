---
title: "explore outliers"
author: "gg"
date: "Nov 28, 2016"
output: pdf_document
---

```{r lib, include=FALSE, message=FALSE, warning=FALSE, echo=TRUE}

library(zCompositions)
library(ALDEx2)
library(CoDaSeq)

plot=FALSE
```

```{r explore1, message=FALSE, warning=FALSE, echo=TRUE, fig.cap='Initial exploratory PCA plot of the SNF1 KO dataset. .', fig.height=4, fig.width=14}

d <- read.table("data/countfinal2.tsv", header=T, row.names=1, sep="\t", stringsAsFactors=F, comment.char="")

m <- read.table("data/ERP004763_sample_mapping.tsv", , header=T, row.names=1, sep="\t", stringsAsFactors=F, comment.char="")

library(CoDaSeq)
library(zCompositions)

min(apply(d,1,sum))
# at least one gene has 0 in every sample

# need to remove genes with 0 reads in all samples
d.n0 <- codaSeq.filter(d, min.reads=0, min.prop=0, min.count=0, samples.by.row=FALSE)

# replace 0 values
d.czm <- cmultRepl(t(d.n0), method="CZM", label=0)

# convert to clr
d.clr <- codaSeq.clr(d.czm)

# do the SVD
d.pcx <- prcomp(d.clr)

# make the biplot
biplot(d.pcx, var.axes=F, scale=0, cex=c(0.5,2))

```

