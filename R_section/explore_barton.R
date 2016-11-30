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

source("script/load_biplot.r")

```

```{r agg, message=FALSE, warning=FALSE, echo=TRUE, fig.cap='Initial exploratory PCA plot of the SNF2 KO dataset. .', fig.height=4, fig.width=14}

nms.agg <- paste(m[,2], m[,3], sep=":")

d.agg <- aggregate(t(d), by=list(nms.agg), FUN=sum)
rownames(d.agg) <- d.agg$Group.1
d.agg$Group.1 <- NULL

# remove rows with 0 counts
# returns data with samples by columns!!!!
d.agg.gt0 <- codaSeq.filter(d.agg, min.reads=0, min.prop=0, min.count=0, samples.by.row=TRUE)

# estimate 0 values (zCompositions)
d.agg.n0 <- cmultRepl(t(d.agg.gt0), method="CZM", label=0)

# clr transform
d.agg.n0.clr <- codaSeq.clr(d.agg.n0)

pcx.agg  <- prcomp(d.agg.n0.clr)
mvar.agg.clr <- sum(pcx.agg$sdev^2)

par(mfrow=c(1,1))
biplot(pcx.agg, var.axes=F, scale=0, cex=c(2,.5))

plot(pcx.agg$x[,1], pcx.agg$x[,2], pch=NA,
    xlim=c(min(pcx.agg$x[,1] )-5, max(pcx.agg$x[,1] + 5)),
    xlab=paste("PC1: ", round(sum(pcx.agg$sdev[1]^2)/mvar.agg.clr, 3)),
    ylab=paste("PC2: ", round(sum(pcx.agg$sdev[2]^2)/mvar.agg.clr, 3)),
    main="PCA score plot")
text(pcx.agg$x[,1], pcx.agg$x[,2], labels=rownames(pcx.agg$x), cex=3)

```

```{r outlier, message=FALSE, warning=FALSE, echo=TRUE, fig.cap='outliers', fig.height=4, fig.width=14}

# get WT indices
WT <- grep("WT", rownames(d.agg))
# subset
WT.agg <- d.agg[WT,]

wt.gt0 <- codaSeq.filter(WT.agg, min.reads=0, min.prop=0, min.count=0, samples.by.row=TRUE)

# estimate 0 values (zCompositions)
wt.agg.n0 <- cmultRepl(t(wt.gt0), method="CZM", label=0)

# clr transform
wt.agg.n0.clr <- codaSeq.clr(wt.agg.n0)

pcx.wt  <- prcomp(wt.agg.n0.clr)
mvar.wt.clr <- sum(pcx.wt$sdev^2)

par(mfrow=c(1,1))
biplot(pcx.wt, var.axes=F, scale=0, cex=c(2,.5))

WT.g <- codaSeq.outlier(wt.agg.n0.clr, plot.me=TRUE)

# now do the same for SNF2
# get SNF indices
SNF <- grep("SNF", rownames(d.agg))
# subset
SNF.agg <- d.agg[SNF,]

SNF.gt0 <- codaSeq.filter(SNF.agg, min.reads=0, min.prop=0, min.count=0, samples.by.row=TRUE)

# estimate 0 values (zCompositions)
SNF.agg.n0 <- cmultRepl(t(SNF.gt0), method="CZM", label=0)

# clr transform
SNF.agg.n0.clr <- codaSeq.clr(SNF.agg.n0)

pcx.SNF  <- prcomp(SNF.agg.n0.clr)
mvar.SNF.clr <- sum(pcx.SNF$sdev^2)

par(mfrow=c(1,1))
biplot(pcx.SNF, var.axes=F, scale=0, cex=c(2,.5))

SNF.g <- codaSeq.outlier(SNF.agg.n0.clr, plot.me=TRUE)


```

```{r good_data_pca, message=FALSE, warning=FALSE, echo=TRUE, fig.cap='outliers', fig.height=4, fig.width=14}


# check names and see if they have . or :
head(SNF.g$good)
# if they are separated by ., then do the following
WT.g$good <- gsub("WT.", "WT:", WT.g$good)
SNF.g$good <- gsub("SNF2.", "SNF2:", SNF.g$good)

d.good <- rbind(d.agg[SNF.g$good,],d.agg[WT.g$good,])

d.good.gt0 <- codaSeq.filter(d.good, min.reads=0, min.prop=0, min.count=0, samples.by.row=TRUE)

# estimate 0 values (zCompositions)
d.good.agg.n0 <- cmultRepl(t(d.good.gt0), method="CZM", label=0)

# clr transform
d.good.agg.n0.clr <- codaSeq.clr(d.good.agg.n0)

pcx.good  <- prcomp(d.good.agg.n0.clr)
mvar.good <- sum(pcx.good$sdev^2)

par(mfrow=c(1,1))
biplot(pcx.good, var.axes=F, scale=0, cex=c(2,.5))

write.table(d.good.gt0, file="data/filtered_table.txt", sep="\t", quote=F, col.names=NA)
```

```{r aldex}

library(ALDEx2)

# make a vector of conditions

# you guys to make the distribution of possible values
x <- aldex.clr(d.good.gt0)

conds <- c(rep("SNF", length(SNF.g$good)), rep("WT", length(WT.g$good)))
# me
x <- aldex.clr(d.good.gt0, conds, mc.samples=16)

x.e <- aldex.effect(x, conds)
x.t <- aldex.ttest(x, conds)
x.all <- data.frame(x.e,x.t)


# now explore1 BLand Altman Plot
aldex.plot(x.all, type="MA")

# volcano plot
plot(x.all$diff.btw, x.all$wi.eBH, log="y"

# effect size plot
aldex.plot(x.all, type="MW")
# we use effect sizes of 2 or more


# subset to find those with large effect sizes
rownames(x.all)[abs(x.all$effect) > 4]


```
```{r edgeR}


# from ALDEx above
group <- factor(conds)

y <- DGEList(counts=d.good.gt0,group=group)

y <- calcNormFactors(y)

# p
design <- model.matrix(~group)


y <- estimateDisp(y)

# exact test using negative binomial
et <- exactTest(y)
```
