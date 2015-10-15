# Assignment 3
# Oct 15, 2015

# For this assignment you will need to install ALDEx2 from Bioconductor
# http://bioconductor.org/packages/release/bioc/html/ALDEx2.html

# there is some documentation on how to run this at:
# https://www.bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.pdf

# lets play with the tongue/saliva dataset
d <- read.table("tongue_saliva.txt", header=T, row.names=1, sep="\t")

# get the taxon list and remove that column from the table
tax <- d$tax.0
d$tax.0 <- NULL

# we will subset it a bit so we only have a few samples to keep the time required within bounds

# take a random sample of 20 each
d.all <- data.frame(d[,c("sa_014403","sa_034170","sa_024546","sa_023380","sa_098739","sa_024190",
 "sa_032226","sa_037796","sa_035377","sa_024339","sa_102523","sa_102663",
 "sa_095217","sa_033676","sa_015068","sa_016696","sa_106985","sa_016972",
 "sa_111578","sa_014965","td_021282","td_111751","td_016782","td_023978",
 "td_035163","td_096596","td_038959","td_016626","td_024341","td_111445",
 "td_024192","td_035791","td_014971","td_024909","td_014515","td_098872",
 "td_106496","td_038744","td_106341","td_024398")])

# remove all rows that have 0 counts total
data <- d.all[which(apply(d.all, 1, sum) > 0),]
tax.data <- tax[which(apply(d.all, 1, sum) > 0)]

library(ALDEx2)
####
# Question 1
# Comment the following lines to demonstrate that you have an understanding of the
# mechanics of the below commands
####

conds <- c(rep("S", 20), rep("T", 20))

# generate the random instances of the data and transform them
# we will use 32 to start
x <- aldex.clr(data,mc.samples=32)
x.e <- aldex.effect(x, conds, include.sample.summary=TRUE)

# this could take anywhere from 30 seconds to 2-3 minutes depending on your computer
x.t <- aldex.ttest(x, conds)

x.all <- data.frame(x.e, x.t)

# positive numbers mean more abundant in T, negative numbers more abundant in S
aldex.plot(x.all, test="wilcox")

###### End commenting here
# this shows you how to subset to recover which taxa are below the significance threshold
x.all[x.all$wi.eBH < 0.1,]
tax.data[x.all$wi.eBH < 0.1]

####
# Question 2
# Explain the results of the ALDEx2 effect size plot
# what is the within/between group variation in this dataset for the bulk of the taxa?
# what taxa are different between samples?
####

####
# Question 3
# Explain how you could get a more stable estimate of the expected values
# returned by ALDEx2 and why that may or may not be desirable
####
# this will make a biplot from the ALDEx expected abundance values
# requires include.sample.summary=TRUE
expected <- x.all[,grep("rab.sample", colnames(x.all))]
colnames(expected) <- gsub("rab.sample.","", colnames(expected))
pcx <- prcomp(t(expected))
biplot(pcx, cex=c(0.5,0.4) )

####
# Question 4
# Generate a biplot from the ALDEx expected values, and using the method you
# used for the last assignment
#
# Given the results of the ALDEx2 output and the two biplots what can you now
# conclude about the distinguishability of these two datasets?
####
