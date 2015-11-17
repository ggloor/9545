##############
# assignment 2
##############
#
# THE PURPOSE OF THIS ASSIGNMENT IS TO FAMILIARIZE YOU WITH EXPLORATORY ANALYSIS OF YOUR DATA
#
# This assignment is due Wed Oct 7 before 9:30 am
#
# I have verified that all the example code works as expected
#
# below is example code that you can use to examine the HMP microbiota dataset
# saliva and tongue samples. I have munged it a bit to remove everything that is
# vanishingly rare. You will need to get the file: tongue_saliva.txt

# question 1
# run the code verbatim
# does it appear that the tongue and saliva microbiomes are distinctive? justify

# question 2
# modify the code to reduce the dataset by occurrence as done in class using the
# example code below. the cutoff of 0.8 is arbitrary and not necessarily suggested
# does this exploration change your assessment from question 1? justify

# question 3
# examine the tongue, and saliva samples independently
# does this exploration change your assessment?
# is there a taxon or group of taxon in common that accounts for the observed split?
# give a few brief questions that you would like to ask the people who did the wet-lab work.

# question 4
# identify the taxon, or group of taxa, that most differentiate tongue and saliva
# identify the taxon, or group of taxa, that are irrelevant to the differentiation of these two
# why did you choose these?


############### EXAMPLE CODE
# required libraries for compositional analysis
library(compositions)
library(zCompositions)

# read the entire data table
# tongue samples start with td_
# saliva samples start with sa_
d <- read.table("tongue_saliva.txt", header=T, row.names=1, sep="\t")

# get the taxon list and remove that column from the table
tax <- d$tax.0
d$tax.0 <- NULL

# samples must be by rows so use t()
d.n0 <- cmultRepl(t(d), label=0, method="CZM")

# the apply function rotates the data when by row, so turn it back to samples by row
# how I hate R
d.n0.clr <- t(apply(d.n0, 1, function(x){log2(x) - mean(log2(x))}))

####
# this set of code replaces row and column names with single values
# and the genus name of each taxon, if known
# replace rownames in d.n0.clr with one character values
# gsub is substitute the every occurence of the pattern with the replacement
# pattern is read as: "sa_ followed by one or more other characters (.+)"
# See Regular_expression on wikipedia if you are keen, this is POSIX grep
rownames(d.n0.clr) <- gsub("sa_.+", "S", rownames(d.n0.clr))
rownames(d.n0.clr) <- gsub("td_.+", "T", rownames(d.n0.clr))

# make sure you are using the correct taxon list if you reduce the dataset by cutoff
# as it is now, it is using the entire dataset
# see class code for example
colnames(d.n0.clr) <- gsub(".+__", "", tax)
palette=palette(c("darkcyan","coral3"))
conds <- data.frame(c(rep(1,length(grep("T", rownames(d.n0.clr)))), rep(2, length(grep("S", rownames(d.n0.clr))))))
colnames(conds) <- "cond"
####

# generate the PCA object
pcx <- prcomp(d.n0.clr)

# plot it
# you control the color of the taxa with the rgb(red,green,blue,transparency) function
# you control the size of the sample and tax labels with the cex() function
coloredBiplot(pcx,col=rgb(0,0,0,0.3),cex=c(0.8,0.4), xlabs.col=conds$cond, var.axes=F, arrow.len=0.05)
# this gives the proportion of variance explained in the first five components
(pcx$sdev^2/mvar(d.n0.clr))[1:5]

#### STOP HERE FOR QUESTION 1, USE THIS FOR OTHER QUESTIONS AS NEEDED
##### methods that you can use for exploration
### method to keep only taxa present in x proportion or more samples the dataset
# the which function returns a list of row numbers where the number of 0 values in the row is greater than the cutoff fraction
cutoff = .8
d.subset <- data.frame(d[which(apply(d, 1, function(x){length(which(x != 0))/length(x)}) > cutoff),])
tax.subset <- tax[which(apply(d, 1, function(x){length(which(x != 0))/length(x)}) > cutoff)]

#### EXAMPLE to make a tongue or saliva only dataset
d.t <- d[,grep("td_", colnames(d))]
d.s <- d[,grep("sa_", colnames(d))]
