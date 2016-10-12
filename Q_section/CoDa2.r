Tiger <- round(runif(100, 1800,2200))
Ladybug <- round(runif(100, 8000,12000))
Alien <- round(runif(100, 450,550))

d <- data.frame(cbind(Tiger, Ladybug, Alien))

d.rare <- as.data.frame(sapply(d, function(x) rdirichlet(1,x)))

plot(d$Tiger, d$Ladybug, main=round(cor(d$Tiger, d$Ladybug), 2), cex.lab=1.8,
    cex.main=2, pch=19,
    col=rgb(1,0,0,0.5), xlab="Tiger", ylab="Ladybug")
plot(d$Tiger, d$Alien, main=round(cor(d$Tiger, d$Alien), 2), cex.lab=1.8,
    cex.main=2, pch=19,
    col=rgb(1,0,0,0.5), xlab="Tiger", ylab="Alien")
plot(d$Ladybug, d$Alien, main=round(cor(d$Ladybug, d$Alien), 2), cex.lab=1.8,
    cex.main=2, pch=19, col=rgb(1,0,0,0.5), xlab="Ladybug", ylab="Alien")
scatterplot3d(d, type="h", lty.hplot=3, angle=40, pch=19, cex.lab=1.5)

par(mfrow=c(2,3),mar=c(5,5,4,1))
plot(d.rare$Tiger, d.rare$Ladybug, main=round(cor(d.rare$Tiger, d.rare$Ladybug), 2),
    cex.main=2, pch=19, col=rgb(1,0,0,0.5), cex.lab=1.8, xlab="Tiger", ylab="Ladybug")
plot(d.rare$Tiger, d.rare$Alien, main=round(cor(d.rare$Tiger, d.rare$Alien), 2),
    cex.main=2, pch=19, col=rgb(1,0,0,0.5), cex.lab=1.8, xlab="Tiger", ylab="Alien")
plot(d.rare$Ladybug, d.rare$Alien, main=round(cor(d.rare$Ladybug, d.rare$Alien), 2),
    cex.main=2, pch=19, col=rgb(1,0,0,0.5), cex.lab=1.8, xlab="Ladybug", ylab="Alien")
scatterplot3d(d.rare, type="h", lty.hplot=3, cex.lab=1.5, angle=40, pch=19)
plot(acomp(d.rare), scale=T, center=T, pch=19, cex=1, axes=FALSE)

par(mfrow=c(1,3))
plot(d$Tiger, d.rare$Tiger)
plot(d$Ladybug, d.rare$Ladybug)
plot(d$Alien, d.rare$Alien)

aitchison.mean <- function( n, log=FALSE ) {

    # Input is a vector of non-negative integer counts.
    # Output is a probability vector of expected frequencies.
    # If log-frequencies are requested, the uninformative subspace is removed.

    n <- round( as.vector( n, mode="numeric" ) )
    if ( any( n < 0 ) ) stop("counts cannot be negative")

    a <- n + 0.5
    sa <- sum(a)

    log.p <- digamma(a) - digamma(sa)
    log.p <- log.p - mean(log.p)

    if ( log ) return(log.p)

    p <- exp( log.p - max(log.p) )
    p <- p / sum(p)
    return(p)
}

 d.aitch <- as.data.frame(sapply(d,aitchison.mean))
scatterplot3d(d.aitch, type="h", lty.hplot=3, cex.lab=1.5, angle=40, pch=19)

par(mfrow=c(3,3), mar=c(4,4,5,1)+0.1)
plot(d$Tiger, d.rare$Tiger)
plot(d$Ladybug, d.rare$Ladybug)
plot(d$Alien, d.rare$Alien)

plot(d$Tiger, d.aitch$Tiger)
plot(d$Ladybug, d.aitch$Ladybug)
plot(d$Alien, d.aitch$Alien)

plot(d.rare$Tiger, d.aitch$Tiger)
abline(0,1)
plot(d.rare$Ladybug, d.aitch$Ladybug)
abline(0,1)
plot(d.rare$Alien, d.aitch$Alien)
abline(0,1)
