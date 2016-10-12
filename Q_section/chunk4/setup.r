# setup for tiger,ladybug,alien example

Tiger <- round(runif(100, 1800,2200))
Ladybug <- round(runif(100, 8000,12000))
Alien <- round(runif(100, 450,550))

d <- data.frame(Tiger,Ladybug,Alien)

myplot <- function(x, y, lab=NULL,
                   mycol=rgb(1,0,0,.25)){
  # x contains the x axis to be plotted
  # y contains the y axis to be plotted
  plot(x, y, pch=19, cex=0.5,
       col=mycol,xlab=lab[1], ylab=lab[2])
}

myplot(d$Ladybug, d$Alien, lab=c("LBug", "Alien"))

