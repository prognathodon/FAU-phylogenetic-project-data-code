library(phytools)
library(geiger)
tree <- read.nexus("finaltree2")
data <- read.csv("Summer 2025 datasets - Mosasaur.csv")
data <- data[data$Species != "", ]
data <- data[is.na(data$Jaw.Skull.length..mm.) == FALSE, ]
setdiff(tree$tip.label, data$Species)
tree <- drop.tip(tree, tip = setdiff(tree$tip.label, data$Species))
setdiff(data$Species, tree$tip.label)
jawlength <- as.vector(data$Jaw.Skull.length..mm.)
names(jawlength) <- data$Species
#traitgram
par(mar=c(7, 7, 7, 7))
phenogram(tree,jawlength, xlim = c(0, 50), ylim = c(0, 1700), fsize=0.6, spread.cost=c(0.2,1))
abline(a=1000, b=0, col = "red")
#heatmap
contMap(tree, jawlength)
#evolutionary modeling
bm <- fitContinuous(tree, jawlength, model = "BM") #Brownian motion
ou <- fitContinuous(tree, jawlength, model = "OU") #Single peak Ornstein-Uhlenbeck model
eb <- fitContinuous(tree, jawlength, model = "EB") # Early burst
mt <- fitContinuous(tree, jawlength, model = "trend") # Trend model
dr <- fitContinuous(tree, jawlength, model = "drift") # Drift model

# Repeated model testing to see if the parameter recovered from trend model is genuine, WARNING: extremely CPU costly
letsee <- function(a, b, d){
  dats <- vector(length=d)
  datai <- vector(length=d)
  for(i in 1:d){
    modi <- fitContinuous(a, b, model = "trend")
    dats[i] <- modi$opt$slope
    datai[i] <- modi$opt$aic
  }
  trenddata <- merge(dats, datai)
  return(trenddata)
}
trendtest <- letsee(tree, jawlength, 1000)
