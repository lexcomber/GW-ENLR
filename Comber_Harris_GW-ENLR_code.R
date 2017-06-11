library(GISTools)
library(rgdal)
library(GWmodel)
library(repmis)
library(knitr)
library(spdep)
library(GISTools)
library(rgdal)
library(GWmodel)
library(repmis)
library(knitr)
library(spdep)
library(spgwr)
library(glmnet)
library(tidyverse)
require(XML)
require(RCurl)
require(XLConnect)
library(maps)
library(ecospat)
library(scales)

# ELN see http://www.onthelambda.com/2015/08/19/kickin-it-with-elastic-net-regression/
# http://r-sig-geo.2731867.n2.nabble.com/Hat-matrix-in-ggwr-models-with-adaptative-kernel-td6812565.html
# ELN http://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html

### Part 1. Data setup
## 1.1 Voting data: 
## 1.1.1 load County data
#setwd("~/Desktop/my_docs_mac/leeds_work/research/COSIT/anal/")
#county2 <- readShapePoly("cb_2015_us_county_5m/cb_2015_us_county_5m.shp")
proj4string(county2) <- CRS("+proj=longlat ")
map.county <- map("county", plot = FALSE, fill = TRUE, res = 0)
## convert 'map' to 'SpatialPolygons'
county <- map2SpatialPolygons(map.county, 
               IDs = map.county$names, proj4string = CRS("+proj=longlat "))
county <- county2[county,]
# do fips
sp.func <- function(x) {
	return(sprintf("%s%s", x[,1], x[,2]))
}
x.i <- county@data[,1:2]
fips <- sp.func(x.i)
county <- SpatialPolygonsDataFrame(county, data = data.frame(county, fips = fips))

## 1.1.2 Census Data: now get some census data to counties
switchlex = "DONE"
if (switchlex != "DONE" ){
	data.mat <- matrix(nrow = 0, ncol = 5)
	for (i in 1:nrow(county)) {
		county.i <- as.vector(county$fips[i]) 
		tit <- paste("https://www.census.gov/quickfacts/download.php?fips=", county.i, ",00", sep = "")
		tmp <- download.file(tit, "tmp.xls")
		wb = loadWorkbook("tmp.xls")
		df = readWorksheet(wb, sheet = "Worksheet", header = TRUE)
		# employed, coled, over65,pop density, white
		data.i <- as.numeric(c(df[59,2], df[54,2], df[14,2], df[87,2], df[19,2]))
		data.mat <- rbind(data.mat, data.i)
		if (i %/% 100 != 0) cat(i, "\t")
		# data.it <- c(df[59,], df[54,], df[14,], df[87,], df[19,])
		
	}
	# Cath the county that had notes - data as text
	myurl <- "https://www.census.gov/quickfacts/download.php?fips=51019,00"
	tmp <- download.file(myurl, "tmp.xls")
	wb = loadWorkbook("tmp.xls")
	df = readWorksheet(wb, sheet = "Worksheet", header = TRUE)
	data.mat[812,1] <- 61
	data.mat[812,2] <- 26.7
	data.mat[812,3] <- 19.9
	data.mat[812,5] <- 89.8
	setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
	save.image(file = "new.pap_v1.RData")
	head(data.mat)
	rownames(data.mat) <- county$fips
	colnames(data.mat) <- c("PCemp", "PCcol", "PCo65", "PopD", "PCwhi")
	save(data.mat, file = "data.mat.RData")

} else {
	setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
	load("data.mat.RData")	
}

## 1.1.3 Voting data: open csv data, determine PA and attach all
# as here: https://simonrogers.net/2016/11/16/us-election-2016-how-to-download-county-level-results-data/
# http://www.personal.psu.edu/users/a/c/acr181/2004_Election.zip
source_data("https://github.com/lexcomber/GW-ENLR/blob/master/US_County_Level_Presidential_Results_08-16.RData?raw=True")
head(voting)
head(county@data)
index <- match(as.numeric(as.character(county$GEOID)), 
	as.numeric(as.character(voting$FIPS2)))
summary(index)
which(is.na(index))
tmp <- SpatialPolygonsDataFrame(county, 
	data = data.frame(county, voting[index,11:14]))
# county fix - see https://en.wikipedia.org/wiki/Oglala_Lakota_County,_South_Dakota
tmp@data[which(is.na(index)),11:14] <- voting[voting$FIPS2 == 46113,11:14]
# determine Trump PA
tmp$Trump = apply(tmp@data[,12:14], 1, which.max) - 1
county <- tmp
county <- SpatialPolygonsDataFrame(county, 
	data = data.frame(county@data, data.mat)) 

## 1.2 Presence Absence Species data : 
## 1.2.1 load Presence Absence data
data(ecospat.testData)
plot(ecospat.testNiche.nat[,1:2], asp = 1, cex = 0.1)
points(ecospat.testNiche.nat[ecospat.testNiche.nat$species_occ == 1,1:2], asp = 1, cex = 0.1, col = "red")
# aetpet: Ratio of actual to potential evapotranspiration.
# gdd: Growing degree-days above 5 degrees C.
# p: Annual amount of precipitations.
# pet: Potential evapotranspiration.
# stdp: Annual variation of precipitations.
# tmax: Maximum temperature of the warmest month.
# tmin: Minimum temperature of the coldest month.
# tmp: Annual mean temperature.

## 1.2.3 Create spdf and clip 
usa <- gUnaryUnion(county)
coords = ecospat.testNiche.nat[, 1:2]
nat.sp <- SpatialPointsDataFrame(coords, 
  data = data.frame(ecospat.testNiche.nat), 
  proj4string = CRS("+proj=longlat "))
nat.sp <- nat.sp[usa,]
plot(nat.sp)
plot(nat.sp[nat.sp$species_occ == 1,], cex = 0.3, col = "red", add = T)
plot(usa, lwd = 2, add = T)

## 1.3 Final Transform for distance etc analyses
new.proj <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96
                +x_0=0.0 +y_0=0.0 +ellps=GRS80 +units=m +datum=NAD83 +no_defs
                +towgs84=0,0,0")
county <- spTransform(county, new.proj)
usa <- spTransform(usa, new.proj)
nat.sp <- spTransform(nat.sp, new.proj)

## 1.4 Some initial maps
## 1.4.1 Voting variables
# voting.name <- c("Percentage working" , "College education", 
#	"Percentage Over 65s", "Population density (persons per sq. mile)", "Percentage white") 
voting.name <- c("PC Employed", "PC College Education", "PC Over 65", 
	"Population Density", "PC White") 
   
county.plot.func <- function(i = 16, pal = brewer.pal(5, "Reds")) {
	val.i <- county@data[, i]
	sh <- auto.shading(val.i, n = 5, cols = pal)
	choropleth(county, v = val.i, shading = sh, border = NA)
	choro.legend("bottomleft", sh = sh, box.lwd = 0, cex = 0.8) 
	tit <- voting.name[i-15]
	title(tit)
}
#quartz(w = 12, h = 5)
setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
png(filename = "F1.png", w = 12, h = 5, units = "in", res = 150)
par(mfrow = c(2,3))
par(mar = c(0,0,2,0))
props <- (county$dem_2016/(county$dem_2016+county$gop_2016)*100)
cols <- colorRampPalette(brewer.pal(9,'RdBu'))(100)
i = cols[findInterval(props, 1:100)]
plot(county, col = i, border = NA, main = "Voting: Democrat (blue) to Republican (red)")
pal.list <- list(brewer.pal(5, "Reds"),
				brewer.pal(5, "Blues"),
				brewer.pal(5, "Greens"),
				brewer.pal(5, "YlOrBr"),
				brewer.pal(5, "OrRd"))
for (i in 16:20) {
	county.plot.func(i, pal = pal.list[[i-15]])
}
dev.off()
## 1.4.2 Species variables  
envvar.plot.func <- function(i = 4, pal = brewer.pal(5, "Reds")) {
	val.i <- nat.sp@data[, i]
	sh <- auto.shading(val.i, n = 5, cols = pal)
	choropleth(nat.sp, v = val.i, shading = sh, cex = 0.6, pch = 19)
	choro.legend("bottomleft", sh = sh, box.lwd = 0, cex = 0.8) 
	tit <- envar.name[i-2]
	title(tit)
}
# gdd+p+pet+stdp+tmp
envar.name <- c("Ratio of actual to potential evapotranspiration", 
	"Growing degree-days above 5C", 
	"Annual amount of precipitation",
	"Potential evapotranspiration", 
	"Annual variation of precipitation", 
	"Maximum temperature of the warmest month", 
	"Minimum temperature of the coldest month",
	"Annual mean temperature")
# quartz(w = 12, h = 5)
setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
png(filename = "F2.png", w = 12, h = 5, units = "in", res = 150)
par(mfrow = c(2,3))
par(mar = c(0,0,2,0))
plot(nat.sp, col = "lightgrey", cex = 0.6, pch = 19)
plot(nat.sp[nat.sp$species_occ == 1,], cex = 0.6, add = T, pch = 19)
title("Species Occurence (presence in black)")
count = 1
for (i in c(4,5,6,7,10)) {
	envvar.plot.func(i, pal = pal.list[[count]])
	count <- count+1
}
dev.off()

### Part 2. Initial exploration of data

## 2.1 Voting Data - Conditional boxplots
d.a <- county[, 10:20]
d.a$PopD <- log(d.a$PopD)
# convert to Factor for ggplots below
d.a$Trump <- as.factor(d.a$Trump)
## 2.1.2 Corellations and conditional boxplots
#plot(d.a@data[,7:11])
cor(d.a@data[,7:11])
### ggplot this for nice boxplots
ylim1 = boxplot.stats(d.a$PCemp)$stats[c(1, 5)]
p1 <- ggplot(d.a@data, aes(Trump, PCemp, fill = Trump)) +
    geom_boxplot(outlier.shape = NA) + 
    coord_cartesian(ylim = ylim1*1.05) + 
    theme(legend.position = "none") +
    labs(subtitle = "PC Employed") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(),
          legend.position = "none") 

#ggplot_build(p1)$data

ylim1 = boxplot.stats(d.a$PCcol)$stats[c(1, 5)]
p2 <- ggplot(d.a@data, aes(Trump, PCcol, fill = Trump)) +
    geom_boxplot(outlier.shape = NA) + 
    coord_cartesian(ylim = ylim1*1.05) + 
    theme(legend.position = "none") +
    labs(subtitle = "PC College Education") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(),
          legend.position = "none") 

ylim1 = boxplot.stats(d.a$PCo65)$stats[c(1, 5)]
p3 <- ggplot(d.a@data, aes(Trump, PCo65, fill = Trump)) +
    geom_boxplot(outlier.shape = NA) + 
    coord_cartesian(ylim = ylim1*1.05) + 
    theme(legend.position = "none") +
    labs(subtitle = "PC Over 65") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(),
          legend.position = "none") 

ylim1 = boxplot.stats(d.a$PopD)$stats[c(1, 5)]
p4 <- ggplot(d.a@data, aes(Trump, PopD, fill = Trump)) +
    geom_boxplot(outlier.shape = NA) + 
    coord_cartesian(ylim = ylim1*1.05) + 
    theme(legend.position = "none") +
    labs(subtitle = "log of Population Density") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(),
          legend.position = "none") 

ylim1 = boxplot.stats(d.a$PCwhi)$stats[c(1, 5)]
p5 <- ggplot(d.a@data, aes(Trump, PCwhi, fill = Trump)) +
    geom_boxplot(outlier.shape = NA) + 
   # coord_cartesian(ylim = ylim1*1.05) + 
    theme(legend.position = "none") + 
    labs(subtitle = "PC White") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(),
          legend.position = "none")     
    
multiplot2 <- function(plot.list, file, cols=3, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  numPlots = length(plot.list)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plot.list[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                          layout.pos.col = matchidx$col))
    }
  }
}
#quartz(w = 12, h = 4)
setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
png(filename = "F3a.png", w = 12, h = 4, units = "in", res = 150)
multiplot2(list(p1,p2,p3,p4,p5), cols = 5)
dev.off()

## 2.2 Species Data - Conditional boxplots
# convert to Factor for ggplots below
nat.sp$species_occ <- as.factor(nat.sp$species_occ)
## 2.1.2 Corellations and conditional boxplots
#plot(nat.sp@data[,c(4,5,6,7,10)], cex = 0.2)
round(cor(nat.sp@data[,c(4,5,6,7,10)]), 3)
### ggplot this for nice boxplots

ylim1 = boxplot.stats(nat.sp$gdd)$stats[c(1, 5)]
p1 <- ggplot(nat.sp@data, aes(species_occ, gdd, fill = species_occ)) +
  geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = ylim1*1.05) + 
  theme(legend.position = "none") + 
    labs(subtitle = envar.name[2]) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(),
          legend.position = "none")   

ylim1 = boxplot.stats(nat.sp$p)$stats[c(1, 5)]
p2 <- ggplot(nat.sp@data, aes(species_occ, p, fill = species_occ)) +
  geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = ylim1*1.05) + 
  theme(legend.position = "none") + 
    labs(subtitle = envar.name[3]) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(),
          legend.position = "none")  

ylim1 = boxplot.stats(nat.sp$pet)$stats[c(1, 5)]
p3 <- ggplot(nat.sp@data, aes(species_occ, pet, fill = species_occ)) +
  geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = ylim1*1.05) + 
  theme(legend.position = "none") + 
    labs(subtitle = envar.name[4]) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(),
          legend.position = "none") 

ylim1 = boxplot.stats(nat.sp$stdp)$stats[c(1, 5)]
p4 <- ggplot(nat.sp@data, aes(species_occ, stdp, fill = species_occ)) +
  geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = ylim1*1.05) + 
  theme(legend.position = "none") + 
    labs(subtitle = envar.name[5]) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(),
          legend.position = "none")  

ylim1 = boxplot.stats(nat.sp$tmp)$stats[c(1, 5)]
p5 <- ggplot(nat.sp@data, aes(species_occ, tmp, fill = species_occ)) +
  geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = ylim1*1.05) + 
  theme(legend.position = "none") + 
    labs(subtitle = envar.name[8]) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(),
          legend.position = "none") 

setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
png(filename = "F3b.png", w = 12, h = 4, units = "in", res = 150)
multiplot2(list(p1,p2,p3,p4,p5), cols = 5)
dev.off()

## 2.3 Cluster and spatial auto-correlations
## morans i of the data

## 2.3.1 Voting data 
rownames(county@data) <- 1:nrow(county@data)
nb = poly2nb(county, queen = T, row.names = rownames(county), snap = 15200)
lw = nb2listw(nb)
mI <- matrix(nrow = 0, ncol = 2)
for (i in 7:11) {
	mi <- moran(d.a@data[,i], lw, length(nb), Szero(lw), zero.policy=T)
	mI <- rbind(mI, unlist(mi))
}
mT <- matrix(nrow = 0, ncol = 2)
for (i in 7:11) {
	mi <- moran.test(d.a@data[,i], lw)
	mT <- rbind(mT, c(mi$estimate[1], mi$p.value))
}
rownames(mI) <- names(d.a[7:11])
rownames(mT) <- names(d.a[7:11])
colnames(mT) <- c("I", "p-value")
moran.tab1 <- round(cbind(mT,mI[,2]), 3) 
colnames(moran.tab1) <- c("I", "p-value", "K")
### suggest spatial dependence in the data of each variable

## 2.3.2 Species data
rownames(nat.sp@data) <- 1:nrow(nat.sp@data)
rm(tmp)
tmp <- spTransform(nat.sp, CRS("+proj=longlat "))
tmp <- as(tmp, "SpatialPixelsDataFrame")
tmp <- as(tmp, "SpatialPolygonsDataFrame")
nb = poly2nb(tmp, queen = T)
lw = nb2listw(nb)
mI <- matrix(nrow = 0, ncol = 2)
for (i in c(4,5,6,7,10)) {
	mi <- moran(tmp@data[,i], lw, length(nb), Szero(lw), zero.policy=T)
	mI <- rbind(mI, unlist(mi))
}
mT <- matrix(nrow = 0, ncol = 2)
for (i in c(4,5,6,7,10)) {
	mi <- moran.test(tmp@data[,i], lw)
	mT <- rbind(mT, c(mi$estimate[1], mi$p.value))
}
rownames(mI) <- names(tmp@data[c(4,5,6,7,10)])
rownames(mT) <- names(tmp@data[c(4,5,6,7,10)])
colnames(mT) <- c("I", "p-value")
moran.tab2 <- round(cbind(mT,mI[,2]), 3) 
colnames(moran.tab2) <- c("I", "p-value", "K")
voting.name <- c("PC Employed", "PC College Education", "PC Over 65", 
	"Population Density", "PC White") 
rownames(moran.tab1) <- voting.name

envar.name <- c("Ratio of actual to potential evapotranspiration", 
	"Growing degree-days above 5C", 
	"Annual amount of precipitation",
	"Potential evapotranspiration", 
	"Annual variation of precipitation", 
	"Maximum temperature of the warmest month", 
	"Minimum temperature of the coldest month",
	"Annual mean temperature")

rownames(moran.tab2) <- envar.name[c(2,3,4,5,8)]
setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')

write.csv(moran.tab1, file = "Tab1a.csv")
write.csv(moran.tab2, file = "Tab1b.csv")

### suggest spatial dependence in the data of each variable


### Part 3. Initial evaluation of GLMs and GGWRs
##### with rescaled data
d.a.voting <- d.a
d.a.species <- nat.sp
d.a.voting@data[,c(7:11)] <- apply(d.a.voting@data[,c(7:11)], 2, function(x) rescale(x, c(0.001, 1)))
d.a.species@data[,c(3:10)] <- apply(d.a.species@data[,c(3:10)], 2, function(x) rescale(x, c(0.001, 1)))

## 3.1 Voting Data 
# convert back from Factor
d.a.voting$Trump <- as.numeric(as.character(d.a.voting$Trump))
#colnames(data.mat) <- c("PCemp", "PCcol", "PCo65", "PopD", "PCwhi")
EUDM <-  gw.dist(coordinates(d.a.voting))
## GLM
glm.tr <- ggwr.basic(Trump~PCemp+PCcol+PCo65+PopD+PCwhi, 
	kernel = "boxcar", bw = max(EUDM+10000), 
	data = d.a.voting, adaptive = F, family = "binomial", dMat = EUDM)
## GGWR
bw.tr <- bw.ggwr(Trump~PCemp+PCcol+PCo65+PopD+PCwhi, data = d.a.voting, adaptive = T, 
	family = "binomial", approach = "AIC", dMat = EUDM)
bw.tr /nrow(d.a) 
ggwr.tr <- ggwr.basic(Trump~PCemp+PCcol+PCo65+PopD+PCwhi, bw = bw.tr, 
	data = d.a.voting, adaptive = T, family = "binomial", dMat = EUDM)
## ELN - voting
x <- as.matrix(d.a.voting@data[c("PCemp", "PCcol", "PCo65", "PopD", "PCwhi")])
y <- as.matrix(d.a.voting@data["Trump"])
eln.tr <- glmnet(x = x, y = y, family = "binomial", standardize = T)

## Predictions - Voting
# GLM
glm.tr.pred <- gwr.predict(Trump~PCemp+PCcol+PCo65+PopD+PCwhi, bw = max(EUDM+10000), 
	data = d.a.voting, predictdata = d.a.voting, adaptive = F, 
	dMat1 = EUDM,
	dMat2 = EUDM)
fitted.results <- glm.tr.pred$SDF$prediction
fitted.results <- as.vector(ifelse(fitted.results > 0.9,1,0))
glm.tr.predCorrect = mean(fitted.results == d.a.voting$Trump)
# GGWR
ggwr.tr.pred <- gwr.predict(Trump~PCemp+PCcol+PCo65+PopD+PCwhi, bw = bw.tr, 
	data = d.a.voting, predictdata = d.a.voting, adaptive = T, 
	dMat1 = EUDM,
	dMat2 = EUDM)
fitted.results <- ggwr.tr.pred$SDF$prediction
fitted.results <- as.vector(ifelse(fitted.results > 0.9,1,0))
ggwr.tr.predCorrect = mean(fitted.results == d.a.voting$Trump)
# ELM
pred.vals <- as.matrix(d.a.voting@data[
	c("PCemp", "PCcol", "PCo65", "PopD", "PCwhi")])
predictions = predict(eln.tr, newx = pred.vals, 
	s=c(0.1), type = "class") 
eln.tr.predCorrect <- mean(predictions == d.a.voting@data["Trump"])

## 3.2 Species Data 
# convert back from Factor
d.a.species$species_occ <- as.numeric(as.character(nat.sp$species_occ))
# actually use absence - too few presence for convegence
d.a.species$PA <- (nat.sp$species_occ == 0) + 0
# gdd, pet, stdp, tmp
EUDM.sp <-  gw.dist(coordinates(d.a.species))
## GLM
glm.sp <- ggwr.basic(PA~gdd+p+pet+stdp+tmp, 
      kernel = "boxcar", bw = max(EUDM.sp+10000), 
      data = d.a.species, adaptive = F, family = "binomial", dMat = EUDM.sp)
## GGWR
bw.sp <- bw.ggwr(PA~gdd+p+pet+stdp+tmp, data = d.a.species, adaptive = F, 
      family = "binomial", approach = "AIC", dMat = EUDM.sp)
bw.sp  /1000 
ggwr.sp <- ggwr.basic(PA~gdd+p+pet+stdp+tmp, bw = bw.sp, 
      data = d.a.species, adaptive = F, family = "binomial", dMat = EUDM.sp)
## ELN - species
x <- as.matrix(d.a.species@data[c("gdd","p", "pet", "stdp", "tmp")])
y <- as.matrix(d.a.species@data["PA"])
eln.sp <- glmnet(x = x, y = y, family = "binomial", standardize = T)

## Predictions - Species
# GLM
glm.sp.pred <- gwr.predict(PA~gdd+p+pet+stdp+tmp, bw = max(EUDM.sp+10000), 
	data = d.a.species, predictdata = d.a.species, adaptive = F, 
	dMat1 = EUDM.sp,
	dMat2 = EUDM.sp)
fitted.results <- glm.sp.pred$SDF$prediction
fitted.results <- as.vector(ifelse(fitted.results > 0.9,1,0))
glm.sp.predCorrect = mean(fitted.results == d.a.species$PA)
# GGWR
ggwr.sp.pred <- gwr.predict(PA~gdd+p+pet+stdp+tmp, bw = bw.sp, 
	data = d.a.species, predictdata = d.a.species, adaptive = F, 
	dMat1 = EUDM.sp,
	dMat2 = EUDM.sp)
fitted.results <- ggwr.sp.pred$SDF$prediction
fitted.results <- as.vector(ifelse(fitted.results > 0.9,1,0))
ggwr.sp.predCorrect = mean(fitted.results == d.a.species$PA)
# ELN
pred.vals <- as.matrix(d.a.species@data[
	c("gdd","p", "pet", "stdp", "tmp")])
predictions = predict(eln.sp, newx = pred.vals, 
	s=c(0.1), type = "class") 
eln.sp.predCorrect <- mean(predictions == d.a.species@data["PA"])

## Combine to a table
tab <- matrix(c( 
	glm.tr$GW.diagnostic$AICc, 
	glm.tr$GW.diagnostic$gw.R2, 
	glm.tr$GW.diagnostic$gwR2.adj, 
	glm.tr.predCorrect,
	eln.tr.predCorrect,
	bw.tr/nrow(d.a.voting)*100,
	ggwr.tr$GW.diagnostic$AICc, 
	ggwr.tr$GW.diagnostic$gw.R2, 
	ggwr.tr$GW.diagnostic$gwR2.adj,
	ggwr.tr.predCorrect,
	glm.sp$GW.diagnostic$AICc, 
	glm.sp$GW.diagnostic$gw.R2, 
	glm.sp$GW.diagnostic$gwR2.adj, 
	glm.sp.predCorrect,
	eln.sp.predCorrect,
	bw.sp/1000,
	ggwr.sp$GW.diagnostic$AICc, 
	ggwr.sp$GW.diagnostic$gw.R2, 
	ggwr.sp$GW.diagnostic$gwR2.adj,
	ggwr.sp.predCorrect), ncol = 2)
rownames(tab) <- c("GLM AICc", "GLM R2", "GLM R2 adj", "GLM prediction", "ENLR prediction",
					"GGWR Bandwith (km or %)", "GGWR AICc", "GGWR R2", "GGWR R2 adj", "GGWR prediction")
colnames(tab) <- c("Voting", "Species P-A")	
write.csv(round(tab, 3), file = "Tab2.csv")

## Coefficient tab2
c1 <- unlist(c(as.vector(glm.tr$SDF@data[1, 1:6]), as.vector(glm.sp$SDF@data[1, 1:6]))) 
c2 <- unlist(c(as.vector(coef(eln.tr,s=0.1)), as.vector(coef(eln.sp,s=0.1)))) 
c3<- rbind(t(apply(ggwr.tr$SDF@data[, 1:6], 2, summary)[c(2,3,5),]),
			t(apply(ggwr.sp$SDF@data[, 1:6], 2, summary)[c(2,3,5),]))
tab3 <- cbind(cbind(c1, c2), c3)
colnames(tab3) <- c("GLM", "ENLR", "GGWR 1stQ", "GGWR median", "GGWR 3rdQ")
rownames(tab3) <- c("Intercept", rownames(moran.tab1), "Intercept", rownames(moran.tab2))
write.csv(round(tab3, 3), file = "Tab3.csv")

## 4. GW elastic net

## 4.1 BW selection for Voting 
# for each bw try to predict each point
# using data minus that point & sum scores for each BW
# Trump with adaptive bandwidth 

switchlex = "DONE"
#switchlex = "NOTDONE"
if (switchlex != "DONE" ){
	d.a <- d.a.voting
	dMat <- as.matrix(dist(coordinates(d.a))) 
	bwd.vec <- round(c(seq(0.01,1,by=0.01))*nrow(d.a), 0)
	cv.errors <- rep(NA, length(bwd.vec))
	na.count <- rep(NA, length(bwd.vec))
	for (i in 1: length(bwd.vec)) {
		bw <- bwd.vec[i]
		cv.error.j = 0
		na.count.j = 0	
		for (j in 1:nrow(d.a)) {
			# point
			pt.j <- d.a[j,]
			# index for not point within bw
			index.j <- order(dMat[j, ])[1:bw]
			index.j <- setdiff(index.j, j)
			# select data
			data.j <- d.a[index.j, ]
			#colnames(data.mat) <- c("PCemp", "PCcol", "PCo65", "PopD", "PCwhi")
			x <- as.matrix(data.j@data[c("PCemp", "PCcol", "PCo65", "PopD", "PCwhi")])
			y <- as.matrix(data.j@data["Trump"])
			pred.vals <- pt.j@data[c("PCemp", "PCcol", "PCo65", "PopD", "PCwhi")]
			#plot(d.a, col = "lightgrey", cex = 0.3)
			#plot(data.j, cex = 0.5, add = T)
			#plot(pt.j, add = T, col = "red")
			# check for 0s and 1s
			if (sum(y) > 4 & sum(y) < (length(y)-2)) {
				# dists
				# weights
				dists.j <- gw.dist(coordinates(pt.j), coordinates(data.j))
		    		weights.j <- as.vector(gw.weight(dMat[j,],bw, 
		    			kernel = "bisquare", adaptive = T))	

		    		weights.j <- as.vector(gw.weight(dists.j,bw, 
		    			kernel = "bisquare", adaptive = T))
		    		x <- x * weights.j		
				#plot(data.j, cex = weights.j, pch = 19)
				#plot(pt.j, col = "red", pch = 19, add = T)
				
				eln.mod <- glmnet(x = x, y = y, 
					family = "binomial", standardize = T)
			    pred.vals <- d.a[c(j, index.j),]
				pred.vals <- as.matrix(pred.vals@data[
					c("PCemp", "PCcol", "PCo65", "PopD", "PCwhi")])
			    predictions = predict(eln.mod, newx = pred.vals, 
			    		s=c(0.1), type = "class") 
				#testing s
			    #p = predict(eln.mod, newx = pred.vals, 
			    	#	s=c(0.1), type = "class") 
			    #mean(as.numeric(predictions)==as.numeric(p))
 				#sum(as.numeric(predictions)) 
 				#sum(as.numeric(p))
 				#(as.numeric(predictions[1]) == pt.j@data["Trump"])+0
			    
			    cv.error.j = cv.error.j + 
			    		(as.numeric(predictions[1]) == pt.j@data["Trump"])+0
			}   else  {
				na.count.j <- na.count.j + 1
			}
			#cat("\t", j)
		}
		cv.errors[i] <- cv.error.j
		na.count[i] <- na.count.j
		cat("\t", i)		
	}
	cv.errors.tr <- cv.errors
	na.count.tr <- na.count
	setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
	#save.image("new_pap_v2_pt2.RData")
	save(list = c("cv.errors.tr", "na.count.tr"), file = "new_pap_v3_pt1_bw.RData")
} else {
	setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
	load("new_pap_v2_pt1_bw.RData")
}


## 4.2 Bandwidth selection for Species P-A
# for each bw try to predict each point
# using data minus that point & sum scores for each BW
# PA with fixed bandwidth 

#switchlex = "DONE"
if (switchlex != "DONE" ){
	d.a <- d.a.species
	EUDM <-  gw.dist(coordinates(d.a))
	bwd.vec <- c(seq(50000,max(EUDM+10000),by=50000))
	dMat <- EUDM 
	cv.errors <- rep(NA, length(bwd.vec))
	na.count <- rep(NA, length(bwd.vec))
	for (i in 1: length(bwd.vec)) {
	  bw <- bwd.vec[i]
	  cv.error.j = 0
	  na.count.j = 0	
	  for (j in 1:nrow(d.a)) {
	    # point
	    pt.j <- d.a[j,]
	    # index for not point within bw
	    # index.j <- order(dMat.i)[ 1:bw.gwda.ab]  
	    index.j <- which(dMat[j,] < bw)
	    index.j <- setdiff(index.j, j)
	    # select data
	    data.j <- d.a[index.j, ]
	    #colnames(data.mat) <- c("PCemp", "PCcol", "PCo65", "PopD", "PCwhi")
	    x <- as.matrix(data.j@data[c("gdd","p", "pet", "stdp", "tmp")])
	    y <- as.matrix(data.j@data["PA"])
	    pred.vals <- pt.j@data[c(c("gdd","p", "pet", "stdp", "tmp"))]
	    #plot(d.a, col = "lightgrey", cex = 0.3)
	    #plot(data.j, cex = 0.5, add = T)
	    #plot(pt.j, add = T, col = "red")
	    # check for 0s and 1s
	    if (sum(y) > 4 & sum(y) < (length(y)-2)) {
	      # dists
	      # weights
	      dists.j <- gw.dist(coordinates(pt.j), coordinates(data.j))
	      weights.j <- as.vector(gw.weight(dists.j,bw, 
				kernel = "bisquare", adaptive = F))
	      x <- x * weights.j
	      #plot(data.j, cex = weights.j, pch = 19)
	      #plot(pt.j, col = "red", pch = 19, add = T)
	     
	      eln.mod <- glmnet(x = x, y = y, 
				family = "binomial", standardize = T)
	      pred.vals <- d.a[c(j, index.j),]
	      pred.vals <- as.matrix(pred.vals@data[
	        	c(c("gdd","p", "pet", "stdp", "tmp"))])
	      predictions = predict(eln.mod, newx = pred.vals, 
	            s=c(0.1), type = "class") 
	      cv.error.j = cv.error.j + 
	        (as.numeric(predictions[1]) == pt.j@data["PA"])+0
	    }   else  {
	      na.count.j <- na.count.j + 1
	    }
	    #cat("\t", j)
	  }
	  cv.errors[i] <- cv.error.j
	  na.count[i] <- na.count.j
	  cat("\t", i)		
	}
	cv.errors.sp <- cv.errors
	na.count.sp <- na.count
	setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
	#save.image("new_pap_v3_pt2.RData")
	save(list = c("cv.errors.sp", "na.count.sp"), file = "new_pap_v3_pt2_bw.RData")
} else {
	setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
	load("new_pap_v2_pt2_bw.RData")
}

### 4.3 GW elastics net logistic regression - Voting
## 4.3.1 visuisation of bw 
# 	save(list = c("cv.errors.tr", "na.count.tr"), file = "new_pap_v2_pt1_bw.RData")
nrow(d.a.voting)
bwd.vec <- round(c(seq(0.01,1,by=0.01))*nrow(d.a.voting), 0)
#bw = bwd.vec[which(cv.errors.tr/(nrow(d.a.voting)-na.count.tr) == 
#	max(cv.errors.tr/(nrow(d.a.voting)-na.count.tr)))]

bw = bwd.vec[which.max(cv.errors.tr)] 
bw/nrow(d.a.voting) # proportion
y = cv.errors.tr/nrow(d.a.voting)
#y = cv.errors.tr/(nrow(d.a.voting)-na.count.tr)
x <- (bwd.vec / nrow(d.a.voting)) * 100
xy <- data.frame(x,y)
tit <- paste("GW ENLR adaptive bandwidth: ", 100*round(bw/nrow(d.a.voting),2), "%",sep = "")
setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
png(filename = "F4a.png", w = 6, h = 3.5, units = "in", res = 150)
ggplot() + 
    geom_point(data = xy, aes(x=x, y=y), size = 1) +
    geom_line(data = xy, aes(x=x, y=y)) +
    labs(
    		subtitle = tit, 
    		x = "Kernel size: % of data points", 
    	y = "Proportion of points correctly predicted")
dev.off()
## 4.3.2 GW ENLR voting
# rescale
d.a <- d.a.voting
# d.a@data[,c(7:11)] <- apply(d.a@data[,c(7:11)], 2, function(x) rescale(x, c(0.001, 1)))
dMat <- as.matrix(dist(coordinates(d.a))) 
coef.mat <- matrix(ncol = 6, nrow = nrow(d.a))
pred.correct.tr <- vector()
for (j in 1:nrow(d.a)) {
	pt.j <- d.a[j,]
	# index for not point within bw
	index.j <- order(dMat[j, ])[1:bw]
	index.j <- setdiff(index.j, j)
	# select data
	data.j <- d.a[index.j, ]
	#colnames(data.mat) <- c("PCemp", "PCcol", "PCo65", "PopD", "PCwhi")
	x <- as.matrix(data.j@data[c("PCemp", "PCcol", "PCo65", "PopD", "PCwhi")])
	y <- as.matrix(data.j@data["Trump"])
	pred.vals <- pt.j@data[c("PCemp", "PCcol", "PCo65", "PopD", "PCwhi")]
	#plot(d.a, col = "lightgrey", cex = 0.3)
	#plot(data.j, cex = 0.5, add = T)
	#plot(pt.j, add = T, col = "red")
	# check for 0s and 1s
	if (sum(y) > 1 & sum(y) < (length(y)-2)) {
		# dists
		# weights
		dists.j <- gw.dist(coordinates(pt.j), coordinates(data.j))
		weights.j <- as.vector(gw.weight(dists.j,bw, 
		kernel = "bisquare", adaptive = T))
		x <- x * weights.j		
		#plot(data.j, cex = weights.j, pch = 19)
		#plot(pt.j, col = "red", pch = 19, add = T)
		eln.mod <- glmnet(x = x, y = y, 
			family = "binomial", standardize = T)
		pred.vals <- d.a[c(j, index.j),]
		pred.vals <- as.matrix(pred.vals@data[
			c("PCemp", "PCcol", "PCo65", "PopD", "PCwhi")])
		predictions = predict(eln.mod, newx = pred.vals, 
			s=c(0.0001), type = "class") 
		pred.correct.tr = append(pred.correct.tr,   
	    		(as.numeric(predictions[1]) == pt.j@data["Trump"])+0)					
	    coef.j <- as.vector(coef(eln.mod,s=0.1))
	}   else  {
		coef.j <- c(NA,NA,NA,NA,NA,NA)
		pred.correct.tr =  append(pred.correct.sp, NA)
	}
	coef.mat[j,] <- coef.j
	if(j %% 100 == 0) cat(j)
}
coef.mat.tr <- coef.mat
summary(coef.mat.tr)
sum(pred.correct.tr)/nrow(d.a.voting)

### 4.4 GW elastics net logistic regression - Species
## 4.4.1 visuisation of bw 
nrow(d.a.species)
EUDM <-  gw.dist(coordinates(nat.sp))
bwd.vec <- c(seq(50000,max(EUDM+10000),by=50000))
#bw = bwd.vec[which(cv.errors.sp/(nrow(d.a.species)-na.count.sp) == 
#	max(cv.errors.sp/(nrow(d.a.species)-na.count.sp), na.rm = T))]

bw = bwd.vec[which.max(cv.errors.sp)]
bw/1000 # km
bw / max(EUDM)
y = cv.errors.sp/nrow(d.a.species)
#y = cv.errors.sp/(nrow(d.a.species)-na.count.sp)
x <- bwd.vec / 1000
xy <- data.frame(x,y)
tit <- paste("GW ENLR fixed bandwidth: ", bw/1000, "km",sep = "")
setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
png(filename = "F4b.png", w = 6, h = 3.5, units = "in", res = 150)
ggplot() + 
    geom_point(data = xy, aes(x=x, y=y), size = 1) +
    geom_line(data = xy, aes(x=x, y=y)) +
    labs(subtitle = tit, 
    		x = "Kernel size: km", 
    	y = "Proportion of points correctly predicted")
dev.off()
## 4.4.2 GW ENLR species
d.a <- d.a.species
dMat <- as.matrix(dist(coordinates(d.a))) 
# rescale
# d.a@data[,c(3:10)] <- apply(d.a@data[,c(3:10)], 2, function(x) rescale(x, c(0.001, 1)))
coef.mat <- matrix(ncol = 6, nrow = nrow(d.a))
pred.correct.sp <- vector()
for (j in 1:nrow(d.a)) {
	pt.j <- d.a[j,]
	# index for not point within bw
	index.j <- which(dMat[j,] < bw)
    index.j <- setdiff(index.j, j)# select data
	data.j <- d.a[index.j, ]
    x <- as.matrix(data.j@data[c("gdd","p","pet", "stdp", "tmp")])
	y <- as.matrix(data.j@data["species_occ"])
	pred.vals <- pt.j@data[c(c("gdd","p","pet", "stdp", "tmp"))]
	#plot(d.a, col = "lightgrey", cex = 0.3)
	#plot(data.j, cex = 0.5, add = T)
	#plot(pt.j, add = T, col = "red")
	# check for 0s and 1s
	if (sum(y) > 1 & sum(y) < (length(y)-2)) {
		# dists
		# weights
		dists.j <- gw.dist(coordinates(pt.j), coordinates(data.j))
		weights.j <- as.vector(gw.weight(dists.j,bw, 
		kernel = "bisquare", adaptive = T))
		x <- x * weights.j		
		#plot(data.j, cex = weights.j, pch = 19)
		#plot(pt.j, col = "red", pch = 19, add = T)
		eln.mod <- glmnet(x = x, y = y, 
			family = "binomial", standardize = T)
		pred.vals <- d.a[c(j, index.j),]
	    pred.vals <- as.matrix(pred.vals@data[
	    		c(c("gdd","p", "pet", "stdp", "tmp"))])
	    predictions = predict(eln.mod, newx = pred.vals, 
	    		s=c(0.0001), type = "class") 
	    pred.correct.sp =  append(pred.correct.sp, 
	        (as.numeric(predictions[1]) == pt.j@data["species_occ"])+0)
	    coef.j <- as.vector(coef(eln.mod,s=0.1))
	}   else  {
		coef.j <- c(NA,NA,NA,NA,NA,NA)
		pred.correct.sp =  append(pred.correct.sp, NA)
	}
	coef.mat[j,] <- coef.j
	if(j %% 100 == 0) cat(j)
}
coef.mat.sp <- coef.mat
summary(coef.mat.sp)
index <- is.na(pred.correct.sp)
sum(pred.correct.sp, na.rm = T)/(nrow(d.a.species))
sum(pred.correct.tr)/nrow(d.a.voting)

setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
save.image("anal_v3_all.RData")

### 4.5 Final Analyses
## coefficient tab
colnames(coef.mat.tr) <- c("Intercept", rownames(moran.tab1))
tab4a <- t(apply(coef.mat.tr, 2, fivenum))
colnames(tab4a) <- c("Min", "1stQ", "Median", "3rdQ", "Max")

colnames(coef.mat.sp) <- c("Intercept", rownames(moran.tab2))
tab4b <- t(apply(coef.mat.sp, 2, fivenum))
colnames(tab4b) <- c("Min", "1stQ", "Median", "3rdQ", "Max")

tab4 <- rbind(tab4a, tab4b)
setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
write.csv(round(tab4, 3), file = "Tab4.csv")

## Maps of GW-EMLR
summary(coef.mat.tr)
voting.name <- c("PC Employed", "PC College Education", "PC Over 65", 
	"Population Density", "PC White") 
#quartz(w = 11, h = 3.5)
setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
png(filename = "F5a.png", w = 11, h = 3, units = "in", res = 150)
par(mfrow = c(1,3))
par(mar = c(0,0,2,0))
pal.list <- c("Reds", "Blues", "Greens", "YlOrBr", "OrRd")
for (i in c(3,5,6)) {
	props <- (coef.mat.tr[,i])
	cols <- colorRampPalette(rev(brewer.pal(9, pal.list[i-1])))(100)
	j = cols[findInterval(props, seq(min(props), max(props), c( (max(props)-min(props)) / 100)))]
	j[is.na(j)] <- "#CCCCCC"
	tit <- paste(colnames(coef.mat.tr)[i], "\n darkest ", 
		format(round(min(props, na.rm = T), 2), nsmall = 2),
		", lightest ", 
		format(round(max(props, na.rm = T), 2), nsmall = 2),
	 	sep = "")
	plot(county, col = j, border = NA, main = tit)
}
dev.off()
#check
#quartz()
#val.i <- (coef.mat.tr[,3])
#	sh <- auto.shading(val.i, n = 9, cols = brewer.pal(9, "Reds"))
#	choropleth(county, v = val.i, shading = sh, border = NA)
#	choro.legend("bottomleft", sh = sh, box.lwd = 0, cex = 0.8) 
	
summary(coef.mat.sp)
#quartz(w = 11, h = 3)
setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
png(filename = "F5b.png", w = 11, h = 3, units = "in", res = 150)
par(mfrow = c(1,3))
par(mar = c(0,0,2,0))
pal.list <- c("Reds", "Blues", "Greens", "YlOrBr", "OrRd")
for (i in c(2, 4, 6)) {
	props <- (coef.mat.sp[,i])
	cols <- colorRampPalette(rev(brewer.pal(9, pal.list[i-1])))(100)
	j = cols[findInterval(props, seq(min(props, na.rm = T), max(props, na.rm = T), 
		c( (max(props, na.rm = T)-min(props, na.rm = T)) / 100)))]
	j[is.na(j)] <- "#CCCCCC"
	tit <- paste(colnames(coef.mat.sp)[i], "\n darkest ", 
		format(round(min(props, na.rm = T), 2), nsmall = 2),
		", lightest ", 
		format(round(max(props, na.rm = T), 2), nsmall = 2),
	 	sep = "")
	plot(nat.sp, col = j, main = tit,  pch = 19, cex = 0.6)

}
dev.off()

sum(coef.mat.sp[,3] < 0, na.rm = T)
sum(coef.mat.sp[,5] > 0, na.rm = T)

setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/')
save.image("anal_v3_all.RData")


#### END ####