
library(repmis)
library(maps)
library(GISTools)
library(ecospat)
library(tidyverse)
library(scales)
library(GWmodel)
library(glmnet)

##### Data
# download the zip file at https://github.com/lexcomber/GW-ENLR
# put into a local directory
# set the working directory
# eg setwd('~/Desktop/')
 
##### FUNCTIONS
## 1. GGWR-ELN BW Function
bw.ggwreln <- function(d.a = d.a.voting,
	x.term = c("PCemp", "PCcol", "PCo65", "PopD", "PCwhi"),
	y.term = "Trump", 
	dMat = as.matrix(dist(coordinates(d.a))),
	lambda = 0.1,
	alpha = 0.75,
	adaptive = T,
	adptbwd.thresh = 0.1,
	verbose = F){	 
		if(adaptive) bwd.vec <- round(c(seq(adptbwd.thresh,1,
			length.out = 100))*nrow(d.a), 0)
		if(!adaptive) bwd.vec <- round(seq(max(dMat) * 0.01, max(dMat)*1.01,length.out = 100))
		cv.correct <- rep(NA, length(bwd.vec))
		na.count <- rep(NA, length(bwd.vec))
		for (i in 1: length(bwd.vec)) {
			bw <- bwd.vec[i]
			cv.correct.j = 0
			na.count.j = 0	
			for (j in 1:nrow(d.a)) {				
				pt.j <- d.a[j,]
				if(adaptive) index.j <- order(dMat[j, ])[1:bw]
				if(!adaptive) index.j <- which(dMat[j,] < bw)	
				data.j <- d.a[index.j, ]
				y <- as.matrix(data.j@data[y.term])
				# NA check for sims
				index.jj <- as.vector(is.na(y))
				if (sum(index.jj) > 0) data.j <- data.j[!index.jj, ]
				y <- as.matrix(data.j@data[y.term])
				x <- as.matrix(data.j@data[x.term])			
								
				if (sum(y) > 4 & sum(y) < (length(y)-2)) {
					dists.j <- gw.dist(coordinates(pt.j), coordinates(data.j))
			    		weights.j <- as.vector(gw.weight(c(dists.j),bw, 
			    			kernel = "bisquare", adaptive = adaptive))
			    		x <- x * weights.j		
					eln.mod <- glmnet(x = x, y = y, 
						family = "binomial", standardize = T, alpha = alpha)
					pred.vals <- as.matrix(data.j@data[x.term])
					predictions = predict(eln.mod, newx = pred.vals, 
				    		s=c(lambda), type = "class")  
				    cv.correct.j = sum(cv.correct.j, 
				    		(as.numeric(predictions[1]) == pt.j@data[y.term])+0,
				    		na.rm = T)
				}   else  {
					na.count.j <- na.count.j + 1
				}
			#cat(cv.correct.j, "\t")
			}
			cv.correct[i] <- cv.correct.j
			na.count[i] <- na.count.j
			if(verbose) cat("\t", i)		
		}
		bw = bwd.vec[which.max(cv.correct)] 	
		return(list(cv.correct, na.count, bw))
}
## 2. GGWR ELN function
ggwreln <- function(d.a = d.a.voting,
	bw = bw,
	x.term = c("PCemp", "PCcol", "PCo65", "PopD", "PCwhi"),
	y.term = "Trump", 
	dMat = as.matrix(dist(coordinates(d.a))),
	lambda = 0.02,
	alpha = 0.75,
	adaptive = T,
	verbose = F) {
	
	coef.mat <- matrix(nrow = nrow(d.a), ncol = length(x.term)+1)
	local_CN <- numeric(nrow(d.a))
	y.hat <- numeric(nrow(d.a))
	for (j in 1:nrow(d.a)) {
		pt.j <- d.a[j,]
		if(adaptive) index.j <- order(dMat[j, ])[1:bw]
		if(!adaptive) index.j <- which(dMat[j,] < bw)	
		data.j <- d.a[index.j, ]
		y <- as.matrix(data.j@data[y.term])
		# NA check for sims
		index.jj <- as.vector(is.na(y))
		if (sum(index.jj) > 0) data.j <- data.j[!index.jj, ]
		y <- as.matrix(data.j@data[y.term])
		x <- as.matrix(data.j@data[x.term])			
		local_CN[j] <- BKWcn(cbind(1,x))
		if (sum(y) > 4 & sum(y) < (length(y)-2)) {
			dists.j <- gw.dist(coordinates(pt.j), coordinates(data.j))
	    	weights.j <- as.vector(gw.weight(c(dists.j),bw, 
	    			kernel = "bisquare", adaptive = adaptive))
	    	x <- x * weights.j		
			eln.mod <- glmnet(x = x, y = y, 
				family = "binomial", standardize = T, alpha = alpha)
			coef.mat[j,] <- as.vector(coef(eln.mod, s = lambda))
			x <- as.matrix(data.j@data[x.term])		
		    y.hat[j] = as.numeric(predict(eln.mod, newx = x, 
				s=lambda, type = "class")[1] )
		}   else  {
			coef.mat[j,] <- NA
			y.hat[j] <- NA
		}
	}
	return(list(coef.mat, y.hat, local_CN))
}

## 3. Helper functions
# convert data to sp
make.sp <- function(NsHc) {
	coords <- NsHc[,c("Coord.X","Coord.Y")]
	sp <- SpatialPointsDataFrame(coords, data = data.frame(NsHc))
	return(sp)
}
# Condition number
BKWcn <- function(X) {
	p <- dim(X)[2]
	Xscale <- sweep(X, 2, sqrt(colSums(X^2)), "/")
	Xsvd <- svd(Xscale)$d
	Xsvd[1] / Xsvd[p]}

##### Simulation 
# assumes data have been downloded, unzipped to local directory and that the working dorectory has been set

## declare some objects
st <- Sys.time()
sims = 100
x.term = c("x1", "x2", "x3")
y.term = "y"		
cn.global.vec <- numeric(sims) # Global CN
bw.sims.vec <- numeric(sims) # BWs for each run from CV
local_CN.mat <- matrix(nrow = 159, ncol = sims) # counts local CN > 30 for each sim, each loc
coef.arr <- array(data = 0, dim = c(159,length(x.term)+1,sims)) #  SDF for each sim
y.arr <- array(data = 0, dim = c(159,2,sims)) # yhat for each sim
bw.cv <- matrix(nrow = 100, ncol = sims) # bw by sims

for (k in 29:sims) {
	## Read Data and Set up	
	tit <- paste0("sim_data_sp6_", k, ".csv")
	d.a <- read.csv(tit) 
	d.a <- make.sp(d.a)
	d.a$y <- d.a$Binom.y
	X <- as.matrix(cbind(1, d.a@data[x.term]))
	cn.global.vec[k] <- BKWcn(X)
	# do adapt BW
	dMat = as.matrix(dist(coordinates(d.a))) 
	res.k <- bw.ggwreln(d.a = d.a,
		x.term = x.term,
		y.term = y.term,
		dMat = dMat,
		lambda = 0.02,
		adaptive = T,
		adptbwd.thresh = 0.1,
		verbose = F)
	bw <- res.k[[3]]
	# do GGWR-ELN
	res2.k <- ggwreln(d.a = d.a,
		bw = bw,
		x.term = x.term,
		y.term = y.term,
		dMat = dMat,
		lambda = 0.02,
		adaptive = T,
		verbose = F)
	bw.sims.vec[k] <- bw
	local_CN.mat[,k] <- (res2.k[[3]] > 30)+0
	coef.arr[,,k] <- res2.k[[1]]
	y.arr[,,k] <- cbind(y = d.a$y, res2.k[[2]])
	bw.cv[,k] = res.k[[1]]/(nrow(d.a)-res.k[[2]])
	cat("\t", k)		
}
Sys.time() - st

setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/resub')
save.image("Comber_Harris_GW-ENLR_code_resub6_sims.RData")

##### Load Code 
setwd('~/Desktop/my_docs_mac/leeds_work/research/conservation/fb_approach/resub')
load("Comber_Harris_GW-ENLR_code_resub6_sims.RData")

library(spgwr)
data(georgia)

###### Simulation assessment
# where local collinearity was found how well did GW-ELN do
# 1. variable shrinkage maps
# for each x how many were shrunk at each location? (1)
# for each x, which were shrink at each location? (3)
# 2. dealling with local collin maps 
# how many times was CN identified at each location (1)
# of those times how manywhat proprtion were dealt with by variable shrinkage (1)

# 1. variable shrinkage maps
# for each x how many were shrunk at each location? (1)
# for each x, which were shrink at each location? (3)
count.vec <- numeric(nrow(d.a))
for( i in 1:sims) {
  coef.i <- coef.arr[,,i]
  count.i <- apply(coef.i[,2:4], 1,function(x) any(x==0))+0
  count.i[is.na(count.i)] <- 0
  count.vec <- count.vec+count.i
}
gSRDF$shrink.prop <- count.vec/sims

x1.vec <- numeric(nrow(d.a))
x2.vec <- numeric(nrow(d.a))
x3.vec <- numeric(nrow(d.a))
for( i in 1:sims) {
  coef.i <- coef.arr[,,i]
  
  count.i <- (coef.i[,2] == 0)+0
  count.i[is.na(count.i)] <- 0
  x1.vec <- x1.vec + count.i 
 
  count.i <- (coef.i[,3] == 0)+0
  count.i[is.na(count.i)] <- 0
  x2.vec <- x2.vec + count.i 
  
  count.i <- (coef.i[,4] == 0)+0
  count.i[is.na(count.i)] <- 0
  x3.vec <- x3.vec + count.i 
}
# x1.vec+x2.vec+x3.vec
gSRDF$x1.prop <- x1.vec/sims
gSRDF$x2.prop <- x2.vec/sims
gSRDF$x3.prop <- x3.vec/sims

# 2. dealling with local collin maps 
# how many times was CN identified at each location (1)
gSRDF$collin.count <- rowSums(local_CN.mat, na.rm = T)
# of those times how manywhat proprtion were dealt with by variable shrinkage (1)
collin.vec <- numeric(nrow(d.a))
collin.and.shrink.vec <- numeric(nrow(d.a))
for (i in 1:sims) {
  local.cn.i <- local_CN.mat[,i]
  coef.i <- coef.arr[,,i]
  for (j in 1:nrow(d.a)) {
    if(local.cn.i[j] == 1 & any(coef.i[j,2:4]==0, na.rm = T)) 
      collin.and.shrink.vec[j] <- collin.and.shrink.vec[j] + 1
    if(local.cn.i[j] == 1) 
      collin.vec[j] <- collin.vec[j] + 1
  }
}  
gSRDF$collin.addressed <- collin.and.shrink.vec/collin.vec

##### Summaries
##### ggplot function
gtwr.coeff.plot <- function(gwr.pred = shape.input, 
	var = "collin.count", 
	name = "NDVI", 
	sh = "Reds") {
 	MyGrad <- scale_fill_distiller(type="seq", direction = 1,
	#	palette = sh, name = "", limits = c(0,1)) 
		palette = sh, name = "") 
	gwr.pred$id <- rownames( gwr.pred@data)
 	gwr.pred.points <- fortify(gwr.pred, region="id")
	gwr.pred.df = full_join(gwr.pred.points, gwr.pred@data, by="id")		
	p.i <- ggplot(gwr.pred.df) + 
			aes(long,lat,group=group,fill=gwr.pred.df[,var]) + 
			geom_polygon() +
			#coord_equal() + 
			MyGrad + 
			labs(title = name) +
			theme(axis.title.x=element_blank(),
		        axis.text.x=element_blank(),
		        axis.ticks.x=element_blank(), 
		        axis.title.y=element_blank(),
		        axis.text.y=element_blank(),
		        #panel.background = element_blank(),
		        axis.ticks.y=element_blank(), 
		        legend.position = "bottom",
				legend.key.width=unit(2,"cm"))
	return(p.i)
}
slot(gSRDF, "polygons") <- lapply(slot(gSRDF, "polygons"), checkPolygonsHoles)

tab <- rbind(summary(cn.global.vec)[1:6], 
				summary(bw.sims.vec/nrow(d.a))[1:6],
				summary(gSRDF$collin.count)[1:6],
				summary(gSRDF$collin.addressed)[1:6])
rownames(tab) <- c("Global Collinearity (100 sims)",
					"GW-ELN Flexible Bandwidth",
					"Counts of local CN > 30 (for 159 counties, 100 sims)",
					"Proportion of local collinearity dealt by shrinkage")
round(tab, 3)
write.csv(tab, "Tab1.csv")

# for each x how many were shrunk at each location? (1)
p1 <- gtwr.coeff.plot(gSRDF, "shrink.prop", "At least one covariate shrunk to 0", "BuPu")
# for each x, which were shrink at each location? (3)
p2 <- gtwr.coeff.plot(gSRDF, "x1.prop", 
	expression(paste("Proportion of ",beta[1], " shrunk")),"BuPu")
p3 <- gtwr.coeff.plot(gSRDF, "x2.prop", 
	expression(paste("Proportion of ",beta[2], " shrunk")),"BuPu")
p4<- gtwr.coeff.plot(gSRDF, "x3.prop", 
	expression(paste("Proportion of ",beta[3], " shrunk")),"BuPu")

#quartz(w = 12, h = 4)
png(filename = "F4.png", w = 9, h = 12, units = "in", res = 150)
multiplot2(list(p1,p3,p2,p4), cols = 2)
dev.off()

# 2. dealling with local collin maps 
# how many times was CN identified at each location (1)
p1 <- gtwr.coeff.plot(gSRDF, "collin.count","Counts of local CN > 30", "YlOrRd")
# of those times how manywhat proprtion were dealt with by variable shrinkage (1)
p2 <- gtwr.coeff.plot(gSRDF, "collin.addressed","Proportion of local collinearity dealt by shrinkage", "YlOrRd")
png(filename = "F5.png", w = 9, h = 6, units = "in", res = 150)
multiplot2(list(p1,p2), cols = 2)
dev.off()

###### END Simulations
