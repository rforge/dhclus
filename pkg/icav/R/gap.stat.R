icav.numeric.tolerance<-1e-20

do.means<-function(data,class){
k<-length(unique(class))
withiness<-vector()


if(sum(class==1)>1) MEAN<-colMeans(drop(data[class==1,]))
else if(sum(class==i)==1) MEAN <-data[class==i,]
else {warning("empty cluster");MEAN<-rep(NA,ncol(data))}
withiness<-withinSum(data,class)

if(k==1) return(list(means = MEAN, withiness = withiness))

for(i in 2:k) {
if(sum(class==i)>1) {MEAN<-rbind(MEAN,colMeans(drop(data[class==i,])))}
else if(sum(class==i)==1) MEAN<-rbind(MEAN, data[class==i,] )
}


return(list(means = MEAN, withiness = withiness))
}

Gap.plot <- function(x, k.max = 20, M = 100,clustering=NULL,mc.cores=8){
  k.max <- trunc(k.max)
  if(k.max < 2) stop("'k.max' has to be >= 2")
  if (is.null(clustering)) stop("NULL clustering method")
  k <- 1
  km.new <- NULL
  gap.new <- gapStat.mc(data = x, class = rep(1, nrow(x)), M = M,clustering=clustering,mc.cores=mc.cores)
  gap.stat <- gap.new
  one.clust<-do.means(x,rep(1, nrow(x)) )
  ONE.TIME<-FALSE
  km.plot<-NULL
  
  repeat{
    km.old <- km.new
    gap.old <- gap.new
    k <- k + 1

    if(k > k.max){
#       warning("'k.max' reached kmeans result for k.max centers returned.")
      break
    }

      km.new <- clustering(x, k)
      if(k==2){clust.obs<-km.new}
      else{clust.obs<-cbind(clust.obs,km.new)}
    gap.new <- gapStat.mc(data = x, class = km.new, M = M,clustering=clustering,mc.cores=mc.cores)
    gap.stat <- rbind(gap.stat, gap.new)
    
    if(gap.old[1] - gap.new[1] - gap.new[2] >= 0){
      if(!ONE.TIME){
	if (!is.null(km.old)){
	    km.plot<-list(cluster= km.old, withinss= withinSum(x,km.old), centers= NULL, size= NULL)
	    class(km.plot) <- "kmeans"
	} 
	ONE.TIME<-TRUE
	}
    }
  }

  rownames(gap.stat) <- NULL
  gap.stat <- cbind(1:nrow(gap.stat), gap.stat)
  colnames(gap.stat) <- c("number of clusters", "Gap statistic", "SE of simulation", "Null Withiness", "Obs. Withiness")


  res <- list(gapStat = gap.stat, clust.obs = clust.obs)

  class(res) <- "clusterGap"
  
  return(res)
}

gapStat.mc <- function (data, class = rep(1, nrow(data)),k.max = 20, M = 100, clustering=NULL,mc.cores=8){
	if (!(length(class) == nrow(data)))
	stop("Length of class vector differs from nrow of data")
	if(M <= 0)
	stop("'M' has to be a positive integer")
	
	if (is.null(clustering)) stop("NULL clustering method in gapStat.mc")
	data <- as.matrix(data)

	ws<-withinSum(data,class)
	ws[ws<icav.numeric.tolerance]<-icav.numeric.tolerance
	temp1 <- log10(sum( ws ))

	data <- scale(data, center = TRUE, scale = F)
	M <- trunc(M)

	veigen <- svd(data)$v
	x1 <- crossprod(t(data), veigen)

	min.x <- apply(x1,2,min)
	max.x <- apply(x1,2,max)
	N<-nrow(x1)
	clust.num<-length(unique(class))
	tots <- vector(length = M)
	
	do.rand.matrix<-function(dim.i,N.data,min.val,max.val){return(runif(N.data,min.val[dim.i],max.val[dim.i]))}
	
	foo.null<-function(itera,NCOL, NData, Min.x, Max.x ,clustering){
# 		for (j in 1:NCOL) {
# 			z1[, j] <- runif(NData, min = Min.x[j], max = Max.x[j])
# 		}
		z1<- matrix( unlist(lapply(1:NCOL, do.rand.matrix, NData, Min.x, Max.x) ) , nc=NCOL, nr=NData )
		if (clust.num>1){
			clust.null<-clustering(z1, clust.num); 
			ws<-withinSum(z1,clust.null)
					ws[ws<icav.numeric.tolerance]<-icav.numeric.tolerance
					tots <- log10( sum(ws) )
		}
		else {
			ws<-withinSum(z1,rep(1, NData))
					ws[ws<icav.numeric.tolerance]<-icav.numeric.tolerance
					tots <- log10( sum(ws) )
		}
	}
	
	tots <- unlist(mclapply(1:M, foo.null , ncol(x1), N, min.x, max.x, clustering, mc.cores=mc.cores))
# 	print(tots)
	out <- c(mean(tots) - temp1, sqrt(1 + 1/M) * sd(tots), mean(tots), temp1)
			names(out) <- c("Gap statistic", "SE of simulation", "Null Withiness", "Obs. Withiness")
			return(out)
}

## plot method for objects of class "clusterGap"

plot.clusterGap <- function(x, ...,title=NULL,ylab="",xlab="",cex=1,cex.lab=1,cex.axis=1){
min.val<-min(x$gapStat[,2]-x$gapStat[,3])
max.val<-max(x$gapStat[,2]+x$gapStat[,3])

  plot(1:nrow(x$gapStat), x$gapStat[,2], ylab = ylab, xlab = xlab, type = "l",ylim=c(min.val,max.val),cex=cex,cex.lab=cex.lab,cex.axis=cex.axis)
  plotrix::plotCI(1:nrow(x$gapStat), x$gapStat[,2], x$gapStat[,3], add = TRUE,yaxt="n",xaxt="n",cex=cex,cex.lab=cex.lab,cex.axis=cex.axis)
  if(is.null(title)) title("Gap statistic for clustering")
  else title(title)
}

plot.GapStat <- function(x, ...){
min.val<-min(x$gapStat[,4],x$gapStat[,5])
max.val<-max(x$gapStat[,4],x$gapStat[,5])

  matplot(1:nrow(x$gapStat), cbind(x$gapStat[,4], x$gapStat[,5] ), ylab = "gap statistic", xlab = "number of clusters", type = "b", col=c(1,3),pch=c(1,3),lty=1 ,ylim=c(min.val,max.val))
  legend(1, (max.val+min.val)/2 , c("Null Withiness", "Obs. Withiness"),col=c(1,3),pch=c(1,3),lty=1,cex=1.25)
  title("Gap statistic for clustering")
}
