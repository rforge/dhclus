library(pknng)
data(Three.Rings.HighNoise)

d.iso.a<-pknng(data[,-3],3,diss=FALSE,mu="mean",method="euclidean",conn="one",penalize=1)
HCmethod<-"average"

a<-hclust(as.dist(d.iso.a),method = HCmethod)
label.a<-cutree(a,k=5)
split.screen(c(1,2))
screen(1)
plot(data[,1],data[,2],col=data[,3],pch=16,cex=0.5)
title("True labels")
screen(2)
plot(data[,1],data[,2],col=label.a,pch=16,cex=0.5)
title("Clustering Result")
close.screen(all = TRUE) 