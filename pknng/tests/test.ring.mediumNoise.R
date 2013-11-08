library(clustmp)
data(Three.Rings.MediumNoise.Rdata)

d.iso.a<-pknng(data[,-3],3,diss=F,mu="mean",method="euclidean",conn="one",penalize=1)
HCmethod<-"average"

a<-hclust(as.dist(d.iso.a),method = HCmethod)
label.a<-cutree_c(a,k=5)
plot(data[,1],data[,2],col=data[,3],pch=16,cex=0.5)
title("True labels")

plot(data[,1],data[,2],col=label.a,pch=16,cex=0.5)
title("Clustering Result")
