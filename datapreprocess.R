jpeg("rle_0.jpeg")
es1=es-apply(es,1,median)
bp1=boxplot(es1,use.cols=T,col=c(rep("red",94),rep("black",94),rep("blue",94)),pch=19,cex=.1,outcol="grey",border=c(rep("red",94),rep("black",94),rep("blue",94)))
dev.off()

source("RUV.R")
k=1:40
x=matrix(0,ncol(es),19)
x[,1:2]=1
cpd=sort(unique(as.vector(as.matrix(hname[2,3:ncol(hname)]))))
cpd=cpd[-c(grep("DMSO",cpd),grep("Media",cpd))] ##remove DMSO and media
j=which(hname[2,]=="Media")
x[j-2,2]=0
for (i in 1:14){
	j=which(hname[2,]==cpd[i])
	x[j-2,i+2]=1
}
j=which(hname[3,]==12)
x[j-2,17]=1
j=which(hname[3,]==24)
x[j-2,18]=1
j=which(hname[4,]=="IC20")
x[j-2,19]=1
esruv=RUVra(t(es),x,k)

iqrmax=rep(0,40)
for (j in 1:40){
	es1=t(esruv[[j]])-apply(t(esruv[[j]]),1,median)
	jpeg(paste("rle_",k[j],".jpeg",sep=""))
	bp1=boxplot(es1,use.cols=T,col=c(rep("red",94),rep("black",94),rep("blue",94)),pch=19,cex=.1,outcol="grey",border=c(rep("red",	94),rep("black",94),rep("blue",94)))
	iqrmax[j]=max(bp1$stats[4,]-bp1$stats[2,])
	dev.off()
}

pdf("iqrmax.pdf")
plot(k,iqrmax,type="b")
dev.off()


uniname=sort(unique(chipname))
gprobe=as.vector(data[,2])
unigene=unique(gprobe)
t6c1=grep("6&1/10 of IC20",uniname)
t6c2=grep("6&IC20",uniname)
t12c1=grep("12&1/10 of IC20",uniname)
t12c2=grep("12&IC20",uniname)
t24c1=grep("24&1/10 of IC20",uniname)
t24c2=grep("24&IC20",uniname)

foldchange=function(ex){
	aex=matrix(0,nrow(ex),length(uniname))
	for (j in 1:length(uniname)){
		i=which(chipname==uniname[j])
		aex[,j]=rowMeans(ex[,i,drop=F])
	}
	fc=list()
	fc[[1]]=aex[,t6c1]-aex[,27]
	fc[[2]]=aex[,t6c2]-aex[,27]
	fc[[3]]=aex[,t12c1]-aex[,25]
	fc[[4]]=aex[,t12c2]-aex[,25]
	fc[[5]]=aex[,t24c1]-aex[,26]
	fc[[6]]=aex[,t24c2]-aex[,26]
	return(fc)
}