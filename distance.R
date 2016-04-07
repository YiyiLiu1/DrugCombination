#####################################################
## DIPS distance
## Merging lists for a drug at different concentrations and time points
# Spearman's Footrule
sdelta=function(p1,p2){#input ranked list of genes
	r1=p1
	r2=p2
	r1[p1]=1:length(p1)
	r2[p2]=1:length(p2) #genes' positions in the ranked list
	sum(abs(r1-r2))
}

# Borda Merging Function
bmerge=function(p1,p2){#input ranked list of genes
	r1=p1
	r2=p2
	r1[p1]=1:length(p1)
	r2[p2]=1:length(p2)
	return(order(r1+r2,decreasing=F))
}

prl=function(lists){
	n=ncol(lists)
	dv=matrix(0,3,n*(n-1)/2)
	k=1
	for (i in 1:(n-1)){
		for (j in (i+1):n){
			dv[c(1,2),k]=c(i,j)
			dv[3,k]=sdelta(lists[,i],lists[,j])
			k=k+1
		}
	}
	#print(dv)
	while (n>1){
		n=n-1
		i=which.min(dv[3,])
		mi=min(dv[1,i[1]],dv[2,i[1]])
		mj=max(dv[1,i[1]],dv[2,i[1]])	
		lists[,mi]=bmerge(lists[,mi],lists[,mj])
		lists=as.matrix(lists[,-mj])
		if (n>1){
			j=c(which(dv[1,]==mj),which(dv[2,]==mj))
			dv=as.matrix(dv[,-j])
			j=which(dv[1,]>mj)
			if (length(j)>0) dv[1,j]=dv[1,j]-1
			j=which(dv[2,]>mj)
			if (length(j)>0) dv[2,j]=dv[2,j]-1
			j=c(which(dv[1,]==mi),which(dv[2,]==mi))
			dv[1,j]=mi
			dv[2,j]=seq(1,n)[-mi]
			dv0=c()
			for (k in seq(1,n)[-mi]){
				dv0=c(dv0,sdelta(lists[,mi],lists[,k]))
			}
			dv[3,j]=dv0
		}
	}
	return(as.vector(lists[,1]))
}

source("GSEA.1.0.R")
dips=function(prl1,prl2,p,q){
	es12p=GSEA.EnrichmentScore2(prl2,prl1[1:p],0)$ES
	es12q=GSEA.EnrichmentScore2(prl2,tail(prl1,q),0)$ES
	tes12=1/2-(es12p-es12q)/4
	es21p=GSEA.EnrichmentScore2(prl1,prl2[1:p],0)$ES
	es21q=GSEA.EnrichmentScore2(prl1, tail(prl2,q), 0)$ES
	tes21=1/2-(es21p-es21q)/4
	return((tes12+tes21)/2)
}


################################################################
## Direction distance
direm=function(fc1,fc2,fcthresh){
	dfc1=dfc2=rep(0,length(fc1))
	dfc1[fc1> fcthresh]=1
	dfc1[fc1< (-fcthresh)]=-1
	dfc2[fc2> fcthresh]=1
	dfc2[fc2< (-fcthresh)]=-1
	return((1-mean(dfc1*dfc2))/2)
}


################################################################
## Correlation distance
dcorp=function(fc1,fc2,gx){
	return((1-cor(fc1[gx],fc2[gx]))/2)
}
dcors=function(fc1,fc2,gx){
	return((1-cor(fc1[gx],fc2[gx],method="spearman"))/2)
}
