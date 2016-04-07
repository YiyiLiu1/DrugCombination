cla=rep(0,length(class))
cla[class=="Synergistic"]=1

#####################################################################
## structural similarity
datass=data.frame(cla=cla,ss=scale(ss))
aucss=rep(0,n)
mfpr=seq(0,1,length.out=100)
for (i in 1:n){
	mtprss=rep(0,100)
	for (v in 1:nfold){
		model=glm(cla~ss,data=datass[-cvdata[[i]][[v]],],family=binomial(link="logit"))
		ypred=predict(model,newdata=datass[cvdata[[i]][[v]],],type="response")
		rocss=roc(ypred,1-cla[cvdata[[i]][[v]]])
		mtprss=mtprss+approx(x=c(0,rocss$fpr),y=c(0,rocss$tpr),xout=mfpr)$y
	}
	mtprss=mtprss/nfold
	aucss[i]=auc(mtprss,mfpr)
}
aucssmean=mean(aucss)


#####################################################################
## direction-based expression similarity
aucdire=list()
for (l in 1:length(fcthresh)){
	aucdire[[l]]=list()
	for (t in 1:6){
		d=-direm0[[l]][[t]][lower.tri(direm0[[l]][[t]],diag=F)]
		aucdire[[l]][[t]]=aucdist(d,ss,cla)
	}
}

aucdiremean=aucdiressmean=paucdire=paucdiress=matrix(0,length(fcthresh),6)
for (l in 1:length(fcthresh)){
	for (t in 1:6){
		aucdiremean[l,t]=mean(aucdire[[l]][[t]]$aucd)
		aucdiressmean[l,t]=mean(aucdire[[l]][[t]]$aucdss)
		paucdire[l,t]=mean(aucdire[[l]][[t]]$maucssp > aucdiressmean[l,t])
		paucdiress[l,t]=mean(aucdire[[l]][[t]]$maucdp > aucdiressmean[l,t])
	}
}

colnames(aucdiremean)=colnames(aucdiressmean)=colnames(paucdire)=colnames(paucdiress)=c("6h, IC20 (48h)","6h, IC20 (24h)","12h, IC20 (48h)","12h, IC20 (24h)","24h, IC20 (48h)","24h, IC20 (24h)")
rownames(aucdiremean)=rownames(aucdiressmean)=rownames(paucdire)=rownames(paucdiress)=c("All","FC<1/1.1 or >1.1","FC<1/1.2 or >1.2","FC<1/1.3 or >1.3","FC<1/1.4 or >1.4")


#####################################################################
## GSEA-based expression similarity
aucdips=list()
for (l in 1:length(q)){
	aucdips[[l]]=list()
	d=-dips0[[l]][lower.tri(dips0[[l]],diag=F)],diag=F)]
	aucdips[[l]][[1]]=aucdist(d,ss,cla)
	for (t in 1:6){
		d=-dips_sep0[[l]][[t]][lower.tri(dips_sep0[[l]][[t]],diag=F)]
		aucdips[[l]][[t+1]]=aucdist(d,ss,cla)
	}
}

aucdipsmean=aucdipsssmean=paucdips=paucdipsss=matrix(0,length(q),7)
for (l in 1:length(q)){
	for (t in 1:7){
		aucdipsmean[l,t]=mean(aucdips[[l]][[t]]$aucd)
		aucdipsssmean[l,t]=mean(aucdips[[l]][[t]]$aucdss)
		paucdips[l,t]=mean(aucdips[[l]][[t]]$maucssp > aucdipsssmean[l,t])
		paucdipsss[l,t]=mean(aucdips[[l]][[t]]$maucdp > aucdipsssmean[l,t])
	}
}

colnames(aucdipsmean)=colnames(aucdipsssmean)=colnames(paucdips)=colnames(paucdipsss)=c("Merged","6h, IC20 (48h)","6h, IC20 (24h)","12h, IC20 (48h)","12h, IC20 (24h)","24h, IC20 (48h)","24h, IC20 (24h)")
rownames(aucdipsmean)=rownames(aucdipsssmean)=rownames(paucdips)=rownames(paucdipsss)=c("p=q=500","p=q=1000","p=q=2000","p=q=5000")


#####################################################################
## Pearson correlation-based expression similarity
aucdcorp=list()
for (l in 1:length(rk_cut)){
	aucdcorp[[l]]=list()
	for (t in 1:6){
		d=-dcorp1[[l]][[t]][lower.tri(dcorp1[[l]][[t]],diag=F)]
		aucdcorp[[l]][[t]]=aucdist(d,ss,cla)
	}
}
aucdcorp[[length(rk_cut)+1]]=list()
for (t in 1:6){
	d=-dcorpa[[t]][lower.tri(dcorpa[[t]],diag=F)]
	aucdcorp[[length(rk_cut)+1]][[t]]=aucdist(d,ss,cla)
}

aucdcorpmean=aucdcorpssmean=paucdcorp=paucdcorpss=matrix(0,length(rk_cut)+1,6)
for (l in 1:(length(rk_cut)+1)){
	for (t in 1:6){
		aucdcorpmean[l,t]=mean(aucdcorp[[l]][[t]]$aucd)
		aucdcorpssmean[l,t]=mean(aucdcorp[[l]][[t]]$aucdss)
		paucdcorp[l,t]=mean(aucdcorp[[l]][[t]]$maucssp > aucdcorpssmean[l,t])
		paucdcorpss[l,t]=mean(aucdcorp[[l]][[t]]$maucdp > aucdcorpssmean[l,t])
	}
}

colnames(aucdcorpmean)=colnames(aucdcorpssmean)=colnames(paucdcorp)=colnames(paucdcorpss)=c("6h, IC20 (48h)","6h, IC20 (24h)","12h, IC20 (48h)","12h, IC20 (24h)","24h, IC20 (48h)","24h, IC20 (24h)")
rownames(aucdcorpmean)=rownames(aucdcorpssmean)=rownames(paucdcorp)=rownames(paucdcorpss)=c("Top 1000","Top 2000","Top 5000","Top 10000","All")


#####################################################################
## Spearman correlation-based expression similarity
aucdcors=list()
for (l in 1:length(rk_cut)){
	aucdcors[[l]]=list()
	for (t in 1:6){
		d=-dcors1[[l]][[t]][lower.tri(dcors1[[l]][[t]],diag=F)]
		aucdcors[[l]][[t]]=aucdist(d,ss,cla)
	}
}
aucdcors[[length(rk_cut)+1]]=list()
for (t in 1:6){
	d=-dcorsa[[t]][lower.tri(dcorsa[[t]],diag=F)]
	aucdcors[[length(rk_cut)+1]][[t]]=aucdist(d,ss,cla)
}


aucdcorsmean=aucdcorsssmean=paucdcors=paucdcorsss=matrix(0,length(rk_cut)+1,6)
for (l in 1:(length(rk_cut)+1)){
	for (t in 1:6){
		aucdcorsmean[l,t]=mean(aucdcors[[l]][[t]]$aucd)
		aucdcorsssmean[l,t]=mean(aucdcors[[l]][[t]]$aucdss)
		paucdcors[l,t]=mean(aucdcors[[l]][[t]]$maucssp > aucdcorsssmean[l,t])
		paucdcorsss[l,t]=mean(aucdcors[[l]][[t]]$maucdp > aucdcorsssmean[l,t])
	}
}

colnames(aucdcorsmean)=colnames(aucdcorsssmean)=colnames(paucdcors)=colnames(paucdcorsss)=c("6h, IC20 (48h)","6h, IC20 (24h)","12h, IC20 (48h)","12h, IC20 (24h)","24h, IC20 (48h)","24h, IC20 (24h)")
rownames(aucdcorsmean)=rownames(aucdcorsssmean)=rownames(paucdcors)=rownames(paucdcorsss)=c("Top 1000","Top 2000","Top 5000","Top 10000","All")
