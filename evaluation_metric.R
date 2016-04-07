###################################################################
#### PC-index and resampled Spearman correlation analysis
## PC-index
pcindex=function(p,r,er){   # p relates to predicted synergistic score, the larger p[i] is the more synergistic pair i is, # r is the experimentally measured synergistic score, the larger r[i] is the more synergistic pair i is, # er is the standard error
	le=length(p)
	cp1=c()
	cp2=c()
	for (i in 1:(le-1)){
		cp1=c(cp1,rep(i,le-i))
		cp2=c(cp2,(i+1):le)
	}
	sp=1/2+(pnorm((r[cp1]-r[cp2])/sqrt(er[cp1]^2+er[cp2]^2),0,1/sqrt(2))-0.5)*sign(p[cp1]-p[cp2])
	return(sum(sp)*2/(le*(le-1)))
}

## Resampled spearman correlation
response=read.table("NCI-DREAMSubchallenge2_GoldStandard.txt",header=T,sep="\t")
response=response[order(response[,1],response[,2]),]
r=-response[,3]
er=response[,4]
newsample=matrix(0,91,10000)
for (i in 1:10000) newsample[,i]=rnorm(91,r,er)
save(newsample,file="newsample.RData")

load("newsample.RData")
rscorr=function(p){
	return(mean(apply(newsample,2,cor,y=p,method="spearman")))
}

## null distributions
randompre=replicate(10000,sample(91,91))
pcinull=apply(randompre,2,pcindex,r=r,er=er)
rsnull=rep(0,10000)
for (i in 1:10000){
	rsnull[i]=rscorr(randompre[,i])
}
save(pcinull,rsnull,file="nullscore.RData")

## pvalue
load("nullscore.RData")
pcipv=function(pci){
	return(mean(pci<=pcinull))
}
rspv=function(rs){
	return(mean(rs<=rsnull))
}

evaluation_metric=function(p,r,er){
	pci0=pcindex(p,r,er)
	rs0=rscorr(p)
	pcipv0=pcipv(pci0)
	rspv0=rspv(rs0)
	return(list(pci=pci0,rs=rs0,pci_p=pcipv0,rs_p=rspv0))
}


###################################################################
#### AUC analysis
roc=function(prob,cla){  ## cla for true class, prob for the predicted probability of belonging to class 1
	up=sort(unique(prob),decreasing=T)
	le=length(up)
	tp=fp=tn=fn=rep(0,le)
	n=length(cla)
	n1=sum(cla==1)
	j=which(prob==up[1])
	k=sum(cla[j]==1)
	tp[1]=k
	fp[1]=length(j)-k
	fn[1]=n1-tp[1]
	tn[1]=n-n1-fp[1]
	if (le>1){
	for (i in 2:le){
		j=which(prob==up[i])
		k=sum(cla[j]==1)
		tp[i]=tp[i-1]+k
		fp[i]=fp[i-1]+length(j)-k
		fn[i]=n1-tp[i]
		tn[i]=n-n1-fp[i]
	}
	}
	tpr=tp/(tp+fn)
	fpr=fp/(fp+tn)
	return(list(tpr=tpr,fpr=fpr))
}

auc=function(tpr,fpr){
	auc=0
	le=length(tpr)
	for (i in 2:le) auc=auc+(fpr[i]-fpr[i-1])*(tpr[i]+tpr[i-1])/2
	return(auc)
}

n=100
lr=length(r)
cla=rep(0,lr)
cla[r-er*2>0]=1
nfold=3

cvFolds=function(Y,V){
	Y0=split(sample(which(Y==0)),rep(1:V,length=sum(Y==0)))
	Y1=split(sample(which(Y==1)),rep(1:V,length=sum(Y==1)))
	folds=vector("list",length=V)
	for (v in 1:V) folds[[v]]=c(Y0[[v]],Y1[[v]])
	return(folds)
}

cvdata=vector("list",length=n)
for (i in 1:n){
	cvdata[[i]]=cvFolds(cla,nfold)
}

save(cvdata,file="auc_cvdata.RData")

nperm=1000
perm=replicate(nperm,sample(lr,lr))
save(perm,file="auc_permutation.RData")

load("auc_permutation.RData")
aucdist=function(d,ss,cla){
	data=data.frame(cla=cla,d=scale(d),ss=scale(ss))

	aucd=aucdss=rep(0,n)
	aucdp=aucssp=matrix(0,nperm,n)
	mfpr=seq(0,1,length.out=100)

	for (i in 1:n){
		print(i)
		mtprd=mtprdss=rep(0,100)
		for (v in 1:nfold){
			model=glm(cla~d,data=data[-cvdata[[i]][[v]],],family=binomial(link="logit"))
			ypred=predict(model,newdata=data[cvdata[[i]][[v]],],type="response")
			rocd=roc(ypred,cla[cvdata[[i]][[v]]])
			mtprd=mtprd+approx(x=c(0,rocd$fpr),y=c(0,rocd$tpr),xout=mfpr)$y

			model=glm(cla~d+ss,data=data[-cvdata[[i]][[v]],],family=binomial(link="logit"))
			ypred=predict(model,newdata=data[cvdata[[i]][[v]],],type="response")
			rocdss=roc(ypred,cla[cvdata[[i]][[v]]])
			mtprdss=mtprdss+approx(x=c(0,rocdss$fpr),y=c(0,rocdss$tpr),xout=mfpr)$y
		}
		mtprd=mtprd/nfold
		mtprdss=mtprdss/nfold
		aucd[i]=auc(mtprd,mfpr)
		aucdss[i]=auc(mtprdss,mfpr)

		for (j in 1:nperm){
			data1=data2=data
			data1$d=data$d[perm[,j]]
			data2$ss=data$ss[perm[,j]]
			mtprdp=mtprssp=rep(0,100)
			for (v in 1:nfold){
				model=glm(cla~d+ss,data=data1[-cvdata[[i]][[v]],],family=binomial(link="logit"))
				ypred=predict(model,newdata=data1[cvdata[[i]][[v]],],type="response")
				rocdp=roc(ypred,1-cla[cvdata[[i]][[v]]])
				mtprdp=mtprdp+approx(x=c(0,rocdp$fpr),y=c(0,rocdp$tpr),xout=mfpr)$y
				model=glm(cla~d+ss,data=data2[-cvdata[[i]][[v]],],family=binomial(link="logit"))
				ypred=predict(model,newdata=data2[cvdata[[i]][[v]],],type="response")
				rocssp=roc(ypred,1-cla[cvdata[[i]][[v]]])
				mtprssp=mtprssp+approx(x=c(0,rocssp$fpr),y=c(0,rocssp$tpr),xout=mfpr)$y
			}
			mtprdp=mtprdp/nfold
			mtprssp=mtprssp/nfold
			aucdp[j,i]=auc(mtprdp,mfpr)
			aucssp[j,i]=auc(mtprssp,mfpr)
		}
	}
	return(list(aucd=aucd,aucdss=aucdss,maucdp=rowMeans(aucdp),maucssp=rowMeans(aucssp)))
}


