#####################################################################
## structural similarity
eva_ss=evaluation_metric(-ss,rs,er)

#####################################################################
## direction-based expression similarity
pcidire=rsdire=ppcidire=prsdire=matrix(0,5,6)
fcthresh=log(c(1,1.1,1.2,1.3,1.4))/log(2)
direm0=list()
edirem0=list()
for (l in 1:length(fcthresh)){
	direm0[[l]]=list()
	for (t in 1:6){
		direm0[[l]][[t]]=matrix(0,14,14)
		for (i in 1:13){
			for (j in (i+1):14){
				direm0[[l]][[t]][i,j]=direm0[[l]][[t]][j,i]=direm(X[[t]][,i],X[[t]][,j],fcthresh[l])
			}
		}
	}
	edirem0[[l]]=list()
	for (t in 1:6){
		d=-direm0[[l]][[t]][lower.tri(direm0[[l]][[t]],diag=F)]
		edirem0[[l]][[t]]=evaluation_metric(d,rs,er)
	}
}
for (t in 1:6){
	for (l in 1:length(fcthresh)){
		pcidire[l,t]=edirem0[[l]][[t]]$pci
		rsdire[l,t]=edirem0[[l]][[t]]$rs
		ppcidire[l,t]=edirem0[[l]][[t]]$pci_p
		prsdire[l,t]=edirem0[[l]][[t]]$rs_p
	}
}

colnames(pcidire)=colnames(rsdire)=colnames(ppcidire)=colnames(prsdire)=c("6h, IC20 (48h)","6h, IC20 (24h)","12h, IC20 (48h)","12h, IC20 (24h)","24h, IC20 (48h)","24h, IC20 (24h)")
rownames(pcidire)=rownames(rsdire)=rownames(ppcidire)=rownames(prsdire)=c("All","FC<1/1.1 or >1.1","FC<1/1.2 or >1.2","FC<1/1.3 or >1.3","FC<1/1.4 or >1.4")


#####################################################################
## GSEA-based expression similarity
p=q=c(500,1000,2000,5000)
Xs=list()
for (t in 1:6){
	Xs[[t]]=matrix(0,nrow(X[[1]]),14)
	for (i in 1:14){
		Xs[[t]][,i]=order(X[[t]][,i],decreasing=T)
	}
}

prl0=matrix(0,nrow(X[[1]]),14)
for (i in 1:14){
	Xi=Xs[[1]][,i]
	for (t in 2:6) Xi=cbind(Xi,Xs[[t]][,i])
	prl0[,i]=prl(Xi)
}
dips0=edips0=list()
for (l in 1:length(q)){
	dips0[[l]]=matrix(0,14,14)
	for (i in 1:13){
		for (j in (i+1):14){
			dips0[[l]][i,j]=dips0[[l]][j,i]=dips(prl0[,i],prl0[,j],p[l],q[l])
		}
	}
	d=-dips0[[l]][lower.tri(dips0[[l]],diag=F)]
	edips0[[l]]=evaluation_metric(d,rs,er)
}

dips_sep0=edips_sep0=list()
for (l in 1:length(q)){
	dips_sep0[[l]]=edips_sep0[[l]]=list()
	for (t in 1:6){
		dips_sep0[[l]][[t]]=matrix(0,14,14)
		for (i in 1:13){
			for (j in (i+1):14){
				dips_sep0[[l]][[t]][i,j]=dips_sep0[[l]][[t]][j,i]=dips(Xs[[t]][,i],Xs[[t]][,j],p[l],q[l])
			}
		}
		d=-dips_sep0[[l]][[t]][lower.tri(dips_sep0[[l]][[t]],diag=F)]
		edips_sep0[[l]][[t]]=evaluation_metric(d,rs,er)
	}
}

pcidips=rsdips=ppcidips=prsdips=matrix(0,4,7)
for (l in 1:length(q)){
	pcidips[l,1]=edips0[[l]]$pci
	rsdips[l,1]=edips0[[l]]$rs
	ppcidips[l,1]=edips0[[l]]$pci_p
	prsdips[l,1]=edips0[[l]]$rs_p
}
for (t in 1:6){
	for (l in 1:length(q)){
		pcidips[l,t+1]=edips_sep0[[l]][[t]]$pci
		rsdips[l,t+1]=edips_sep0[[l]][[t]]$rs
		ppcidips[l,t+1]=edips_sep0[[l]][[t]]$pci_p
		prsdips[l,t+1]=edips_sep0[[l]][[t]]$rs_p
	}
}

colnames(pcidips)=colnames(rsdips)=colnames(ppcidips)=colnames(prsdips)=c("Merged","6h, IC20 (48h)","6h, IC20 (24h)","12h, IC20 (48h)","12h, IC20 (24h)","24h, IC20 (48h)","24h, IC20 (24h)")
rownames(pcidips)=rownames(rsdips)=rownames(ppcidips)=rownames(prsdips)=c("p=q=500","p=q=1000","p=q=2000","p=q=5000")


#####################################################################
## correlation-based expression similarity
rk_cut=c(1000,2000,5000,10000)
dcorpa=dcorsa=list()
for (t in 1:6){
	dcorpa[[t]]=dcorsa[[t]]=matrix(0,14,14)
	for (i in 1:13){
		for (j in (i+1):14){
			dcorpa[[t]][i,j]=dcorpa[[t]][j,i]=dcorp(X[[t]][,i],X[[t]][,j],1:nrow(X[[t]]))
			dcorsa[[t]][i,j]=dcorsa[[t]][j,i]=dcors(X[[t]][,i],X[[t]][,j],1:nrow(X[[t]]))
		}
	}
}
edcorpa=edcorsa=list()
for (t in 1:6){
	d=-dcorpa[[t]][lower.tri(dcorpa[[t]],diag=F)]
	edcorpa[[t]]=evaluation_metric(d,rs,er)
	d=-dcorsa[[t]][lower.tri(dcorsa[[t]],diag=F)]
	edcorsa[[t]]=evaluation_metric(d,rs,er)
}

dcorp1=dcors1=list()
edcorp1=edcors1=list()
Xas=list()
for (t in 1:6){
	Xas[[t]]=matrix(0,nrow(X[[1]]),14)
	for (i in 1:14){
		Xas[[t]][,i]=order(abs(X[[t]][,i]),decreasing=T)
	}
}

for (l in 1:length(rk_cut)){
	dcorp1[[l]]=dcors1[[l]]=list()
	for (t in 1:6){
		dcorp1[[l]][[t]]=dcors1[[l]][[t]]=matrix(0,14,14)
		for (i in 1:13){
			for (j in (i+1):14){
				dcorp1[[l]][[t]][i,j]=dcorp1[[l]][[t]][j,i]=dcorp(X[[t]][,i],X[[t]][,j],unique(c(Xas[[t]][1:rk_cut[l],i],Xas[[t]][1:rk_cut[l],j])))
				dcors1[[l]][[t]][i,j]=dcors1[[l]][[t]][j,i]=dcors(X[[t]][,i],X[[t]][,j],unique(c(Xas[[t]][1:rk_cut[l],i],Xas[[t]][1:rk_cut[l],j])))
			}
		}
	}
	edcorp1[[l]]=edcors1[[l]]=list()
	for (t in 1:6){
		d=-dcorp1[[l]][[t]][lower.tri(dcorp1[[l]][[t]],diag=F)]
		edcorp1[[l]][[t]]=evaluation_metric(d,rs,er)
		d=-dcors1[[l]][[t]][lower.tri(dcors1[[l]][[t]],diag=F)]
		edcors1[[l]][[t]]=evaluation_metric(d,rs,er)
	}
}

pcidcorp=rsdcorp=ppcidcorp=prsdcorp=matrix(0,5,6)
for (t in 1:6){
	pcidcorp[5,t]=edcorpa[[t]]$pci
	rsdcorp[5,t]=edcorpa[[t]]$rs
	ppcidcorp[5,t]=edcorpa[[t]]$pci_p
	prsdcorp[5,t]=edcorpa[[t]]$rs_p
	for (l in 1:length(rk_cut)){
		pcidcorp[l,t]=edcorp1[[l]][[t]]$pci
		rsdcorp[l,t]=edcorp1[[l]][[t]]$rs
		ppcidcorp[l,t]=edcorp1[[l]][[t]]$pci_p
		prsdcorp[l,t]=edcorp1[[l]][[t]]$rs_p
	}
}

colnames(pcidcorp)=colnames(rsdcorp)=colnames(ppcidcorp)=colnames(prsdcorp)=c("6h, IC20 (48h)","6h, IC20 (24h)","12h, IC20 (48h)","12h, IC20 (24h)","24h, IC20 (48h)","24h, IC20 (24h)")
rownames(pcidcorp)=rownames(rsdcorp)=rownames(ppcidcorp)=rownames(prsdcorp)=c("Top 1000","Top 2000","Top 5000","Top 10000","All")


pcidcors=rsdcors=ppcidcors=prsdcors=matrix(0,5,6)
for (t in 1:6){
	pcidcors[5,t]=edcorsa[[t]]$pci
	rsdcors[5,t]=edcorsa[[t]]$rs
	ppcidcors[5,t]=edcorsa[[t]]$pci_p
	prsdcors[5,t]=edcorsa[[t]]$rs_p
	for (l in 1:length(rk_cut)){
		pcidcors[l,t]=edcors1[[l]][[t]]$pci
		rsdcors[l,t]=edcors1[[l]][[t]]$rs
		ppcidcors[l,t]=edcors1[[l]][[t]]$pci_p
		prsdcors[l,t]=edcors1[[l]][[t]]$rs_p
	}
}

colnames(pcidcors)=colnames(rsdcors)=colnames(ppcidcors)=colnames(prsdcors)=c("6h, IC20 (48h)","6h, IC20 (24h)","12h, IC20 (48h)","12h, IC20 (24h)","24h, IC20 (48h)","24h, IC20 (24h)")
rownames(pcidcors)=rownames(rsdcors)=rownames(ppcidcors)=rownames(prsdcors)=c("Top 1000","Top 2000","Top 5000","Top 10000","All")
