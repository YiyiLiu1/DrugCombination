RUVra=function(z,x,k){ ##z gene expression n samples by J genes, x covariates of interest n samples by p covariates, k number of factors (can be a vector to save computation compared with doing them one by one)
	r=z-x%*%solve(t(x)%*%x)%*%t(x)%*%z
	rsvd=svd(r,max(k),0)
	rl=list()
	for (j in 1:length(k)){
		w=rsvd$u[,1:k[j],drop=F]
		rl[[j]]=z-w%*%solve(t(w)%*%w)%*%t(w)%*%z
	}
	return(rl)
}

