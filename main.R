#####################################################################
## read gene expression data
data=read.table("ACalifano_DLBCL_Ly3_14Comp_treatment.txt",skip=4,header=F,sep="\t")
es=as.matrix(data[,3:ncol(data)])
hname=read.table("ACalifano_DLBCL_Ly3_14Comp_treatment.txt",nrow=4,header=F,sep="\t")
chipname=c()
for (j in 3:ncol(hname)){
	chipname=c(chipname,paste(hname[2,j],hname[3,j],hname[4,j],sep="&"))
}

#####################################################################
## read experimentally measured EOB
response=read.table("NCI-DREAMSubchallenge2_GoldStandard.txt",header=T,sep="\t")
response=response[order(response[,1],response[,2]),]
rs=-response[,3]
er=response[,4]
class=rep("Additive",length(rs))
i=which(rs-er*2>0)
class[i]="Synergistic"
i=which(rs+er*2<0)
class[i]="Antagonistic"

#####################################################################
## preprocess gene expression data with RUV
source("datapreprocess.R")
k=26
X=foldchange(t(esruv[[k]])) ## expression data after removing batch effects (k=26)

#####################################################################
## read structure similarity data
ssim=read.table("pubchem_2d_sorted.txt",header=T,sep="\t")
ss=ssim[,3]/100

#####################################################################
## PC-index and resampled Spearman correlation analysis
source("distance.R")
source("evaluation_metric.R")
source("evaluation1.R")

# Negative structural similarity PC-index
print(eva_ss$pci)

# Negative structural similarity PC-index p-value
print(eva_ss$pci_p)

# Negative structural similarity resampled Spearman correlation
print(eva_ss$rs)

# Negative structural similarity resampled Spearman correlation p-value
print(eva_ss$rs_p)

# Direction-based gene expression similarity PC-index
print(pcidire)

# Direction-based gene expression similarity PC-index p-value
print(ppcidire)

# Direction-based gene expression similarity resampled Spearman correlation
print(rsdire)

# Direction-based gene expression similarity resampled Spearman correlation p-value
print(prsdire)

# GSEA-based gene expression similarity PC-index
print(pcidips)

# GSEA-based gene expression similarity PC-index p-value
print(ppcidips)

# GSEA-based gene expression similarity resampled Spearman correlation
print(rsdips)

# GSEA-based gene expression similarity resampled Spearman correlation p-value
print(prsdips)

# Pearson correlation-based gene expression similarity PC-index
print(pcidcorp)

# Pearson correlation-based gene expression similarity PC-index p-value
print(ppcidcorp)

# Pearson correlation-based gene expression similarity resampled Spearman correlation
print(rsdcorp)

# Pearson correlation-based gene expression similarity resampled Spearman correlation p-value
print(prsdcorp)

# Spearman correlation-based gene expression similarity PC-index
print(pcidcors)

# Spearman correlation-based gene expression similarity PC-index p-value
print(ppcidcors)

# Spearman correlation-based gene expression similarity resampled Spearman correlation
print(rsdcors)

# Spearman correlation-based gene expression similarity resampled Spearman correlation p-value
print(prsdcors)


#####################################################################
## AUC analysis
source("evaluation2.R")

# Only structural similarity
print(aucssmean)

# Only direction-based expression similarity
print(aucdiremean)

# Direction-based expression similarity + structural similarity
print(aucdiressmean)

# P-value, improvement by incorporating structural similarity to direction-based expression similarity
print(paucdiress)

# P-value, improvement by incorporating direction-based expression similarity to structural similarity
print(paucdire)

# Only GSEA-based expression similarity
print(aucdipsmean)

# GSEA-based expression similarity + structural similarity
print(aucdipsssmean)

# P-value, improvement by incorporating structural similarity to GSEA-based expression similarity
print(paucdipsss)

# P-value, improvement by incorporating GSEA-based expression similarity to structural similarity
print(paucdips)

# Only Pearson correlation-based expression similarity
print(aucdcorpmean)

# Pearson correlation-based expression similarity + structural similarity
print(aucdcorpssmean)

# P-value, improvement by incorporating structural similarity to Pearson correlation-based expression similarity
print(paucdcorpss)

# P-value, improvement by incorporating Pearson correlation-based expression similarity to structural similarity
print(paucdcorp)

# Only Pearson correlation-based expression similarity
print(aucdcorsmean)

# Pearson correlation-based expression similarity + structural similarity
print(aucdcorsssmean)

# P-value, improvement by incorporating structural similarity to Pearson correlation-based expression similarity
print(paucdcorsss)

# P-value, improvement by incorporating Pearson correlation-based expression similarity to structural similarity
print(paucdcors)
