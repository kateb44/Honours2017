rosmap= read.csv('rosmap_xsect_slopes_n3118_021916.csv')
rownames(rosmap)= rosmap$projid

Seq = read.delim("RSEM_gene_FPKM_Quantile_Combat_log2_normalized.txt")
rownames(Seq) = Seq$gene_id
Seq = Seq[,-c(1,2)]
rownames(Seq) = unlist(lapply(strsplit(rownames(Seq),split = '\\.'),function(x)x[1]))


X = read.delim("RNA-Seq_projid.txt")
projID = X[,3]
names(projID) = X[,1]

colnames(Seq) = projID[gsub('X','',colnames(Seq))]
Seq = Seq[rowSums(Seq)!=0,]

rnaData = Seq

load('eQTLcisKate.RData')
PrevalFX = do.call('rbind',lapply(eQTL,function(x)x[,'fxPreval']))
OLSFX = do.call('rbind',lapply(eQTL,function(x)x[,'fxOLS']))
LassoFX = do.call('rbind',lapply(eQTL,function(x)x[,'fxLasso']))

library(limma)
design = model.matrix(~pathoAD,rosmap)
u = intersect(colnames(PrevalFX),rownames(design))

fit = lmFit(rnaData[rownames(OLSFX),u],design[u,])
fit = ebayes(fit)
hist(fit$p.value[,'pathoAD'])


fitPV = lmFit(PrevalFX[,u],design[u,])
fitPV = ebayes(fitPV)
hist(fitPV$p.value[,'pathoAD'],xlab="p-values",ylab="Frequency",main="p-values of the TSPL Method")

fitOLS = lmFit(OLSFX[,u],design[u,])
fitOLS = ebayes(fitOLS)
hist(fitOLS$p.value[,'pathoAD'],xlab="p-values",ylab="Frequency",main="p-values of the TSLS Method")

fitLasso = lmFit(LassoFX[,u],design[u,])
fitLasso = ebayes(fitLasso)
hist(fitLasso$p.value[,'pathoAD'],xlab="p-values",ylab="Frequency",main="p-values of the TSL Method")

plot(fitLasso$t[,'pathoAD'],fitOLS$t[,'pathoAD'],xlab="TSL",ylab="TSLS",main="Test Statistics of TSL vs TSLS")
plot(fitPV$t[,'pathoAD'],fitOLS$t[,'pathoAD'],xlab="TSLS",ylab="Frequency",main="Test Statistics of TSPL vs TSLS")
plot(fit$t[,'pathoAD'],fitOLS$t[,'pathoAD'])

GWAS = read.delim('gwas_catalog_v1.0-associations_e84_r2016-03-13.tsv',header = TRUE)

x = GWAS[,c('SNP_GENE_IDS_ENSEMBL',"ENSEMBL_DOWNSTREAM_GENE_ID","ENSEMBL_UPSTREAM_GENE_ID")]
Disease = split(x ,GWAS$DISEASE.TRAIT)
Disease = lapply(Disease,function(x)as.character(unlist(x)))
Disease = lapply(Disease,function(x)unlist(strsplit(as.character(x),', ')))
Disease = Disease[unlist(lapply(Disease,length))>50]

snps = intersect(unlist(Disease[grep('Alz',names(Disease))]),rownames(OLSFX))

plot(fitLasso$t[,'pathoAD'],fitOLS$t[,'pathoAD'],xlab="TSL",ylab="TSLS",main="Test Statistics of TSL vs TSLS")
points(fitLasso$t[snps,'pathoAD'],fitOLS$t[snps,'pathoAD'],col=2,pch = 19)

plot(fitPV$t[,'pathoAD'],fitOLS$t[,'pathoAD'],xlab="TSPL",ylab="TSLS",main="Test Statistics of TSPL vs TSLS")
points(fitPV$t[snps,'pathoAD'],fitOLS$t[snps,'pathoAD'],col=2,pch = 19)

plot(fit$t[,'pathoAD'],fitPV$t[,'pathoAD'])
points(fit$t[snps,'pathoAD'],fitPV$t[snps,'pathoAD'],col=2,pch = 19)

library(goseq)
library(GO.db)
#
BP = getgo(rownames(PrevalFX),genome = 'hg19',id = 'ensGene')
BP2ens = Biobase::reverseSplit(BP)
xx <- as.list(GOTERM)
gomap = unlist(lapply(xx,Term))
names(BP2ens) = gomap[names(BP2ens)]
BPlen = unlist(lapply(BP2ens,length))



pathOR = NULL
pathOR2 = list()
mod2ens = list()
mod2ens[['DE']] = names(which(fit$p.value[,'pathoAD']<0.05))
mod2ens[['OLS']] = names(which(fitOLS$p.value[,'pathoAD']<0.05))
mod2ens[['Lasso']] = names(which(fitLasso$p.value[,'pathoAD']<0.05))
mod2ens[['PV']] = names(which(fitPV$p.value[,'pathoAD']<0.05))



for(j in names(mod2ens)){
  library(limma)
  P = K = Q = G = NULL
  u = intersect(unique(unlist(BP2ens)),unique(unlist(mod2ens)))
  topG = intersect(mod2ens[[j]],u)
  for(i in names(which(BPlen>20&BPlen<500))){
    q = length(intersect(BP2ens[[i]],topG))  
    m = sum(topG%in%u)
    n = sum(!u%in%topG)
    k = sum(u%in%BP2ens[[i]])
    K[i] = k
    Q[i] = q
    #G[i] = paste(geneSymbol[intersect(BP2ens[[i]],topG)],collapse = ', ')
    P[i] = phyper(q,m,n,k,lower.tail = FALSE)+dhyper(q,m,n,k)
  }
  pathOR[j] = paste(names(head(sort(P),6)),collapse = ', ')
  pathOR2[[j]] = data.frame(mod = j,Pathway = names(P),P,K,Q)
}

pathOR


