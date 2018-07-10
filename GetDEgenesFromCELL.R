library(limma)
library(affy)

# read .CEL files
Data <- ReadAffy()

# perform RMA normalization (log2)
eset <- rma(Data)
rma_eset = exprs(eset)

pData(eset)
strain <- c("lrp+","lrp+","lrp+","lrp+","lrp+","lrp-","lrp-","lrp-","lrp-")
design <- model.matrix(~factor(strain))
colnames(design) <- c("lrp+","lrp-vs+")
design
fit <- lmFit(eset, design)
fit <- eBayes(fit)
options(digits=2)
topTable(fit, coef=2, n=40, adjust="BH")


res<-topTable(fit, number=Inf, adjust.method="none", coef=1)
write.table(res,"oct4_dif_exp.txt",sep="\t")
