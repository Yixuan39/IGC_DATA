taucountmatrix <- MG[44:57,]+MG[60:73,]
scale0matrix <- taucountmatrix/MG[28:41,]
omegamatrix <- rbind(MG[8,], MG[8,], MG[8,], MG[8,],MG[8,], MG[8,], MG[8,], MG[8,], MG[8,], MG[8,], MG[8,], MG[8,], MG[8,], MG[8,])
rownames(omegamatrix) <- rownames(MG[12:25,])
rownames(scale0matrix) <- rownames(MG[12:25,])
rownames(taucountmatrix) <- rownames(MG[12:25,])

newmatrix=matrix(rep(1,14*45),nrow=14)
for(i in 1:14){
  newmatrix[i,]=scale0matrix[i,]/MG[1,]
}
newmatrix=newmatrix/MG[12:25,]

genenamesmatrix <- rbind(colnames(taucountmatrix),colnames(taucountmatrix),colnames(taucountmatrix),colnames(taucountmatrix), colnames(taucountmatrix),colnames(taucountmatrix),colnames(taucountmatrix),colnames(taucountmatrix), colnames(taucountmatrix),colnames(taucountmatrix),colnames(taucountmatrix),colnames(taucountmatrix), colnames(taucountmatrix),colnames(taucountmatrix) )
branchnamesmatrix <- cbind(rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix), rownames(taucountmatrix))
taucountvector <- c(taucountmatrix)
scale0vector <- c(scale0matrix)
omegavector <- c(omegamatrix)
genenamesvector <- c(genenamesmatrix)
branchnamesvector <- c(branchnamesmatrix)


data=data.frame(IGCevent=taucountvector,scale0=log(scale0vector),omega=omegavector,genename=genenamesvector,
               branchname=branchnamesvector, seqdiff=c(newmatrix ))

library(lme4)
glm(round(IGCevent)~scale0+genename+seqdiff+seqdiff*genename,data=data,family = "poisson")
glm(IGCevent~scale0+branchname+genename,data=data,family = poisson)


