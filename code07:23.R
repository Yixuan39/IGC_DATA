# Sensitivity test
# Forcing tau from 0.1 to 0.6
library(gridExtra)

Sensitivity <- ListFunction("TauValues", "Matrix")
IGC_estimate <-as.vector(IGC1_Full_noforce_noclock$MG211)
Sensitivity <- cbind(Sensitivity, IGC_estimate)
save(Sensitivity, file = "sensitivity.Rdata")


pdf(file = "sensitivityTest.pdf", height = 30, width = 12)
grid.table(round(Sensitivity,4))
dev.off()

# branch specific tau v.s. overall tau
# 26:41 is row of branch specific tau, 9 is row of overall tau
a <- c()
b <- c()
for(i in IGC1_Full_noforce_noclock){
  a <- append(a, i[9])
  b <- append(b, i[26:41])
}
b <- matrix(b, nrow = 16)
rownames(b) <- names(IGC1_Full_noforce_noclock[[1]][26:41])
names(a) <- names(IGC1_Full_noforce_noclock)
colnames(b) <- names(IGC1_Full_noforce_noclock)

pdf("branch.specific.vs.overall.pdf",height=6,width=8)
for(j in 1:16){
  n <- names(IGC1_Full_noforce_noclock[[1]][26:41])[j]
  plot(a, b[j,], type = "n", xlab = "overallTau", ylab = n)
  text(a,b[j,], substr(names(b[j,]),3,length(names(b[j,]))))
  abline(a=0, b=1)
}
dev.off()