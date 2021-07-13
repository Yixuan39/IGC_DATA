outputs <- IGC2_Full_noforce_noclock
gene.set <- names(outputs)

### edge lists and taus on the esges ###
edge.list.set <- list()

for (gene in gene.set){
	estimate0 <- outputs[[gene]]
	names(estimate0)
	n.par0 <- length(estimate0)
	n.branch0 <- (n.par0 - 9)/5
	branches0 <- names(estimate0)[9+1:n.branch0]
	branches1 <- chartr(")"," ",chartr("("," ",branches0))
	branches2 <- gsub(pattern=" ",replacement="",branches1)
	write.csv(branches2,"temp.csv",quote=F)
	edge.list0 <- read.csv("temp.csv")
	edge.list.set[[gene]] <- edge.list0
}

### defining branches as the subsets of taxa  ####
### unique labeling by sorting the elements   ####
branch.set <- list()
for (gene in gene.set){
  node1 <- edge.list.set[[gene]][,1]
  node2 <- edge.list.set[[gene]][,2]

  branch <- c()
  for (i in 1:length(node2)){
    node.set <- node2[i]
    x <- substr(node.set,1,1)
    while(sum(x=="N")>0){
      node.set0 <- node.set[x=="N"]
      for (element in node.set0){
        node.set1 <- setdiff(node.set,element)
        node.set <- union(node.set1,node2[node1==element])
      } 
        x <- substr(node.set,1,1)
    }
   node.set <- sort(node.set)
   species <- node2[substr(node2,1,1)!="N"]
   node.setc <- sort(setdiff(species,node.set))
   node.set.option <- list(node.set,node.setc)
   short <- order(c(length(node.set),length(node.setc)))[1]
   node.set1 <- node.set.option[[short]]
   branch0 <- NULL
   for (node in node.set1){
     branch0 <- paste(branch0,node,sep="_")
   }
   branch[i] <- branch0
  }
  
  branch.set[[gene]] <- branch
}

##### analysis of genes without missing #####
n.branch_par <- c()
for (gene in gene.set){
   n.branch_par <- c(n.branch_par,length(outputs[[gene]])-9)
}
names(n.branch_par) <- gene.set
table(n.branch_par)

n.branch.max <- max(n.branch_par)
gene.set0 <- gene.set[n.branch_par==n.branch.max]
n.branch0 <- n.branch.max/5

nsites.set0 <- c(); branch.set0 <- c(); bl.set0 <- c();tau_branch.set0 <- c()
igc12_branch.set0 <- c(); igc21_branch.set0 <- c(); mut_branch.set0 <- c()
for (gene in gene.set0){
   res.gene <- outputs[[gene]]
   nsites.set0 <- cbind(nsites.set0,rep(res.gene[1],n.branch0))
   branch.set0 <- cbind(branch.set0,branch.set[[gene]])
   bl.set0 <- cbind(bl.set0,res.gene[9+1:n.branch0])
   tau_branch.set0 <- cbind(tau_branch.set0,
				res.gene[9+n.branch0+1:n.branch0])
   igc12_branch.set0 <- cbind(igc12_branch.set0,
				res.gene[9+2*n.branch0+1:n.branch0])
   igc21_branch.set0 <- cbind(igc21_branch.set0,
				res.gene[9+3*n.branch0+1:n.branch0])
   mut_branch.set0 <- cbind(mut_branch.set0,
				res.gene[9+4*n.branch0+1:n.branch0])
}

data0 <- data.frame(gene=rep(gene.set0,each=n.branch0),
			  nsites=as.numeric(nsites.set0),
			  branch=as.character(branch.set0),
                    bl=as.numeric(bl.set0),
			  tau=as.numeric(tau_branch.set0),
			  igc12=as.numeric(igc12_branch.set0),
			  igc21=as.numeric(igc21_branch.set0),
			  mut=as.numeric(mut_branch.set0))

data0$igc <- data0$igc12 + data0$igc21
data0$all <- data0$igc + data0$mut
head(data0)

### gene effect and branch effect ###
igc.logistic.gene <- glm(cbind(igc,mut) ~ gene,
                   family=binomial, data=data0)
summary(igc.logistic.gene)

igc.logistic.branch <- glm(cbind(igc,mut) ~ branch,
                   family=binomial, data=data0)
summary(igc.logistic.branch)

igc.logistic.gene_branch <- glm(cbind(igc,mut) ~ gene + branch,
                   family=binomial, data=data0)
summary(igc.logistic.gene_branch)

AIC(igc.logistic.gene)
AIC(igc.logistic.branch)
AIC(igc.logistic.gene_branch)

coef <- coef(igc.logistic.gene_branch)
n.gene <- length(unique(data0$gene))
gene.effect <- c(0,coef[2:n.gene])
names(gene.effect) <- sort(unique(data0$gene))
branch.effect <- c(0,coef[(n.gene+1):length(coef)])
names(branch.effect) <- sort(unique(data0$branch))

pdf("gene_branch_effect.pdf",height=6,width=15)
par(mar=c(5,5,4,2))
barplot(gene.effect,las=3,cex.names=1,xlab="gene",ylab="log odds ratio",
	col="blue",main="gene effect")
par(mar=c(15,5,4,2))
barplot(branch.effect,las=3,cex.names=1,xlab="branch",ylab="log odds ratio",
	col="blue",main="branch effect")
dev.off()



