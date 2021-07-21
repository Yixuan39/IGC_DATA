library(ape)
files <- list.files("Hiro/TGD_CDS_new")
trees <- list()
gene.set <- c()
n.tip <- c()
for (file in files){
  gene <- paste("MG",file,sep="")
  tree0 <- scan(paste("Hiro/TGD_CDS_new/",file,"/Pillar",file,".newick",sep=""),"")
  write(paste(tree0,";",sep=""),"temp.newick")
  trees[[gene]] <- read.tree("temp.newick")
  gene.set <- c(gene.set,gene)
  n.tip <- c(n.tip,length(trees[[gene]]$tip.label))
}
names(n.tip) <- gene.set
table(n.tip)

pdf("Hiro/tree_node_label.pdf")
for (gene in gene.set){
plot(trees[[gene]],show.node.label=TRUE,main=gene)
}
dev.off()

pdf("Hiro/tree_node_label0.pdf")
for (gene in gene.set){
tree1 <- trees[[gene]]
if (n.tip[gene]==8){
plot(tree1,show.node.label=TRUE,main=gene)}
}
dev.off()

