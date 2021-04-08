IGC0 <- IGC_Full_noforce_noclock
dim(IGC0)
dimnames(IGC0)
branch <- IGC0[10:25,]
tau0 <- IGC0[9,]
tau <- IGC0[26:41,]
dimnames(tau)
apply(tau,1,mean)
tau1 <- tau[-2,]
branch1 <- branch[-2,]

tau.data <- data.frame(tau=as.numeric(tau1),
				branch=rep(rownames(branch1),dim(tau1)[2]),
				paralog=rep(colnames(tau1),each=dim(tau1)[1]))

hist(tau.data[,1],nclass=50)

#############
tau.gamma <- glm(tau~-1+branch+paralog,family=Gamma(link = "log"),
			 data=tau.data)
tau.gamma.pred <- predict(tau.gamma,type="response")
tau.res <- tau.data[,1]/tau.gamma.pred
hist(tau.res)
summary(tau.gamma)
################

matplot(tau1,type="b",lty=1,pch=1)
matplot(log(tau1),type="b",lty=1,pch=1)

###########
logtau.lm <- lm(tau~-1+branch+paralog,data=tau.data)
summary(logtau.lm)
##### Adjusted R-squared:  0.9162 #######
names(coef(logtau.lm))
pdf("log_tau_lm.pdf",width=12,height=6)
par(mfrow=c(1,2))
barplot(coef(logtau.lm)[1:15],las=3,main="branch effect",
	names=substr(names(coef(logtau.lm))[1:15],7,20),
	cex.names=0.7,col=2)
barplot(coef(logtau.lm)[16:59],las=3,main="paralog effect",
	names=substr(names(coef(logtau.lm))[16:59],8,20),
	cex.names=0.6,col=2)
par(mfrow=c(1,1))
dev.off()

