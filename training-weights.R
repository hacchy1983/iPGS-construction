############################################################################################
## Load data
## File paths in this section should be modified for each research institute
############################################################################################

## Load case-control flags
flag <- read.table("/path/to/case-control-flag.txt")
colnames(flag) <- c("ID", "BINARY_OUTCOME")
data <- flag

## Load polygenic scores
score_001 <- read.table("/path/to/best-score-file-001.txt")[,c(1,4)]
colnames(score_001) <- c("ID", "model-name-001")
data <- merge(data, score_001, by="ID")

score_002 <- read.table("/path/to/best-score-file-002.txt")[,c(1,4)]
colnames(score_002) <- c("ID", "model-name-002")
data <- merge(data, score_002, by="ID")

... 

score_XXX <- read.table("/path/to/best-score-file-XXX.txt")[,c(1,4)]
colnames(score_XXX) <- c("ID", "model-name-XXX")
data <- merge(data, score_XXX, by="ID")

## Load covariates (such as age and sex)
cov <- read.table("/path/to/covariates.txt", header=TRUE)
data <- merge(data, cov, by="ID")

## Load principal components
pc <- read.table("/path/to/PCs.txt")
colnames(pc) <- c("ID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", 
              "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")
data <- merge(data, pc, by="ID")


############################################################################################
## Determine weights for input PGS models using Elastic Net Logistic Regression
## `glmnet` package should be installed before executing this section
############################################################################################
library(glmnet)
alpha <- seq(0.01, 0.99, by=0.01)
dev.df <- NULL
auc.df <- NULL

scaling.factors <- apply(data[,3:ncol(data)], 2, sd)
write.table(scaling.factors, "scaling-factors.txt", quote=F, row.names=T, col.names=F, sep="\t")

X <- scale(as.matrix(data[,3:ncol(data)]))
Y <- data[,2]
for( i in 1:length(alpha) ) {
  set.seed(88778807)
  m      <- cv.glmnet( x=X, y=Y, family="binomial", alpha=alpha[i], type.measure="auc" )
  tmp.df <- data.frame( ALPHA=alpha[i], CVM=max(m$cvm), LAMBDA=m$lambda.min )
  auc.df <- rbind( auc.df, tmp.df )
}
auc.best_df <- auc.df[auc.df$CVM == max(auc.df$CVM),][1,]

auc.best.model <- glmnet(x=X, y=Y, family="binomial", lambda=auc.best_df$LAMBDA, alpha=auc.best_df$ALPHA)
write.table(as.data.frame(auc.best.model$beta[,1]), file="iPGS.auc-best.model.beta", quote=F, row.names=T, col.names=F, sep="\t")
write.table(auc.best_df, file="iPGS.auc-best.model.param", quote=F, row.names=F, col.names=T, sep="\t")


m <- glm(as.formula(paste("Y ~ ", paste(paste0("X[,\"", colnames(X), "\"]"), collapse=" + ") ) ), family="binomial" )
s <- summary(m)
write.table(s$coefficients, "glm.beta", quote=F, row.names=T, col.names=T, sep="\t")
