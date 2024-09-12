rm(list=ls())
library(glmnet)
otus <- readRDS("OTU_data_spike.rds")#[,1:3]
otus[1:5,1:5]
indoxyl <- read.csv(file = "indoxylSulfate.csv")
rownames(indoxyl) <- indoxyl[,1]
head(indoxyl)
indoxyl <- indoxyl[, -1, drop=FALSE]
head(indoxyl)
all(rownames(otus)==rownames(indoxyl))

otus <- data.frame(otus)
fit <- lm(formula = indoxyl[,1]~. , data=otus)

pred <- predict(object = fit, newdata = otus)
plot(pred, indoxyl[,1])
cor(pred, indoxyl[,1])

dim(otus)

preds <- NULL
for(i in 1:nrow(otus)){
  otus.test <- otus[i, , drop=FALSE]
  otus.train <- otus[-i, , drop=FALSE]
  indoxyl.test <- indoxyl[i,1]
  indoxyl.train <- indoxyl[-i,1]
  fit <- lm(formula = indoxyl.train~. , data=otus.train)
  pred <- predict(object = fit, newdata = otus.test)
  preds <- c(preds, pred)
}
preds
plot(x= indoxyl[,1], y=preds)
cor(x= indoxyl[,1], y=preds)

fit <- cv.glmnet(x=as.matrix(otus), y = indoxyl[,1], nfolds = 10)
plot(fit)
pred <- predict(object = fit, newx = as.matrix(otus), s = "lambda.min")
plot(pred, indoxyl[,1])
cor(pred, indoxyl[,1])

preds <- NULL
for(i in 1:nrow(otus)){
  otus.test <- otus[i, , drop=FALSE]
  otus.train <- otus[-i, , drop=FALSE]
  indoxyl.test <- indoxyl[i,1]
  indoxyl.train <- indoxyl[-i,1]
  fit <- cv.glmnet(x=as.matrix(otus.train), y = indoxyl.train, nfolds = 3)
  pred <- predict(object = fit, newx = as.matrix(otus.test), s = "lambda.min")
  preds <- c(preds, pred)
}
coefs <- coef(fit, s = "lambda.min")
coefs[coefs[,1]!=0,]

plot(fit)
plot(preds, indoxyl[,1])
cor(preds, indoxyl[,1])
fit0 <- glmnet(x=as.matrix(otus), y = indoxyl[,1])
par(mfrow=c(1,2))
plot(fit)
plot(fit0)

rm(list=ls())
load("mmml_vsn.rda")
ls()
expr <- exprs(mmml.vsn)
expr[1:5,1:5]
#str(mmml.vsn@assayData)
#mmml.vsn@assayData$exprs
#mmml.vsn@assayData[["exprs"]]

pheno <- pData(mmml.vsn)
pheno
GCBscore <- pheno[,"GCBABC"]
names(GCBscore) <- rownames(pheno)
expr[1:5,1:5]
colnames(expr)==names(GCBscore)
all(colnames(expr)==names(GCBscore))
Ids.train <- sample(1:ncol(expr), size = 2*ncol(expr)/3, replace = FALSE)
Ids.test <- c(1:ncol(expr))[!c(1:ncol(expr)) %in% Ids.train]
intersect(Ids.train, Ids.test)
x.training <- t(expr[,Ids.train])
x.test <- t(expr[,Ids.test])
y <- c(1,0)[GCBscore]
#y[1:20]
#GCBscore[1:20]
y.training <- y[Ids.train]
y.test <- y[Ids.test]

dim(x.training)
fit <- cv.glmnet(x = x.training, y = y.training, family="binomial", alpha=0)
predictions <- predict(object = fit, newx = x.test, s = "lambda.min", type="response")
par(mfrow=c(1,1))
hist(predictions)

library(ROCR)
?prediction
pred <- prediction(predictions = predictions, labels = y.test)
perf <- performance(pred, 'tpr', 'fpr')
plot(perf)
abline(a = 0, b = 1, col='red')

perf <- performance(pred, 'auc')
str(perf)
perf@y.values

?performance
perf <- performance(pred, 'rec', 'prec')
plot(perf)

fit.ridge <- glmnet(x = x.training[,1:1000], y = y.training, family="binomial", alpha=0)
fit.lasso <- glmnet(x = x.training[,1:1000], y = y.training, family="binomial", alpha=1)
par(mfrow=c(1,2))

fit.ridge$beta[1:5,50:55]
fit.ridge$lambda[50:55]
plot(fit.ridge)
plot(fit.lasso)


sum(fit.lasso$beta[,1]!=0)
sum(fit.lasso$beta[,50]!=0)


fit.lasso$lambda[50:55]


plot(fit)
