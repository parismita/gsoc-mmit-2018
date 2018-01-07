library(trtf)
library(survival)

data(neuroblastomaProcessed, package="penaltyLearning")

finite.targets <- with(neuroblastomaProcessed, {
  data.frame(log.penalty=target.mat[is.finite(target.mat)])
})

m <- ctm(as.basis(~log.penalty, data=finite.targets), todistr="Normal")

train.Surv <- with(neuroblastomaProcessed, { Surv(target.mat[, 1], target.mat[,2], type="interval2")})

## Train on n=50 observations 
## model (predict all 1).
train.feature.mat <- neuroblastomaProcessed$feature.mat
train.df <- data.frame(log.penalty=train.Surv, train.feature.mat)[1:50,]
mlt.fit <- mlt(m, data=train.df)
tree.fit <- trafotree(m, formula = log.penalty ~ ., data=train.df, mltargs=list(theta=coef(mlt.fit)))

##prediction
pred.vec <- predict(tree.fit)
pred.log.penalty <- tree.fit$coef[names(pred.vec), "(Intercept)"]
is.lo <- pred.log.penalty < neuroblastomaProcessed$target.mat[1:50, 1]
is.hi <- neuroblastomaProcessed$target.mat[1:50, 2] < pred.log.penalty
is.error <- is.lo | is.hi
mean(is.error)

plot(tree.fit)


## Train on all 3418 observations and 2 features
train.feature.mat <- neuroblastomaProcessed$feature.mat[, c("log.n", "log.mad")]
train.df <- data.frame(log.penalty=train.Surv, train.feature.mat)
mlt.fit <- mlt(m, data=train.df)
tree.fit <- trafotree(m, formula = log.penalty ~ ., data=train.df, mltargs=list(theta=coef(mlt.fit)))
node_party(tree.fit) # 15-node tree.
tree.fit$coef # 8 leaf nodes.

##prediction
pred.vec <- predict(tree.fit)
pred.log.penalty <- tree.fit$coef[names(pred.vec), "(Intercept)"]
is.lo <- pred.log.penalty < neuroblastomaProcessed$target.mat[, 1]
is.hi <- neuroblastomaProcessed$target.mat[, 2] < pred.log.penalty
is.error <- is.lo | is.hi
mean(is.error)

## train on all 3418 observations and 117 features.
train.df <- data.frame(log.penalty=train.Surv, neuroblastomaProcessed$feature.mat)
mlt.fit <- mlt(m, data=train.df)
tree.fit <- trafotree(m, formula = log.penalty ~ ., data=train.df, mltargs=list(theta=coef(mlt.fit)))
node_party(tree.fit) # 19-node tree.
tree.fit$coef 

plot(tree.fit)

##prediction
pred.vec <- predict(tree.fit)
pred.log.penalty <- tree.fit$coef[names(pred.vec), "(Intercept)"]
is.lo <- pred.log.penalty < neuroblastomaProcessed$target.mat[, 1]
is.hi <- neuroblastomaProcessed$target.mat[, 2] < pred.log.penalty
is.error <- is.lo | is.hi
mn<- mean(is.error)
print (mn)
