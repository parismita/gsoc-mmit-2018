library(rpart)
library(rpart.plot)
library(partykit)
library(caret)
set.seed(290875)

#model how well cells in an image are segmented
data("segmentationData")
segmentationData$Cell<-NULL
train<-subset(segmentationData,Case=="Train")
test<-subset(segmentationData,Case=="Test")
train$Case<-NULL
test$Case<-NULL

#formula->all expt y,time,status
form <- Class ~ . - Class
ctrl <- rpart.control(minsplit = 15, minbucket = 5)
tree <- rpart(formula = form, data = train, na.action = na.rpart, control = ctrl, method = "class")
rpart.plot(tree)
#summary(tree)
plotcp(tree)
printcp(tree)

#pruned tree to prevent overfitting, where inappropriate nodes are removed
pdtree<- prune(tree, cp=tree$cptable[which.min(tree$cptable[,"xerror"]),"CP"])
rpart.plot(pdtree)

#Model Testing
out<-predict(pdtree,test)
status_predicted<- colnames(out)[max.col(out, ties.method = c("first"))] # predicted
status_input<- as.character (test$Class) # actuals
mean (status_input != status_predicted) # misclassification %

par(mfrow=c(1,2))
rsq.rpart(pdtree)
par(mfrow=c(1,1))



