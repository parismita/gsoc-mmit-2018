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
ctrl <- ctree_control(teststat = "quad",minsplit = 40, minbucket = 20)
tree <- ctree(formula = form, data = train, control = ctrl, method = "class")
plot(tree)

#Model Testing
out<-predict(tree,test,type = "prob")
status_predicted<- colnames(out)[max.col(out, ties.method = c("first"))] # predicted
status_input<- as.character (test$Class) # actuals
m <- mean (status_input != status_predicted) # misclassification %
print(m)

#Estimated conditional class probabilities depending on the first split variable.
prob <- predict(tree,test, type = "prob")[,1] + runif(nrow(test), min = -0.01, max = 0.01)
splitvar <- character_split(split_node(node_party(tree)), data = data_party(tree))$name
plot(test[[splitvar]], prob, pch = as.numeric(test$Class), ylab = "Conditional Class Prob.",xlab = splitvar)
abline(v = split_node(node_party(tree))$breaks, lty = 2)
legend(0.15, 0.7, pch = 1:2, legend = levels(test$Class), bty = "n")

