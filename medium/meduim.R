library(partykit)
library(rpart)
library(rpart.plot)

#assuming all features are numerical vectors
#############################################gini calculation##########################################
gini_process<-function(classes,splitvar = NULL){
  if (is.null(splitvar)){
    base_prob <-table(classes)/length(classes)
    return(1-sum(base_prob**2))
  }
  
  base_prob <-table(splitvar)/length(splitvar)
  crosstab <- table(classes,splitvar)
  crossprob <- prop.table(crosstab,2)
  Gini <- c()
  for(i in 1:length(crossprob[1,])){
    Node_Gini <- 1-sum(crossprob[,i]**2)
    Gini <- c(Gini,Node_Gini)
  }
  return(sum(base_prob * Gini))
}

#########################################split index#########################################################
get_split <- function(dataset,target,i){
  b <- c(0,0)
  b_score <- 2
  for (index in 1:(length(dataset[0,]))){ #4
    if(names(dataset[index])!=names(dataset[i])){
      if((max(dataset[,index])==Inf)||(min(dataset[,index])==-Inf)||(is.factor(dataset[names(dataset[index])]))){
        if(is.factor(dataset[names(dataset[index])])){
          spv <- t(dataset[names(dataset[index])])
        }
        else{
          for (row in 1:length(dataset[index,])){ #150
            #print(c(row,index))
            br<-dataset[row,index]
            spv <- t(dataset[names(dataset[index])] < as.numeric(br))[1,]
            if(any(spv)){
              gini <- gini_process(target, spv)
              #print(gini)
              
              if (gini < b_score){
                b <- c(index,br)
                b_score <- gini
              }
            }
          }}
      }
      else{
        for (k in min(dataset[,index]):max(dataset[,index])){ #150
          spv <- t(dataset[names(dataset[index])] < as.numeric(k))[1,]
          if(any(spv)){
            gini <- gini_process(target, spv)
            #print(gini)
            #print(c(index,k))
            
            if (gini < b_score){
              b <- c(index,k)
              b_score <- gini
            }
          }
        }
      }
    }
  }
  #print(b_score)
  #print(b)
  #print(gini)
  #print(c(index,k))
  if(b_score==0){ return(c(0,0))}
  return(b)
}

#######################################right-left split###################################################
test_split <- function(dataset, index, value){
  left <- c()
  right <- c()
  for (row in 1:length(dataset[,1])){
    if (as.numeric(dataset[row,index]) < value){
      left <- cbind(left,row)
    }
    else{
      right <- cbind(right,row)
    }
  }
  group <-c()
  group$left <- left
  group$right <- right
  return (group)
}

#####################################3tree############################################################3
# Create child splits for a node or make terminal
split <- function(group, depth, maxdepth,min_size,target,i,dataset,b_index){
  rsize<- length(group$left)
  lsize<- length(group$right)
  #b_index
  l<- length(b_index[,1])
  id <- as.numeric(b_index[l,1])+1
  term <- 0
  varid<-0
  kidi<-0
  kidf<-0
  br<-0
  info <- ''
  tag<-'no'
  
  # check for a no split
  if (!((any(group$left))&&(any(group$right)))){
    return(b_index)
  }
  
  #gini for right
  tag<-'right'
  #print(length(dataset[group$right]))
  b <- get_split(dataset[group$right,],target[group$right],i)
  if((b[1]==0)||(as.numeric(depth)>as.numeric(maxdepth))||(as.numeric(rsize)<as.numeric(min_size))){
    term <- 1
    info <- to_terminal(group$right,target)
    group$right <- NULL
    node <- c(id,varid,br,term,kidi,kidf,info,tag)
    b_index <- rbind(b_index,node)
  }
  else{
    varid <- b[1]
    br<- b[2]
    info <- to_terminal(group$right,target)
    node <- c(id,varid,br,term,kidi,kidf,info,tag)
    b_index <- rbind(b_index,node)
    gr <- test_split(dataset[group$right,], b[1], b[2])
    b_index <- split(gr,depth+1,maxdepth,min_size,target[group$right],i,dataset[group$right,],b_index)
  }
  
  #b_index
  tag<-'left'
  l<- length(b_index[,1])
  id <- as.numeric(b_index[l,1])+1
  term <- 0
  varid<-0
  br<-0
  info <- ''
  
  #gini for left
  b <- get_split(dataset[group$left,],target[group$left],i)
  if((b[1]==0)||(as.numeric(depth)>as.numeric(maxdepth))||(as.numeric(lsize)<as.numeric(min_size))){
    term <- 1
    info <- to_terminal(group$left,target)
    group$left <- FALSE
    node <- c(id,varid,br,term,kidi,kidf,info,tag)
    b_index <- rbind(b_index,node)
  }
  else{
    varid <- b[1]
    br<- b[2]
    info <- to_terminal(group$left,target)
    node <- c(id,varid,br,term,kidi,kidf,info,tag)
    b_index <- rbind(b_index,node)
    gr <- test_split(dataset[group$left,], b[1], br)
    b_index <- split(gr,depth+1,maxdepth,min_size,target[group$left],i,dataset[group$left,],b_index)
  }
  
  return(b_index)
}



#########################################terminal output#############################################
to_terminal <- function(group,target){
  outcomes = target[group]
  #print(sort(-table(outcomes)))
  return (names(sort(-table(outcomes)))[1])
}


#########################################terminal output#############################################
kid <- function(b_index){
  n<- 0
  k<-c()
  
  for(i in 2:(length(b_index[,1])-1)){
    if(b_index[i,8]=="right"){
      b_index[i-1,6] = i
      if(b_index[i+1,8]=="left"){
        b_index[i-1,5] = i+1      
        }
    }
    if((b_index[i,8]=="left")&&(b_index[i-2,5]==0)){
      k<-c(k,i)
    }
  }
  c<-1
  for(i in length(b_index[,1]):1){
    if((b_index[i,5]==0)&&(b_index[i,4]==0)){
      b_index[i,5]=k[c]
      c=c+1
    }
  }
  return (b_index)
}


#############################################cart algo################################################
cart <- function(dataset,target,i,maxdepth,min_size){
  b_index <- c(1)
  b <- get_split(dataset,target,i)
  b_index <- data.frame(b_index,b[1],b[2],0,0,0,'',"root",stringsAsFactors=FALSE)
  colnames(b_index) <- c("id","varid","break","terminal","kidi","kidf","info","tag")
  group <- test_split(dataset, b[1], b[2])
  b_index <- split(group,1,maxdepth,min_size,target,i,dataset,b_index)
  b_index <- kid(b_index)
  #View(b_index)
  
  nodelist <- list(
    #root node
    list(id = 1L, split = partysplit(varid = as.integer(b_index[1,2]), breaks = as.numeric(b_index[1,3])), 
         kids = c(as.integer(b_index[1,5]),as.integer(b_index[1,6]))))
  
  for(i in 2:length(b_index[,1])){
    if(as.numeric(b_index[i,4])){
      nodelist[i] <- list(list(id = as.integer(b_index[i,1]), info = b_index[i,7]))
    }
    else{
      nodelist[i] <- list(
        list(id = as.integer(b_index[i,1]), split = partysplit(varid = as.integer(b_index[i,2]), 
                                                               breaks = as.numeric(b_index[i,3])), 
             kids = c(as.integer(b_index[i,5]),as.integer(b_index[i,6])),info = b_index[i,7]))
    }
  }
  ## convert to a recursive structure
  node <- as.partynode(nodelist)
  
  ## set up party object
  dataset <- model.frame(Kyphosis ~ ., data = kyphosis)
  tree <- party(node, data = dataset, fitted = data.frame("(fitted)" = fitted_node(node, data = dataset),
                                                          "(response)" = model.response(dataset),
                                                          check.names = FALSE), terms = terms(dataset))
  tree <- as.constparty(tree)
  return(tree)
}


#############################################3main program################################################
data(iris)
data(kyphosis)
kyphosis

tree<- cart(kyphosis,kyphosis$Kyphosis,1L,4,10)
plot(tree)

plot(as.party(rpart(kyphosis,formula = Kyphosis~.)))

