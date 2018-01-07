library(partykit)
library(survival)
set.seed(290875)


#########################################split index#########################################################
get_split <- function(y, dataset,i,min_size){
  split<-NULL
  split$log<- -Inf
  split$group$right<-NULL
  split$group$left<-NULL
  split$index <- 0
  split$row <-0
  split$br <- 0
  
  for (index in 1:(length(dataset[0,]))){ #4
    if(i!=index){
    form <-as.formula(paste(y, paste(names(dataset[index]), collapse="+"),collapse = "+"))
    for (row in 1:length(dataset[,index])){ #150
      br<-dataset[row,index]
      group <- test_split(dataset, index, br)
      ctrl<-survreg.control(maxiter=30, rel.tolerance=1e-06,toler.chol=1e-7)
      if((length(group$left)>min_size)&&(length(group$right)>min_size))
      {
      fitl <- survreg(form, dataset[group$left,] ,dist = "t",control = ctrl)
      l<-logLik(fitl)
      fitr <- survreg(form, dataset[group$right,] ,dist = "t",control = ctrl)
      r<-logLik(fitr)
      log<- l+r

      if(is.nan(l)){
        log<-r
      }
      if(is.nan(r)){
        log<-l
      }
      if(is.nan(log)){
        log<- split$log
      }
      #print(log)
       if (split$log < as.numeric(log)){
        split$index <- index
        split$row <-row
        split$br <- br
        split$log <- log
        split$group<-group
        
       }
      }
     }
    }
   }
  print(split$log)
  print(split$row)
  print(split$index)
  print("end")
  return(split)
}

#######################################right-left split###################################################
test_split <- function(dataset, index, value){
  left <- c()
  right <- c()
  for (row in 1:length(dataset[,2])){
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
split <- function(y,sp, depth, maxdepth,min_size,i,dataset,b_index){
  View(b_index)
  rsize<- length(sp$group$right)
  lsize<- length(sp$group$left)
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
  if (!((any(sp$group$left))&&(any(sp$group$right)))){
    return(b_index)
  }
  
  if (length(dataset[,2])<=min_size){
    term<-1
    tag<-b_index[id-1,8]
    node <- c(id,varid,br,term,kidi,kidf,info,tag)
    b_index <- rbind(b_index,node)
    return(b_index)
  }
  
  #gini for right
  tag<-'right'
  print("right")
  sp<- get_split(y,dataset[sp$group$right,],i,min_size)
  rsize<- length(sp$group$right)
  lsize<- length(sp$group$left)
  if((sp$log==-Inf)||(sp$log==0)||(as.numeric(depth)>as.numeric(maxdepth))||(as.numeric(rsize)<as.numeric(min_size))){
    term <- 1
    sp$group$right <- NULL
    node <- c(id,varid,br,term,kidi,kidf,info,tag)
    b_index <- rbind(b_index,node)
  }
  else{
    varid <- sp$index
    node <- c(id,varid,sp$br,term,kidi,kidf,info,tag)
    b_index <- rbind(b_index,node)
    print(sp$group$right)
    print(length(dataset[,2]))
    print(length(dataset[,sp$group$right]))
    b_index <- split(y,sp,depth+1,maxdepth,min_size,i,dataset[sp$group$right,],b_index)
  }
  
  #b_index
  tag<-'left'
  print("left")
  l<- length(b_index[,1])
  id <- as.numeric(b_index[l,1])+1
  term <- 0
  varid<-0
  br<-0
  info <- ''
  
  #gini for left
  sp <- get_split(y,dataset[sp$group$left,],i,min_size)
  rsize<- length(sp$group$right)
  lsize<- length(sp$group$left)
  if((sp$log==-Inf)||(sp$log==0)||(as.numeric(depth)>as.numeric(maxdepth))||(as.numeric(lsize)<as.numeric(min_size))){
    term <- 1
    sp$group$left <- FALSE
    node <- c(id,varid,br,term,kidi,kidf,info,tag)
    b_index <- rbind(b_index,node)
  }
  else{
    varid <- sp$index
    br<- sp$br
    node <- c(id,varid,br,term,kidi,kidf,info,tag)
    b_index <- rbind(b_index,node)
    gr <- sp$group
    b_index <- split(y,gr,depth+1,maxdepth,min_size,i,dataset[gr$left,],b_index)
  }
  
  return(b_index)
}



#########################################terminal output#############################################
#not required
to_terminal <- function(group){
  if(length(group$left)>length(group$right)){
    outcomes = length(group$left)
  }
  else{
    outcomes = length(group$right)
  }
  #print(sort(-table(outcomes)))
  return (outcomes)
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
cart <- function(y,dataset,i,maxdepth,min_size){
  b_index <- c(1)
  b <- get_split(y,dataset,i,min_size)
  b_index <- data.frame(b_index,b$index,b$br,0,0,0,'',"root",stringsAsFactors=FALSE)
  colnames(b_index) <- c("id","varid","break","terminal","kidi","kidf","info","tag")
  b_index <- split(y,b,1,maxdepth,min_size,i,dataset,b_index)
  b_index <- kid(b_index)
  View(b_index)
  
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
  n<- c(names(dataset[i]),".")
  form <-as.formula(paste(y, paste(n, collapse="+"),collapse = "+"))
  dataset <- model.frame(form, data = dataset)
  tree <- party(node, data = dataset, fitted = data.frame("(fitted)" = fitted_node(node, data = dataset),
                                                          "(response)" = model.response(dataset),
                                                          check.names = FALSE), terms = terms(dataset))
  tree <- as.constparty(tree)
  
  return(tree)
}






#########################main code##################################################################
data(neuroblastomaProcessed, package="penaltyLearning")

survTarget <- with(neuroblastomaProcessed, {Surv(target.mat[, 1], target.mat[,2], type="interval2")})
trainData <- data.frame(log.penalty=survTarget, neuroblastomaProcessed$feature.mat)


y <- "log.penalty~"
form <-as.formula(paste(y, paste(names(dataset[2]), collapse="+"),collapse = "+"))
fit <- survreg(form, trainData,dist = "t")
l<-logLik(fit)

fit<- cart(y,trainData[1:45,],1,4,1)
plot(fit)

#doubts
#1) loglik = 0 wht to do with tht
#2) prob in dataset[sp$group$right,] for length =1 
#3) using 't' is it correct


