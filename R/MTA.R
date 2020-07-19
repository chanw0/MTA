### The final function
MTA=function(Y1,Y_new=NULL,phy.tree=NULL,timevec=NULL, transf="None",
             M=NULL,proportion.explained=0.85,
             k=5,lambda1.set=c(0.01,0.1,1,10),
             lambda2.set=seq(0,2,0.1),lambda3.set=c(0),num.sam=10,alpha=0.05)
{

############### phylogentic tree to laplacian penalty
if (is.null(phy.tree)) Laplacian.matrix=NULL else {

  distance.tree=distTips(phy.tree, tips = "all", method ="patristic")
  distance.tree=(distance.tree-min(distance.tree))/(max(distance.tree)-min(distance.tree))
  adjacency.tree=as.matrix(distance.tree)

  adjacency.tree=1-adjacency.tree
  diag(adjacency.tree)=1

  diag.L=rowSums(adjacency.tree)
  Laplacian.matrix=diag(diag.L)-adjacency.tree
}



########## data transformation before the analysis

  if(transf=="Log") {

    if(!(is.list(Y1))) {

      small.value=min(c(Y1)[c(Y1)!=0])/2
      for (ii in 1:dim(Y1)[1]) for(tt in 1:dim(Y1)[3]) {x=Y1[ii,,tt];
      if(min(x)==0) x=(x+small.value)/(sum((x+small.value)))
      Y1[ii,,tt]=log(x) }} else {

      small.value=min(unlist(Y1)[unlist(Y1)!=0])/2
      for(j in 1:lenght(Y1)) for (ii in 1:dim(Y1[[j]])[1]) for(tt in 1:dim(Y1[[j]])[3]) {x=Y1[[j]][ii,,tt];
      if(min(x)==0) x=(x+small.value)/(sum((x+small.value)))
      Y1[[j]][ii,,tt]=log(x) }}

    if(!(is.null(Y_new))) {

      small.value=min(c(Y_new)[c(Y_new)!=0])/2
      for (ii in 1:dim(Y_new)[1]) for(tt in 1:dim(Y_new)[3]) {x=Y_new[ii,,tt];
      if(min(x)==0) x=(x+small.value)/(sum((x+small.value)))
      Y_new[ii,,tt]=log(x) }}
  }


  if(transf=="ALR") {

    if(!(is.list(Y1))) {
      Y22=Y1[,-dim(Y1)[2],]
      small.value=min(c(Y1)[c(Y1)!=0])/2
      for (ii in 1:dim(Y1)[1]) for(tt in 1:dim(Y1)[3]) {x=Y1[ii,,tt];
      if(min(x)==0) x=(x+small.value)/(sum((x+small.value)))
      Y22[ii,,tt]=log(x[-length(x)]/x[length(x)]) }

      Y1=Y22

      } else { YY2=list()

        small.value=min(unlist(Y1)[unlist(Y1)!=0])/2
        for(j in 1:lenght(Y1)) {

          Y22=Y1[[j]][,-dim(Y1[[j]])[2],]

          for (ii in 1:dim(Y1[[j]])[1]) for(tt in 1:dim(Y1[[j]])[3]) {x=Y1[[j]][ii,,tt];
        if(min(x)==0) x=(x+small.value)/(sum((x+small.value)))
        Y22[ii,,tt]=log(x[-length(x)]/x[length(x)]) }

          YY2[[j]]=Y22
          }}

    if(!(is.null(Y_new))) {

      Y22=Y_new[,-dim(Y_new)[2],]
      small.value=min(c(Y_new)[c(Y_new)!=0])/2
      for (ii in 1:dim(Y_new)[1]) for(tt in 1:dim(Y_new)[3]) {x=Y_new[ii,,tt];
      if(min(x)==0) x=(x+small.value)/(sum((x+small.value)))
      Y22[ii,,tt]=log(x[-length(x)]/x[length(x)]) }

      Y_new=Y22
    }
}



    if(transf=="CLR") {

      if(!(is.list(Y1))) {

        small.value=min(c(Y1)[c(Y1)!=0])/2
        for (ii in 1:dim(Y1)[1]) for(tt in 1:dim(Y1)[3]) {x=Y1[ii,,tt];
        if(min(x)==0) x=(x+small.value)/(sum((x+small.value)))
        Y1[ii,,tt]=log(x/(cumprod(x)[length(x)])^(1/length(x))) }} else {

          small.value=min(unlist(Y1)[unlist(Y1)!=0])/2
          for(j in 1:lenght(Y1)) for (ii in 1:dim(Y1[[j]])[1]) for(tt in 1:dim(Y1[[j]])[3]) {x=Y1[[j]][ii,,tt];
          if(min(x)==0) x=(x+small.value)/(sum((x+small.value)))
          Y1[[j]][ii,,tt]=log(x/(cumprod(x)[length(x)])^(1/length(x))) }}

      if(!(is.null(Y_new))) {

        small.value=min(c(Y_new)[c(Y_new)!=0])/2
        for (ii in 1:dim(Y_new)[1]) for(tt in 1:dim(Y_new)[3]) {x=Y_new[ii,,tt];
        if(min(x)==0) x=(x+small.value)/(sum((x+small.value)))
        Y_new[ii,,tt]=log(x/(cumprod(x)[length(x)])^(1/length(x))) }}
    }


######## analysis involving three parts

  if(!(is.list(Y1))) {trend.result=MTA_pattern(Y1,M=M,proportion.explained=proportion.explained, k=k,
                                           Laplacian.matrix=Laplacian.matrix,timevec=timevec,
                                           lambda1.set=lambda1.set,lambda2.set=lambda2.set,lambda3.set=lambda3.set)

                  ff=trend.result[[2]]; colnames(ff)=paste("contributing to", 1:ncol(ff), "common trend")

              result.sum=list(trend.result[[1]],ff)
              names(result.sum)=c("Trend Plot","Factor Score")

              }


  if(is.list(Y1) & is.null(Y_new)) {

   pairwise=combn(length(Y1),2)

   if(is.null(names(Y1))) label.training=1:length(Y1) else label.training=names(Y1)

    comparsion.res=list()
    for(i in 1:ncol(pairwise)) {

    tem=pairwise[,i]
    result.sum=MTA_comparison(Controls=Y1[[tem[1]]], Cases=Y1[[tem[2]]],k=k,Laplacian.matrix=Laplacian.matrix,timevec=timevec,lambda1.set=lambda1.set,
                            lambda2.set=lambda2.set,lambda3.set=lambda3.set,num.sam=num.sam,alpha=alpha)

    comparsion.res[[i]]=result.sum
    }
    names(comparsion.res)=apply(pairwise,2,function(x){return(paste("Comparison between group", x[1], "and group", x[2]))})

    result.sum=comparsion.res
   }


  if(is.list(Y1) & !(is.null(Y_new))) {

    if(is.null(names(Y1))) label.training=1:length(Y1) else label.training=names(Y1)

    result.sum=MTA_classification(group_new=Y_new, training.set=Y1, label.training=label.training, proportion.explained=proportion.explained, k=k,
                                  Laplacian.matrix=Laplacian.matrix,timevec=timevec,lambda1.set=lambda1.set,
                                  lambda2.set=lambda2.set,lambda3.set=lambda3.set)



    }


result.sum

}
