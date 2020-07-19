############# Group comparison

MTA_classification=function(group_new, training.set, label.training, proportion.explained, k,Laplacian.matrix,timevec,lambda1.set,
                            lambda2.set,lambda3.set)
{

  if(is.null(timevec)) timevec=1:dim(training.set[[1]])[3]

  ######## cubic spline basis
  BS = create.bspline.basis(timevec, norder=4)
  B = getbasismatrix(timevec, BS)

  pred.matrix=list()

  for(i in 1:length(label.training)) {

  tem1=MTA_pattern(training.set[[i]],M=NULL,proportion.explained=proportion.explained, k=k,
                   Laplacian.matrix=Laplacian.matrix,timevec=timevec,lambda1.set=lambda1.set,
                  lambda2.set=lambda2.set,lambda3.set=lambda3.set)

  ff=tem1[[2]];cc=tem1[[3]]
  pred.matrix[[i]]=ff%*%cc%*%(t(B))

    }

  aa=sapply(1:dim(group_new)[1], function(x,group_new,pred.matrix,label.training){

    dist=sapply(1:length(label.training), function(y,pred.matrix,subj.new) {
                         return(sqrt(sum((subj.new-pred.matrix[[y]])^2)))}, pred.matrix=pred.matrix,subj.new=group_new[x,,])

    return(label.training[which.min(dist)])},
    group_new=group_new,pred.matrix=pred.matrix,label.training=label.training)

res=cbind(1:dim(group_new)[1],aa)
colnames(res)=c("New subject","Group label")
return(data.frame(res,check.names =FALSE))
}
