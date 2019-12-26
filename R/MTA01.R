############Extracted the dynamic trends from a group of subjects

MTA01=function(All.x,k,Laplacian.matrix,timevec,lambda1.set,lambda2.set,lambda3.set)
{

  N=dim(All.x)[1]; P=dim(All.x)[2];T=dim(All.x)[3];
  if(is.null(timevec)) timevec=timevec=1:T

  BS = create.bspline.basis(timevec, norder=4)
  B = getbasismatrix(timevec, BS)  #basis function
  Omega = getbasispenalty(BS)


  if(N==1) {

    Ex.pro=NULL
    for(lambda1 in lambda1.set) for(lambda3 in lambda3.set) for(lambda2 in lambda2.set) {

          res1=MTA00(All.x,Laplacian.matrix,timevec,lambda1,lambda2,lambda3)
          c.updated=res1[[1]];f.updated=res1[[2]];

          Ex.pro=rbind(Ex.pro,c(lambda1,lambda2,lambda3,
                                sum((as.matrix(f.updated)%*%c.updated%*%(t(B)))^2)/sum(All.x^2)))
        }

    tem.in=which.max(Ex.pro[,4])
    lambda1.fin=Ex.pro[tem.in,1];lambda2.fin=Ex.pro[tem.in,2];lambda3.fin=Ex.pro[tem.in,3];

    res1=MTA00(All.x,Laplacian.matrix,timevec,lambda1.fin,lambda2.fin,lambda3.fin)

    } else {res1=MTA.cv(All.x,k,Laplacian.matrix,timevec,lambda1.set,lambda2.set,lambda3.set)}


  c.updated=res1[[1]];f.updated=res1[[2]];
  return(list(c.updated,f.updated))

}
