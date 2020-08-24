### Cross Validation

MTA.cv=function(All.x,k,Laplacian.matrix,timevec,lambda1.set,lambda2.set,lambda3.set)
{

  if(is.null(timevec)) timevec=1:dim(All.x)[3]

  ######## cubic spline basis
  BS = create.bspline.basis(timevec, norder=4)
  B = getbasismatrix(timevec, BS)
  Omega = getbasispenalty(BS)

for (jj in 1:100) {
  
  cv.group=sample(rep(seq_len(k), length.out = dim(All.x)[1]))

  allmse=NULL
  for(lambda1 in lambda1.set) for(lambda3 in lambda3.set) for(lambda2 in lambda2.set) {

        mse=0
        for(k_in in 1:k)
        {
          tem.x=All.x[cv.group!=k_in,,]

          #       if(length(dim(tem.x))<3) {tem.x=array(NA, dim = c(1,dim(tem.x)[1],dim(tem.x)[2]))
          #        tem.x[1,,]=All.x[-rr,,]}

          tem.res=MTA00(x=tem.x,Laplacian.matrix=Laplacian.matrix,timevec=timevec,
                      lambda1=lambda1,lambda2=lambda2,lambda3=lambda3)

          c.updated=tem.res[[1]]; f.updated=tem.res[[2]]

          if(sum(abs(f.updated))!=0) {

            res=sapply(which(cv.group==k_in), function(x,aa,bb) {return(sum((aa[x,,]-bb)^2))},
                   aa=All.x,bb=as.matrix(f.updated)%*%c.updated%*%(t(B)))

            mse=mse+mean(res)}  else {mse=0; break;}
        }

        allmse=rbind(allmse,c(lambda1,lambda2,lambda3,mse))
      }

  large.in=which(allmse[,4]>0)
  
  if(length(large.in)>0) break;}

  tem.in=which.min(allmse[large.in,4]); tem.in=large.in[tem.in]


  lambda1.fin=allmse[tem.in,1]; lambda2.fin=allmse[tem.in,2];lambda3.fin=allmse[tem.in,3];

  #print(allmse[tem.in,])

  tem.res=MTA00(x=All.x,Laplacian.matrix=Laplacian.matrix,timevec=timevec,
              lambda1=lambda1.fin,lambda2=lambda2.fin,lambda3=lambda3.fin)

  c.updated=tem.res[[1]];f.updated=tem.res[[2]]

  return(list(c.updated,f.updated))
}

