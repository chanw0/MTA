
########### Spline coefficients

### (a) in step 2
trend.coef.func = function(f,x,B, Om, lambda1){

  N=dim(x)[1]

  ## Ginv instead of solve
  theta = Ginv(N*t(B)%*%B+lambda1*Om)%*%(t(B))%*%(t(apply(x, c(2,3), sum)))%*%f

  return(as.vector(theta))
}

##### (b) in step 2

####### score

soft <- function(vec,lam){
  return(sign(vec)*pmax(0,abs(vec)-lam))
}

score.func=function(c,x,B, Laplacian.matrix, lambda2,lambda3) {

  N=dim(x)[1]; P=dim(x)[2]

  tem.eq1=t(c)%*%(t(B))%*%(t(apply(x, c(2,3), sum)))

  if(is.null(Laplacian.matrix)) tem.eq2=Ginv(diag(c(N*t(c)%*%(t(B))%*%B%*%c),P)) else tem.eq2=
    Ginv(diag(c(N*t(c)%*%(t(B))%*%B%*%c),P)+lambda3*Laplacian.matrix)

  # tem.eq2=solve(N*kronecker(t(c)%*%(t(B))%*%B%*%c, diag(1,P))+lambda3*Laplacian.matrix)
  # lambda2=BinarySearch(tem.eq1,sumabs)

  if(lambda2==0) tem.eq1=tem.eq1 else tem.eq1=soft(tem.eq1,lambda2/2)

  ff=tem.eq1 %*% tem.eq2

  return(as.vector(ff))

}

################
MTA00=function(x,Laplacian.matrix,timevec,lambda1,lambda2,lambda3,num.iter=100)
{

  N=dim(x)[1]; P=dim(x)[2]; T=dim(x)[3]

  if(is.null(timevec)) timevec=1:T

  BS = create.bspline.basis(timevec, norder=4)
  B = getbasismatrix(timevec, BS)
  Omega = getbasispenalty(BS)


  qq=apply(x,c(2,3),mean)
  qq=svd(cor(t(qq)))

  f.updated=qq$u[,1]

  ff=f.updated; cc=rep(0,T+2)
  for(i in 1:num.iter)
  {
   # print(i);
    c.updated=trend.coef.func(f.updated,x,B, Omega, lambda1)
    f.updated=score.func(c.updated,x,B, Laplacian.matrix,lambda2,lambda3)

    if(sum(abs(f.updated))!=0) {

      f.updated=f.updated/sqrt(sum(f.updated^2))

      c.updated=round(c.updated,4)
      f.updated=round(f.updated,4)

      cc=rbind(cc,c.updated)
      ff=rbind(ff,f.updated)

     #print(c.updated);print(f.updated)

    threshold=abs(c(ff[i+1,]-ff[i,],cc[i+1,]-cc[i,],(ff[i+1,ff[i,]!=0]-ff[i,ff[i,]!=0])/ff[i,ff[i,]!=0],
                    (cc[i+1,cc[i,]!=0]-cc[i,cc[i,]!=0])/cc[i,cc[i,]!=0]))

    if(max(threshold)<10^(-4)) break;

    } else break;}


  return(list(c.updated,f.updated))
}

